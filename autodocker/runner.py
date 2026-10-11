#!/usr/bin/env python3
"""AutoDocker: structure-based virtual screening with AutoDock Vina.

The pipeline runs five stages:

  1. Ligand library preparation: local SDF/MOL2/PDB/PDBQT files, or a fixed
     list of approved drugs downloaded from PubChem; optional Lipinski
     rule-of-five filter and MMFF94 minimization (SDF inputs), conversion to
     PDBQT with Gasteiger charges (Open Babel).
  2. Receptor preparation: chain selection, cleaning, hydrogens at a chosen
     pH, PDBQT conversion, and a docking box from fpocket pockets, a box
     around the whole receptor (--no-fpocket) or explicit coordinates.
     Optionally keeps waters near the site, metal ions and cofactors, and
     builds flexible side chains.
  3. Docking with Vina (or SMINA with --smina-only), seeded and resumable;
     optional SMINA (Vinardo) rescoring of Vina's top poses for a consensus
     ranking.
  4. Analysis: pose-consistency metrics, optional RMSD clustering of poses.
  5. Outputs: ranking.csv, text reports, an HTML report and
     run_manifest.json (versions, settings and checksums of the tools used).

Run ``runner.py --help`` for the options.
"""
import argparse
import hashlib
import json
import logging
import math
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
from datetime import datetime, timezone
from typing import Dict, List, Optional

from vspipeline import defaults
from vspipeline import outputs

_FALLBACK_VERSION = "2.0.0"


def _read_version(base_dir: Optional[str] = None) -> str:
    """Version string: autodocker/VERSION if present, else $AUTODOCKER_VERSION,
    else the version this file was released with."""
    base_dir = base_dir or os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base_dir, "VERSION")
    try:
        with open(path, encoding="utf-8") as f:
            version = f.read().strip()
        if version:
            return version
    except OSError:
        pass
    return os.environ.get("AUTODOCKER_VERSION", "").strip() or _FALLBACK_VERSION


__version__ = _read_version()

# When executed as a script (__name__ == "__main__"), alias this module as
# "runner" so the vspipeline submodules' ``import runner`` binds THIS running
# instance instead of importing a second copy of the file. That keeps the
# tool-path globals below (and main()'s mutations of them) shared everywhere.
sys.modules.setdefault("runner", sys.modules[__name__])

logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)

# Exit codes.
EXIT_OK = 0
EXIT_FAILURE = 1            # the pipeline failed
EXIT_USAGE = 2              # invalid options/environment or missing tools
EXIT_TOO_MANY_FAILURES = 3  # outputs written, but too many ligands failed to dock
EXIT_INTERRUPTED = 130      # Ctrl-C (128 + SIGINT)
EXIT_TERMINATED = 143       # SIGTERM, e.g. docker stop (128 + SIGTERM)


class TerminatedBySignal(KeyboardInterrupt):
    """Raised in the main process when it receives SIGTERM. A subclass of
    KeyboardInterrupt, so a terminated run takes the interrupted path:
    workers and engine processes are stopped, the checkpoint is flushed and
    run_manifest.json records status "interrupted"."""

DEFAULT_SEED = defaults.SEED
DEFAULT_COMMAND_TIMEOUT = defaults.COMMAND_TIMEOUT  # seconds
# Residues within this distance (Angstrom) of the box center are candidates
# for --auto-flexible.
AUTO_FLEX_RADIUS = defaults.AUTO_FLEX_RADIUS
# Box edge lengths (Angstrom) outside this range get a warning.
TYPICAL_GRID_SIZE_RANGE = (15.0, 40.0)
# Timeout (seconds) for the "--version" probes of the external tools.
VERSION_PROBE_TIMEOUT = 30


class ConfigError(ValueError):
    """Invalid option or environment value (reported without a traceback)."""


def find_tool(*names: str) -> Optional[str]:
    """Return the path of the first *name* found on PATH, or None.

    Returns None when nothing is found instead of a placeholder name so
    callers can reliably distinguish 'found' from 'not found'.
    """
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    return None


def default_seed() -> int:
    """Seed used when --seed is not given: $VS_SEED if set, else 42.

    Raises ConfigError for a non-integer $VS_SEED. It is evaluated after
    argument parsing, so a bad value never breaks --help or --version.
    """
    value = os.environ.get("VS_SEED", "").strip()
    if not value:
        return DEFAULT_SEED
    try:
        return int(value)
    except ValueError:
        raise ConfigError(f"VS_SEED must be an integer, got {value!r}") from None


# Problems with environment variables read at import time; reported by
# main() after argument parsing (so --version still works).
_ENV_ERRORS: List[str] = []


def _env_command_timeout() -> int:
    raw = os.environ.get("VS_COMMAND_TIMEOUT", "").strip()
    if not raw:
        return DEFAULT_COMMAND_TIMEOUT
    try:
        value = int(raw)
        if value <= 0:
            raise ValueError
        return value
    except ValueError:
        _ENV_ERRORS.append(
            f"VS_COMMAND_TIMEOUT must be a positive integer (seconds), got {raw!r}")
        return DEFAULT_COMMAND_TIMEOUT


# Tool paths (prefer modern Vina for consistency).
OBABEL = find_tool("obabel", "OpenBabel", "OpenBabel.exe")
VINA_PRIMARY = find_tool("vina", "vina.exe", "qvina02", "qvina2", "qvina")
# The fallback engine equals the primary (no fallback) unless --qvina-fallback
# is given; main() then sets it to QuickVina 2.
VINA_FALLBACK = VINA_PRIMARY
QVINA = find_tool("qvina02", "qvina2")
# Deprecated environment switch for --qvina-fallback (versions <= 1.1.x).
_QVINA_FALLBACK_ENV = "VS_ENABLE_QVINA_FALLBACK"
# Optional engines and helpers.
SMINA = find_tool("smina", "smina.exe")
FPOCKET = find_tool("fpocket")
# reduce and AutoDock4/AutoGrid4 are not used by the pipeline (receptor
# hydrogens come from Open Babel); the paths are kept for the real-tool tests.
REDUCE = find_tool("reduce")
AUTODOCK = find_tool("autodock4", "autodock")
AUTOGRID = find_tool("autogrid4", "autogrid")
COMMAND_TIMEOUT = _env_command_timeout()


def _require_tools(smina_only: bool = False, consensus: bool = False,
                   fpocket: bool = False) -> None:
    """Fail when a tool the requested run needs is missing.

    Called from main() after argument parsing, so --help/--version work and
    --vina-bin can supply a docking binary that is not on PATH. ``fpocket``:
    the run detects pockets (see _uses_fpocket).
    """
    needed = [("obabel", OBABEL)]
    if not smina_only:
        needed.append(("vina", VINA_PRIMARY))
    if smina_only or consensus:
        needed.append(("smina", SMINA))
    if fpocket:
        needed.append(("fpocket", FPOCKET))
    missing = [name for name, path in needed if path is None]
    if missing:
        hint = "Install them or put them on PATH (Vina can also be given with --vina-bin)."
        if "fpocket" in missing:
            hint += (" Without fpocket, pass --no-fpocket (box around the whole receptor) "
                     "or give the box with --grid-center and --grid-size.")
        raise RuntimeError(
            f"Required tool(s) not found in PATH: {', '.join(missing)}. {hint}")


def _uses_fpocket(args) -> bool:
    """True when the run detects pockets with fpocket: not with --no-fpocket,
    and not when --grid-center and --grid-size give the whole box."""
    return not (args.no_fpocket
                or (args.grid_center is not None and args.grid_size is not None))


def _log_tool_availability() -> None:
    logger.debug(f"SMINA: {SMINA or 'not found (--consensus and --smina-only unavailable)'}")
    logger.debug(f"fpocket: {FPOCKET or 'not found'}")

# ---------------------------------------------------------------------------
# Pipeline internals (vspipeline/*).
# These imports re-export the names the CLI, tests and scripts use through
# ``runner.<name>`` (listed in __all__); the vspipeline modules resolve shared
# state (tool paths, ``run``/``run_obabel``) through THIS module at call time,
# so monkeypatching runner.<name> reaches them. Docking workers receive the
# engine paths and timeout explicitly (docking.EngineConfig, built from these
# globals when docking starts).
# ---------------------------------------------------------------------------

from vspipeline.net import _http_get_bytes  # noqa: E402
from vspipeline.report import (  # noqa: E402
    NGL_CDN_URL,
    _PDBQT_KEEP_PREFIXES,
    _PDBQT_STRIP_PREFIXES,
    _build_viewer_section,
    _fetch_ngljs,
    _pdbqt_to_pdb_string,
    calculate_rmsd,
    cluster_poses,
    extract_coordinates_from_pdbqt,
    generate_html_report,
)
from vspipeline.pocket import (  # noqa: E402
    COFACTOR_RESNAMES,
    METAL_ELEMENTS,
    VINA_METAL_TYPE_MAP,
    WATER_RESNAMES,
    _append_hetatm_to_receptor,
    _hetatm_to_pdbqt_line,
    _iter_hetatm_residues,
    _max_pdbqt_serial,
    _parse_hetatm,
    detect_cofactors,
    detect_metal_ions,
    detect_water_molecules,
)
from vspipeline.flexible import (  # noqa: E402
    FLEX_BACKBONE_NAMES,
    _ATOMIC_RADII,
    _BOND_TOLERANCE,
    _build_flex_residue_block,
    _parse_residue_key,
    _pdbqt_atom_element,
    _remove_residues_from_pdbqt,
    _residue_atoms_from_pdbqt,
    _resnum_matches,
    build_flexible_residue_pdbqt,
    detect_flexible_residues,
    missing_residues,
)
from vspipeline.consensus import (  # noqa: E402
    _run_smina_scoring,
    rescore_with_smina,
    run_smina_docking,
    scorer_agreement_spearman,
    smina_dock_all,
    write_consensus_outputs,
)
from vspipeline.charges import (  # noqa: E402
    _ensure_pdbqt_has_charges,
    _extract_pdbqt_charge,
    _fix_pdbqt_charges,
    _pdbqt_has_atoms,
    _pdbqt_has_nonzero_charges,
    _sanitize_receptor_pdbqt,
)
from vspipeline.admet import (  # noqa: E402
    ADMETFilter,
    LipinskiFilter,
    _is_valid_sdf,
    _obabel_descriptors,
)
from vspipeline.checkpoint import (  # noqa: E402
    DockingCheckpoint,
    OutputDirLockedError,
    acquire_lock,
    release_lock,
)
from vspipeline.library import LibraryManager  # noqa: E402
from vspipeline.utils import (  # noqa: E402
    _parse_chain_selection,
    available_cpus,
    ligand_name,
    _safe_ligand_id,
    get_chains,
    parse_pdb_coords,
    run,
    run_obabel,
)
from vspipeline.protein import ProteinPreparation  # noqa: E402
from vspipeline import protein as protein_module  # noqa: E402
from vspipeline.docking import (  # noqa: E402
    AFFINITY_FLOOR,
    DockJob,
    EngineConfig,
    Status,
    _VINA_CPU_FLAG_CACHE,
    _build_vina_command,
    _calculate_metrics,
    _calculate_protocol_simscore,
    _detect_vina_type,
    _engine_identity,
    _extract_score,
    _log_tool_version,
    _mode_table_lines,
    _parse_grid_config,
    _parse_grid_triplet,
    _parse_vina_modes,
    _tool_label,
    _vina_cpu_flag,
    affinity_in_window,
    dock_all,
    dock_ligand,
    parse_pdbqt_results,
)
from vspipeline.results import (  # noqa: E402
    ResultsAnalyzer,
    _score_csv,
    _score_sort_key,
    _score_text,
)

# Re-exported names (see above); listed so linters treat them as used.
__all__ = [
    # net, report
    "_http_get_bytes", "NGL_CDN_URL", "_PDBQT_KEEP_PREFIXES", "_PDBQT_STRIP_PREFIXES",
    "_build_viewer_section", "_fetch_ngljs", "_pdbqt_to_pdb_string", "calculate_rmsd",
    "cluster_poses", "extract_coordinates_from_pdbqt", "generate_html_report",
    # pocket
    "COFACTOR_RESNAMES", "METAL_ELEMENTS", "VINA_METAL_TYPE_MAP", "WATER_RESNAMES",
    "_append_hetatm_to_receptor", "_hetatm_to_pdbqt_line", "_iter_hetatm_residues",
    "_max_pdbqt_serial", "_parse_hetatm", "detect_cofactors", "detect_metal_ions",
    "detect_water_molecules",
    # flexible
    "FLEX_BACKBONE_NAMES", "_ATOMIC_RADII", "_BOND_TOLERANCE", "_build_flex_residue_block",
    "_parse_residue_key", "_pdbqt_atom_element", "_remove_residues_from_pdbqt",
    "_residue_atoms_from_pdbqt", "_resnum_matches", "build_flexible_residue_pdbqt",
    "detect_flexible_residues",
    # consensus
    "_run_smina_scoring", "rescore_with_smina", "run_smina_docking",
    "scorer_agreement_spearman", "smina_dock_all", "write_consensus_outputs",
    # charges
    "_ensure_pdbqt_has_charges", "_extract_pdbqt_charge", "_fix_pdbqt_charges",
    "_pdbqt_has_atoms", "_pdbqt_has_nonzero_charges", "_sanitize_receptor_pdbqt",
    # admet (Lipinski filter; ADMETFilter is the deprecated alias)
    "ADMETFilter", "LipinskiFilter", "_is_valid_sdf", "_obabel_descriptors",
    # checkpoint, library, protein
    "DockingCheckpoint", "LibraryManager", "ProteinPreparation", "OutputDirLockedError",
    "acquire_lock", "release_lock",
    # utils
    "_parse_chain_selection", "available_cpus", "ligand_name", "_safe_ligand_id",
    "get_chains", "parse_pdb_coords", "run", "run_obabel",
    # docking
    "AFFINITY_FLOOR", "DockJob", "EngineConfig", "Status", "_VINA_CPU_FLAG_CACHE",
    "_build_vina_command", "_calculate_metrics", "_calculate_protocol_simscore",
    "_detect_vina_type", "_engine_identity", "_extract_score", "_log_tool_version",
    "_mode_table_lines", "_parse_grid_config", "_parse_grid_triplet", "_parse_vina_modes",
    "_tool_label", "_vina_cpu_flag", "affinity_in_window", "dock_all", "dock_ligand",
    "parse_pdbqt_results",
    # results
    "ResultsAnalyzer", "_score_csv", "_score_sort_key", "_score_text",
]

# =============================
# COMMAND LINE
# =============================

DEFAULT_BINDING_MODES = defaults.BINDING_MODES
DEFAULT_ENERGY_RANGE = defaults.ENERGY_RANGE
DEFAULT_MAX_FAILURE_FRACTION = defaults.MAX_FAILURE_FRACTION

# Options whose explicit use _validate_args checks: the --library custom
# filters, --max-compounds (fda/custom) and --grid-center-frame (needs
# --grid-center); _parse_args detects whether they were given.
_UNSET = object()
_CUSTOM_ONLY_DESTS = ("mw_min", "mw_max", "logp_max")
_EXPLICIT_DESTS = _CUSTOM_ONLY_DESTS + ("max_compounds", "grid_center_frame")

# Seeds accepted by --seed/$VS_SEED. Vina, QuickVina 2 and SMINA read the
# seed as a 32-bit signed int (they reject larger values); Vina takes 0 to
# mean "pick a random seed" and RDKit's ETKDG (3D coordinates of 2D ligands)
# takes -1 to mean random and rejects other negative seeds. So only positive
# 32-bit seeds give a reproducible run.
SEED_RANGE = (1, 2**31 - 1)
# pH range accepted by --ph.
PH_RANGE = (0.0, 14.0)


class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                     argparse.RawDescriptionHelpFormatter):
    """Keeps the description's line breaks and appends option defaults
    (except None/False, which mean "not set")."""

    def _get_help_string(self, action):
        if (action.default is None or action.default is False or action.default == []
                or isinstance(action, (argparse._StoreFalseAction, argparse._StoreTrueAction))):
            return action.help
        return super()._get_help_string(action)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="runner.py",
        description=(
            f"AutoDocker {__version__}: structure-based virtual screening with AutoDock Vina.\n"
            "Stages: [1/5] ligand preparation, [2/5] receptor and box preparation,\n"
            "[3/5] docking (Vina, or SMINA with --smina-only; optional SMINA rescoring),\n"
            "[4/5] pose analysis, [5/5] reports and run_manifest.json.\n"
            "Exit codes: 0 success, 1 pipeline error, 2 invalid options/environment,\n"
            "missing tools or an output directory in use by another run, 3 more ligands\n"
            "failed than --max-failure-fraction allows (outputs are still written),\n"
            "130 interrupted (Ctrl-C), 143 terminated (SIGTERM)."),
        formatter_class=_HelpFormatter,
        # Abbreviated long options (--lig for --ligands) are not accepted:
        # the container entry point (entry.sh) must recognise -r/-l/-o and
        # --library exactly as the runner does.
        allow_abbrev=False,
    )

    io = parser.add_argument_group("input/output")
    io.add_argument("-r", "--receptor", required=True, help="Receptor (protein) PDB file")
    io.add_argument("-l", "--ligands", default=None,
                    help="Ligand file or directory (searched recursively) of .sdf, .mol2, "
                    ".pdb and/or .pdbqt files; one molecule per file. Files sharing a name "
                    "are resolved by precedence SDF > PDBQT > MOL2 > PDB (the others are "
                    "reported as superseded). Required with --library local; not used by "
                    "the PubChem libraries")
    io.add_argument("-o", "--output", required=True,
                    help="Output directory (if it lies inside -l it is skipped when "
                    "collecting ligands)")

    lib = parser.add_argument_group("ligand library")
    lib.add_argument("--library", choices=["fda", "custom", "local"], default="local",
                     help="local: the files of -l; fda: a fixed list of "
                     f"{len(LibraryManager.FDA_DRUG_NAMES)} well-known "
                     "approved drugs downloaded from PubChem (not the full FDA list); "
                     "custom: the same list filtered by --mw-min/--mw-max/--logp-max")
    lib.add_argument("--mw-min", type=float, default=defaults.CUSTOM_MW_MIN,
                     help="--library custom only: minimum molecular weight (Da)")
    lib.add_argument("--mw-max", type=float, default=defaults.CUSTOM_MW_MAX,
                     help="--library custom only: maximum molecular weight (Da); "
                     "replaces Lipinski's MW <= 500")
    lib.add_argument("--logp-max", type=float, default=defaults.CUSTOM_LOGP_MAX,
                     help="--library custom only: maximum logP (Open Babel's "
                     "Wildman-Crippen logP, the model used everywhere); replaces "
                     "Lipinski's logP <= 5")
    lib.add_argument("--max-compounds", type=int, default=LibraryManager.MAX_COMPOUNDS,
                     help="--library fda/custom only: number of PubChem compounds to "
                     "prepare (the first ones of the list that pass)")
    lib.add_argument("--no-lipinski", action="store_true",
                     help="Skip the Lipinski rule-of-five filter. The filter (SDF and "
                     "PubChem inputs only) rejects a ligand with ANY violation of MW <= "
                     "500, logP <= 5, H-bond donors <= 5, acceptors <= 10 (Open Babel "
                     "descriptors HBD/HBA2); it is not an ADMET model")
    lib.add_argument("--no-admet", action="store_true",
                     help="Deprecated alias of --no-lipinski")
    lib.add_argument("--no-minimize", action="store_true",
                     help="Skip MMFF94 ligand minimization (applies to SDF and PubChem "
                     "ligands only; MOL2/PDB/PDBQT keep their input coordinates, except "
                     "2D MOL2/PDB, which get seeded RDKit 3D coordinates)")
    lib.add_argument("--minimize-steps", type=int, default=defaults.MINIMIZE_STEPS,
                     help="MMFF94 minimization steps (SDF and PubChem ligands only)")

    rec = parser.add_argument_group("receptor and docking box")
    rec.add_argument("--chain", default=None,
                     help="Protein chain(s), e.g. A, A,B, or all (IDs are case-sensitive). "
                     "Default: the first chain in sorted order (a warning lists the others)")
    rec.add_argument("--ph", type=float, default=defaults.PH,
                     help="pH at which receptor hydrogens are added (Open Babel -p); "
                     f"{PH_RANGE[0]:g} to {PH_RANGE[1]:g}")
    rec.add_argument("--keep-hetero", action="store_true",
                     help="Keep HETATM records of the selected chains (single-atom metal "
                     "ions are added with Vina metal types; other elements without Open "
                     "Babel Gasteiger parameters stop the run)")
    rec.add_argument("--pockets",
                     help="fpocket pocket number(s), e.g. 1 or 2,3. Default: best fpocket "
                     "score (druggability breaks ties). Not combinable with --no-fpocket or "
                     "with --grid-center and --grid-size together (no pocket is used then)")
    rec.add_argument("--padding", type=float, default=defaults.PADDING,
                     help="Box padding around the selected pocket(s), Angstrom")
    rec.add_argument("--no-fpocket", action="store_true",
                     help="Skip fpocket; dock in a box around the whole receptor (heavy-atom "
                     f"bounding box + padding, each edge at least "
                     f"{protein_module.MIN_RECEPTOR_BOX_SIZE:g} Angstrom)")
    rec.add_argument("--grid-center", default=None,
                     help="Explicit box center 'X,Y,Z' (Angstrom), e.g. 12.3,4.5,6.7; "
                     "overrides auto-detection (frame: see --grid-center-frame). With "
                     "--grid-size as well, fpocket is not run")
    rec.add_argument("--grid-center-frame", choices=["input", "docking"], default="input",
                     help="Frame of --grid-center: 'input' (receptor PDB coordinates) or "
                     "'docking' (the centered frame of receptor.pdbqt and grid.conf; "
                     "versions <= 1.1.4 used this frame). Needs --grid-center")
    rec.add_argument("--grid-size", default=None,
                     help="Explicit box size 'SX,SY,SZ' (Angstrom), e.g. 22.5,22.5,22.5; "
                     "overrides auto-detection")
    rec.add_argument("--keep-waters", action="store_true",
                     help="Keep water molecules near the binding site in the receptor")
    rec.add_argument("--water-distance", type=float, default=defaults.WATER_DISTANCE,
                     help="Distance from the box center within which waters are kept (Angstrom)")
    rec.add_argument("--detect-metals", action="store_true",
                     help="Keep metal ions in the receptor (Vina metal atom types)")
    rec.add_argument("--detect-cofactors", action="store_true",
                     help="Keep known cofactors in the receptor")

    flex = parser.add_argument_group("flexible receptor")
    flex.add_argument("--flexibility", type=int, choices=range(0, 11), metavar="{0..10}",
                      help="On/off switch kept for compatibility: any value > 0 enables "
                      "flexible side chains (the number itself has no effect), 0 forces "
                      "rigid docking. Implied by --flexible-residues or --auto-flexible")
    flex.add_argument("--flexible-residues",
                      help="Flexible residues, e.g. 'A45,A78,B102' (chain + residue "
                      "number [+ insertion code]; chain IDs are case-sensitive; "
                      "'chain:number' for any chain, e.g. 1:45; a blank chain is A: A45 "
                      "or 45). Every residue must exist in the docked chain(s) and be "
                      "buildable as a flexible side chain, else the run stops. Not "
                      "combinable with --auto-flexible")
    flex.add_argument("--auto-flexible", type=int,
                      help=f"Make the N residues closest to the box center (within "
                      f"{AUTO_FLEX_RADIUS:g} A) flexible; if none qualifies, docking is "
                      "rigid (with a warning; run_manifest.json records the residues used)")

    dock = parser.add_argument_group("docking")
    dock.add_argument("-p", "--processes", type=int, default=defaults.PROCESSES,
                      help="Parallel docking processes: N > 0 uses N; 0 or -1 uses all "
                      "available CPUs; -N (N > 1) uses available CPUs / N")
    dock.add_argument("--exhaustiveness", type=int, default=defaults.EXHAUSTIVENESS,
                      help="Vina exhaustiveness (>= 1; 8 is Vina's default, higher is "
                      "more thorough and slower)")
    dock.add_argument("--binding-modes", type=int, default=DEFAULT_BINDING_MODES,
                      help="Maximum number of binding modes Vina reports")
    dock.add_argument("--energy-range", type=float, default=DEFAULT_ENERGY_RANGE,
                      help="Vina energy range of the reported modes (kcal/mol)")
    dock.add_argument("--seed", type=int, default=None,
                      help=f"Random seed, {SEED_RANGE[0]} to {SEED_RANGE[1]} (Vina treats 0 "
                      f"as 'random'; default: $VS_SEED, else {DEFAULT_SEED})")
    dock.add_argument("--threads", type=int, default=None,
                      help="CPU threads per Vina process, passed as --cpu (default: all "
                      "cores, or cores divided by -p when -p > 1)")
    dock.add_argument("--min-valid-affinity", type=float, default=defaults.MIN_VALID_AFFINITY,
                      help="Upper bound of the accepted affinity window (kcal/mol): scores "
                      f">= this, or < {AFFINITY_FLOOR:g}, are reported as OUT_OF_WINDOW; "
                      f"must be > {AFFINITY_FLOOR:g}")
    dock.add_argument("--vina-bin",
                      help="Docking binary to use instead of the auto-detected Vina/QuickVina")
    dock.add_argument("--qvina-fallback", action="store_true",
                      help="Re-dock a ligand with QuickVina 2 when Vina fails on it (crash "
                      "or timeout). The two engines' scores then share one ranking; "
                      "ranking.csv's Engine column shows which engine scored each ligand. "
                      f"Replaces the deprecated ${_QVINA_FALLBACK_ENV}=1")
    dock.add_argument("--vina-extra", action="append", default=[],
                      help="Extra argument string passed to the docking binary; repeat as needed")
    dock.add_argument("--timeout", type=int, default=COMMAND_TIMEOUT,
                      help="Per-command timeout in seconds; the default is "
                      f"$VS_COMMAND_TIMEOUT if set, else {DEFAULT_COMMAND_TIMEOUT}")
    dock.add_argument("--no-resume", action="store_true",
                      help="Start fresh (ignore the checkpoint)")
    dock.add_argument("--max-failure-fraction", type=float,
                      default=DEFAULT_MAX_FAILURE_FRACTION,
                      help="Exit with code 3 (after writing all outputs) when more than this "
                      "fraction of the ligands failed: (ligands that could not be prepared "
                      "+ ligands whose docking failed) / (ligands prepared or not, without "
                      "superseded files). Scores outside the accepted window are not "
                      "failures; 1.0 never fails, not even when every ligand failed")
    dock.add_argument("--consensus", action="store_true",
                      help="Rescore Vina's top pose of each ligand with SMINA (Vinardo "
                      "scoring) and write consensus_ranking.csv and consensus_summary.json "
                      "(Spearman rho between the two scores)")
    dock.add_argument("--smina-only", action="store_true",
                      help="Dock with SMINA (default scoring) instead of Vina (1 mode per "
                      "ligand, no checkpoint/resume; not combinable with --consensus, "
                      "--vina-bin, --vina-extra or non-default --binding-modes/--energy-range)")

    out = parser.add_argument_group("analysis and reports")
    out.add_argument("--top-n", type=int, default=defaults.TOP_N,
                     help="Number of top-ranked compounds highlighted in Top_hits.txt and in the HTML report's top-hits table and 3D viewer (ranking.csv always lists every ligand)")
    out.add_argument("--html-report", dest="html_report", action="store_true",
                     default=True,
                     help="Write results_report.html (on by default; the flag is accepted "
                     "for compatibility)")
    out.add_argument("--no-html-report", dest="html_report", action="store_false",
                     help="Do not write the HTML report")
    out.add_argument("--cluster-poses", action="store_true",
                     help="Write pose_clusters.json (greedy RMSD clustering of each "
                     "ligand's poses). The HTML report clusters poses with "
                     "--rmsd-threshold whether or not this is given")
    out.add_argument("--rmsd-threshold", type=float, default=defaults.RMSD_THRESHOLD,
                     help="RMSD clustering threshold (Angstrom)")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Debug logging (all pipeline modules)")
    parser.add_argument("--version", action="version",
                        version=f"AutoDocker {__version__}")
    return parser


def _parse_args(parser: argparse.ArgumentParser, argv: List[str]):
    """Parse ``argv``; return (args, set of _EXPLICIT_DESTS given explicitly).

    argparse only fills defaults for attributes missing from the namespace,
    so a sentinel there shows whether an option was given.
    """
    namespace = argparse.Namespace(**{dest: _UNSET for dest in _EXPLICIT_DESTS})
    args = parser.parse_args(argv, namespace=namespace)
    explicit = set()
    for dest in _EXPLICIT_DESTS:
        if getattr(args, dest) is _UNSET:
            setattr(args, dest, parser.get_default(dest))
        else:
            explicit.add(dest)
    return args, explicit


# Float options; argparse's float() also accepts nan and inf.
_FLOAT_DESTS = ("ph", "padding", "water_distance", "energy_range", "min_valid_affinity",
                "max_failure_fraction", "mw_min", "mw_max", "logp_max", "rmsd_threshold")


def _option(dest: str) -> str:
    return "--" + dest.replace("_", "-")


def _flexible_residue_keys(value: str) -> List[str]:
    """The residue keys of a --flexible-residues value ('A45, A78' -> [A45, A78])."""
    return [r.strip() for r in value.split(",") if r.strip()]


def _check_output_dir(path: str) -> None:
    """-o must be an existing writable directory or creatable as one (main()
    creates it once every other check has passed)."""
    if not path or not path.strip():
        raise ConfigError("-o/--output must not be empty")
    if os.path.exists(path):
        if not os.path.isdir(path):
            raise ConfigError(f"-o/--output {path!r} exists and is not a directory")
        if not os.access(path, os.W_OK | os.X_OK):
            raise ConfigError(f"-o/--output {path!r} is not writable")
        return
    parent = os.path.dirname(os.path.abspath(path))
    while not os.path.exists(parent):
        parent = os.path.dirname(parent)
    if not os.path.isdir(parent):
        raise ConfigError(f"-o/--output {path!r} cannot be created: {parent} is not a "
                          "directory")
    if not os.access(parent, os.W_OK | os.X_OK):
        raise ConfigError(f"-o/--output {path!r} cannot be created: {parent} is not "
                          "writable")


def _check_receptor_selection(args) -> None:
    """--chain and --flexible-residues against the receptor PDB (a cheap read,
    before any ligand is prepared): the chains must exist, and every flexible
    residue must have atoms in the docked chain(s). Uses the same rules as
    receptor preparation (case-sensitive chain IDs, first chain by default)
    and flexible-residue building (chain + residue number and insertion code)."""
    try:
        available = get_chains(args.receptor)
    except OSError as e:
        raise ConfigError(f"-r/--receptor: cannot read {args.receptor}: {e}") from None
    except ValueError as e:  # utils.ChainAmbiguityError: blank chain next to chain A
        raise ConfigError(f"-r/--receptor: {e}") from None
    if args.chain is not None:
        try:
            selected = _parse_chain_selection(args.chain, available)
        except ValueError as e:
            raise ConfigError(f"--chain {args.chain!r}: {e}") from None
    else:
        selected = available[:1]

    if not args.flexible_residues or (args.flexibility is not None and args.flexibility <= 0):
        return
    try:
        lines, _ = protein_module.read_pdb_records(args.receptor)
    except (OSError, UnicodeDecodeError) as e:
        raise ConfigError(f"-r/--receptor: cannot read {args.receptor}: {e}") from None
    missing = missing_residues(lines, selected,
                               _flexible_residue_keys(args.flexible_residues))
    if missing:
        raise ConfigError(
            f"--flexible-residues: no such residue in the docked chain(s) "
            f"{', '.join(selected)} of {args.receptor}: {', '.join(missing)}")


def _validate_args(args, explicit) -> None:
    """Check option values and combinations; raise ConfigError on the first
    problem. Resolves --seed, --processes and the deprecated --no-admet.
    Reads the receptor's chain list, but writes nothing."""
    if _ENV_ERRORS:
        raise ConfigError("; ".join(_ENV_ERRORS))
    for dest in _FLOAT_DESTS:
        value = getattr(args, dest)
        if not math.isfinite(value):
            raise ConfigError(f"{_option(dest)} must be a finite number, got {value}")
    seed_source = "--seed"
    if args.seed is None:
        args.seed = default_seed()
        seed_source = "$VS_SEED"
    if not SEED_RANGE[0] <= args.seed <= SEED_RANGE[1]:
        raise ConfigError(
            f"{seed_source} must be between {SEED_RANGE[0]} and {SEED_RANGE[1]}, got "
            f"{args.seed} (the docking engines take a 32-bit seed, Vina treats 0 as "
            "'random' and RDKit rejects negative seeds)")
    if args.no_admet:
        logger.warning("[!] --no-admet is deprecated; use --no-lipinski")
        args.no_lipinski = True

    if not os.path.isfile(args.receptor):
        raise ConfigError(f"-r/--receptor: no such file: {args.receptor}")
    if args.library == "local":
        if not args.ligands:
            raise ConfigError("-l/--ligands is required with --library local")
        if not os.path.exists(args.ligands):
            raise ConfigError(f"-l/--ligands: no such file or directory: {args.ligands}")
    custom_only = sorted(d for d in explicit if d in _CUSTOM_ONLY_DESTS)
    if args.library != "custom" and custom_only:
        options = ", ".join(_option(d) for d in custom_only)
        raise ConfigError(
            f"{options} only apply to --library custom (with --library {args.library} "
            "the Lipinski cut-offs are used; see --no-lipinski)")
    if args.library == "local" and "max_compounds" in explicit:
        raise ConfigError("--max-compounds only applies to --library fda/custom")
    if args.mw_min < 0 or args.mw_min > args.mw_max:
        raise ConfigError(f"--mw-min ({args.mw_min:g}) must be >= 0 and <= --mw-max "
                          f"({args.mw_max:g})")

    positive_ints = [("--exhaustiveness", args.exhaustiveness),
                     ("--binding-modes", args.binding_modes),
                     ("--timeout", args.timeout),
                     ("--minimize-steps", args.minimize_steps),
                     ("--top-n", args.top_n),
                     ("--max-compounds", args.max_compounds)]
    if args.threads is not None:
        positive_ints.append(("--threads", args.threads))
    if args.auto_flexible is not None:
        positive_ints.append(("--auto-flexible", args.auto_flexible))
    for option, value in positive_ints:
        if value < 1:
            raise ConfigError(f"{option} must be >= 1, got {value}")
    for option, value in (("--energy-range", args.energy_range),
                          ("--rmsd-threshold", args.rmsd_threshold),
                          ("--water-distance", args.water_distance)):
        if not value > 0:
            raise ConfigError(f"{option} must be > 0, got {value:g}")
    if args.padding < 0:
        raise ConfigError(f"--padding must be >= 0, got {args.padding:g}")
    if not PH_RANGE[0] <= args.ph <= PH_RANGE[1]:
        raise ConfigError(f"--ph must be between {PH_RANGE[0]:g} and {PH_RANGE[1]:g}, "
                          f"got {args.ph:g}")
    if not args.min_valid_affinity > AFFINITY_FLOOR:
        raise ConfigError(
            f"--min-valid-affinity must be > {AFFINITY_FLOOR:g} kcal/mol (the lower end "
            f"of the accepted window), got {args.min_valid_affinity:g}: no score could "
            "be accepted")
    if not 0.0 <= args.max_failure_fraction <= 1.0:
        raise ConfigError("--max-failure-fraction must be between 0 and 1, got "
                          f"{args.max_failure_fraction:g}")
    for option, value in (("--grid-center", args.grid_center),
                          ("--grid-size", args.grid_size)):
        if value is not None:
            try:
                triplet = _parse_grid_triplet(value, option.lstrip("-"))
            except ValueError as e:
                raise ConfigError(str(e)) from None
            if option == "--grid-size" and min(triplet) <= 0:
                raise ConfigError(f"--grid-size must be positive, got {value}")

    # Options that would be accepted but have no effect.
    if "grid_center_frame" in explicit and args.grid_center is None:
        raise ConfigError("--grid-center-frame only applies to --grid-center")
    if args.pockets is not None:
        spec = args.pockets.strip()
        tokens = [t for t in re.split(r"[,\s]+", spec) if t]
        if spec.lower() != "auto" and (
                not tokens or not all(re.fullmatch(r"[0-9]+", t) and int(t) >= 1
                                      for t in tokens)):
            raise ConfigError(f"--pockets must list fpocket pocket numbers such as 1 or "
                              f"2,3, got {args.pockets!r}")
        if args.no_fpocket:
            raise ConfigError("--pockets cannot be combined with --no-fpocket (fpocket "
                              "is not run)")
        if args.grid_center is not None and args.grid_size is not None:
            raise ConfigError("--pockets cannot be combined with --grid-center and "
                              "--grid-size together (the box is given, fpocket is not run)")
    if args.flexible_residues is not None:
        if args.auto_flexible is not None:
            raise ConfigError("--flexible-residues and --auto-flexible cannot be combined "
                              "(name the residues, or let N of them be chosen)")
        keys = _flexible_residue_keys(args.flexible_residues)
        if not keys:
            raise ConfigError(f"--flexible-residues lists no residue, got "
                              f"{args.flexible_residues!r} (e.g. A45,A78)")
        invalid = [k for k in keys if _parse_residue_key(k) is None]
        if invalid:
            raise ConfigError(
                f"--flexible-residues: invalid residue key(s) {', '.join(invalid)} "
                "(expected chain + residue number, e.g. A45 or A100B; chain:number "
                "for any chain, e.g. 1:45; a blank chain is A: A45 or 45)")

    if args.vina_bin and not shutil.which(args.vina_bin):
        raise ConfigError(f"--vina-bin {args.vina_bin!r} does not exist or is not executable")

    if args.qvina_fallback:
        fallback = "--qvina-fallback"
        if os.environ.get(_QVINA_FALLBACK_ENV, "0") == "1":
            fallback += (f" (turned on by ${_QVINA_FALLBACK_ENV}=1 in the environment; "
                         "unset it)")
        if args.vina_bin:
            raise ConfigError(f"{fallback} cannot be combined with --vina-bin")
        if args.smina_only:
            raise ConfigError(f"{fallback} cannot be combined with --smina-only")
        if not QVINA:
            raise ConfigError(f"{fallback} needs QuickVina 2 (qvina02 or qvina2) on PATH")

    if args.smina_only:
        conflicts = []
        if args.consensus:
            conflicts.append("--consensus (SMINA would be compared with itself)")
        if args.vina_bin:
            conflicts.append("--vina-bin (SMINA is taken from PATH)")
        if args.vina_extra:
            conflicts.append("--vina-extra")
        if args.binding_modes != DEFAULT_BINDING_MODES:
            conflicts.append("--binding-modes (SMINA mode reports 1 mode)")
        if args.energy_range != DEFAULT_ENERGY_RANGE:
            conflicts.append("--energy-range")
        if conflicts:
            raise ConfigError("--smina-only cannot be combined with " + ", ".join(conflicts))

    _check_receptor_selection(args)
    _check_output_dir(args.output)

    if args.processes in (0, -1):
        args.processes = available_cpus()
    elif args.processes < -1:
        args.processes = max(1, available_cpus() // (-args.processes))


def _configure_logging(verbose: bool) -> None:
    """Set the level of the CLI logger and of every vspipeline.* logger."""
    level = logging.DEBUG if verbose else logging.INFO
    for name in (__name__, "runner", "vspipeline"):
        logging.getLogger(name).setLevel(level)


# =============================
# RUN MANIFEST
# =============================

_ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")
# (version flag, regular expression picking the version line or None).
_VERSION_PROBES = {
    "obabel": ("-V", None),
    "fpocket": ("-h", re.compile(r"fpocket\s+\d[\w.]*")),
}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _sha256(path: str) -> Optional[str]:
    try:
        digest = hashlib.sha256()
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                digest.update(block)
        return digest.hexdigest()
    except OSError:
        return None


def _probe_tool(kind: str, path: Optional[str]) -> Dict:
    """Path, resolved path, version line and sha256 of an external tool."""
    info = {"path": path, "resolved_path": None, "version": None, "sha256": None}
    if not path:
        return info
    resolved = shutil.which(path)
    resolved = os.path.realpath(resolved) if resolved else None
    info["resolved_path"] = resolved
    info["sha256"] = _sha256(resolved) if resolved else None
    flag, pattern = _VERSION_PROBES.get(kind, ("--version", None))
    try:
        proc = subprocess.run([path, flag], capture_output=True, text=True,
                              timeout=VERSION_PROBE_TIMEOUT)
        text = _ANSI_ESCAPE.sub("", (proc.stdout or "") + "\n" + (proc.stderr or ""))
    except (OSError, subprocess.SubprocessError) as e:
        info["version"] = f"unknown ({type(e).__name__}: {e})"
        return info
    if pattern:
        match = pattern.search(text)
        info["version"] = match.group(0) if match else None
    else:
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        info["version"] = lines[0] if lines else None
    return info


def _tools_for_run(args) -> Dict[str, Optional[str]]:
    """External tools the run invokes, by role."""
    tools = {"obabel": OBABEL}
    if not args.smina_only:
        tools["vina"] = VINA_PRIMARY
        if VINA_FALLBACK != VINA_PRIMARY:
            tools["vina_fallback"] = VINA_FALLBACK
    if args.smina_only or args.consensus:
        tools["smina"] = SMINA
    if _uses_fpocket(args):
        tools["fpocket"] = FPOCKET
    return tools


def _probe_tools(args) -> Dict[str, Dict]:
    probes = {}
    for role, path in _tools_for_run(args).items():
        kind = "vina" if role.startswith("vina") else role
        probes[role] = _probe_tool(kind, path)
        version = probes[role]["version"] or "(no version output)"
        logger.info(f"[*] {role}: {path} - {version}")
    return probes


def _environment_snapshot() -> Dict[str, str]:
    return {k: v for k, v in sorted(os.environ.items())
            if k.startswith("VS_") or k in ("OMP_NUM_THREADS", "AUTODOCKER_VERSION",
                                            "PUBCHEM_REST")}


def _python_packages() -> Dict[str, Optional[str]]:
    """Versions of the optional Python packages that change results (None =
    not installed). Without RDKit, 2D ligands get 3D coordinates from obabel
    --gen3d, which is not reproducible."""
    packages: Dict[str, Optional[str]] = {}
    for name in ("rdkit", "numpy"):
        try:
            module = __import__(name)
            packages[name] = str(getattr(module, "__version__", "unknown"))
        except ImportError:
            packages[name] = None
    return packages


def _json_safe(value):
    """``value`` with every non-finite float replaced by its string ("nan",
    "inf"), so run_manifest.json is strict JSON (no NaN/Infinity tokens)."""
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


def write_run_manifest(path: str, manifest: Dict) -> None:
    """Write run_manifest.json with sorted keys (deterministic layout) as
    strict JSON, atomically: a reader (or a run killed mid-write) sees the
    previous file or the new one, never a truncated one."""
    tmp = f"{path}.tmp.{os.getpid()}"
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(_json_safe(manifest), f, indent=2, sort_keys=True, default=str,
                      allow_nan=False)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def _remove_stale_outputs(outdir: str) -> List[str]:
    """Remove the derived outputs (outputs.DERIVED_OUTPUTS) an earlier run left
    in ``outdir``; this run writes again the ones its options produce. The
    resumable docking results (docked/, checkpoint) are kept. Returns the
    names removed."""
    removed = []
    for name in outputs.DERIVED_OUTPUTS:
        path = os.path.join(outdir, name)
        if os.path.isdir(path) and not os.path.islink(path):
            shutil.rmtree(path)
        elif os.path.lexists(path):
            os.remove(path)
        else:
            continue
        removed.append(name)
    if removed:
        logger.info("[*] Removed outputs of an earlier run in this directory (they are "
                    "rewritten if this run produces them): " + ", ".join(removed))
    return removed


def _build_manifest(args, argv: List[str], tools: Dict, started_at: str,
                    status: str, exit_code: int, error: Optional[str],
                    extra: Dict) -> Dict:
    resolved = {k: v for k, v in sorted(vars(args).items())}
    return {
        "autodocker_version": __version__,
        "argv": [sys.argv[0]] + list(argv),
        "args": resolved,
        "seed": getattr(args, "seed", None),
        "status": status,
        "exit_code": exit_code,
        "error": error,
        "python": {"version": platform.python_version(),
                   "implementation": platform.python_implementation(),
                   "executable": sys.executable},
        "platform": platform.platform(),
        "cpu_count": {"os": os.cpu_count(), "available": available_cpus()},
        "tools": tools,
        "python_packages": _python_packages(),
        "environment": _environment_snapshot(),
        "command_timeout_s": COMMAND_TIMEOUT,
        "started_at": started_at,
        "finished_at": _utc_now(),
        **extra,
    }


# =============================
# MAIN PIPELINE
# =============================


def _keep_residue_key(residue: Dict):
    """(chain, resnum, resname) of a detected water/metal/cofactor record."""
    return (residue["chain"], residue["resnum"],
            residue.get("resname") or residue.get("name"))


def _failed_docking(results, metrics_dict) -> List[str]:
    """Ligands whose docking failed: Status FAILED (an OUT_OF_WINDOW result
    is a valid engine score, not a failure)."""
    failed = []
    for name, score in results:
        if score is not None:
            continue
        status = Status.of((metrics_dict or {}).get(name, {}).get("status"), score)
        if status is not Status.OUT_OF_WINDOW:
            failed.append(name)
    return failed


def _failure_exit_code(n_failed: int, n_ligands: int, max_fraction: float) -> int:
    """EXIT_TOO_MANY_FAILURES if more than ``max_fraction`` of the ligands
    failed, else EXIT_OK. ``n_failed`` and ``n_ligands`` include the ligands
    that could not be prepared (see _failure_counts)."""
    if n_ligands and n_failed / n_ligands > max_fraction:
        return EXIT_TOO_MANY_FAILURES
    return EXIT_OK


def _failure_counts(n_prep_failed: int, n_dock_failed: int, n_docked: int,
                    max_fraction: float) -> Dict:
    """The --max-failure-fraction bookkeeping recorded in run_manifest.json.
    A ligand fails when it cannot be prepared or its docking fails; files
    superseded by a same-name file of higher precedence are not ligands."""
    n_ligands = n_docked + n_prep_failed
    n_failed = n_dock_failed + n_prep_failed
    return {"preparation_failed": n_prep_failed, "docking_failed": n_dock_failed,
            "failed": n_failed, "ligands": n_ligands,
            "fraction": n_failed / n_ligands if n_ligands else 0.0,
            "max_failure_fraction": max_fraction}


def _threads_per_process(args) -> Optional[int]:
    """Threads per engine process (Vina/QuickVina --cpu, SMINA --cpu):
    --threads, else the available cores split between the -p workers when
    -p > 1, else None (the engine uses every core)."""
    if args.threads:
        return args.threads
    if args.processes > 1:
        return max(1, available_cpus() // args.processes)
    return None


def _run_pipeline(args, state: Dict) -> int:
    """Stages 1-5. Returns the exit code; raises on pipeline errors."""
    logger.info("=" * 70)
    logger.info(f"AUTODOCKER VIRTUAL SCREENING PIPELINE v{__version__}")
    logger.info("=" * 70)

    # ===== [1/5] Ligand library =====
    logger.info("[1/5] Ligand library preparation")
    if args.library != "local" and args.ligands:
        logger.info(f"[*] -l {args.ligands} is not used with --library {args.library}")
    lib_manager = LibraryManager(
        args.output,
        args.ligands if args.library == "local" else None,
        minimize=not args.no_minimize,
        minimize_steps=args.minimize_steps,
        seed=args.seed,
    )
    state["library"] = lib_manager
    apply_lipinski = not args.no_lipinski
    if args.library == "fda":
        ligands = lib_manager.create_fda_library(
            apply_lipinski=apply_lipinski, max_compounds=args.max_compounds)
    elif args.library == "custom":
        ligands = lib_manager.create_custom_library(
            mw_min=args.mw_min, mw_max=args.mw_max, logp_max=args.logp_max,
            apply_lipinski=apply_lipinski, max_compounds=args.max_compounds)
    else:
        ligands = lib_manager.prepare_local(apply_lipinski=apply_lipinski)
    if not ligands:
        raise ValueError("No ligands found to dock")

    # ===== [2/5] Receptor and box =====
    logger.info("[2/5] Receptor preparation and docking box")
    protein_prep = ProteinPreparation(args.receptor, args.output)
    state["protein"] = protein_prep
    protein_prep.validate()

    chain = args.chain or protein_prep.select_chain()
    protein_prep.prepare_receptor(chain, keep_hetero=args.keep_hetero, ph=args.ph)

    if args.no_fpocket:
        cx, cy, cz, sx, sy, sz = protein_prep.centroid_grid(padding=args.padding)
        logger.info("[*] Docking box: whole-receptor box (--no-fpocket)")
    elif not _uses_fpocket(args):
        # Both values are set below; no pocket is needed.
        cx = cy = cz = sx = sy = sz = 0.0
        logger.info("[*] Docking box: --grid-center and --grid-size (fpocket not run)")
    else:
        cx, cy, cz, sx, sy, sz = protein_prep.detect_pocket(
            pocket_spec=args.pockets, padding=args.padding)

    if _uses_fpocket(args) or args.no_fpocket:
        logger.info(f"[*] Grid center: {cx:.2f}, {cy:.2f}, {cz:.2f}")
        logger.info(f"[*] Grid size: {sx:.2f}, {sy:.2f}, {sz:.2f}")

    if args.grid_center:
        cx, cy, cz = _parse_grid_triplet(args.grid_center, "grid-center")
        if args.grid_center_frame == "input":
            # receptor.pdbqt and grid.conf are centered at the origin;
            # move the user's input-PDB coordinates into that frame.
            dx, dy, dz = protein_prep.frame_offset()
            logger.info(f"[*] Grid center given in input-PDB coordinates: "
                        f"{cx:.2f}, {cy:.2f}, {cz:.2f}")
            cx, cy, cz = cx - dx, cy - dy, cz - dz
        logger.info(f"[*] Grid center overridden (docking frame): {cx:.2f}, {cy:.2f}, {cz:.2f}")
    if args.grid_size:
        sx, sy, sz = _parse_grid_triplet(args.grid_size, "grid-size")
        logger.info(f"[*] Grid size overridden: {sx:.2f}, {sy:.2f}, {sz:.2f}")

    low, high = TYPICAL_GRID_SIZE_RANGE
    for axis, size in (("x", sx), ("y", sy), ("z", sz)):
        if size < low or size > high:
            logger.warning(f"[!] Grid size_{axis}={size:.2f}A outside typical range "
                           f"({low:g}-{high:g}A)")

    protein_prep.write_grid(cx, cy, cz, sx, sy, sz)

    if args.keep_waters or args.detect_metals or args.detect_cofactors:
        logger.info("[2/5] Receptor waters, metals and cofactors")
        # Detect on the ORIGINAL PDB: HETATM records (waters, metals,
        # cofactors) are stripped from pdb_clean during receptor cleaning.
        source = protein_prep.pdb_file
        keep_residues = []

        if args.keep_waters:
            logger.info("[*] Detecting water molecules near binding site...")
            # (cx, cy, cz) is in the centered receptor frame; waters are
            # read from the original PDB, so compare in that frame.
            dx, dy, dz = protein_prep.frame_offset()
            waters = detect_water_molecules(source, (cx + dx, cy + dy, cz + dz),
                                            distance_threshold=args.water_distance)
            if waters:
                logger.info(f"[*] Found {len(waters)} water residue(s) within "
                            f"{args.water_distance}A of the binding site; keeping them "
                            "in receptor")
                keep_residues += [_keep_residue_key(w) for w in waters]
            else:
                logger.info("[*] No waters detected near binding site")

        if args.detect_metals:
            logger.info("[*] Detecting metal ions...")
            metals = detect_metal_ions(source)
            if metals:
                logger.info(f"[*] Found {len(metals)} metal ion(s); keeping them in receptor:")
                for metal in metals:
                    logger.info(f"   - {metal['name']} ({metal['element']}, atom "
                                f"{metal['atom']}) chain {metal['chain']} "
                                f"residue {metal['resnum']}")
                keep_residues += [_keep_residue_key(m) for m in metals]
            else:
                logger.info("[*] No metal ions detected")

        if args.detect_cofactors:
            logger.info("[*] Detecting cofactors...")
            cofactors = detect_cofactors(source)
            if cofactors:
                logger.info(f"[*] Found {len(cofactors)} cofactor(s); keeping them in receptor:")
                for cof in cofactors:
                    logger.info(f"   - {cof['name']} chain {cof['chain']} residue {cof['resnum']}")
                keep_residues += [_keep_residue_key(c) for c in cofactors]
            else:
                logger.info("[*] No cofactors detected")

        # Detection scans every chain; keep only residues in the chains
        # that were docked (the receptor contains nothing else).
        docked_chains = set(getattr(protein_prep, "selected_chains", None) or [])
        if docked_chains:
            outside = [r for r in keep_residues if (r[0] or "A") not in docked_chains]
            if outside:
                logger.warning(
                    f"[!] Not keeping {len(outside)} residue(s) outside the docked "
                    f"chain(s) {', '.join(sorted(docked_chains))}: "
                    + ", ".join(f"{n} {c}{r}" for c, r, n in outside))
            keep_residues = [r for r in keep_residues if r not in outside]
        keep_residues = list(dict.fromkeys(keep_residues))
        if keep_residues:
            protein_prep.rebuild_receptor_keep_hetatm(keep_residues)

    # The residues actually built as flexible side chains ([] = rigid docking);
    # recorded in run_manifest.json and the closing summary.
    flexible_residues: List[str] = []
    state["flexible_residues"] = flexible_residues
    flex_file = None
    rigid_receptor = None
    if args.flexibility is None and (args.flexible_residues or args.auto_flexible):
        args.flexibility = 1
    if args.flexibility and args.flexibility > 0:
        logger.info("[2/5] Flexible receptor setup")
        requested: List[str] = []
        if args.flexible_residues:
            requested = _flexible_residue_keys(args.flexible_residues)
            logger.info(f"[*] Using specified flexible residues: {', '.join(requested)}")
        elif args.auto_flexible:
            requested = detect_flexible_residues(
                protein_prep.receptor_pdbqt, (cx, cy, cz),
                radius=AUTO_FLEX_RADIUS, max_residues=args.auto_flexible)
            logger.info(f"[*] Auto-selected {len(requested)} flexible residues "
                        "near binding site")
        else:
            logger.warning("[!] Flexibility enabled but no residues specified. "
                           "Use --flexible-residues or --auto-flexible")

        built_residues: List[str] = []
        if requested:
            logger.info(f"[*] Flexible residues: {', '.join(requested)}")
            flex_file = build_flexible_residue_pdbqt(
                protein_prep.receptor_pdbqt, requested, args.output,
                built=built_residues)
        if args.flexible_residues:
            # Residues named by the user must all be built: docking with
            # fewer (or none) would not be the requested experiment.
            built_keys = {_parse_residue_key(k) for k in built_residues}
            not_built = [k for k in requested if _parse_residue_key(k) not in built_keys]
            if not_built:
                raise RuntimeError(
                    f"--flexible-residues: could not build {', '.join(not_built)} as "
                    f"flexible side chain(s) in {protein_prep.receptor_pdbqt} (see the "
                    "warnings above: missing atoms or no rotatable side chain, e.g. "
                    "GLY, ALA, PRO); remove them from --flexible-residues")
        elif args.auto_flexible and not flex_file:
            logger.warning(
                f"[!] --auto-flexible {args.auto_flexible}: no residue within "
                f"{AUTO_FLEX_RADIUS:g} A of the box center could be made flexible; "
                "docking with a RIGID receptor")
        if flex_file:
            # Vina requires flexible residues to be excluded from the
            # rigid receptor passed to --receptor.
            flexible_residues.extend(built_residues)
            rigid_receptor = _remove_residues_from_pdbqt(
                protein_prep.receptor_pdbqt, built_residues,
                os.path.join(args.output, outputs.RIGID_RECEPTOR_PDBQT), flex_file=flex_file)

    # ===== [3/5] Docking =====
    logger.info("[3/5] Docking")
    resume = not args.no_resume
    vina_params = {
        'exhaustiveness': args.exhaustiveness,
        'binding_modes': args.binding_modes,
        'energy_range': args.energy_range,
        'seed': args.seed,
        'min_valid_affinity': args.min_valid_affinity,
        'extra_args': args.vina_extra,
        'threads': _threads_per_process(args),
    }
    state["threads_per_process"] = {
        "threads": vina_params['threads'],
        "source": ("--threads" if args.threads else
                   "available CPUs / --processes" if vina_params['threads'] else
                   "engine default (all cores)"),
        "processes": args.processes,
    }
    dock_receptor = rigid_receptor if flex_file else protein_prep.receptor_pdbqt
    grid_box = (cx, cy, cz, sx, sy, sz)

    state["docking_started"] = True
    if args.smina_only:
        # SMINA replaces Vina: dock with SMINA alone, and report SMINA's
        # own poses and failures (not Vina metadata next to SMINA scores).
        # Ligand names are unique: LibraryManager resolves same-name files.
        logger.info("[*] SMINA-only mode: docking with SMINA instead of Vina")
        results, metrics_dict = smina_dock_all(
            dock_receptor, ligands, args.output, grid_box, vina_params,
            num_processes=args.processes, flex_file=flex_file)
        if not any(score is not None or Status.of(metrics_dict[name].get("status"))
                   is Status.OUT_OF_WINDOW for name, score in results):
            # Reported like any other failure fraction (all outputs, then
            # exit 3 unless --max-failure-fraction 1.0).
            logger.error("[!] SMINA produced no score for any ligand; see "
                         f"{outputs.DOCKED_DIR}/*{outputs.SMINA_LOG_SUFFIX}")
        checkpoint = None
    else:
        # The engine paths and timeout are passed to the workers explicitly
        # (after --vina-bin/--timeout were applied in main()).
        engine = EngineConfig.from_runner()
        state["engine"] = engine
        results, checkpoint, metrics_dict = dock_all(
            dock_receptor, ligands, protein_prep.grid_conf,
            args.output, num_processes=args.processes, resume=resume,
            vina_params=vina_params, flex_file=flex_file, engine=engine)
        state["run_fingerprint"] = getattr(checkpoint, "fingerprint", None)

    scorer_agreement = None
    if args.consensus:
        logger.info("[3/5] Consensus: SMINA rescoring of the Vina top poses")
        smina_scores, scoring = rescore_with_smina(
            dock_receptor, results, args.output, metrics_dict=metrics_dict,
            num_processes=args.processes, flex_file=flex_file)
        scorer_agreement = {
            **write_consensus_outputs(args.output, results, smina_scores, scoring,
                                      metrics_dict=metrics_dict),
            "scoring": scoring}
        state["consensus"] = scorer_agreement

    # ===== [4/5] Analysis =====
    logger.info("[4/5] Pose analysis")
    if resume and checkpoint is not None and checkpoint.completed:
        logger.info(f"[*] Checkpoint updated: {len(checkpoint.completed)} successful dockings")
    html_report_file = os.path.join(args.output, outputs.HTML_REPORT)
    clusters_file = os.path.join(args.output, outputs.POSE_CLUSTERS_JSON)
    if args.cluster_poses:
        logger.info("[*] Performing pose clustering...")
        clusters = cluster_poses(
            os.path.join(args.output, outputs.DOCKED_DIR), rmsd_threshold=args.rmsd_threshold,
            ligands=[name for name, _ in results])
        clusters_data = {
            ligand: [{'id': c['id'], 'pose_count': len(c['poses'])} for c in ligand_clusters]
            for ligand, ligand_clusters in clusters.items()}
        try:
            with open(clusters_file, 'w', encoding="utf-8") as f:
                json.dump(clusters_data, f, indent=2, sort_keys=True)
            logger.info(f"[*] Clustering saved to: {clusters_file}")
        except OSError as e:
            logger.warning(f"Could not save clustering data: {e}")

    # ===== [5/5] Outputs =====
    logger.info("[5/5] Saving results")
    analyzer = ResultsAnalyzer(args.output, top_n=args.top_n)
    analyzer.save_ranking(results, mode="extended", metrics_dict=metrics_dict)
    analyzer.save_metrics_report(results, metrics_dict)

    # The HTML report is generated AFTER ranking.csv exists.
    html_ok = False
    if args.html_report:
        logger.info("[*] Generating HTML report...")
        report_meta = {
            "Protein (receptor)": os.path.basename(args.receptor),
            "Ligand library": args.library,
            "Ligands screened": len(results),
            "Top N reported": args.top_n,
            "Exhaustiveness": args.exhaustiveness,
            "Seed": args.seed,
        }
        html_ok = generate_html_report(
            os.path.join(args.output, outputs.RANKING_CSV), os.path.join(args.output, outputs.DOCKED_DIR),
            html_report_file,
            receptor_pdbqt=protein_prep.receptor_pdbqt,
            meta=report_meta,
            grid_file=protein_prep.grid_conf,
            scorer_agreement=scorer_agreement,
            rmsd_threshold=args.rmsd_threshold,
            top_n=args.top_n)

    failed = _failed_docking(results, metrics_dict)
    state["docking_failures"] = {"failed": len(failed), "docked": len(results)}
    counts = _failure_counts(len(lib_manager.failed), len(failed), len(results),
                             args.max_failure_fraction)
    state["failures"] = counts
    exit_code = _failure_exit_code(counts["failed"], counts["ligands"],
                                   args.max_failure_fraction)

    logger.info("=" * 70)
    if exit_code == EXIT_OK:
        logger.info(f"PIPELINE COMPLETED SUCCESSFULLY (v{__version__})")
    else:
        logger.error(f"PIPELINE COMPLETED, BUT {counts['failed']}/{counts['ligands']} LIGANDS "
                     f"FAILED ({counts['preparation_failed']} in preparation, "
                     f"{counts['docking_failed']} in docking; more than "
                     f"--max-failure-fraction {args.max_failure_fraction:g})")
    logger.info("=" * 70)
    logger.info(f"Results saved to: {args.output}")
    written = [("Ranking", analyzer.csv_file), ("Full results", analyzer.results_file),
               ("Top hits", analyzer.top_hits_file), ("Metrics", analyzer.metrics_file),
               ("Pocket summary", protein_prep.pocket_summary_file),
               ("Ligand metadata", lib_manager.metadata_file),
               ("Grid visualization", protein_prep.grid_box_script),
               ("Docked structures", os.path.join(args.output, outputs.DOCKED_DIR))]
    if args.consensus:
        written += [("Consensus ranking", os.path.join(args.output, outputs.CONSENSUS_RANKING_CSV)),
                    ("Consensus summary", os.path.join(args.output, outputs.CONSENSUS_SUMMARY_JSON))]
    if args.cluster_poses:
        written.append(("Pose clusters", clusters_file))
    if args.html_report:
        if html_ok:
            written.append(("HTML report", html_report_file))
        else:
            logger.warning("[!] HTML report could not be written (see warnings above)")
    written.append(("Run manifest", os.path.join(args.output, outputs.RUN_MANIFEST_JSON)))
    for label, path in written:
        if os.path.exists(path) or label == "Run manifest":
            logger.info(f"  - {label}: {path}")
    if args.flexibility and args.flexibility > 0:
        logger.info(f"  - Flexible residues: "
                    f"{', '.join(flexible_residues) if flexible_residues else 'none (rigid docking)'}")
    logger.info(f"Total compounds docked: {len(results)} ({len(failed)} failed)")
    if lib_manager.failed:
        logger.info(f"Ligands that could not be prepared: {len(lib_manager.failed)}")
    logger.info("=" * 70)
    return exit_code


def main(argv: Optional[List[str]] = None) -> None:
    global VINA_PRIMARY, VINA_FALLBACK, COMMAND_TIMEOUT

    argv = sys.argv[1:] if argv is None else list(argv)
    started_at = _utc_now()
    parser = build_parser()
    args, explicit = _parse_args(parser, argv)
    _configure_logging(args.verbose)
    _log_tool_availability()
    if os.environ.get(_QVINA_FALLBACK_ENV, "0") == "1":
        logger.warning(f"[!] ${_QVINA_FALLBACK_ENV}=1 is deprecated; use --qvina-fallback")
        args.qvina_fallback = True

    try:
        _validate_args(args, explicit)
    except ConfigError as e:
        parser.error(str(e))  # exits with EXIT_USAGE

    if args.vina_bin:
        VINA_PRIMARY = args.vina_bin
        VINA_FALLBACK = args.vina_bin
    elif args.qvina_fallback:
        VINA_FALLBACK = QVINA
    COMMAND_TIMEOUT = args.timeout
    try:
        _require_tools(smina_only=args.smina_only, consensus=args.consensus,
                       fpocket=_uses_fpocket(args))
    except RuntimeError as e:
        logger.error(f"[!] {e}")
        sys.exit(EXIT_USAGE)
    try:
        os.makedirs(args.output, exist_ok=True)
    except OSError as e:
        logger.error(f"[!] -o/--output: cannot create {args.output!r}: {e}")
        sys.exit(EXIT_USAGE)
    # One run per output directory: take its lock before anything is written
    # there (stale-output removal, ligands, receptor, grid, manifest), in
    # every mode, and hold it until the run ends. A run that finds the
    # directory in use writes nothing, not even a manifest.
    try:
        acquire_lock(args.output)
    except OutputDirLockedError as e:
        logger.error(f"[!] {e}")
        sys.exit(EXIT_USAGE)
    except OSError as e:
        logger.error(f"[!] -o/--output: cannot lock {args.output!r}: {e}")
        sys.exit(EXIT_USAGE)
    try:
        exit_code = _main_locked(args, argv, started_at)
    finally:
        release_lock(args.output)
    if exit_code != EXIT_OK:
        sys.exit(exit_code)


def _on_sigterm(signum, frame):
    # Later SIGTERMs must not interrupt the cleanup of the first one.
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    raise TerminatedBySignal()


def _install_sigterm_handler():
    """SIGTERM -> TerminatedBySignal in the main thread; returns the previous
    handler (None when signals cannot be set here, e.g. in a thread)."""
    try:
        return signal.signal(signal.SIGTERM, _on_sigterm)
    except ValueError:
        return None


def _main_locked(args, argv: List[str], started_at: str) -> int:
    """The run itself, with the output directory locked. Returns the exit
    code; always writes run_manifest.json."""
    previous_sigterm = _install_sigterm_handler()
    state: Dict = {"docking_started": False}
    tools: Dict = {}
    status, exit_code, error = "failed", EXIT_FAILURE, None
    try:
        try:
            tools = _probe_tools(args)
            _remove_stale_outputs(args.output)
            exit_code = _run_pipeline(args, state)
            status = "completed" if exit_code == EXIT_OK else "too_many_failures"
        except KeyboardInterrupt as e:
            terminated = isinstance(e, TerminatedBySignal)
            status = "interrupted"
            exit_code = EXIT_TERMINATED if terminated else EXIT_INTERRUPTED
            error = "SIGTERM" if terminated else "KeyboardInterrupt"
            what = "terminated (SIGTERM)" if terminated else "interrupted"
            if not state["docking_started"]:
                logger.warning(f"[!] Pipeline {what} before docking started.")
            elif args.smina_only:
                logger.warning(f"[!] Pipeline {what} during SMINA docking. --smina-only "
                               "keeps no checkpoint; rerun the command to dock again.")
            else:
                logger.warning(f"[!] Pipeline {what}. Docking progress is saved in the "
                               "checkpoint; rerun the same command to resume.")
        except Exception as e:
            error = f"{type(e).__name__}: {e}"
            logger.error(f"Pipeline failed: {error}")
            if args.verbose:
                import traceback
                logger.error(traceback.format_exc())
        finally:
            # The manifest is written whatever happened; a second SIGTERM
            # must not cut it short.
            try:
                signal.signal(signal.SIGTERM, signal.SIG_IGN)
            except ValueError:
                pass
            _write_final_manifest(args, argv, tools, started_at, status, exit_code,
                                  error, state)
    finally:
        if previous_sigterm is not None:
            signal.signal(signal.SIGTERM, previous_sigterm)
    return exit_code


def _write_final_manifest(args, argv, tools, started_at, status, exit_code, error,
                          state) -> None:
    extra = {}
    library = state.get("library")
    if library is not None:
        extra["ligand_preparation"] = {"failed": len(library.failed),
                                       "superseded": len(library.superseded)}
    if "docking_failures" in state:
        extra["docking_failures"] = state["docking_failures"]
    if "failures" in state:
        # What --max-failure-fraction was applied to (preparation + docking).
        extra["failures"] = state["failures"]
    if "threads_per_process" in state:
        # Threads of each Vina/QuickVina/SMINA process (--cpu); null = all cores.
        extra["threads_per_process"] = state["threads_per_process"]
    if state.get("run_fingerprint"):
        # The checkpoint fingerprint: receptor, grid, flex file, docking
        # parameters and the primary engine binary (see docking._run_fingerprint).
        extra["run_fingerprint"] = state["run_fingerprint"]
    if state.get("engine") is not None:
        extra["docking_engines"] = [
            _engine_identity(path) for path, _ in state["engine"].docking_engines()]
    if "consensus" in state:
        extra["consensus"] = state["consensus"]
    if "flexible_residues" in state:
        # The residues docked as flexible side chains ([] = rigid).
        extra["flexible_residues"] = list(state["flexible_residues"])
    manifest_file = os.path.join(args.output, outputs.RUN_MANIFEST_JSON)
    try:
        write_run_manifest(manifest_file, _build_manifest(
            args, argv, tools, started_at, status, exit_code, error, extra))
    except OSError as e:
        logger.warning(f"Could not write {manifest_file}: {e}")


if __name__ == "__main__":
    main()
# vim: ts=4:et
