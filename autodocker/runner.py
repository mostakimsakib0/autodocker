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
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from typing import Dict, List, Optional

_FALLBACK_VERSION = "2.0.0"


def _read_version(base_dir: Optional[str] = None) -> str:
    """Version string: autodocker/VERSION if present, else $AUTODOCKER_VERSION,
    else the version this file was released with."""
    base_dir = base_dir or os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base_dir, "VERSION")
    try:
        with open(path) as f:
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

DEFAULT_SEED = 42
DEFAULT_COMMAND_TIMEOUT = 900  # seconds
# Residues within this distance (Angstrom) of the box center are candidates
# for --auto-flexible.
AUTO_FLEX_RADIUS = 8.0
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
if os.environ.get("VS_ENABLE_QVINA_FALLBACK", "0") == "1":
    VINA_FALLBACK = find_tool("qvina02", "qvina2") or VINA_PRIMARY
else:
    VINA_FALLBACK = VINA_PRIMARY
# Optional engines and helpers.
SMINA = find_tool("smina", "smina.exe")
FPOCKET = find_tool("fpocket")
# reduce and AutoDock4/AutoGrid4 are not used by the pipeline (receptor
# hydrogens come from Open Babel); the paths are kept for the real-tool tests.
REDUCE = find_tool("reduce")
AUTODOCK = find_tool("autodock4", "autodock")
AUTOGRID = find_tool("autogrid4", "autogrid")
COMMAND_TIMEOUT = _env_command_timeout()


def _require_tools(smina_only: bool = False, consensus: bool = False) -> None:
    """Fail when a tool the requested run needs is missing.

    Called from main() after argument parsing, so --help/--version work and
    --vina-bin can supply a docking binary that is not on PATH.
    """
    needed = [("obabel", OBABEL)]
    if not smina_only:
        needed.append(("vina", VINA_PRIMARY))
    if smina_only or consensus:
        needed.append(("smina", SMINA))
    missing = [name for name, path in needed if path is None]
    if missing:
        raise RuntimeError(
            f"Required tool(s) not found in PATH: {', '.join(missing)}. "
            "Install them or put them on PATH (Vina can also be given with --vina-bin).")


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
from vspipeline.checkpoint import DockingCheckpoint  # noqa: E402
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
    "DockingCheckpoint", "LibraryManager", "ProteinPreparation",
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

DEFAULT_BINDING_MODES = 9
DEFAULT_ENERGY_RANGE = 3.0
DEFAULT_MAX_FAILURE_FRACTION = 0.5

# Options that only apply to --library custom (or fda for --max-compounds);
# main() detects whether they were given explicitly.
_UNSET = object()
_EXPLICIT_DESTS = ("mw_min", "mw_max", "logp_max", "max_compounds")


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
            "Exit codes: 0 success, 1 pipeline error, 2 invalid options/environment or\n"
            "missing tools, 3 more ligands failed to dock than --max-failure-fraction\n"
            "allows (outputs are still written), 130 interrupted."),
        formatter_class=_HelpFormatter,
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
                     help="local: the files of -l; fda: a fixed list of 74 well-known "
                     "approved drugs downloaded from PubChem (not the full FDA list); "
                     "custom: the same list filtered by --mw-min/--mw-max/--logp-max")
    lib.add_argument("--mw-min", type=float, default=200.0,
                     help="--library custom only: minimum molecular weight (Da)")
    lib.add_argument("--mw-max", type=float, default=500.0,
                     help="--library custom only: maximum molecular weight (Da); "
                     "replaces Lipinski's MW <= 500")
    lib.add_argument("--logp-max", type=float, default=5.0,
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
                     "ligands only; MOL2/PDB/PDBQT keep their input coordinates)")
    lib.add_argument("--minimize-steps", type=int, default=250,
                     help="MMFF94 minimization steps (SDF and PubChem ligands only)")

    rec = parser.add_argument_group("receptor and docking box")
    rec.add_argument("--chain", default=None,
                     help="Protein chain(s), e.g. A, A,B, or all (IDs are case-sensitive). "
                     "Default: the first chain in sorted order (a warning lists the others)")
    rec.add_argument("--ph", type=float, default=7.4,
                     help="pH at which receptor hydrogens are added (Open Babel -p)")
    rec.add_argument("--keep-hetero", action="store_true",
                     help="Keep matching HETATM records during receptor cleaning")
    rec.add_argument("--pockets",
                     help="fpocket pocket number(s), e.g. 1 or 2,3. Default: best fpocket "
                     "score (druggability breaks ties)")
    rec.add_argument("--padding", type=float, default=6.0,
                     help="Box padding around the selected pocket(s), Angstrom")
    rec.add_argument("--no-fpocket", action="store_true",
                     help="Skip fpocket; dock in a box around the whole receptor (heavy-atom "
                     "bounding box + padding, each edge at least 24 Angstrom)")
    rec.add_argument("--grid-center", default=None,
                     help="Explicit box center 'X,Y,Z' (Angstrom), e.g. 12.3,4.5,6.7; "
                     "overrides auto-detection (frame: see --grid-center-frame)")
    rec.add_argument("--grid-center-frame", choices=["input", "docking"], default="input",
                     help="Frame of --grid-center: 'input' (receptor PDB coordinates) or "
                     "'docking' (the centered frame of receptor.pdbqt and grid.conf; "
                     "versions <= 1.1.4 used this frame)")
    rec.add_argument("--grid-size", default=None,
                     help="Explicit box size 'SX,SY,SZ' (Angstrom), e.g. 22.5,22.5,22.5; "
                     "overrides auto-detection")
    rec.add_argument("--keep-waters", action="store_true",
                     help="Keep water molecules near the binding site in the receptor")
    rec.add_argument("--water-distance", type=float, default=4.0,
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
                      help="Flexible residues, e.g. 'A45,A78,B102'")
    flex.add_argument("--auto-flexible", type=int,
                      help=f"Make the N residues closest to the box center (within "
                      f"{AUTO_FLEX_RADIUS:g} A) flexible")

    dock = parser.add_argument_group("docking")
    dock.add_argument("-p", "--processes", type=int, default=1,
                      help="Parallel docking processes: N > 0 uses N; 0 or -1 uses all "
                      "available CPUs; -N (N > 1) uses available CPUs / N")
    dock.add_argument("--exhaustiveness", type=int, default=8,
                      help="Vina exhaustiveness (>= 1; 8 is Vina's default, higher is "
                      "more thorough and slower)")
    dock.add_argument("--binding-modes", type=int, default=DEFAULT_BINDING_MODES,
                      help="Maximum number of binding modes Vina reports")
    dock.add_argument("--energy-range", type=float, default=DEFAULT_ENERGY_RANGE,
                      help="Vina energy range of the reported modes (kcal/mol)")
    dock.add_argument("--seed", type=int, default=None,
                      help=f"Random seed (default: $VS_SEED, else {DEFAULT_SEED})")
    dock.add_argument("--threads", type=int, default=None,
                      help="CPU threads per Vina process, passed as --cpu (default: all "
                      "cores, or cores divided by -p when -p > 1)")
    dock.add_argument("--min-valid-affinity", type=float, default=-1.0,
                      help="Upper bound of the accepted affinity window (kcal/mol): scores "
                      f">= this, or < {AFFINITY_FLOOR:g}, are reported as OUT_OF_WINDOW")
    dock.add_argument("--vina-bin",
                      help="Docking binary to use instead of the auto-detected Vina/QuickVina")
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
                      "fraction of the docked ligands failed (scores outside the "
                      "accepted window are not failures); 1.0 never fails")
    dock.add_argument("--consensus", action="store_true",
                      help="Rescore Vina's top pose of each ligand with SMINA (Vinardo "
                      "scoring) and write consensus_ranking.csv and consensus_summary.json "
                      "(Spearman rho between the two scores)")
    dock.add_argument("--smina-only", action="store_true",
                      help="Dock with SMINA (default scoring) instead of Vina (1 mode per "
                      "ligand, no checkpoint/resume; not combinable with --consensus, "
                      "--vina-bin, --vina-extra or non-default --binding-modes/--energy-range)")

    out = parser.add_argument_group("analysis and reports")
    out.add_argument("--top-n", type=int, default=20,
                     help="Number of compounds in Top_hits.txt")
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
    out.add_argument("--rmsd-threshold", type=float, default=2.0,
                     help="RMSD clustering threshold (Angstrom)")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Debug logging (all pipeline modules)")
    parser.add_argument("--version", action="version",
                        version=f"AutoDocker {__version__}")
    return parser


def _parse_args(parser: argparse.ArgumentParser, argv: List[str]):
    """Parse ``argv``; return (args, set of _EXPLICIT_DESTS given explicitly).

    argparse only fills defaults for attributes missing from the namespace,
    so a sentinel there shows whether an option was given (abbreviations
    such as --mw-ma included).
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


def _validate_args(args, explicit) -> None:
    """Check option values and combinations; raise ConfigError on the first
    problem. Resolves --seed, --processes and the deprecated --no-admet."""
    if _ENV_ERRORS:
        raise ConfigError("; ".join(_ENV_ERRORS))
    if args.seed is None:
        args.seed = default_seed()
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
    custom_only = sorted(d for d in explicit if d != "max_compounds")
    if args.library != "custom" and custom_only:
        options = ", ".join("--" + d.replace("_", "-") for d in custom_only)
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
    if not 0.0 <= args.max_failure_fraction <= 1.0:
        raise ConfigError("--max-failure-fraction must be between 0 and 1, got "
                          f"{args.max_failure_fraction:g}")
    for option, value in (("--grid-center", args.grid_center),
                          ("--grid-size", args.grid_size)):
        if value:
            try:
                triplet = _parse_grid_triplet(value, option.lstrip("-"))
            except ValueError as e:
                raise ConfigError(str(e)) from None
            if option == "--grid-size" and min(triplet) <= 0:
                raise ConfigError(f"--grid-size must be positive, got {value}")

    if args.vina_bin and not shutil.which(args.vina_bin):
        raise ConfigError(f"--vina-bin {args.vina_bin!r} does not exist or is not executable")

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
    if not args.no_fpocket:
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


def write_run_manifest(path: str, manifest: Dict) -> None:
    """Write run_manifest.json with sorted keys (deterministic layout)."""
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True, default=str)
        f.write("\n")


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


def _failure_exit_code(n_failed: int, n_docked: int, max_fraction: float) -> int:
    """EXIT_TOO_MANY_FAILURES if more than ``max_fraction`` of the docked
    ligands failed, else EXIT_OK."""
    if n_docked and n_failed / n_docked > max_fraction:
        return EXIT_TOO_MANY_FAILURES
    return EXIT_OK


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
    else:
        cx, cy, cz, sx, sy, sz = protein_prep.detect_pocket(
            pocket_spec=args.pockets, padding=args.padding)

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

    flexible_residues: List[str] = []
    flex_file = None
    rigid_receptor = None
    if args.flexibility is None and (args.flexible_residues or args.auto_flexible):
        args.flexibility = 1
    if args.flexibility and args.flexibility > 0:
        logger.info("[2/5] Flexible receptor setup")
        if args.flexible_residues:
            flexible_residues = [
                r.strip() for r in args.flexible_residues.split(',') if r.strip()]
            logger.info(f"[*] Using specified flexible residues: {', '.join(flexible_residues)}")
        elif args.auto_flexible:
            flexible_residues = detect_flexible_residues(
                protein_prep.receptor_pdbqt, (cx, cy, cz),
                radius=AUTO_FLEX_RADIUS, max_residues=args.auto_flexible)
            logger.info(f"[*] Auto-selected {len(flexible_residues)} flexible residues "
                        "near binding site")
        else:
            logger.warning("[!] Flexibility enabled but no residues specified. "
                           "Use --flexible-residues or --auto-flexible")

        if flexible_residues:
            logger.info(f"[*] Flexible residues: {', '.join(flexible_residues)}")
            built_residues: List[str] = []
            flex_file = build_flexible_residue_pdbqt(
                protein_prep.receptor_pdbqt, flexible_residues, args.output,
                built=built_residues)
            if not flex_file:
                logger.warning("[!] Flexible receptor could not be built; running rigid docking")
            else:
                # Vina requires flexible residues to be excluded from the
                # rigid receptor passed to --receptor.
                flexible_residues = built_residues
                rigid_receptor = _remove_residues_from_pdbqt(
                    protein_prep.receptor_pdbqt, built_residues,
                    os.path.join(args.output, "rigid.pdbqt"), flex_file=flex_file)

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
        'threads': args.threads,
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
            raise RuntimeError("SMINA produced no score for any ligand")
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
            **write_consensus_outputs(args.output, results, smina_scores, scoring),
            "scoring": scoring}
        state["consensus"] = scorer_agreement

    # ===== [4/5] Analysis =====
    logger.info("[4/5] Pose analysis")
    if resume and checkpoint is not None and checkpoint.completed:
        logger.info(f"[*] Checkpoint updated: {len(checkpoint.completed)} successful dockings")
    html_report_file = os.path.join(args.output, "results_report.html")
    clusters_file = os.path.join(args.output, "pose_clusters.json")
    if args.cluster_poses:
        logger.info("[*] Performing pose clustering...")
        clusters = cluster_poses(
            os.path.join(args.output, "docked"), rmsd_threshold=args.rmsd_threshold,
            ligands=[name for name, _ in results])
        clusters_data = {
            ligand: [{'id': c['id'], 'pose_count': len(c['poses'])} for c in ligand_clusters]
            for ligand, ligand_clusters in clusters.items()}
        try:
            with open(clusters_file, 'w') as f:
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
            os.path.join(args.output, 'ranking.csv'), os.path.join(args.output, "docked"),
            html_report_file,
            receptor_pdbqt=protein_prep.receptor_pdbqt,
            meta=report_meta,
            grid_file=protein_prep.grid_conf,
            scorer_agreement=scorer_agreement,
            rmsd_threshold=args.rmsd_threshold)

    failed = _failed_docking(results, metrics_dict)
    state["docking_failures"] = {"failed": len(failed), "docked": len(results)}
    exit_code = _failure_exit_code(len(failed), len(results), args.max_failure_fraction)

    logger.info("=" * 70)
    if exit_code == EXIT_OK:
        logger.info(f"PIPELINE COMPLETED SUCCESSFULLY (v{__version__})")
    else:
        logger.error(f"PIPELINE COMPLETED, BUT {len(failed)}/{len(results)} LIGANDS FAILED "
                     f"TO DOCK (more than --max-failure-fraction "
                     f"{args.max_failure_fraction:g})")
    logger.info("=" * 70)
    logger.info(f"Results saved to: {args.output}")
    outputs = [("Ranking", analyzer.csv_file), ("Full results", analyzer.results_file),
               ("Top hits", analyzer.top_hits_file), ("Metrics", analyzer.metrics_file),
               ("Pocket summary", protein_prep.pocket_summary_file),
               ("Ligand metadata", lib_manager.metadata_file),
               ("Grid visualization", protein_prep.grid_box_script),
               ("Docked structures", os.path.join(args.output, "docked"))]
    if args.consensus:
        outputs += [("Consensus ranking", os.path.join(args.output, "consensus_ranking.csv")),
                    ("Consensus summary", os.path.join(args.output, "consensus_summary.json"))]
    if args.cluster_poses:
        outputs.append(("Pose clusters", clusters_file))
    if args.html_report:
        if html_ok:
            outputs.append(("HTML report", html_report_file))
        else:
            logger.warning("[!] HTML report could not be written (see warnings above)")
    outputs.append(("Run manifest", os.path.join(args.output, "run_manifest.json")))
    for label, path in outputs:
        if os.path.exists(path) or label == "Run manifest":
            logger.info(f"  - {label}: {path}")
    if args.flexibility and args.flexibility > 0:
        logger.info(f"  - Flexible residues: "
                    f"{', '.join(flexible_residues) if flexible_residues else 'none'}")
    logger.info(f"Total compounds docked: {len(results)} ({len(failed)} failed)")
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

    try:
        _validate_args(args, explicit)
    except ConfigError as e:
        parser.error(str(e))  # exits with EXIT_USAGE

    if args.vina_bin:
        VINA_PRIMARY = args.vina_bin
        VINA_FALLBACK = args.vina_bin
    COMMAND_TIMEOUT = args.timeout
    try:
        _require_tools(smina_only=args.smina_only, consensus=args.consensus)
    except RuntimeError as e:
        logger.error(f"[!] {e}")
        sys.exit(EXIT_USAGE)

    os.makedirs(args.output, exist_ok=True)
    tools = _probe_tools(args)

    state: Dict = {"docking_started": False}
    status, exit_code, error = "failed", EXIT_FAILURE, None
    try:
        exit_code = _run_pipeline(args, state)
        status = "completed" if exit_code == EXIT_OK else "too_many_failures"
    except KeyboardInterrupt:
        status, exit_code, error = "interrupted", EXIT_INTERRUPTED, "KeyboardInterrupt"
        if state["docking_started"] and not args.smina_only:
            logger.warning("[!] Pipeline interrupted. Docking progress is saved in the "
                           "checkpoint; rerun the same command to resume.")
        else:
            logger.warning("[!] Pipeline interrupted before docking started.")
    except Exception as e:
        error = f"{type(e).__name__}: {e}"
        logger.error(f"Pipeline failed: {error}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
    finally:
        extra = {}
        library = state.get("library")
        if library is not None:
            extra["ligand_preparation"] = {"failed": len(library.failed),
                                           "superseded": len(library.superseded)}
        if "docking_failures" in state:
            extra["docking_failures"] = state["docking_failures"]
        if state.get("run_fingerprint"):
            # The checkpoint fingerprint: receptor, grid, flex file, docking
            # parameters and engine binaries (see docking._run_fingerprint).
            extra["run_fingerprint"] = state["run_fingerprint"]
        if state.get("engine") is not None:
            extra["docking_engines"] = [
                _engine_identity(path) for path, _ in state["engine"].docking_engines()]
        if "consensus" in state:
            extra["consensus"] = state["consensus"]
        manifest_file = os.path.join(args.output, "run_manifest.json")
        try:
            write_run_manifest(manifest_file, _build_manifest(
                args, argv, tools, started_at, status, exit_code, error, extra))
        except OSError as e:
            logger.warning(f"Could not write {manifest_file}: {e}")
    if exit_code != EXIT_OK:
        sys.exit(exit_code)


if __name__ == "__main__":
    main()
# vim: ts=4:et
