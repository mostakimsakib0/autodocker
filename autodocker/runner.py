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
import logging
import math
import os
import re
import shutil
import signal
import sys
from typing import Dict, List, Optional

from vspipeline import (
    checkpoint,
    defaults,
    docking,
    flexible,
    library,
    manifest,
    outputs,
    pipeline,
    protein,
    tools,
    utils,
)
from vspipeline.version import __version__

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
# Residues within this distance (Angstrom) of the box center are candidates
# for --auto-flexible.
AUTO_FLEX_RADIUS = pipeline.AUTO_FLEX_RADIUS


class ConfigError(ValueError):
    """Invalid option or environment value (reported without a traceback)."""


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


# Deprecated environment switch for --qvina-fallback (versions <= 1.1.x).
_QVINA_FALLBACK_ENV = "VS_ENABLE_QVINA_FALLBACK"


def _require_tools(smina_only: bool = False, consensus: bool = False,
                   fpocket: bool = False) -> None:
    """Fail when a tool the requested run needs is missing.

    Called from main() after argument parsing, so --help/--version work and
    --vina-bin can supply a docking binary that is not on PATH. ``fpocket``:
    the run detects pockets (see pipeline.uses_fpocket).
    """
    found = tools.settings()
    needed = [("obabel", found.obabel)]
    if not smina_only:
        needed.append(("vina", found.vina_primary))
    if smina_only or consensus:
        needed.append(("smina", found.smina))
    if fpocket:
        needed.append(("fpocket", found.fpocket))
    missing = [name for name, path in needed if path is None]
    if missing:
        hint = "Install them or put them on PATH (Vina can also be given with --vina-bin)."
        if "fpocket" in missing:
            hint += (" Without fpocket, pass --no-fpocket (box around the whole receptor) "
                     "or give the box with --grid-center and --grid-size.")
        raise RuntimeError(
            f"Required tool(s) not found in PATH: {', '.join(missing)}. {hint}")


def _log_tool_availability() -> None:
    found = tools.settings()
    logger.debug(f"SMINA: {found.smina or 'not found (--consensus and --smina-only unavailable)'}")
    logger.debug(f"fpocket: {found.fpocket or 'not found'}")


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

    _add_input_output_options(parser)
    _add_library_options(parser)
    _add_receptor_options(parser)
    _add_flexible_options(parser)
    _add_docking_options(parser)
    _add_report_options(parser)
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Debug logging (all pipeline modules)")
    parser.add_argument("--version", action="version",
                        version=f"AutoDocker {__version__}")
    return parser


def _add_input_output_options(parser: argparse.ArgumentParser) -> None:
    """-r/-l/-o."""
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


def _add_library_options(parser: argparse.ArgumentParser) -> None:
    """--library and the ligand preparation options."""
    lib = parser.add_argument_group("ligand library")
    lib.add_argument("--library", choices=["fda", "custom", "local"], default="local",
                     help="local: the files of -l; fda: a fixed list of "
                     f"{len(library.LibraryManager.FDA_DRUG_NAMES)} well-known "
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
    lib.add_argument("--max-compounds", type=int, default=library.LibraryManager.MAX_COMPOUNDS,
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


def _add_receptor_options(parser: argparse.ArgumentParser) -> None:
    """Receptor preparation and docking box options."""
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
                     f"{protein.MIN_RECEPTOR_BOX_SIZE:g} Angstrom)")
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


def _add_flexible_options(parser: argparse.ArgumentParser) -> None:
    """Flexible side-chain options."""
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


def _add_docking_options(parser: argparse.ArgumentParser) -> None:
    """Docking engine and run-control options."""
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
                      f">= this, or < {docking.AFFINITY_FLOOR:g}, are reported as OUT_OF_WINDOW; "
                      f"must be > {docking.AFFINITY_FLOOR:g}")
    dock.add_argument("--vina-bin",
                      help="Docking binary to use instead of the auto-detected Vina/QuickVina")
    dock.add_argument("--qvina-fallback", action="store_true",
                      help="Re-dock a ligand with QuickVina 2 when Vina fails on it (crash "
                      "or timeout). The two engines' scores then share one ranking; "
                      "ranking.csv's Engine column shows which engine scored each ligand. "
                      f"Replaces the deprecated ${_QVINA_FALLBACK_ENV}=1")
    dock.add_argument("--vina-extra", action="append", default=[],
                      help="Extra argument string passed to the docking binary; repeat as needed")
    dock.add_argument("--timeout", type=int, default=tools.settings().command_timeout,
                      help="Per-command timeout in seconds; the default is "
                      f"$VS_COMMAND_TIMEOUT if set, else {tools.DEFAULT_COMMAND_TIMEOUT}")
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


def _add_report_options(parser: argparse.ArgumentParser) -> None:
    """Analysis and report options."""
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
        available = utils.get_chains(args.receptor)
    except OSError as e:
        raise ConfigError(f"-r/--receptor: cannot read {args.receptor}: {e}") from None
    except ValueError as e:  # utils.ChainAmbiguityError: blank chain next to chain A
        raise ConfigError(f"-r/--receptor: {e}") from None
    if args.chain is not None:
        try:
            selected = utils._parse_chain_selection(args.chain, available)
        except ValueError as e:
            raise ConfigError(f"--chain {args.chain!r}: {e}") from None
    else:
        selected = available[:1]

    if not args.flexible_residues or (args.flexibility is not None and args.flexibility <= 0):
        return
    try:
        lines, _ = protein.read_pdb_records(args.receptor)
    except (OSError, UnicodeDecodeError) as e:
        raise ConfigError(f"-r/--receptor: cannot read {args.receptor}: {e}") from None
    missing = flexible.missing_residues(lines, selected,
                                        pipeline.flexible_residue_keys(args.flexible_residues))
    if missing:
        raise ConfigError(
            f"--flexible-residues: no such residue in the docked chain(s) "
            f"{', '.join(selected)} of {args.receptor}: {', '.join(missing)}")


def _validate_args(args, explicit) -> None:
    """Check option values and combinations; raise ConfigError on the first
    problem (the checks run in this order). Resolves --seed, --processes and
    the deprecated --no-admet. Reads the receptor's chain list, but writes
    nothing."""
    _check_environment_and_seed(args)
    _check_inputs(args, explicit)
    _check_values(args)
    _check_box_options(args, explicit)
    _check_flexible_options(args)
    _check_engine_options(args)
    _check_receptor_selection(args)
    _check_output_dir(args.output)
    _resolve_processes(args)


def _check_environment_and_seed(args) -> None:
    """$VS_* problems, non-finite floats and the seed range; resolves --seed
    ($VS_SEED, else the default) and the deprecated --no-admet."""
    if tools.ENV_ERRORS:
        raise ConfigError("; ".join(tools.ENV_ERRORS))
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


def _check_inputs(args, explicit) -> None:
    """-r, -l and the --library specific options."""
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


def _check_values(args) -> None:
    """Value ranges of the numeric options and the --grid-center/--grid-size
    triplets."""
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
    if not args.min_valid_affinity > docking.AFFINITY_FLOOR:
        raise ConfigError(
            f"--min-valid-affinity must be > {docking.AFFINITY_FLOOR:g} kcal/mol (the lower end "
            f"of the accepted window), got {args.min_valid_affinity:g}: no score could "
            "be accepted")
    if not 0.0 <= args.max_failure_fraction <= 1.0:
        raise ConfigError("--max-failure-fraction must be between 0 and 1, got "
                          f"{args.max_failure_fraction:g}")
    for option, value in (("--grid-center", args.grid_center),
                          ("--grid-size", args.grid_size)):
        if value is not None:
            try:
                triplet = docking._parse_grid_triplet(value, option.lstrip("-"))
            except ValueError as e:
                raise ConfigError(str(e)) from None
            if option == "--grid-size" and min(triplet) <= 0:
                raise ConfigError(f"--grid-size must be positive, got {value}")


def _check_box_options(args, explicit) -> None:
    """Box options that would be accepted but have no effect, and --pockets."""
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


def _check_flexible_options(args) -> None:
    """--flexible-residues: not with --auto-flexible, and valid residue keys."""
    if args.flexible_residues is not None:
        if args.auto_flexible is not None:
            raise ConfigError("--flexible-residues and --auto-flexible cannot be combined "
                              "(name the residues, or let N of them be chosen)")
        keys = pipeline.flexible_residue_keys(args.flexible_residues)
        if not keys:
            raise ConfigError(f"--flexible-residues lists no residue, got "
                              f"{args.flexible_residues!r} (e.g. A45,A78)")
        invalid = [k for k in keys if flexible._parse_residue_key(k) is None]
        if invalid:
            raise ConfigError(
                f"--flexible-residues: invalid residue key(s) {', '.join(invalid)} "
                "(expected chain + residue number, e.g. A45 or A100B; chain:number "
                "for any chain, e.g. 1:45; a blank chain is A: A45 or 45)")


def _check_engine_options(args) -> None:
    """--vina-bin, --qvina-fallback and the --smina-only conflicts."""
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
        if not tools.settings().qvina:
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


def _resolve_processes(args) -> None:
    """-p 0/-1: all available CPUs; -p -N: available CPUs / N."""
    if args.processes in (0, -1):
        args.processes = utils.available_cpus()
    elif args.processes < -1:
        args.processes = max(1, utils.available_cpus() // (-args.processes))


def _configure_logging(verbose: bool) -> None:
    """Set the level of the CLI logger and of every vspipeline.* logger."""
    level = logging.DEBUG if verbose else logging.INFO
    for name in (__name__, "runner", "vspipeline"):
        logging.getLogger(name).setLevel(level)


# =============================
# RUN
# =============================

def _tools_for_run(args) -> Dict[str, Optional[str]]:
    """External tools the run invokes, by role."""
    found = tools.settings()
    used = {"obabel": found.obabel}
    if not args.smina_only:
        used["vina"] = found.vina_primary
        if found.vina_fallback != found.vina_primary:
            used["vina_fallback"] = found.vina_fallback
    if args.smina_only or args.consensus:
        used["smina"] = found.smina
    if pipeline.uses_fpocket(args):
        used["fpocket"] = found.fpocket
    return used


def _probe_tools(args) -> Dict[str, Dict]:
    """Path, version and checksum of each tool of the run (run_manifest.json)."""
    probes = {}
    for role, path in _tools_for_run(args).items():
        kind = "vina" if role.startswith("vina") else role
        probes[role] = manifest.probe_tool(kind, path)
        version = probes[role]["version"] or "(no version output)"
        logger.info(f"[*] {role}: {path} - {version}")
    return probes


def _failure_exit_code(n_failed: int, n_ligands: int, max_fraction: float) -> int:
    """EXIT_TOO_MANY_FAILURES if more than ``max_fraction`` of the ligands
    failed, else EXIT_OK. ``n_failed`` and ``n_ligands`` include the ligands
    that could not be prepared (see pipeline.failure_counts)."""
    if n_ligands and n_failed / n_ligands > max_fraction:
        return EXIT_TOO_MANY_FAILURES
    return EXIT_OK


def _run_pipeline(args, state: pipeline.RunState) -> int:
    """Stages 1-5 (vspipeline.pipeline). Returns the exit code; raises on
    pipeline errors."""
    logger.info("=" * 70)
    logger.info(f"AUTODOCKER VIRTUAL SCREENING PIPELINE v{__version__}")
    logger.info("=" * 70)

    ligands = pipeline.prepare_ligands(args, state)                       # [1/5]
    receptor = pipeline.prepare_receptor(args, state)                     # [2/5]
    box = pipeline.docking_box(args, receptor)
    pipeline.keep_hetatm_residues(args, receptor, box)
    flex = pipeline.flexible_receptor(args, state, receptor, box)
    docked = pipeline.dock(args, state, ligands, receptor, box, flex)     # [3/5]
    pipeline.consensus_rescoring(args, state, docked)
    clusters_file = pipeline.cluster_poses(args, docked)                  # [4/5]
    analyzer, html_ok = pipeline.write_results(args, receptor, docked)    # [5/5]

    counts = pipeline.count_failures(args, state, docked)
    exit_code = _failure_exit_code(counts["failed"], counts["ligands"],
                                   args.max_failure_fraction)
    pipeline.log_summary(args, state, receptor, docked, analyzer, html_ok, clusters_file,
                         too_many_failures=exit_code != EXIT_OK)
    return exit_code


# =============================
# MAIN
# =============================

def main(argv: Optional[List[str]] = None) -> None:
    argv = sys.argv[1:] if argv is None else list(argv)
    started_at = manifest.utc_now()
    args = _parse_and_validate(argv)      # invalid options: exits with EXIT_USAGE
    _apply_tool_options(args)
    _claim_output_dir(args)               # missing tools, -o unusable: EXIT_USAGE
    try:
        exit_code = _main_locked(args, argv, started_at)
    finally:
        checkpoint.release_lock(args.output)
    if exit_code != EXIT_OK:
        sys.exit(exit_code)


def _parse_and_validate(argv: List[str]):
    """The validated options; exits with EXIT_USAGE (argparse's
    parser.error) on an invalid option or environment value."""
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
    return args


def _apply_tool_options(args) -> None:
    """--vina-bin, --qvina-fallback and --timeout -> tools.settings()."""
    if args.vina_bin:
        tools.configure(vina_primary=args.vina_bin, vina_fallback=args.vina_bin)
    elif args.qvina_fallback:
        tools.configure(vina_fallback=tools.settings().qvina)
    tools.configure(command_timeout=args.timeout)


def _claim_output_dir(args) -> None:
    """Check the tools, create -o and lock it for this run; exits with
    EXIT_USAGE when that fails."""
    try:
        _require_tools(smina_only=args.smina_only, consensus=args.consensus,
                       fpocket=pipeline.uses_fpocket(args))
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
        checkpoint.acquire_lock(args.output)
    except checkpoint.OutputDirLockedError as e:
        logger.error(f"[!] {e}")
        sys.exit(EXIT_USAGE)
    except OSError as e:
        logger.error(f"[!] -o/--output: cannot lock {args.output!r}: {e}")
        sys.exit(EXIT_USAGE)


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
    state = pipeline.RunState()
    tool_probes: Dict = {}
    status, exit_code, error = "failed", EXIT_FAILURE, None
    try:
        try:
            tool_probes = _probe_tools(args)
            state.openbabel = manifest.check_openbabel(args.seed)
            pipeline.remove_stale_outputs(args.output)
            exit_code = _run_pipeline(args, state)
            status = "completed" if exit_code == EXIT_OK else "too_many_failures"
        except KeyboardInterrupt as e:
            terminated = isinstance(e, TerminatedBySignal)
            status = "interrupted"
            exit_code = EXIT_TERMINATED if terminated else EXIT_INTERRUPTED
            error = "SIGTERM" if terminated else "KeyboardInterrupt"
            _log_interrupted(args, state, "terminated (SIGTERM)" if terminated else "interrupted")
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
            _write_final_manifest(args, argv, tool_probes, started_at, status, exit_code,
                                  error, state)
    finally:
        if previous_sigterm is not None:
            signal.signal(signal.SIGTERM, previous_sigterm)
    return exit_code


def _log_interrupted(args, state: pipeline.RunState, what: str) -> None:
    """How to continue after Ctrl-C or SIGTERM."""
    if not state.docking_started:
        logger.warning(f"[!] Pipeline {what} before docking started.")
    elif args.smina_only:
        logger.warning(f"[!] Pipeline {what} during SMINA docking. --smina-only "
                       "keeps no checkpoint; rerun the command to dock again.")
    else:
        logger.warning(f"[!] Pipeline {what}. Docking progress is saved in the "
                       "checkpoint; rerun the same command to resume.")


def _write_final_manifest(args, argv, tool_probes, started_at, status, exit_code, error,
                          state: pipeline.RunState) -> None:
    manifest_file = os.path.join(args.output, outputs.RUN_MANIFEST_JSON)
    try:
        manifest.write_run_manifest(manifest_file, manifest.build_manifest(
            args, argv, tool_probes, started_at, status, exit_code, error,
            state.manifest_entries()))
    except OSError as e:
        logger.warning(f"Could not write {manifest_file}: {e}")


if __name__ == "__main__":
    main()
# vim: ts=4:et
