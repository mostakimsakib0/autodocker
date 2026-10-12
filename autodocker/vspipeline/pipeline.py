#!/usr/bin/env python3
"""The pipeline stages run by runner.py, one function per stage.

Each stage takes the validated command-line options (``args``, the
argparse namespace runner.py builds) and the outputs of the earlier stages,
and records what run_manifest.json reports in the shared :class:`RunState`;
runner.py writes the manifest from it whatever stage the run ended in.

  0. remove_stale_outputs
  1. prepare_ligands                 [1/5]
  2. prepare_receptor, docking_box,  [2/5]
     keep_hetatm_residues, flexible_receptor
  3. dock, consensus_rescoring       [3/5]
  4. cluster_poses                   [4/5]
  5. write_results, count_failures,  [5/5]
     log_summary

Pipeline modules are referenced through their module (``docking.dock_all``)
so a test can replace a stage's helper with ``monkeypatch.setattr``.
"""
import json
import logging
import os
import shutil
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from vspipeline import (
    consensus,
    defaults,
    docking,
    flexible,
    fpocket,
    library,
    ligand_box,
    outputs,
    pdb_records,
    pocket,
    protein,
    replicates,
    report,
    results,
    utils,
)
from vspipeline.version import __version__

logger = logging.getLogger(__name__)

# Residues within this distance (Angstrom) of the box center are candidates
# for --auto-flexible.
AUTO_FLEX_RADIUS = defaults.AUTO_FLEX_RADIUS
# Box edge lengths (Angstrom) outside this range get a warning.
TYPICAL_GRID_SIZE_RANGE = (15.0, 40.0)

# (center x, y, z, size x, y, z) of the docking box, Angstrom, docking frame.
Box = Tuple[float, float, float, float, float, float]


@dataclass
class RunState:
    """What the run has produced so far, for run_manifest.json (see
    ``manifest_entries``). Stages fill it in as they go, so a run that fails
    or is interrupted still records how far it got."""
    # Set when docking starts: an interrupted run then tells how to resume.
    docking_started: bool = False
    library: Optional["library.LibraryManager"] = None
    protein: Optional["protein.ProteinPreparation"] = None
    # Open Babel version and OB_RANDOM_SEED policy (manifest.check_openbabel).
    openbabel: Optional[Dict] = None
    # The residues actually built as flexible side chains ([] = rigid
    # docking); None until the flexible-receptor stage runs.
    flexible_residues: Optional[List[str]] = None
    threads_per_process: Optional[Dict] = None
    engine: Optional[docking.EngineConfig] = None
    run_fingerprint: Optional[str] = None
    consensus: Optional[Dict] = None
    # --box-ligand/--box-from-het: the reference ligand (set with the
    # receptor) and the manifest's "box" entry (set with the box).
    reference: Optional["ligand_box.ReferenceBox"] = None
    box: Optional[Dict] = None
    # --replicates N > 1: replicate_summary.json (seeds, spread, Spearman rho).
    replicates: Optional[Dict] = None
    # run_manifest.json "inputs": receptor and --box-ligand (manifest.run_inputs,
    # hashed when the run starts); the ligands come from ``library``.
    inputs: Optional[Dict] = None
    docking_failures: Optional[Dict] = None
    failures: Optional[Dict] = None

    def manifest_entries(self) -> Dict:
        """The run_manifest.json entries describing this run's progress."""
        extra = {}
        if self.inputs is not None:
            ligands = dict(sorted((self.library.inputs if self.library is not None
                                   else {}).items()))
            extra["inputs"] = {**self.inputs, "ligands": ligands}
        if self.library is not None:
            extra["ligand_preparation"] = {"failed": len(self.library.failed),
                                           "superseded": len(self.library.superseded)}
        if self.docking_failures is not None:
            extra["docking_failures"] = self.docking_failures
        if self.failures is not None:
            # What --max-failure-fraction was applied to (preparation + docking).
            extra["failures"] = self.failures
        if self.threads_per_process is not None:
            # Threads of each Vina/QuickVina/SMINA process (--cpu); null = all cores.
            extra["threads_per_process"] = self.threads_per_process
        if self.run_fingerprint:
            # The checkpoint fingerprint: receptor, grid, flex file, docking
            # parameters and the primary engine binary (see docking._run_fingerprint).
            extra["run_fingerprint"] = self.run_fingerprint
        if self.engine is not None:
            extra["docking_engines"] = [
                docking._engine_identity(path) for path, _ in self.engine.docking_engines()]
        if self.consensus is not None:
            extra["consensus"] = self.consensus
        if self.box is not None:
            # --box-ligand/--box-from-het: reference, rule, center in both frames.
            extra["box"] = self.box
        if self.replicates is not None:
            # --replicates N > 1: seeds, median SD, Spearman rho of replicates 0/1.
            extra["replicates"] = self.replicates
        if self.flexible_residues is not None:
            # The residues docked as flexible side chains ([] = rigid).
            extra["flexible_residues"] = list(self.flexible_residues)
        templates = getattr(self.protein, "residue_templates", None)
        if templates is not None:
            # Trp/His ring hydrogens of the receptor (residue_templates.py).
            extra["receptor_residue_templates"] = templates
        # Open Babel version and OB_RANDOM_SEED policy (manifest.check_openbabel).
        extra.update(self.openbabel or {})
        return extra


@dataclass
class FlexibleReceptor:
    """Output of the flexible-receptor stage. Both paths are None for rigid
    docking."""
    flex_file: Optional[str] = None       # flexible side chains (flex.pdbqt)
    rigid_receptor: Optional[str] = None  # receptor without them (rigid.pdbqt)


@dataclass
class DockingResults:
    """Output of the docking stage."""
    results: List[Tuple[str, Optional[float]]]  # (ligand, score or None)
    metrics: Dict[str, Dict]
    checkpoint: Optional[docking.DockingCheckpoint] = None  # not with --smina-only
    receptor: str = ""                    # the receptor PDBQT docked against
    flex_file: Optional[str] = None
    scorer_agreement: Optional[Dict] = None  # --consensus summary


def uses_fpocket(args) -> bool:
    """True when the run detects pockets with fpocket: not with --no-fpocket,
    and not when --grid-center and --grid-size give the whole box. With
    --box-ligand/--box-from-het fpocket only reports which pocket holds the
    reference ligand (--no-fpocket skips that)."""
    return not (args.no_fpocket
                or (args.grid_center is not None and args.grid_size is not None))


def flexible_residue_keys(value: str) -> List[str]:
    """The residue keys of a --flexible-residues value ('A45, A78' -> [A45, A78])."""
    return [r.strip() for r in value.split(",") if r.strip()]


def threads_per_process(args) -> Optional[int]:
    """Threads per engine process (Vina/QuickVina --cpu, SMINA --cpu):
    --threads, else the available cores split between the -p workers when
    -p > 1, else None (the engine uses every core)."""
    if args.threads:
        return args.threads
    if args.processes > 1:
        return max(1, utils.available_cpus() // args.processes)
    return None


def keep_residue_key(residue: Dict):
    """(chain, resnum, resname) of a detected water/metal/cofactor record."""
    return (residue["chain"], residue["resnum"],
            residue.get("resname") or residue.get("name"))


def failed_docking(scores, metrics_dict) -> List[str]:
    """Ligands whose docking failed: Status FAILED (an OUT_OF_WINDOW result
    is a valid engine score, not a failure)."""
    failed = []
    for name, score in scores:
        if score is not None:
            continue
        status = docking.Status.of((metrics_dict or {}).get(name, {}).get("status"), score)
        if status is not docking.Status.OUT_OF_WINDOW:
            failed.append(name)
    return failed


def failure_counts(n_prep_failed: int, n_dock_failed: int, n_docked: int,
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


# ---------------------------------------------------------------------------
# Stage 0: the output directory
# ---------------------------------------------------------------------------

def remove_stale_outputs(outdir: str) -> List[str]:
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


# ---------------------------------------------------------------------------
# Stage 1: ligands
# ---------------------------------------------------------------------------

def prepare_ligands(args, state: RunState) -> List[str]:
    """[1/5] The ligand library (local files or PubChem), prepared as PDBQT.
    Returns the ligand PDBQT paths; raises when there is none."""
    logger.info("[1/5] Ligand library preparation")
    if args.library != "local" and args.ligands:
        logger.info(f"[*] -l {args.ligands} is not used with --library {args.library}")
    lib_manager = library.LibraryManager(
        args.output,
        args.ligands if args.library == "local" else None,
        minimize=not args.no_minimize,
        minimize_steps=args.minimize_steps,
        seed=args.seed,
    )
    state.library = lib_manager
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
    return ligands


# ---------------------------------------------------------------------------
# Stage 2: receptor, box, kept HETATM residues, flexible side chains
# ---------------------------------------------------------------------------

def reference_box(args) -> Tuple[Optional["ligand_box.ReferenceBox"], str]:
    """(ReferenceBox or None, receptor PDB to prepare). --box-from-het: the
    group is written alone to reference_ligand.pdb and the receptor is
    prepared from a copy without it (receptor_no_reference.pdb), so no later
    step (--keep-hetero, --detect-cofactors, ...) can bring it back."""
    if args.box_ligand:
        return ligand_box.box_from_ligand_file(
            args.box_ligand, args.padding, defaults.LIGAND_BOX_MIN_SIZE), args.receptor
    if not args.box_from_het:
        return None, args.receptor
    reference = ligand_box.box_from_het(args.receptor, args.box_from_het, args.padding,
                                        defaults.LIGAND_BOX_MIN_SIZE)
    stripped = os.path.join(args.output, outputs.RECEPTOR_NO_REFERENCE_PDB)
    removed = ligand_box.write_receptor_without_group(args.receptor, stripped,
                                                      reference.group)
    ligand_box.write_reference_pdb(os.path.join(args.output, outputs.REFERENCE_LIGAND_PDB),
                                   args.receptor, reference.group)
    logger.info(f"[*] --box-from-het {reference.reference}: removed its {removed} records "
                f"from the receptor ({outputs.RECEPTOR_NO_REFERENCE_PDB}); it is not part "
                "of the receptor, also with --keep-hetero")
    return reference, stripped


def prepare_receptor(args, state: RunState) -> "protein.ProteinPreparation":
    """[2/5] Chain selection, cleaning and protonation of the receptor (with
    --box-from-het, of the receptor without the reference group; the
    reference ligand is recorded in ``state.reference``)."""
    logger.info("[2/5] Receptor preparation and docking box")
    reference, receptor_input = reference_box(args)
    state.reference = reference
    protein_prep = protein.ProteinPreparation(receptor_input, args.output, seed=args.seed)
    state.protein = protein_prep
    protein_prep.validate()

    chain = args.chain or protein_prep.select_chain()
    protein_prep.prepare_receptor(chain, keep_hetero=args.keep_hetero, ph=args.ph)
    return protein_prep


def _reference_pockets(protein_prep: "protein.ProteinPreparation",
                       reference: "ligand_box.ReferenceBox") -> Dict:
    """Run fpocket on the prepared receptor only to report which pocket
    holds the reference ligand (pocket_summary.csv, manifest). A failing
    fpocket does not stop the run: the box does not depend on it."""
    try:
        target = protein_prep.pdb_clean
        pockets = fpocket.rank_pockets(fpocket.collect_pockets(
            fpocket.run_fpocket(target), target))
        fpocket.write_pocket_summary(protein_prep.pocket_summary_file, pockets)
    except RuntimeError as e:
        logger.warning(f"[!] fpocket (run to compare with the reference ligand) failed: {e}")
        return {"error": str(e)}
    rows = ligand_box.pocket_overlaps(list(reference.coords), pockets, reference.center)
    summary = ligand_box.fpocket_summary(rows)
    best, top = summary["best_overlapping_pocket"], summary["top_ranked_pocket"]
    if best is None:
        logger.warning(f"[!] No fpocket pocket overlaps the reference ligand (fpocket's "
                       f"top pocket: pocket{top['pocket']})")
    elif summary["top_pocket_contains_reference"]:
        logger.info(f"[*] fpocket's top pocket (pocket{top['pocket']}) holds "
                    f"{best['overlap_fraction']:.0%} of the reference ligand")
    else:
        logger.warning(
            f"[!] The reference ligand is in pocket{best['pocket']} (fpocket rank "
            f"{best['rank']}, {best['overlap_fraction']:.0%} of its heavy atoms), not in "
            f"fpocket's top pocket{top['pocket']}: an fpocket-selected box would have "
            "missed this site")
    return summary


def _reference_ligand_box(args, reference: "ligand_box.ReferenceBox",
                          protein_prep: "protein.ProteinPreparation") -> Tuple[Box, Dict]:
    """--box-ligand/--box-from-het: the box (docking frame) around the
    reference ligand, with its checks and log lines. Returns (box, details
    for the manifest's "box" entry: receptor contact and fpocket report)."""
    offset = protein_prep.frame_offset()
    for line in reference.log_lines(offset):
        logger.info(line)
    xs, ys, zs = pdb_records._heavy_atom_coords(protein_prep.pdb_clean)
    receptor_atoms = list(zip(xs, ys, zs))
    if not ligand_box.points_in_box(receptor_atoms, reference.center, reference.size):
        raise ValueError(
            f"The box around the reference ligand ({reference.describe()}) contains no "
            "receptor atom of the docked chain(s): is the ligand in the receptor PDB's "
            "frame, and is its site in the chains selected with --chain?")
    contact = ligand_box.min_distance(list(reference.coords), receptor_atoms)
    if contact > ligand_box.CONTACT_WARN_DISTANCE:
        logger.warning(
            f"[!] The reference ligand is {contact:.1f} A from the nearest receptor heavy "
            f"atom (> {ligand_box.CONTACT_WARN_DISTANCE:g} A): check that it is in the "
            "receptor PDB's frame and that --chain selects the chains of its site")
    if args.box_ligand and args.keep_hetero:
        same = ligand_box.matching_het_groups(list(reference.coords),
                                              ligand_box.het_groups(args.receptor))
        if same:
            logger.warning(
                f"[!] The reference ligand occupies the receptor's HETATM group "
                f"{same[0].spec}, which --keep-hetero keeps in the receptor; use "
                f"--box-from-het {same[0].spec} to remove it")
    details = {"receptor_contact": contact,
               "fpocket": (_reference_pockets(protein_prep, reference)
                           if uses_fpocket(args) else None)}
    cx, cy, cz = reference.docking_center(offset)
    sx, sy, sz = reference.size
    return (cx, cy, cz, sx, sy, sz), details


def docking_box(args, state: RunState, protein_prep: "protein.ProteinPreparation") -> Box:
    """The docking box (docking frame): from fpocket pockets, around the
    whole receptor (--no-fpocket), around a reference ligand
    (--box-ligand/--box-from-het, recorded in ``state.box``), and/or
    --grid-center/--grid-size. Writes grid.conf, grid_box.pml and
    frame_offset.txt."""
    reference = state.reference
    reference_details: Dict = {}
    if reference is not None:
        (cx, cy, cz, sx, sy, sz), reference_details = _reference_ligand_box(
            args, reference, protein_prep)
    elif args.no_fpocket:
        cx, cy, cz, sx, sy, sz = protein_prep.centroid_grid(padding=args.padding)
        logger.info("[*] Docking box: whole-receptor box (--no-fpocket)")
    elif not uses_fpocket(args):
        # Both values are set below; no pocket is needed.
        cx = cy = cz = sx = sy = sz = 0.0
        logger.info("[*] Docking box: --grid-center and --grid-size (fpocket not run)")
    else:
        cx, cy, cz, sx, sy, sz = protein_prep.detect_pocket(
            pocket_spec=args.pockets, padding=args.padding)

    if reference is None and (uses_fpocket(args) or args.no_fpocket):
        logger.info(f"[*] Grid center: {cx:.2f}, {cy:.2f}, {cz:.2f}")
        logger.info(f"[*] Grid size: {sx:.2f}, {sy:.2f}, {sz:.2f}")

    if args.grid_center:
        cx, cy, cz = docking._parse_grid_triplet(args.grid_center, "grid-center")
        if args.grid_center_frame == "input":
            # receptor.pdbqt and grid.conf are centered at the origin;
            # move the user's input-PDB coordinates into that frame.
            dx, dy, dz = protein_prep.frame_offset()
            logger.info(f"[*] Grid center given in input-PDB coordinates: "
                        f"{cx:.2f}, {cy:.2f}, {cz:.2f}")
            cx, cy, cz = cx - dx, cy - dy, cz - dz
        logger.info(f"[*] Grid center overridden (docking frame): {cx:.2f}, {cy:.2f}, {cz:.2f}")
    if args.grid_size:
        sx, sy, sz = docking._parse_grid_triplet(args.grid_size, "grid-size")
        logger.info(f"[*] Grid size overridden: {sx:.2f}, {sy:.2f}, {sz:.2f}")
    if reference is not None:
        state.box = reference.manifest(
            protein_prep.frame_offset(),
            grid_size_override=(sx, sy, sz) if args.grid_size else None,
            fpocket=reference_details.get("fpocket"),
            receptor_contact=reference_details.get("receptor_contact"))

    low, high = TYPICAL_GRID_SIZE_RANGE
    for axis, size in (("x", sx), ("y", sy), ("z", sz)):
        if size < low or size > high:
            logger.warning(f"[!] Grid size_{axis}={size:.2f}A outside typical range "
                           f"({low:g}-{high:g}A)")

    protein_prep.write_grid(cx, cy, cz, sx, sy, sz)
    return cx, cy, cz, sx, sy, sz


def keep_hetatm_residues(args, protein_prep: "protein.ProteinPreparation",
                         box: Box) -> None:
    """--keep-waters, --detect-metals, --detect-cofactors: append the
    detected residues of the docked chain(s) to the receptor PDBQT."""
    if not (args.keep_waters or args.detect_metals or args.detect_cofactors):
        return
    logger.info("[2/5] Receptor waters, metals and cofactors")
    cx, cy, cz = box[:3]
    # Detect on the ORIGINAL PDB: HETATM records (waters, metals,
    # cofactors) are stripped from pdb_clean during receptor cleaning.
    source = protein_prep.pdb_file
    keep_residues = []

    if args.keep_waters:
        logger.info("[*] Detecting water molecules near binding site...")
        # (cx, cy, cz) is in the centered receptor frame; waters are
        # read from the original PDB, so compare in that frame.
        dx, dy, dz = protein_prep.frame_offset()
        waters = pocket.detect_water_molecules(source, (cx + dx, cy + dy, cz + dz),
                                               distance_threshold=args.water_distance)
        if waters:
            logger.info(f"[*] Found {len(waters)} water residue(s) within "
                        f"{args.water_distance}A of the binding site; keeping them "
                        "in receptor")
            keep_residues += [keep_residue_key(w) for w in waters]
        else:
            logger.info("[*] No waters detected near binding site")

    if args.detect_metals:
        logger.info("[*] Detecting metal ions...")
        metals = pocket.detect_metal_ions(source)
        if metals:
            logger.info(f"[*] Found {len(metals)} metal ion(s); keeping them in receptor:")
            for metal in metals:
                logger.info(f"   - {metal['name']} ({metal['element']}, atom "
                            f"{metal['atom']}) chain {metal['chain']} "
                            f"residue {metal['resnum']}")
            keep_residues += [keep_residue_key(m) for m in metals]
        else:
            logger.info("[*] No metal ions detected")

    if args.detect_cofactors:
        logger.info("[*] Detecting cofactors...")
        cofactors = pocket.detect_cofactors(source)
        if cofactors:
            logger.info(f"[*] Found {len(cofactors)} cofactor(s); keeping them in receptor:")
            for cof in cofactors:
                logger.info(f"   - {cof['name']} chain {cof['chain']} residue {cof['resnum']}")
            keep_residues += [keep_residue_key(c) for c in cofactors]
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


def flexible_receptor(args, state: RunState, protein_prep: "protein.ProteinPreparation",
                      box: Box) -> FlexibleReceptor:
    """--flexible-residues / --auto-flexible: build the flexible side chains
    and the rigid receptor without them. Records the residues built in
    ``state.flexible_residues``. Turns --flexibility on when either option
    names residues."""
    # The residues actually built as flexible side chains ([] = rigid docking);
    # recorded in run_manifest.json and the closing summary.
    flexible_residues: List[str] = []
    state.flexible_residues = flexible_residues
    flex = FlexibleReceptor()
    if args.flexibility is None and (args.flexible_residues or args.auto_flexible):
        args.flexibility = 1
    if not (args.flexibility and args.flexibility > 0):
        return flex
    logger.info("[2/5] Flexible receptor setup")
    requested: List[str] = []
    if args.flexible_residues:
        requested = flexible_residue_keys(args.flexible_residues)
        logger.info(f"[*] Using specified flexible residues: {', '.join(requested)}")
    elif args.auto_flexible:
        requested = flexible.detect_flexible_residues(
            protein_prep.receptor_pdbqt, box[:3],
            radius=AUTO_FLEX_RADIUS, max_residues=args.auto_flexible)
        logger.info(f"[*] Auto-selected {len(requested)} flexible residues "
                    "near binding site")
    else:
        logger.warning("[!] Flexibility enabled but no residues specified. "
                       "Use --flexible-residues or --auto-flexible")

    built_residues: List[str] = []
    if requested:
        logger.info(f"[*] Flexible residues: {', '.join(requested)}")
        flex.flex_file = flexible.build_flexible_residue_pdbqt(
            protein_prep.receptor_pdbqt, requested, args.output,
            built=built_residues)
    if args.flexible_residues:
        # Residues named by the user must all be built: docking with
        # fewer (or none) would not be the requested experiment.
        built_keys = {flexible.residue_identity(k) for k in built_residues}
        not_built = [k for k in requested if flexible.residue_identity(k) not in built_keys]
        if not_built:
            raise RuntimeError(
                f"--flexible-residues: could not build {', '.join(not_built)} as "
                f"flexible side chain(s) in {protein_prep.receptor_pdbqt} (see the "
                "warnings above: missing atoms or no rotatable side chain, e.g. "
                "GLY, ALA, PRO); remove them from --flexible-residues")
    elif args.auto_flexible and not flex.flex_file:
        logger.warning(
            f"[!] --auto-flexible {args.auto_flexible}: no residue within "
            f"{AUTO_FLEX_RADIUS:g} A of the box center could be made flexible; "
            "docking with a RIGID receptor")
    if flex.flex_file:
        # Vina requires flexible residues to be excluded from the
        # rigid receptor passed to --receptor.
        flexible_residues.extend(built_residues)
        flex.rigid_receptor = flexible._remove_residues_from_pdbqt(
            protein_prep.receptor_pdbqt, built_residues,
            os.path.join(args.output, outputs.RIGID_RECEPTOR_PDBQT),
            flex_file=flex.flex_file)
    return flex


# ---------------------------------------------------------------------------
# Stage 3: docking and consensus rescoring
# ---------------------------------------------------------------------------

def dock(args, state: RunState, ligands: List[str],
         protein_prep: "protein.ProteinPreparation", box: Box,
         flex: FlexibleReceptor) -> DockingResults:
    """[3/5] Dock every ligand with Vina (resumable, see docking.dock_all) or,
    with --smina-only, with SMINA."""
    logger.info("[3/5] Docking")
    vina_params = {
        'exhaustiveness': args.exhaustiveness,
        'binding_modes': args.binding_modes,
        'energy_range': args.energy_range,
        'seed': args.seed,
        'min_valid_affinity': args.min_valid_affinity,
        'extra_args': args.vina_extra,
        'threads': threads_per_process(args),
    }
    state.threads_per_process = {
        "threads": vina_params['threads'],
        "source": ("--threads" if args.threads else
                   "available CPUs / --processes" if vina_params['threads'] else
                   "engine default (all cores)"),
        "processes": args.processes,
    }
    receptor = flex.rigid_receptor if flex.flex_file else protein_prep.receptor_pdbqt

    state.docking_started = True
    if args.smina_only:
        # SMINA replaces Vina: dock with SMINA alone, and report SMINA's
        # own poses and failures (not Vina metadata next to SMINA scores).
        # Ligand names are unique: LibraryManager resolves same-name files.
        logger.info("[*] SMINA-only mode: docking with SMINA instead of Vina")
        scores, metrics_dict = consensus.smina_dock_all(
            receptor, ligands, args.output, box, vina_params,
            num_processes=args.processes, flex_file=flex.flex_file)
        if not any(score is not None or docking.Status.of(metrics_dict[name].get("status"))
                   is docking.Status.OUT_OF_WINDOW for name, score in scores):
            # Reported like any other failure fraction (all outputs, then
            # exit 3 unless --max-failure-fraction 1.0).
            logger.error("[!] SMINA produced no score for any ligand; see "
                         f"{outputs.DOCKED_DIR}/*{outputs.SMINA_LOG_SUFFIX}")
        checkpoint = None
    else:
        # The engine paths and timeout are passed to the workers explicitly
        # (after --vina-bin/--timeout were applied in runner.main()).
        engine = docking.EngineConfig.from_settings()
        state.engine = engine
        scores, checkpoint, metrics_dict = docking.dock_all(
            receptor, ligands, protein_prep.grid_conf,
            args.output, num_processes=args.processes, resume=not args.no_resume,
            vina_params=vina_params, flex_file=flex.flex_file, engine=engine,
            replicates=args.replicates)
        state.run_fingerprint = getattr(checkpoint, "fingerprint", None)
        if args.replicates > 1:
            state.replicates = replicate_outputs(args, scores, metrics_dict)
    return DockingResults(results=scores, metrics=metrics_dict, checkpoint=checkpoint,
                          receptor=receptor, flex_file=flex.flex_file)


def replicate_outputs(args, scores, metrics_dict) -> Dict:
    """replicate_scores.csv and replicate_summary.json of a run with
    --replicates N > 1; returns the summary (run_manifest.json)."""
    summary = replicates.summary(metrics_dict, args.replicates, args.seed)
    replicates.write_replicate_scores(
        os.path.join(args.output, outputs.REPLICATE_SCORES_CSV), args.output,
        sorted(name for name, _ in scores), metrics_dict)
    replicates.write_summary(os.path.join(args.output, outputs.REPLICATE_SUMMARY_JSON),
                             summary)
    for line in replicates.log_lines(summary):
        (logger.warning if line.startswith("[!]") else logger.info)(line)
    return summary


def consensus_rescoring(args, state: RunState, docked: DockingResults) -> None:
    """--consensus: rescore Vina's top poses with SMINA and write
    consensus_ranking.csv and consensus_summary.json."""
    if not args.consensus:
        return
    logger.info("[3/5] Consensus: SMINA rescoring of the Vina top poses")
    smina_scores, scoring = consensus.rescore_with_smina(
        docked.receptor, docked.results, args.output, metrics_dict=docked.metrics,
        num_processes=args.processes, flex_file=docked.flex_file)
    docked.scorer_agreement = {
        **consensus.write_consensus_outputs(args.output, docked.results, smina_scores,
                                            scoring, metrics_dict=docked.metrics),
        "scoring": scoring}
    state.consensus = docked.scorer_agreement


# ---------------------------------------------------------------------------
# Stage 4: pose analysis
# ---------------------------------------------------------------------------

def cluster_poses(args, docked: DockingResults) -> str:
    """[4/5] --cluster-poses: write pose_clusters.json. Returns its path
    (written or not)."""
    logger.info("[4/5] Pose analysis")
    checkpoint = docked.checkpoint
    if not args.no_resume and checkpoint is not None and checkpoint.completed:
        logger.info(f"[*] Checkpoint updated: {len(checkpoint.completed)} successful dockings")
    clusters_file = os.path.join(args.output, outputs.POSE_CLUSTERS_JSON)
    if args.cluster_poses:
        logger.info("[*] Performing pose clustering...")
        clusters = report.cluster_poses(
            os.path.join(args.output, outputs.DOCKED_DIR), rmsd_threshold=args.rmsd_threshold,
            ligands=[name for name, _ in docked.results],
            # --replicates: each ligand's reported (best-replicate) pose file.
            pose_files=({name: (docked.metrics.get(name) or {}).get("dock_file")
                         for name, _ in docked.results} if args.replicates > 1 else None))
        clusters_data = {
            ligand: [{'id': c['id'], 'pose_count': len(c['poses'])} for c in ligand_clusters]
            for ligand, ligand_clusters in clusters.items()}
        try:
            with open(clusters_file, 'w', encoding="utf-8") as f:
                json.dump(clusters_data, f, indent=2, sort_keys=True)
            logger.info(f"[*] Clustering saved to: {clusters_file}")
        except OSError as e:
            logger.warning(f"Could not save clustering data: {e}")
    return clusters_file


# ---------------------------------------------------------------------------
# Stage 5: outputs
# ---------------------------------------------------------------------------

def write_results(args, state: RunState, protein_prep: "protein.ProteinPreparation",
                  docked: DockingResults) -> Tuple["results.ResultsAnalyzer", bool]:
    """[5/5] ranking.csv, the text reports and (unless --no-html-report) the
    HTML report. Returns (the analyzer that wrote them, HTML report written)."""
    logger.info("[5/5] Saving results")
    analyzer = results.ResultsAnalyzer(args.output, top_n=args.top_n,
                                       replicates=args.replicates)
    analyzer.save_ranking(docked.results, mode="extended", metrics_dict=docked.metrics)
    analyzer.save_metrics_report(docked.results, docked.metrics)

    # The HTML report is generated AFTER ranking.csv exists.
    html_ok = False
    if args.html_report:
        logger.info("[*] Generating HTML report...")
        report_meta = {
            "Protein (receptor)": os.path.basename(args.receptor),
            "Ligand library": args.library,
            "Ligands screened": len(docked.results),
            "Top N reported": args.top_n,
            "Exhaustiveness": args.exhaustiveness,
            "Seed": args.seed,
        }
        if args.replicates > 1:
            report_meta["Replicates"] = (
                f"{args.replicates} (seeds {args.seed}..{args.seed + args.replicates - 1}; "
                "best replicate reported)")
        if state.box is not None:
            report_meta["Docking box"] = f"reference ligand {state.box['reference']}"
        html_ok = report.generate_html_report(
            os.path.join(args.output, outputs.RANKING_CSV),
            os.path.join(args.output, outputs.DOCKED_DIR),
            os.path.join(args.output, outputs.HTML_REPORT),
            receptor_pdbqt=protein_prep.receptor_pdbqt,
            meta=report_meta,
            grid_file=protein_prep.grid_conf,
            scorer_agreement=docked.scorer_agreement,
            rmsd_threshold=args.rmsd_threshold,
            top_n=args.top_n)
    return analyzer, html_ok


def count_failures(args, state: RunState, docked: DockingResults) -> Dict:
    """The failed ligands (preparation + docking) that --max-failure-fraction
    is applied to; recorded in ``state``. Returns ``failure_counts``."""
    failed = failed_docking(docked.results, docked.metrics)
    state.docking_failures = {"failed": len(failed), "docked": len(docked.results)}
    state.failures = failure_counts(len(state.library.failed), len(failed),
                                    len(docked.results), args.max_failure_fraction)
    return state.failures


def log_summary(args, state: RunState, protein_prep: "protein.ProteinPreparation",
                docked: DockingResults, analyzer: "results.ResultsAnalyzer",
                html_ok: bool, clusters_file: str, too_many_failures: bool) -> None:
    """The closing summary: outcome and the files written."""
    counts = state.failures
    lib_manager = state.library
    logger.info("=" * 70)
    if not too_many_failures:
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
    if args.replicates > 1:
        written += [("Replicate scores",
                     os.path.join(args.output, outputs.REPLICATE_SCORES_CSV)),
                    ("Replicate summary",
                     os.path.join(args.output, outputs.REPLICATE_SUMMARY_JSON))]
    if args.box_from_het:
        written.append(("Reference ligand",
                        os.path.join(args.output, outputs.REFERENCE_LIGAND_PDB)))
    if args.html_report:
        if html_ok:
            written.append(("HTML report", os.path.join(args.output, outputs.HTML_REPORT)))
        else:
            logger.warning("[!] HTML report could not be written (see warnings above)")
    written.append(("Run manifest", os.path.join(args.output, outputs.RUN_MANIFEST_JSON)))
    for label, path in written:
        if os.path.exists(path) or label == "Run manifest":
            logger.info(f"  - {label}: {path}")
    if args.flexibility and args.flexibility > 0:
        flexible_residues = state.flexible_residues
        logger.info(f"  - Flexible residues: "
                    f"{', '.join(flexible_residues) if flexible_residues else 'none (rigid docking)'}")
    logger.info(f"Total compounds docked: {len(docked.results)} "
                f"({counts['docking_failed']} failed)")
    if lib_manager.failed:
        logger.info(f"Ligands that could not be prepared: {len(lib_manager.failed)}")
    logger.info("=" * 70)
