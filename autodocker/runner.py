#!/usr/bin/env python3
"""
Virtual Screening Pipeline

PHASE 1: Advanced Result Analysis
  - HTML professional reports with charts
  - Pose clustering (RMSD-based)
  - Pose diversity analysis
  - Binding mode statistics

PHASE 2: Smart Preprocessing
  - Water molecule handling (keep strategic waters)
  - Metal ion detection & parameterization
  - Cofactor handling
  - pH-based protonation (pKa calculation)

PHASE 3: Flexible Receptor Docking
  - Select rotatable residues in binding site
  - Flexible side-chain docking
  - Residue flexibility control (1-10)

PHASE 4: Consensus Scoring
  - Multi-engine docking (Vina + SMINA)
  - Consensus ranking
  - Score agreement analysis

Original Features:
  - PubChem library download (FDA & custom libraries)
  - SimScore calculation (RMSD-based pose consistency)
  - ADMET filtering (Lipinski's rule of five)
  - Advanced Vina parameters

Workflow:
  1. Library preparation (PubChem download with ADMET filtering)
  2. Protein preparation (chain selection + fpocket)
  3. Smart preprocessing (optional: water/metal/cofactor)
  4. Docking (flexible or rigid, single or multi-engine)
  5. Advanced results analysis (HTML reports + clustering)
"""
import argparse
import csv
import json
import logging
import os
import shutil
import sys
from multiprocessing import cpu_count as mp_cpu_count
from typing import Optional

__version__ = "1.0.2"

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


# Tool paths (prefer modern Vina for consistency)
OBABEL = find_tool("obabel", "OpenBabel", "OpenBabel.exe")
# Prefer modern Vina for reproducibility/consistency
VINA_PRIMARY = find_tool("vina", "vina.exe", "qvina02", "qvina2", "qvina")
if os.environ.get("VS_ENABLE_QVINA_FALLBACK", "0") == "1":
    VINA_FALLBACK = find_tool("qvina02", "qvina2") or VINA_PRIMARY
else:
    VINA_FALLBACK = VINA_PRIMARY
VINA = VINA_PRIMARY  # Start with primary tool
# Optional for consensus scoring
SMINA = find_tool("smina", "smina.exe")
COMMAND_TIMEOUT = int(os.environ.get("VS_COMMAND_TIMEOUT", "900"))

# Verify tools exist - STRICT validation
missing = [name for name, path in (("obabel", OBABEL), ("vina", VINA_PRIMARY))
           if path is None]
if missing:
    raise RuntimeError(
        f"Critical tool(s) not found in PATH: {', '.join(missing)}. "
        "Ensure they are installed and in PATH.")

if SMINA is None:
    logger.info("SMINA not found - consensus scoring will be skipped")

# ---------------------------------------------------------------------------
# Pipeline internals (split of the former monolith into vspipeline/*).
# These imports re-export every name the CLI, tests, and scripts rely on;
# the moved code itself resolves shared state through THIS module at call
# time, so monkeypatching runner.<name> behaves exactly as before.
# ---------------------------------------------------------------------------

from vspipeline.net import _http_get_bytes
from vspipeline.report import (
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
from vspipeline.pocket import (
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
from vspipeline.flexible import (
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
from vspipeline.consensus import (
    _run_smina_scoring,
    consensus_rank,
    run_smina_docking,
    scorer_agreement_spearman,
)
from vspipeline.charges import (
    _assign_simple_charges,
    _ensure_pdbqt_has_charges,
    _extract_pdbqt_charge,
    _fix_pdbqt_charges,
    _pdbqt_has_atoms,
    _sanitize_receptor_pdbqt,
)
from vspipeline.admet import (
    ADMETFilter,
    _is_valid_sdf,
    _obabel_descriptors,
)
from vspipeline.checkpoint import DockingCheckpoint
from vspipeline.library import LibraryManager
from vspipeline.utils import (
    _parse_chain_selection,
    _safe_ligand_id,
    get_chains,
    parse_pdb_coords,
    run,
)
from vspipeline.protein import ProteinPreparation
from vspipeline.docking import (
    _VINA_THREADS_CACHE,
    _build_vina_command,
    _calculate_coordinate_simscore,
    _calculate_metrics,
    _calculate_protocol_simscore,
    _calculate_rmsd,
    _detect_vina_type,
    _extract_score,
    _log_tool_version,
    _parse_grid_config,
    _parse_grid_triplet,
    _parse_vina_modes,
    _tool_label,
    _vina_supports_threads,
    dock_all,
    dock_ligand,
)
from vspipeline.results import (
    ResultsAnalyzer,
    _score_csv,
    _score_sort_key,
    _score_text,
)

# =============================
# MAIN PIPELINE
# =============================


def main():
    global VINA, VINA_PRIMARY, VINA_FALLBACK, COMMAND_TIMEOUT

    parser = argparse.ArgumentParser(
        description="Virtual Screening Pipeline v2.0 - ENHANCED with 4 Phases:\n"
        "  PHASE 1: Advanced result analysis (HTML reports, clustering)\n"
        "  PHASE 2: Smart preprocessing (water/metal/cofactor handling)\n"
        "  PHASE 3: Flexible receptor docking\n"
        "  PHASE 4: Consensus scoring (Vina + SMINA)\n"
        "Original: PubChem library integration, SimScore, ADMET filtering",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-r", "--receptor", required=True,
                        help="Protein PDB file")
    parser.add_argument("-l", "--ligands", required=True,
                        help="Ligands .sdf file or directory with .sdf files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory")

    parser.add_argument("--library", choices=["fda", "custom", "local"], default="local",
                        help="Library: fda (PubChem approved drugs), custom (PubChem filtered), local (manual SDF)")
    parser.add_argument("--mw-min", type=int, default=200, help="Min MW")
    parser.add_argument("--mw-max", type=int, default=500, help="Max MW")
    parser.add_argument("--logp-max", type=float, default=5.0,
                        help="Max LogP (drug-likeness)")
    parser.add_argument("--no-admet", action="store_true",
                        help="Skip ADMET (Lipinski) filtering")
    parser.add_argument("--no-minimize", action="store_true",
                        help="Skip MMFF94 ligand minimization")
    parser.add_argument("--minimize-steps", type=int,
                        default=250, help="MMFF94 minimization steps")

    parser.add_argument(
        "--chain", help="Protein chain(s), e.g. A, A,B, or all (auto-detect if not specified)")
    parser.add_argument("--keep-hetero", action="store_true",
                        help="Keep matching HETATM records during receptor cleaning")
    parser.add_argument(
        "--pockets", help="Pocket number(s) from fpocket, e.g. 1 or 2,3. Default: best druggability")
    parser.add_argument("--padding", type=float, default=6.0,
                        help="Grid padding around selected pocket(s), Angstrom")
    parser.add_argument("--no-fpocket", action="store_true",
                        help="Skip fpocket, use default grid")
    parser.add_argument("--grid-center", default=None,
                        help="Explicit grid center 'X,Y,Z' (Angstrom), e.g. 12.3,4.5,6.7. Overrides auto-detection")
    parser.add_argument("--grid-size", default=None,
                        help="Explicit grid box size 'SX,SY,SZ' (Angstrom), e.g. 22.5,22.5,22.5. Overrides auto-detection")

    parser.add_argument("-p", "--processes", type=int,
                        default=1, help="Parallel processes")
    parser.add_argument("--exhaustiveness", type=int, default=8,
                        help="Vina exhaustiveness (1-32, higher=thorough)")
    parser.add_argument("--binding-modes", type=int,
                        default=9, help="Vina binding modes")
    parser.add_argument("--energy-range", type=float,
                        default=3.0, help="Vina energy range (kcal/mol)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility")
    parser.add_argument("--threads", type=int, default=None,
                        help="Vina threads (default: all cores; e.g. 2 to bound CPU when docking in parallel)")
    parser.add_argument("--min-valid-affinity", type=float, default=-1.0,
                        help="Reject docking scores >= this threshold (kcal/mol)")
    parser.add_argument(
        "--vina-bin", help="Docking binary to use instead of auto-detected Vina/QuickVina")
    parser.add_argument("--vina-extra", action="append", default=[],
                        help="Extra argument string passed to the docking binary; repeat as needed")
    parser.add_argument("--timeout", type=int, default=COMMAND_TIMEOUT,
                        help="Per-command timeout in seconds")
    parser.add_argument("--no-resume", action="store_true",
                        help="Start fresh (ignore checkpoint)")
    parser.add_argument("--top-n", type=int, default=20,
                        help="Number of compounds in Top_hits.txt")

    # ===== PHASE 1: Advanced Result Analysis =====
    parser.add_argument("--html-report", dest="html_report", action="store_true",
                        default=True,
                        help="Generate professional HTML report with charts (default: on)")
    parser.add_argument("--no-html-report", dest="html_report", action="store_false",
                        help="Disable the HTML report")
    parser.add_argument("--cluster-poses", action="store_true",
                        help="Cluster similar poses by RMSD")
    parser.add_argument("--rmsd-threshold", type=float,
                        default=2.0, help="RMSD clustering threshold (Angstrom)")

    # ===== PHASE 2: Smart Preprocessing =====
    parser.add_argument("--keep-waters", action="store_true",
                        help="Keep water molecules near binding site")
    parser.add_argument("--detect-metals", action="store_true",
                        help="Detect and parameterize metal ions")
    parser.add_argument("--detect-cofactors", action="store_true",
                        help="Detect cofactors in protein")
    parser.add_argument("--water-distance", type=float, default=4.0,
                        help="Distance threshold for water detection (Å)")

    # ===== PHASE 3: Flexible Receptor Docking =====
    parser.add_argument("--flexibility", type=int, choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                        help="Flexibility level (0=rigid, 1=slight, 10=very flexible). Enables flexible side-chains")
    parser.add_argument("--flexible-residues",
                        help="Specify flexible residues, e.g. 'A45,A78,B102'")
    parser.add_argument("--auto-flexible", type=int,
                        help="Auto-select N closest residues to binding site")

    # ===== PHASE 4: Consensus Scoring =====
    parser.add_argument("--consensus", action="store_true",
                        help="Run consensus scoring (Vina + SMINA)")
    parser.add_argument("--smina-only", action="store_true",
                        help="Use SMINA scoring instead of Vina")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="Debug logging")
    parser.add_argument("--version", action="version",
                        version=f"AutoDocker {__version__}")

    args = parser.parse_args()

    if args.vina_bin:
        VINA = args.vina_bin
        VINA_PRIMARY = args.vina_bin
        VINA_FALLBACK = args.vina_bin
    COMMAND_TIMEOUT = args.timeout

    if args.processes in [0, -1]:
        args.processes = mp_cpu_count()
    elif args.processes < -1:
        args.processes = mp_cpu_count() // (- args.processes)

    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    os.makedirs(args.output, exist_ok=True)

    _log_tool_version(VINA_PRIMARY, "Primary docking tool")
    if VINA_FALLBACK != VINA_PRIMARY:
        _log_tool_version(VINA_FALLBACK, "Fallback docking tool")
    if SMINA:
        _log_tool_version(SMINA, "SMINA")

    try:
        logger.info("=" * 70)
        logger.info("VIRTUAL SCREENING PIPELINE (COMPLETE)")
        logger.info("=" * 70)

        # ===== STEP 1: Library Preparation =====
        logger.info("[1/5] Library Preparation")
        lib_manager = LibraryManager(
            args.output,
            args.ligands,
            minimize=not args.no_minimize,
            minimize_steps=args.minimize_steps,
        )
        apply_admet = not args.no_admet

        if args.library == "fda":
            ligands = lib_manager.create_fda_library(apply_admet=apply_admet)
        elif args.library == "custom":
            ligands = lib_manager.create_custom_library(
                mw_min=args.mw_min, mw_max=args.mw_max, logp_max=args.logp_max,
                apply_admet=apply_admet
            )
        else:
            ligands = lib_manager._prepare_local_sdf(apply_admet=apply_admet)

        if not ligands:
            raise ValueError("No ligands found to dock")

        # ===== STEP 2: Protein Preparation =====
        logger.info("[2/5] Protein Preparation + Grid")
        protein_prep = ProteinPreparation(args.receptor, args.output)
        protein_prep.validate()

        chain = args.chain or protein_prep.select_chain()
        protein_prep.prepare_receptor(chain, keep_hetero=args.keep_hetero)

        if args.no_fpocket:
            cx, cy, cz, sx, sy, sz = protein_prep._protein_centroid_grid(
                padding=args.padding)
            logger.warning(
                "[!] Using centroid fallback grid (fpocket skipped)")
        else:
            cx, cy, cz, sx, sy, sz = protein_prep.detect_pocket(
                pocket_spec=args.pockets,
                padding=args.padding,
            )

        logger.info(f"[*] Grid center: {cx:.2f}, {cy:.2f}, {cz:.2f}")
        logger.info(f"[*] Grid size: {sx:.2f}, {sy:.2f}, {sz:.2f}")

        if args.grid_center:
            cx, cy, cz = _parse_grid_triplet(args.grid_center, "grid-center")
            logger.info(
                f"[*] Grid center overridden: {cx:.2f}, {cy:.2f}, {cz:.2f}")
        if args.grid_size:
            sx, sy, sz = _parse_grid_triplet(args.grid_size, "grid-size")
            logger.info(
                f"[*] Grid size overridden: {sx:.2f}, {sy:.2f}, {sz:.2f}")

        for axis, size in (("x", sx), ("y", sy), ("z", sz)):
            if size < 15 or size > 40:
                logger.warning(
                    f"[!] Grid size_{axis}={size:.2f}A outside typical range (15-40A)")

        protein_prep.write_grid(cx, cy, cz, sx, sy, sz)

        # ===== PHASE 2: Smart Preprocessing (Optional) =====
        if args.keep_waters or args.detect_metals or args.detect_cofactors:
            logger.info("[PHASE 2] Smart Preprocessing")
            # Detect on the ORIGINAL PDB: HETATM records (waters, metals,
            # cofactors) are stripped from pdb_clean during receptor cleaning.
            source = protein_prep.pdb_file

            keep_residues = []
            waters = metals = cofactors = []

            if args.keep_waters:
                logger.info(
                    "[*] Detecting water molecules near binding site...")
                waters = detect_water_molecules(source, (cx, cy, cz),
                                                distance_threshold=args.water_distance)
                if waters:
                    logger.info(
                        f"[✔] Found {len(waters)} water residue(s) within "
                        f"{args.water_distance}A of the binding site; keeping them in receptor")
                    keep_residues += [
                        (w["chain"], w["resnum"], w["resname"]) for w in waters]
                else:
                    logger.info("[*] No waters detected near binding site")

            if args.detect_metals:
                logger.info("[*] Detecting metal ions...")
                metals = detect_metal_ions(source)
                if metals:
                    logger.info(
                        f"[✔] Found {len(metals)} metal ion(s); keeping them in receptor:")
                    for metal in metals:
                        logger.info(
                            f"   - {metal['name']} ({metal['element']}) chain {metal['chain']} residue {metal['resnum']}")
                    keep_residues += [
                        (m["chain"], m["resnum"], m["resname"]) for m in metals]
                else:
                    logger.info("[*] No metal ions detected")

            if args.detect_cofactors:
                logger.info("[*] Detecting cofactors...")
                cofactors = detect_cofactors(source)
                if cofactors:
                    logger.info(
                        f"[✔] Found {len(cofactors)} cofactor(s); keeping them in receptor:")
                    for cof in cofactors:
                        logger.info(
                            f"   - {cof['name']} chain {cof['chain']} residue {cof['resnum']}")
                    keep_residues += [
                        (c["chain"], c["resnum"], c["name"]) for c in cofactors]
                else:
                    logger.info("[*] No cofactors detected")

            keep_residues = list(dict.fromkeys(keep_residues))
            if keep_residues:
                protein_prep.rebuild_receptor_keep_hetatm(keep_residues)

        # ===== PHASE 3: Flexible Receptor Setup (Optional) =====
        flexible_residues = []
        flex_file = None
        rigid_receptor = None
        if args.flexibility and args.flexibility > 0:
            logger.info(
                f"[PHASE 3] Flexible Receptor Setup (Level {args.flexibility}/10)")

            if args.flexible_residues:
                flexible_residues = [
                    r.strip() for r in args.flexible_residues.split(',') if r.strip()]
                logger.info(
                    f"[*] Using specified flexible residues: {flexible_residues}")
            elif args.auto_flexible:
                flexible_residues = detect_flexible_residues(
                    protein_prep.receptor_pdbqt, (cx, cy, cz),
                    radius=8.0, max_residues=args.auto_flexible)
                logger.info(
                    f"[*] Auto-selected {len(flexible_residues)} flexible residues near binding site")
            else:
                logger.warning(
                    "[!] Flexibility enabled but no residues specified. Use --flexible-residues or --auto-flexible")

            if flexible_residues:
                logger.info(
                    f"[✔] Flexible residues: {', '.join(flexible_residues)}")
                flex_file = build_flexible_residue_pdbqt(
                    protein_prep.receptor_pdbqt, flexible_residues, args.output)
                if not flex_file:
                    logger.warning(
                        "[!] Flexible receptor could not be built; running rigid docking")
                else:
                    # Vina requires flexible residues to be excluded from the
                    # rigid receptor passed to --receptor.
                    rigid_receptor = _remove_residues_from_pdbqt(
                        protein_prep.receptor_pdbqt, flexible_residues,
                        os.path.join(args.output, "rigid.pdbqt"))

        # ===== STEP 3: Docking =====
        logger.info("[3/6] Docking with Advanced Parameters")
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

        results, checkpoint, metrics_dict = dock_all(
            dock_receptor, ligands, protein_prep.grid_conf,
            args.output, num_processes=args.processes, resume=resume,
            vina_params=vina_params, flex_file=flex_file
        )

        # ===== PHASE 4: Consensus Scoring (Optional) =====
        scorer_agreement = None
        if args.consensus or args.smina_only:
            logger.info("[PHASE 4] Consensus Scoring")
            if SMINA is None:
                if args.smina_only:
                    raise RuntimeError(
                        "--smina-only requested but SMINA binary not found in PATH")
                logger.warning(
                    "[!] SMINA not found - continuing with Vina-only ranking")
            else:
                grid_box = (cx, cy, cz, sx, sy, sz)
                smina_scores = _run_smina_scoring(
                    dock_receptor, ligands, args.output, grid_box, vina_params,
                    num_processes=args.processes)

                if args.smina_only:
                    smina_file = os.path.join(args.output, "smina_ranking.csv")
                    smina_rows = sorted(
                        smina_scores.items(), key=lambda kv: kv[1])
                    with open(smina_file, "w", newline="") as f:
                        writer = csv.writer(f)
                        writer.writerow(["Ligand", "SMINA_Affinity"])
                        for name, score in smina_rows:
                            writer.writerow([name, _score_csv(score)])
                    logger.info(
                        f"[✔] SMINA-only ranking saved: {smina_file} "
                        f"({len(smina_rows)} ligands)")
                    # Rewire the pipeline's primary results so ranking.csv,
                    # top hits and the HTML report reflect SMINA scores.
                    if smina_rows:
                        results = [(name, score) for name, score in smina_rows]
                    else:
                        logger.warning(
                            "[!] SMINA scored 0 ligands — falling back to "
                            "Vina ranking for downstream outputs")
                else:
                    consensus_rows = []
                    for name, vina_score in results:
                        smina_score = smina_scores.get(name)
                        if vina_score is None or smina_score is None:
                            continue
                        row = consensus_rank(vina_score, smina_score)
                        consensus_rows.append((name, row))
                    consensus_file = os.path.join(
                        args.output, "consensus_ranking.csv")
                    with open(consensus_file, "w", newline="") as f:
                        writer = csv.writer(f)
                        writer.writerow(
                            ["Ligand", "Vina_Affinity", "SMINA_Affinity",
                             "Consensus", "Agreement (heuristic)"])
                        for name, row in sorted(
                                consensus_rows,
                                key=lambda kv: kv[1]["consensus"]):
                            writer.writerow([
                                name, _score_csv(row["vina"]),
                                _score_csv(row["smina"]),
                                _score_csv(row["consensus"]),
                                f"{row['agreement']:.1f}%"])
                    # Global, citable scorer concordance (Spearman rank corr).
                    rho, n_agree = scorer_agreement_spearman(
                        [row["vina"] for _, row in consensus_rows],
                        [row["smina"] for _, row in consensus_rows])
                    scorer_agreement = {"rho": rho, "n": n_agree}
                    try:
                        with open(consensus_file, "a", newline="") as f:
                            f.write("\n")
                            if rho is not None:
                                f.write(
                                    f"# Global scorer agreement: Spearman rho "
                                    f"= {rho:.3f} (n={n_agree}); Vina vs SMINA "
                                    f"rank correlation (citable metric).\n")
                            else:
                                f.write(
                                    f"# Global scorer agreement: undefined "
                                    f"(n={n_agree}).\n")
                    except Exception as e:
                        logger.warning(f"Could not append agreement summary: {e}")
                    logger.info(
                        f"[✔] Consensus ranking saved: {consensus_file} "
                        f"({len(consensus_rows)} ligands)")

        # ===== STEP 4: Results Analysis =====
        logger.info("[4/6] Results Analysis + Metrics")
        if resume and checkpoint.completed:
            logger.info(
                f"[✔] Checkpoint updated: {len(checkpoint.completed)} successful dockings")

        # ===== PHASE 1: Advanced Result Analysis =====
        html_report_file = os.path.join(args.output, "results_report.html")
        if args.html_report or args.cluster_poses:
            logger.info("[PHASE 1] Advanced Result Analysis")

            # Generate clustering analysis
            if args.cluster_poses:
                logger.info("[*] Performing pose clustering...")
                poses_dir = os.path.join(args.output, "docked")
                clusters = cluster_poses(
                    poses_dir, rmsd_threshold=args.rmsd_threshold)

                # Save clustering data
                clusters_file = os.path.join(args.output, "pose_clusters.json")
                try:
                    with open(clusters_file, 'w') as f:
                        # Convert cluster data to serializable format
                        clusters_data = {}
                        for ligand, ligand_clusters in clusters.items():
                            clusters_data[ligand] = [
                                {'id': c['id'], 'pose_count': len(c['poses'])}
                                for c in ligand_clusters
                            ]
                        json.dump(clusters_data, f, indent=2)
                    logger.info(f"[✔] Clustering saved to: {clusters_file}")
                except Exception as e:
                    logger.warning(f"Could not save clustering data: {e}")

        # ===== STEP 5: Save Results =====
        logger.info("[5/6] Saving Results")
        analyzer = ResultsAnalyzer(args.output, top_n=args.top_n)
        analyzer.save_ranking(results, mode="extended",
                              metrics_dict=metrics_dict)
        analyzer.save_metrics_report(results, metrics_dict)

        # Generate the HTML report AFTER ranking.csv exists.
        if args.html_report:
            logger.info("[*] Generating professional HTML report...")
            ranking_csv = os.path.join(args.output, 'ranking.csv')
            report_meta = {
                "Protein (receptor)": os.path.basename(args.receptor),
                "Ligand library": args.library,
                "Ligands screened": len(results),
                "Top N reported": args.top_n,
                "Exhaustiveness": args.exhaustiveness,
                "Seed": args.seed,
            }
            generate_html_report(ranking_csv, os.path.join(
                args.output, "docked"), html_report_file,
                receptor_pdbqt=protein_prep.receptor_pdbqt,
                meta=report_meta,
                grid_file=protein_prep.grid_conf,
                scorer_agreement=scorer_agreement)

        logger.info("=" * 70)
        logger.info("PIPELINE COMPLETED SUCCESSFULLY (v2.0)")
        logger.info("=" * 70)
        logger.info(f"Results saved to: {args.output}")
        logger.info(f"\n📊 STANDARD OUTPUT:")
        logger.info(f"  ✓ Ranking: {analyzer.csv_file}")
        logger.info(f"  ✓ Full results: {analyzer.results_file}")
        logger.info(f"  ✓ Top hits: {analyzer.top_hits_file}")
        logger.info(f"  ✓ Metrics: {analyzer.metrics_file}")
        if os.path.exists(protein_prep.pocket_summary_file):
            logger.info(
                f"  ✓ Pocket summary: {protein_prep.pocket_summary_file}")
        if os.path.exists(lib_manager.metadata_file):
            logger.info(f"  ✓ Ligand metadata: {lib_manager.metadata_file}")
        if os.path.exists(protein_prep.grid_box_script):
            logger.info(
                f"  ✓ Grid visualization: {protein_prep.grid_box_script}")
        logger.info(
            f"  ✓ Docked structures: {os.path.join(args.output, 'docked')}")

        # Show new Phase 1 outputs
        if args.html_report:
            logger.info(f"\n📈 PHASE 1 (Advanced Analysis):")
            logger.info(f"  ✓ HTML Report: {html_report_file}")
        if args.cluster_poses:
            logger.info(
                f"  ✓ Pose Clustering: {os.path.join(args.output, 'pose_clusters.json')}")

        # Show Phase 2 outputs
        if args.keep_waters or args.detect_metals or args.detect_cofactors:
            logger.info(f"\n🔧 PHASE 2 (Smart Preprocessing):")
            if args.detect_metals or args.detect_cofactors:
                logger.info(f"  ✓ See above for detected metals/cofactors")

        # Show Phase 3 outputs
        if args.flexibility and args.flexibility > 0:
            logger.info(f"\n🎯 PHASE 3 (Flexible Docking):")
            logger.info(
                f"  ✓ Flexible residues: {', '.join(flexible_residues) if flexible_residues else 'None'}")

        logger.info(f"\n📊 Total compounds docked: {len(results)}")
        logger.info("=" * 70)

    except KeyboardInterrupt:
        logger.warning(
            "\n[!] Pipeline interrupted. Progress saved to checkpoint.")
        logger.info("Resume next time with same command.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
# vim: ts=4:et
