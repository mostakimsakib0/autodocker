# AutoDocker 2.0.0: reference-ligand box and seed replicates

Two features answer limitations of the 1.x screen that the paper names:

* **`--box-ligand FILE` / `--box-from-het RESNAME[:CHAIN[:RESNUM]]`**: the
  docking box is built around a ligand known to bind the site. fpocket's
  top-ranked pocket is not always the ligand site (case study), and smina has
  `--autobox_ligand` for the same reason.
* **`--replicates N`**: every ligand is docked N times with seeds
  `--seed + k` (k = 0..N-1), so the screen reports how much a score depends
  on the seed instead of one stochastic sample ("single seed" limitation).

The core logic is in two new modules with their own tests:
`autodocker/vspipeline/ligand_box.py` (`tests/test_ligand_box.py`) and
`autodocker/vspipeline/replicates.py` (`tests/test_replicates.py`). The
command-line wiring was written as the patch at the end of this page
(against `e02da2d`, before the runner refactor) and has been ported by hand
to the refactored code: the options and their validation are in `runner.py`
(`_add_receptor_options`, `_add_docking_options`, `_check_values`,
`_check_box_options`, `_check_engine_options`), the rest in the stages of
`vspipeline/pipeline.py` (`reference_box`/`prepare_receptor`,
`docking_box`, `dock`/`replicate_outputs`, `cluster_poses`,
`write_results`, `log_summary`, `RunState.manifest_entries`). The patch is
kept as the record of the design; the runner-level tests are in
`tests/test_runner.py` (usage errors) and `tests/test_runner_e2e.py`
(`--box-from-het MK1`, `--box-ligand`, `--replicates 2` and its resume).
One difference: an unreadable, empty or 2D `--box-ligand` file is also a
usage error (exit 2), checked with the other options.

## Reference-ligand box

### Box rule

* **Center**: midpoint of the bounding box of the reference ligand's heavy
  atoms, the center rule of the fpocket-pocket and whole-receptor boxes
  (`protein.bounding_box`).
* **Size, per axis**: `max(extent + 2 x padding, 22.5)` Angstrom with
  padding = `--padding` (`defaults.PADDING`, 6 A); by default
  `max(22.5, extent + 12)`.

Why this rule and not smina's `--autobox_add 4` (extent + 8):

* The library ligands are not the reference: a screen needs room for
  ligands larger than, or oriented differently from, the reference. 6 A per
  side is the padding AutoDocker already uses around fpocket pockets, so
  `--padding` keeps one meaning (the space added around whatever defines the
  site) and one default.
* The minimum edge of 22.5 A is the rule of the re-docking and enrichment
  benchmarks (`scripts/_benchlib.ligand_box`); a fragment or a single ion as
  reference still gives a box in which a drug-like ligand (up to ~20 A long)
  can turn. Proposed as `defaults.LIGAND_BOX_MIN_SIZE`.
* For indinavir in 1HSG (MK1, extent 14.9 x 10.2 x 8.6 A) the box is
  26.9 x 22.5 x 22.5 A, inside the 15-40 A range the runner checks.

The benchmarks use a cubic box centered on the heavy-atom centroid (largest
extent on all axes). The pipeline keeps the per-axis bounding-box rule of
its other boxes; the benchmarks keep passing their own box with
`--grid-center/--grid-size`, so their results do not change.

### Reference sources

* `--box-ligand FILE`: SDF/MOL (V2000 and V3000), MOL2, PDB or PDBQT; the
  first molecule/model is used (the log and manifest record how many the
  file holds). Hydrogens are ignored (element H/D; PDBQT types H, HD, HS). A
  2D file (every z = 0) is refused. Pure Python: the box depends only on the
  coordinates in the file, no Open Babel run.
* `--box-from-het RESNAME[:CHAIN[:RESNUM]]`: a HETATM group of the input
  receptor PDB, e.g. `MK1`, `MK1:B`, `MK1:B:902`, `MK1:B:902A` (insertion
  code), `MK1::902` (any chain). Residue names are matched
  case-insensitively, chain IDs exactly (as `--chain`); a blank chain is A.
  Records are read like receptor cleaning reads them
  (`protein.read_pdb_records`: first model, the highest-occupancy altloc of
  each residue). Water is refused. Not found, ambiguous, or empty: the error
  lists the groups present (`MK1:B:902 (45 heavy atoms); and 127 water
  residue(s)`) and the exact spec to use.
* The group is **removed from the receptor**, also with `--keep-hetero`,
  `--detect-cofactors` or `--keep-waters`: the runner writes
  `receptor_no_reference.pdb` (the input without the group's HETATM/ANISOU
  records, in every model and altloc; every other byte kept) and prepares
  the receptor from it, so no later step can bring it back. The group alone
  is written to `reference_ligand.pdb` (input frame) for viewing.

### Frames

The reference coordinates are in the input receptor PDB frame. The docking
frame (receptor.pdbqt, grid.conf, poses) is the input frame minus
`ProteinPreparation.frame_offset()`, the mean of the cleaned receptor's heavy
atoms: exactly the conversion of `--grid-center X --grid-center-frame input`
(`to_docking_frame`, inverse `to_input_frame`). The offset is taken after
receptor preparation, so with `--box-from-het` it is the offset of the
receptor without the reference.

### Option rules

| Combination | Rule |
|---|---|
| `--box-ligand` + `--box-from-het` | error (one reference) |
| + `--grid-center` (and so `--grid-center-frame`) | error (the reference gives the center) |
| + `--pockets` | error (no pocket selects the box) |
| + `--grid-size` | accepted: the center from the reference, the size given (e.g. a fixed 22.5 A cube); the manifest records both sizes |
| + `--padding` | the padding of the rule |
| + `--no-fpocket` | accepted: skips the fpocket comparison below |
| + `--keep-hetero` with `--box-ligand` | accepted; a warning names a receptor HETATM group that occupies the reference's place (>= 50 % of its heavy atoms within 1 A) and suggests `--box-from-het` |
| + `--chain` | the box must contain receptor atoms of the docked chain(s), else the run stops; a warning when no docked heavy atom is within 5 A of the reference (other frame, or its site in another chain; 1HSG's MK1 is chain B and its site spans A and B: use `--chain all`) |

`--box-from-het` is resolved during option validation (a cheap read of the
receptor), so a wrong group is a usage error (exit 2) before any ligand is
prepared.

### fpocket transparency

Unless `--no-fpocket` is given, fpocket still runs, only to report where the
reference ligand is: `pocket_summary.csv` as usual, and in the manifest the
top-ranked pocket, the pocket holding most of the reference (fraction of its
heavy atoms inside the pocket's alpha spheres, `pocketN_vert.pqr`), and
whether fpocket's default choice would have been that pocket. A warning says
so when it would not. The box never depends on fpocket; a failing fpocket is
a warning.

### Log lines

From `--box-from-het MK1 --chain all` on the 1HSG site fixture:

```
[*] --box-from-het MK1:B:902: removed its 45 records from the receptor (receptor_no_reference.pdb); it is not part of the receptor, also with --keep-hetero
[*] Docking box: reference ligand --box-from-het MK1:B:902 (HETATM group of the receptor), 45 heavy atoms
[*] Reference ligand extent: 14.85, 10.21, 8.63 A; box = max(extent + 2 x 6, 22.5) A per axis
[*] Grid center (input-PDB frame): 11.63, 22.47, 5.86
[*] Grid center (docking frame): -2.15, -0.07, 0.69
[*] Grid size: 26.85, 22.50, 22.50
[*] fpocket's top pocket (pocket1) holds 100% of the reference ligand
```

When the reference is elsewhere: `[!] The reference ligand is in pocketN
(fpocket rank R, P% of its heavy atoms), not in fpocket's top pocketM: an
fpocket-selected box would have missed this site`.

### run_manifest.json `box`

```json
"box": {
  "source": "het_group",                 // or "ligand_file"
  "reference": "MK1:B:902",              // or the file path
  "reference_sha256": null,              // sha256 of the --box-ligand file
  "het_group": {"resname": "MK1", "chain": "B", "resid": "902", "spec": "MK1:B:902"},
  "reference_heavy_atoms": 45, "reference_records_in_file": 1,
  "reference_extent": [14.854, 10.214, 8.635],
  "rule": "...", "padding": 6.0, "min_edge": 22.5,
  "center_input_frame": [11.634, 22.47, 5.859],
  "center_docking_frame": [...], "frame_offset": [...],
  "size_computed": [26.854, 22.5, 22.5], "size": [...], "size_source": "reference ligand",
  "reference_atoms_in_box": 45, "closest_receptor_heavy_atom": 2.6,
  "fpocket": {"pockets_found": ..., "top_ranked_pocket": {...},
              "best_overlapping_pocket": {...}, "top_pocket_contains_reference": false,
              "overlapping_pockets": [...], "overlap_definition": "..."}
}
```

## Seed replicates

### What is reported

* Replicate k of every ligand is docked with seed `--seed + k`; replicate 0
  is the run's seed. `seed + N - 1` must stay <= 2^31 - 1; N <= 100 (a typo
  guard). Vina/QuickVina only (`--smina-only` refuses it); `--consensus`
  rescores the reported pose.
* **Binding_Affinity = the best (lowest) in-window score of the
  replicates** (ties: the lowest k), and every per-ligand metric (pose file,
  log, SimScore, modes, engine, status) is that replicate's: the reported
  score and the reported pose belong together. This is what more search
  would give (N seeds at exhaustiveness E sample like one longer search);
  it is optimistic by construction (the minimum of N draws), so the mean is
  reported next to it.
* ranking.csv, only when N > 1, appends `n_replicates` (replicates with an
  in-window score), `score_mean`, `score_sd` (sample SD, n - 1; blank below
  two scores), `score_min`, `score_max` after the existing columns (mean and
  SD rounded to 3 decimals, min/max as Vina wrote them). Statistics cover
  in-window replicates only; an OUT_OF_WINDOW or FAILED replicate is counted
  in `replicate_scores.csv` and the summary, not in the statistics.
* A ligand is OUT_OF_WINDOW when no replicate is in-window and one is out
  of the window, FAILED when all replicates failed (the error lists each).
  A ligand with some failed replicates keeps its score from the others; the
  failed replicates are not checkpointed, so a resumed run docks them again.
* Top_hits.txt adds `Replicates: 3/3 scored; mean -7.25, SD 0.35, range
  -7.50 to -7.00 (reported: replicate 1, seed 43)`.
* `replicate_scores.csv`: one row per (ligand, replicate): seed, status,
  score, raw score, whether it is the reported one, pose and log files.
* `replicate_summary.json` and run_manifest.json `replicates`: seeds,
  ligands, median and max per-ligand SD, median range, **Spearman rho of
  replicate 0 vs replicate 1 scores** over the ligands scored in both
  (`consensus.scorer_agreement_spearman`, tie-aware; null when undefined),
  how often each replicate gave the reported score, and the ligands with
  failed replicates.

### Files, checkpoint, fingerprint

* Replicate k of `<name>` writes `docked/<name>_rep<k>_out.pdbqt` and
  `docked/<name>_rep<k>_vina.log`, and is the checkpoint entry
  `<name>_rep<k>`. `dock_ligand` derives every output name from the ligand
  file name, so the replicate docks a byte-identical copy
  `docked/replicate_inputs/<name>_rep<k>.pdbqt` (rewritten only when the
  source changes; same SHA-256 as the source, so the checkpoint's ligand
  hash still detects a changed ligand). The docking worker, its log, the
  fallback-engine bookkeeping and the checkpoint need no replicate logic;
  `-p` workers run (ligand, replicate) jobs in parallel. The log's
  `--ligand` path names the copy.
* The stem `<name>_rep<k>` is injective (`parse_replicate_stem` splits at
  the last `_rep<digits>`), so ligand `a_rep1` cannot collide with replicate
  1 of `a`.
* The run fingerprint includes `"replicates": N` only when N > 1. Changing
  N re-docks everything (the N > 1 files are not the N = 1 files, and the
  reported score depends on N); a resumed run with the same N docks only the
  missing (ligand, replicate) pairs.

### N = 1 is byte-identical

With N = 1 (the default) `dock_all` takes none of the replicate branches:
no staged copies, the same job list, log lines, checkpoint keys and
fingerprint (`fingerprint_params` returns the parameters unchanged), the
same ranking.csv columns, Top_hits.txt, Results_full.txt and HTML report
(clustering scans `docked/` as before). Checked by running the unpatched
and the patched runner on the same input (1HSG site, 2 ligands, explicit
box): `diff -r` of the two output directories shows only the output
directory's name inside the Vina logs, the checkpoint timestamps, the lock
PID and run_manifest.json (timestamps; its `args` gains `box_ligand`,
`box_from_het` and `replicates`). The checkpoint fingerprint is identical.

## Wiring patch

Historical: against `e02da2d` (runner.py before the refactor); ported as
described at the top of this page, do not apply it. It touches `runner.py`,
`docking.py`, `results.py`, `report.py`, `outputs.py`, `defaults.py`, and
the two new modules (their proposed constants become `defaults`/`outputs`
names). Extract and apply with

```
sed -n '/^```diff$/,/^```$/p' docs/features_2_0.md | sed '1d;$d' > /tmp/wiring.diff
git apply --3way /tmp/wiring.diff
```

What each part does:

* `defaults.py`: `LIGAND_BOX_MIN_SIZE = 22.5`, `REPLICATES = 1`.
* `outputs.py`: `RECEPTOR_NO_REFERENCE_PDB`, `REFERENCE_LIGAND_PDB`,
  `REPLICATE_TAG`, `REPLICATE_INPUTS_DIR`, `REPLICATE_SCORES_CSV`,
  `REPLICATE_SUMMARY_JSON` (the four files are added to `DERIVED_OUTPUTS`),
  and `RANKING_REPLICATE_COLUMNS`.
* `docking.py`: `dock_all(..., replicates=1)` and
  `_run_fingerprint(..., replicates=1)`; every new line runs only when
  N > 1 (units from `replicates.expand_ligands`, per-unit seed, aggregation
  with `replicates.aggregate` before returning).
* `results.py`: `ResultsAnalyzer(..., replicates=1)`: the replicate columns
  and the Top_hits.txt line only when N > 1; a FAILED ligand of an N > 1 run
  points to `docked/<name>_rep<k>_vina.log`.
* `report.py`: `cluster_poses(..., pose_files=None)`; the HTML report
  passes the reported pose files when ranking.csv has `n_replicates`.
* `runner.py`: options `--box-ligand`, `--box-from-het`, `--replicates`;
  validation (`ligand_box.check_options`, the HET group looked up in the
  receptor, `replicates.check_replicates`, `--smina-only` conflict); the
  receptor prepared from `receptor_no_reference.pdb` with `--box-from-het`;
  `_apply_reference_box` (log lines, receptor-in-box error, contact and
  `--keep-hetero` warnings, fpocket comparison); `--grid-size` override
  recorded; `dock_all(replicates=N)`; `_replicate_outputs` (CSV, JSON, log
  lines); clustering of the reported poses; HTML meta lines; the list of
  written files; manifest entries `box` and `replicates`.

Verified with the patch applied: full suite (830 passed, 2 skipped), ruff
clean, the N = 1 byte-identity run above, and CLI runs on the 1HSG fixture:
`--box-from-het MK1 --replicates 2 --keep-hetero --cluster-poses` (no MK1
atom in receptor.pdbqt, 45/45 MK1 atoms in the box, fpocket pocket 1 holds
100 % of MK1, 4 (ligand, replicate) jobs, rerun: 4 cached, pose clusters of
the best replicates, HTML report with the replicate line), `--box-ligand
mk1.sdf --keep-hetero --no-fpocket` (warns that MK1:B:902 stays in the
receptor), and usage errors (exit 2) for an unknown HET group and for
`--box-ligand` with `--grid-center`. Runner-level tests for these paths
(tests/test_runner.py, test_runner_e2e.py) belong to the refactored runner.

```diff
diff --git a/autodocker/runner.py b/autodocker/runner.py
index 0474acd..25317d7 100755
--- a/autodocker/runner.py
+++ b/autodocker/runner.py
@@ -44,9 +44,11 @@ from vspipeline import (
     flexible,
     hydrogens,
     library,
+    ligand_box,
     outputs,
     pocket,
     protein,
+    replicates,
     report,
     results,
     tools,
@@ -136,7 +138,9 @@ def _require_tools(smina_only: bool = False, consensus: bool = False,
 
 def _uses_fpocket(args) -> bool:
     """True when the run detects pockets with fpocket: not with --no-fpocket,
-    and not when --grid-center and --grid-size give the whole box."""
+    and not when --grid-center and --grid-size give the whole box. With
+    --box-ligand/--box-from-het fpocket only reports which pocket holds the
+    reference ligand (--no-fpocket skips that)."""
     return not (args.no_fpocket
                 or (args.grid_center is not None and args.grid_size is not None))
 
@@ -262,7 +266,8 @@ def build_parser() -> argparse.ArgumentParser:
                      "score (druggability breaks ties). Not combinable with --no-fpocket or "
                      "with --grid-center and --grid-size together (no pocket is used then)")
     rec.add_argument("--padding", type=float, default=defaults.PADDING,
-                     help="Box padding around the selected pocket(s), Angstrom")
+                     help="Box padding around the selected pocket(s) or reference ligand, "
+                     "Angstrom")
     rec.add_argument("--no-fpocket", action="store_true",
                      help="Skip fpocket; dock in a box around the whole receptor (heavy-atom "
                      f"bounding box + padding, each edge at least "
@@ -278,6 +283,18 @@ def build_parser() -> argparse.ArgumentParser:
     rec.add_argument("--grid-size", default=None,
                      help="Explicit box size 'SX,SY,SZ' (Angstrom), e.g. 22.5,22.5,22.5; "
                      "overrides auto-detection")
+    rec.add_argument("--box-ligand", default=None, metavar="FILE",
+                     help="Box around a reference ligand (SDF/MOL/MOL2/PDB/PDBQT, first "
+                     "molecule; coordinates in the receptor PDB's frame): centered on its "
+                     "heavy atoms' bounding box, each edge max(extent + 2 x --padding, "
+                     f"{defaults.LIGAND_BOX_MIN_SIZE:g}) Angstrom. fpocket then only reports "
+                     "which pocket holds the ligand. Not combinable with --grid-center, "
+                     "--pockets or --box-from-het; --grid-size replaces the size")
+    rec.add_argument("--box-from-het", default=None, metavar="RESNAME[:CHAIN[:RESNUM]]",
+                     help="As --box-ligand, with the reference ligand taken from a HETATM "
+                     "group of the receptor PDB (e.g. MK1, MK1:B or MK1:B:902; residue "
+                     "number with insertion code, e.g. 902A); the group is removed from "
+                     "the receptor, also with --keep-hetero")
     rec.add_argument("--keep-waters", action="store_true",
                      help="Keep water molecules near the binding site in the receptor")
     rec.add_argument("--water-distance", type=float, default=defaults.WATER_DISTANCE,
@@ -318,6 +335,10 @@ def build_parser() -> argparse.ArgumentParser:
     dock.add_argument("--seed", type=int, default=None,
                       help=f"Random seed, {SEED_RANGE[0]} to {SEED_RANGE[1]} (Vina treats 0 "
                       f"as 'random'; default: $VS_SEED, else {DEFAULT_SEED})")
+    dock.add_argument("--replicates", type=int, default=defaults.REPLICATES,
+                      help="Dock every ligand N times with seeds --seed, --seed + 1, ...; "
+                      "the reported score is the best replicate's (ranking.csv adds "
+                      "n_replicates and the score mean/SD/min/max when N > 1). Vina only")
     dock.add_argument("--threads", type=int, default=None,
                       help="CPU threads per Vina process, passed as --cpu (default: all "
                       "cores, or cores divided by -p when -p > 1)")
@@ -520,6 +541,10 @@ def _validate_args(args, explicit) -> None:
     for option, value in positive_ints:
         if value < 1:
             raise ConfigError(f"{option} must be >= 1, got {value}")
+    try:
+        replicates.check_replicates(args.replicates, args.seed, SEED_RANGE[1])
+    except ValueError as e:
+        raise ConfigError(str(e)) from None
     for option, value in (("--energy-range", args.energy_range),
                           ("--rmsd-threshold", args.rmsd_threshold),
                           ("--water-distance", args.water_distance)):
@@ -548,6 +573,15 @@ def _validate_args(args, explicit) -> None:
             if option == "--grid-size" and min(triplet) <= 0:
                 raise ConfigError(f"--grid-size must be positive, got {value}")
 
+    try:
+        ligand_box.check_options(args.box_ligand, args.box_from_het, args.pockets,
+                                 args.grid_center)
+        if args.box_from_het is not None:
+            # Cheap read of the receptor: the group must exist and be unique.
+            ligand_box.select_het_group(args.receptor, args.box_from_het)
+    except ligand_box.LigandBoxError as e:
+        raise ConfigError(str(e)) from None
+
     # Options that would be accepted but have no effect.
     if "grid_center_frame" in explicit and args.grid_center is None:
         raise ConfigError("--grid-center-frame only applies to --grid-center")
@@ -607,6 +641,8 @@ def _validate_args(args, explicit) -> None:
             conflicts.append("--binding-modes (SMINA mode reports 1 mode)")
         if args.energy_range != DEFAULT_ENERGY_RANGE:
             conflicts.append("--energy-range")
+        if args.replicates > 1:
+            conflicts.append("--replicates (replicates are Vina dockings)")
         if conflicts:
             raise ConfigError("--smina-only cannot be combined with " + ", ".join(conflicts))
 
@@ -912,6 +948,107 @@ def _threads_per_process(args) -> Optional[int]:
     return None
 
 
+def _reference_box(args):
+    """(ReferenceBox or None, receptor PDB to prepare). --box-from-het: the
+    group is written alone to reference_ligand.pdb and the receptor is
+    prepared from a copy without it (receptor_no_reference.pdb)."""
+    if args.box_ligand:
+        return ligand_box.box_from_ligand_file(
+            args.box_ligand, args.padding, defaults.LIGAND_BOX_MIN_SIZE), args.receptor
+    if not args.box_from_het:
+        return None, args.receptor
+    reference = ligand_box.box_from_het(args.receptor, args.box_from_het, args.padding,
+                                        defaults.LIGAND_BOX_MIN_SIZE)
+    stripped = os.path.join(args.output, outputs.RECEPTOR_NO_REFERENCE_PDB)
+    removed = ligand_box.write_receptor_without_group(args.receptor, stripped,
+                                                      reference.group)
+    ligand_box.write_reference_pdb(os.path.join(args.output, outputs.REFERENCE_LIGAND_PDB),
+                                   args.receptor, reference.group)
+    logger.info(f"[*] --box-from-het {reference.reference}: removed its {removed} records "
+                f"from the receptor ({outputs.RECEPTOR_NO_REFERENCE_PDB}); it is not part "
+                "of the receptor, also with --keep-hetero")
+    return reference, stripped
+
+
+def _reference_pockets(protein_prep, reference) -> Dict:
+    """Run fpocket on the prepared receptor only to report which pocket
+    holds the reference ligand (pocket_summary.csv, manifest). A failing
+    fpocket does not stop the run: the box does not depend on it."""
+    try:
+        target = protein_prep.pdb_clean
+        pockets = protein.rank_pockets(protein_prep._collect_pockets(
+            protein_prep._run_fpocket(target), target))
+        protein_prep._write_pocket_summary(pockets)
+    except RuntimeError as e:
+        logger.warning(f"[!] fpocket (run to compare with the reference ligand) failed: {e}")
+        return {"error": str(e)}
+    rows = ligand_box.pocket_overlaps(list(reference.coords), pockets, reference.center)
+    summary = ligand_box.fpocket_summary(rows)
+    best, top = summary["best_overlapping_pocket"], summary["top_ranked_pocket"]
+    if best is None:
+        logger.warning(f"[!] No fpocket pocket overlaps the reference ligand (fpocket's "
+                       f"top pocket: pocket{top['pocket']})")
+    elif summary["top_pocket_contains_reference"]:
+        logger.info(f"[*] fpocket's top pocket (pocket{top['pocket']}) holds "
+                    f"{best['overlap_fraction']:.0%} of the reference ligand")
+    else:
+        logger.warning(
+            f"[!] The reference ligand is in pocket{best['pocket']} (fpocket rank "
+            f"{best['rank']}, {best['overlap_fraction']:.0%} of its heavy atoms), not in "
+            f"fpocket's top pocket{top['pocket']}: an fpocket-selected box would have "
+            "missed this site")
+    return summary
+
+
+def _apply_reference_box(args, reference, protein_prep, state: Dict):
+    """Box (docking frame) around the reference ligand, with the checks and
+    log lines; fills state["reference_*"] for the manifest entry."""
+    offset = protein_prep.frame_offset()
+    for line in reference.log_lines(offset):
+        logger.info(line)
+    xs, ys, zs = protein._heavy_atom_coords(protein_prep.pdb_clean)
+    receptor_atoms = list(zip(xs, ys, zs))
+    if not ligand_box.points_in_box(receptor_atoms, reference.center, reference.size):
+        raise ValueError(
+            f"The box around the reference ligand ({reference.describe()}) contains no "
+            "receptor atom of the docked chain(s): is the ligand in the receptor PDB's "
+            "frame, and is its site in the chains selected with --chain?")
+    contact = ligand_box.min_distance(list(reference.coords), receptor_atoms)
+    state["reference_contact"] = contact
+    if contact > ligand_box.CONTACT_WARN_DISTANCE:
+        logger.warning(
+            f"[!] The reference ligand is {contact:.1f} A from the nearest receptor heavy "
+            f"atom (> {ligand_box.CONTACT_WARN_DISTANCE:g} A): check that it is in the "
+            "receptor PDB's frame and that --chain selects the chains of its site")
+    if args.box_ligand and args.keep_hetero:
+        same = ligand_box.matching_het_groups(list(reference.coords),
+                                              ligand_box.het_groups(args.receptor))
+        if same:
+            logger.warning(
+                f"[!] The reference ligand occupies the receptor's HETATM group "
+                f"{same[0].spec}, which --keep-hetero keeps in the receptor; use "
+                f"--box-from-het {same[0].spec} to remove it")
+    if _uses_fpocket(args):
+        state["reference_fpocket"] = _reference_pockets(protein_prep, reference)
+    cx, cy, cz = reference.docking_center(offset)
+    sx, sy, sz = reference.size
+    return cx, cy, cz, sx, sy, sz
+
+
+def _replicate_outputs(args, scores, metrics_dict) -> Dict:
+    """replicate_scores.csv and replicate_summary.json of a run with
+    --replicates N > 1; returns the summary (run_manifest.json)."""
+    summary = replicates.summary(metrics_dict, args.replicates, args.seed)
+    replicates.write_replicate_scores(
+        os.path.join(args.output, outputs.REPLICATE_SCORES_CSV), args.output,
+        sorted(name for name, _ in scores), metrics_dict)
+    replicates.write_summary(os.path.join(args.output, outputs.REPLICATE_SUMMARY_JSON),
+                             summary)
+    for line in replicates.log_lines(summary):
+        (logger.warning if line.startswith("[!]") else logger.info)(line)
+    return summary
+
+
 def _run_pipeline(args, state: Dict) -> int:
     """Stages 1-5. Returns the exit code; raises on pipeline errors."""
     logger.info("=" * 70)
@@ -945,14 +1082,17 @@ def _run_pipeline(args, state: Dict) -> int:
 
     # ===== [2/5] Receptor and box =====
     logger.info("[2/5] Receptor preparation and docking box")
-    protein_prep = protein.ProteinPreparation(args.receptor, args.output, seed=args.seed)
+    reference, receptor_input = _reference_box(args)
+    protein_prep = protein.ProteinPreparation(receptor_input, args.output, seed=args.seed)
     state["protein"] = protein_prep
     protein_prep.validate()
 
     chain = args.chain or protein_prep.select_chain()
     protein_prep.prepare_receptor(chain, keep_hetero=args.keep_hetero, ph=args.ph)
 
-    if args.no_fpocket:
+    if reference is not None:
+        cx, cy, cz, sx, sy, sz = _apply_reference_box(args, reference, protein_prep, state)
+    elif args.no_fpocket:
         cx, cy, cz, sx, sy, sz = protein_prep.centroid_grid(padding=args.padding)
         logger.info("[*] Docking box: whole-receptor box (--no-fpocket)")
     elif not _uses_fpocket(args):
@@ -963,7 +1103,7 @@ def _run_pipeline(args, state: Dict) -> int:
         cx, cy, cz, sx, sy, sz = protein_prep.detect_pocket(
             pocket_spec=args.pockets, padding=args.padding)
 
-    if _uses_fpocket(args) or args.no_fpocket:
+    if reference is None and (_uses_fpocket(args) or args.no_fpocket):
         logger.info(f"[*] Grid center: {cx:.2f}, {cy:.2f}, {cz:.2f}")
         logger.info(f"[*] Grid size: {sx:.2f}, {sy:.2f}, {sz:.2f}")
 
@@ -980,6 +1120,12 @@ def _run_pipeline(args, state: Dict) -> int:
     if args.grid_size:
         sx, sy, sz = docking._parse_grid_triplet(args.grid_size, "grid-size")
         logger.info(f"[*] Grid size overridden: {sx:.2f}, {sy:.2f}, {sz:.2f}")
+    if reference is not None:
+        state["box"] = reference.manifest(
+            protein_prep.frame_offset(),
+            grid_size_override=(sx, sy, sz) if args.grid_size else None,
+            fpocket=state.pop("reference_fpocket", None),
+            receptor_contact=state.pop("reference_contact", None))
 
     low, high = TYPICAL_GRID_SIZE_RANGE
     for axis, size in (("x", sx), ("y", sy), ("z", sz)):
@@ -1150,8 +1296,11 @@ def _run_pipeline(args, state: Dict) -> int:
         scores, ckpt, metrics_dict = docking.dock_all(
             dock_receptor, ligands, protein_prep.grid_conf,
             args.output, num_processes=args.processes, resume=resume,
-            vina_params=vina_params, flex_file=flex_file, engine=engine)
+            vina_params=vina_params, flex_file=flex_file, engine=engine,
+            replicates=args.replicates)
         state["run_fingerprint"] = getattr(ckpt, "fingerprint", None)
+        if args.replicates > 1:
+            state["replicates"] = _replicate_outputs(args, scores, metrics_dict)
 
     scorer_agreement = None
     if args.consensus:
@@ -1175,7 +1324,9 @@ def _run_pipeline(args, state: Dict) -> int:
         logger.info("[*] Performing pose clustering...")
         clusters = report.cluster_poses(
             os.path.join(args.output, outputs.DOCKED_DIR), rmsd_threshold=args.rmsd_threshold,
-            ligands=[name for name, _ in scores])
+            ligands=[name for name, _ in scores],
+            pose_files=({name: (metrics_dict.get(name) or {}).get("dock_file")
+                         for name, _ in scores} if args.replicates > 1 else None))
         clusters_data = {
             ligand: [{'id': c['id'], 'pose_count': len(c['poses'])} for c in ligand_clusters]
             for ligand, ligand_clusters in clusters.items()}
@@ -1188,7 +1339,8 @@ def _run_pipeline(args, state: Dict) -> int:
 
     # ===== [5/5] Outputs =====
     logger.info("[5/5] Saving results")
-    analyzer = results.ResultsAnalyzer(args.output, top_n=args.top_n)
+    analyzer = results.ResultsAnalyzer(args.output, top_n=args.top_n,
+                                       replicates=args.replicates)
     analyzer.save_ranking(scores, mode="extended", metrics_dict=metrics_dict)
     analyzer.save_metrics_report(scores, metrics_dict)
 
@@ -1204,6 +1356,12 @@ def _run_pipeline(args, state: Dict) -> int:
             "Exhaustiveness": args.exhaustiveness,
             "Seed": args.seed,
         }
+        if args.replicates > 1:
+            report_meta["Replicates"] = (
+                f"{args.replicates} (seeds {args.seed}..{args.seed + args.replicates - 1}; "
+                "best replicate reported)")
+        if "box" in state:
+            report_meta["Docking box"] = f"reference ligand {state['box']['reference']}"
         html_ok = report.generate_html_report(
             os.path.join(args.output, outputs.RANKING_CSV), os.path.join(args.output, outputs.DOCKED_DIR),
             html_report_file,
@@ -1243,6 +1401,14 @@ def _run_pipeline(args, state: Dict) -> int:
                     ("Consensus summary", os.path.join(args.output, outputs.CONSENSUS_SUMMARY_JSON))]
     if args.cluster_poses:
         written.append(("Pose clusters", clusters_file))
+    if args.replicates > 1:
+        written += [("Replicate scores",
+                     os.path.join(args.output, outputs.REPLICATE_SCORES_CSV)),
+                    ("Replicate summary",
+                     os.path.join(args.output, outputs.REPLICATE_SUMMARY_JSON))]
+    if args.box_from_het:
+        written.append(("Reference ligand",
+                        os.path.join(args.output, outputs.REFERENCE_LIGAND_PDB)))
     if args.html_report:
         if html_ok:
             written.append(("HTML report", html_report_file))
@@ -1402,6 +1568,12 @@ def _write_final_manifest(args, argv, tool_probes, started_at, status, exit_code
             docking._engine_identity(path) for path, _ in state["engine"].docking_engines()]
     if "consensus" in state:
         extra["consensus"] = state["consensus"]
+    if "box" in state:
+        # --box-ligand/--box-from-het: reference, rule, center in both frames.
+        extra["box"] = state["box"]
+    if "replicates" in state:
+        # --replicates N > 1: seeds, median SD, Spearman rho of replicates 0/1.
+        extra["replicates"] = state["replicates"]
     if "flexible_residues" in state:
         # The residues docked as flexible side chains ([] = rigid).
         extra["flexible_residues"] = list(state["flexible_residues"])
diff --git a/autodocker/vspipeline/defaults.py b/autodocker/vspipeline/defaults.py
index c7eb91e..9509598 100644
--- a/autodocker/vspipeline/defaults.py
+++ b/autodocker/vspipeline/defaults.py
@@ -20,10 +20,12 @@ COMMAND_TIMEOUT = 900           # seconds per external command
 # in docking.AFFINITY_FLOOR) are reported as OUT_OF_WINDOW.
 MIN_VALID_AFFINITY = -1.0       # kcal/mol
 MAX_FAILURE_FRACTION = 0.5      # exit 3 above this fraction of failed ligands
+REPLICATES = 1                  # dockings per ligand, seeds --seed + 0..N-1
 
 # Receptor and docking box
 PH = 7.4                        # receptor protonation pH
-PADDING = 6.0                   # Angstrom around the selected pocket(s)
+PADDING = 6.0                   # Angstrom around the selected pocket(s) or reference ligand
+LIGAND_BOX_MIN_SIZE = 22.5      # Angstrom; smallest edge of a --box-ligand/--box-from-het box
 WATER_DISTANCE = 4.0            # Angstrom; --keep-waters radius around the box center
 AUTO_FLEX_RADIUS = 8.0          # Angstrom; --auto-flexible candidates around the box center
 AUTO_FLEX_MAX_RESIDUES = 10
diff --git a/autodocker/vspipeline/docking.py b/autodocker/vspipeline/docking.py
index 8c0dfe8..407007b 100644
--- a/autodocker/vspipeline/docking.py
+++ b/autodocker/vspipeline/docking.py
@@ -33,6 +33,7 @@ from enum import Enum
 from typing import Callable, Dict, List, Optional, Sequence, Tuple, Union
 from vspipeline import defaults  # noqa: E402
 from vspipeline import outputs  # noqa: E402
+from vspipeline import replicates as replicates_mod  # noqa: E402
 from vspipeline import tools, utils  # noqa: E402
 from vspipeline.charges import _pdbqt_has_atoms, _pdbqt_has_nonzero_charges  # noqa: E402
 from vspipeline.checkpoint import DockingCheckpoint  # noqa: E402
@@ -924,7 +925,8 @@ def _receptor_sha256(path: Optional[str]) -> str:
 
 
 def _run_fingerprint(receptor: str, grid_file: str, flex_file: Optional[str],
-                     vina_params: Dict, engine: Optional[EngineConfig] = None) -> str:
+                     vina_params: Dict, engine: Optional[EngineConfig] = None,
+                     replicates: int = 1) -> str:
     """Hash of everything that determines a ligand's docking result apart from
     the ligand itself: receptor and flex files (byte-exact), grid file,
     docking parameters and the primary engine binary (path, SHA-256,
@@ -935,9 +937,13 @@ def _run_fingerprint(receptor: str, grid_file: str, flex_file: Optional[str],
     primary engine produced stays valid whether or not a fallback is
     configured. A result the fallback produced is tagged with the fallback's
     SHA-256 (metrics ``fallback_sha256``) and re-docked on its own when the
-    fallback changes or is turned off (see dock_all)."""
+    fallback changes or is turned off (see dock_all).
+
+    ``replicates`` (--replicates N) is part of the parameters only when N > 1,
+    so the fingerprint of a plain run is unchanged."""
     engine = engine or EngineConfig.from_settings()
-    params = {k: v for k, v in vina_params.items() if k != "threads"}
+    params = replicates_mod.fingerprint_params(
+        {k: v for k, v in vina_params.items() if k != "threads"}, replicates)
     payload = json.dumps({
         "receptor": _file_sha256(receptor),
         "grid": _file_sha256(grid_file),
@@ -1175,7 +1181,8 @@ def dock_all(receptor: str, ligands: List[str], grid_file: str,
              resume: bool = True,
              vina_params: Optional[Dict] = None,
              flex_file: Optional[str] = None,
-             engine: Optional[EngineConfig] = None
+             engine: Optional[EngineConfig] = None,
+             replicates: int = 1
              ) -> Tuple[List[Tuple[str, Optional[float]]], DockingCheckpoint, Dict]:
     """Dock all ligands with resume and metrics.
 
@@ -1190,6 +1197,12 @@ def dock_all(receptor: str, ligands: List[str], grid_file: str,
     When no ligand produced any score, an error is logged and every ligand
     is returned as FAILED (the caller writes the outputs and decides the exit
     code).
+
+    ``replicates`` N > 1 docks every ligand N times, replicate k with seed
+    ``vina_params["seed"] + k`` (vspipeline.replicates): every (ligand,
+    replicate) is a job and a checkpoint entry of its own (<name>_rep<k>),
+    and the results and metrics returned are per ligand (best replicate,
+    with the replicate statistics in its metrics).
     """
     engine = engine or EngineConfig.from_settings()
     dock_dir = os.path.join(outdir, outputs.DOCKED_DIR)
@@ -1206,6 +1219,14 @@ def dock_all(receptor: str, ligands: List[str], grid_file: str,
         )
 
     vina_params = {**_default_vina_params(), **(vina_params or {})}
+    library_names = ligand_names
+    if replicates > 1:
+        # One docking unit per (ligand, replicate): a copy of the ligand named
+        # <name>_rep<k>, so its outputs and checkpoint key carry the replicate.
+        units = replicates_mod.expand_ligands(ligands, replicates, dock_dir)
+        unit_replicate = {utils.ligand_name(path): k for path, k in units}
+        ligands = [path for path, _ in units]
+        ligand_names = [utils.ligand_name(lig) for lig in ligands]
 
     # Each Vina process uses every core by default, so N workers would
     # oversubscribe the CPU N-fold. Split the cores between workers unless
@@ -1216,7 +1237,7 @@ def dock_all(receptor: str, ligands: List[str], grid_file: str,
     ligand_hashes = {utils.ligand_name(lig): _file_sha256(lig) for lig in ligands}
     checkpoint = DockingCheckpoint(
         outdir, fingerprint=_run_fingerprint(receptor, grid_file, flex_file,
-                                             vina_params, engine),
+                                             vina_params, engine, replicates),
         resume=resume)
 
     # A result of the fallback engine is only valid with the same fallback
@@ -1245,9 +1266,16 @@ def dock_all(receptor: str, ligands: List[str], grid_file: str,
         checkpoint.completed.pop(utils.ligand_name(lig), None)
 
     min_valid_affinity = float(vina_params['min_valid_affinity'])
-    logger.info(
-        f"[*] Docking {len(ligands_todo)} new ligands "
-        f"({len(ligands) - len(ligands_todo)} cached)")
+    if replicates > 1:
+        logger.info(
+            f"[*] Docking {len(ligands_todo)} new (ligand, replicate) pairs "
+            f"({len(ligands) - len(ligands_todo)} cached): {len(library_names)} ligands x "
+            f"{replicates} replicates, seeds {vina_params['seed']}.."
+            f"{vina_params['seed'] + replicates - 1}")
+    else:
+        logger.info(
+            f"[*] Docking {len(ligands_todo)} new ligands "
+            f"({len(ligands) - len(ligands_todo)} cached)")
     logger.info(f"[*] Vina parameters: exhaustiveness={vina_params['exhaustiveness']}, "
                 f"modes={vina_params['binding_modes']}, "
                 f"range={vina_params['energy_range']} kcal/mol, "
@@ -1296,7 +1324,10 @@ def dock_all(receptor: str, ligands: List[str], grid_file: str,
     if not ligands_todo:
         logger.info("[*] All ligands already docked (checkpoint)")
     else:
-        jobs = [DockJob(receptor, lig, grid_file, dock_dir, vina_params, flex_file, engine)
+        jobs = [DockJob(receptor, lig, grid_file, dock_dir,
+                        vina_params if replicates <= 1 else replicates_mod.replicate_params(
+                            vina_params, unit_replicate[utils.ligand_name(lig)]),
+                        flex_file, engine)
                 for lig in ligands_todo]
         start_time = time.time()
         total = len(jobs)
@@ -1352,4 +1383,8 @@ def dock_all(receptor: str, ligands: List[str], grid_file: str,
             "[!] All docking attempts failed; no ligand produced a score. Inspect "
             f"{os.path.join(dock_dir, '*' + outputs.VINA_LOG_SUFFIX)} for root causes.")
 
+    if replicates > 1:
+        all_results, metrics_dict = replicates_mod.aggregate(
+            library_names, dict(all_results), metrics_dict, replicates,
+            int(vina_params['seed']))
     return all_results, checkpoint, metrics_dict
diff --git a/autodocker/vspipeline/ligand_box.py b/autodocker/vspipeline/ligand_box.py
index d3c0472..f761fc7 100644
--- a/autodocker/vspipeline/ligand_box.py
+++ b/autodocker/vspipeline/ligand_box.py
@@ -41,11 +41,11 @@ from vspipeline.protein import (
 
 Coord = Tuple[float, float, float]
 
-# Proposed defaults.LIGAND_BOX_MIN_SIZE: smallest box edge (Angstrom) around
-# a reference ligand; the benchmarks' rule (scripts/_benchlib.py). A small
-# reference (a fragment, a metal ion) still gives a box in which a drug-like
-# library ligand (up to ~20 A long) can turn.
-MIN_LIGAND_BOX_EDGE = 22.5
+# Smallest box edge (Angstrom) around a reference ligand; the benchmarks'
+# rule (scripts/_benchlib.py). A small reference (a fragment, a metal ion)
+# still gives a box in which a drug-like library ligand (up to ~20 A long)
+# can turn.
+MIN_LIGAND_BOX_EDGE = defaults.LIGAND_BOX_MIN_SIZE
 # Padding on each side: --padding (defaults.PADDING, the fpocket-pocket padding).
 LIGAND_BOX_PADDING = defaults.PADDING
 # Warn when no receptor heavy atom is this close (Angstrom) to the reference
diff --git a/autodocker/vspipeline/outputs.py b/autodocker/vspipeline/outputs.py
index d424e42..5b2de15 100644
--- a/autodocker/vspipeline/outputs.py
+++ b/autodocker/vspipeline/outputs.py
@@ -29,6 +29,10 @@ GRID_CONF = "grid.conf"
 GRID_BOX_PML = "grid_box.pml"
 POCKET_SUMMARY_CSV = "pocket_summary.csv"
 FRAME_OFFSET_TXT = "frame_offset.txt"
+# --box-from-het: the input receptor without the reference HETATM group (the
+# file receptor preparation reads) and that group alone (input frame).
+RECEPTOR_NO_REFERENCE_PDB = "receptor_no_reference.pdb"
+REFERENCE_LIGAND_PDB = "reference_ligand.pdb"
 
 # fpocket's own output directory (next to PROTEIN_CLEAN_PDB, named by fpocket)
 FPOCKET_OUT_DIR = PROTEIN_CLEAN_PDB[:-len(".pdb")] + "_out"
@@ -36,6 +40,14 @@ FPOCKET_OUT_DIR = PROTEIN_CLEAN_PDB[:-len(".pdb")] + "_out"
 # Docking
 DOCKED_DIR = "docked"   # <name>_out.pdbqt poses and <name>_vina.log per ligand
 
+# --replicates N > 1: replicate k of a ligand is <name>_rep<k> in DOCKED_DIR
+# (poses, log, checkpoint key), docked from a copy of the ligand in
+# DOCKED_DIR/REPLICATE_INPUTS_DIR; per-replicate scores and the run summary.
+REPLICATE_TAG = "_rep"
+REPLICATE_INPUTS_DIR = "replicate_inputs"
+REPLICATE_SCORES_CSV = "replicate_scores.csv"
+REPLICATE_SUMMARY_JSON = "replicate_summary.json"
+
 # Outputs derived from a run's options and results. Every run removes the ones
 # left by an earlier run in the same output directory before it starts, so a
 # file that this run does not write (--no-html-report, no --consensus, no
@@ -46,7 +58,8 @@ DERIVED_OUTPUTS = (
     RANKING_CSV, TOP_HITS_TXT, METRICS_TXT, RESULTS_FULL_TXT, HTML_REPORT,
     POSE_CLUSTERS_JSON, CONSENSUS_RANKING_CSV, CONSENSUS_SUMMARY_JSON,
     POCKET_SUMMARY_CSV, FPOCKET_OUT_DIR, GRID_CONF, GRID_BOX_PML, FRAME_OFFSET_TXT,
-    RIGID_RECEPTOR_PDBQT, FLEX_PDBQT,
+    RIGID_RECEPTOR_PDBQT, FLEX_PDBQT, RECEPTOR_NO_REFERENCE_PDB, REFERENCE_LIGAND_PDB,
+    REPLICATE_SCORES_CSV, REPLICATE_SUMMARY_JSON,
 )
 
 # Per-ligand files in DOCKED_DIR: <ligand name> + suffix
@@ -80,3 +93,13 @@ RANKING_COLUMNS = (
     COL_MW, COL_LOGP, COL_HBA, COL_HBD, COL_TPSA, COL_ROTORS, COL_SOURCE_URL,
     COL_DOCKED_FILE, COL_STATUS, COL_ENGINE, COL_MINIMIZED,
 )
+# Appended to RANKING_COLUMNS only with --replicates N > 1 (Binding_Affinity
+# is then the best replicate's score; these summarise the in-window replicates).
+COL_N_REPLICATES = "n_replicates"
+COL_SCORE_MEAN = "score_mean"
+COL_SCORE_SD = "score_sd"
+COL_SCORE_MIN = "score_min"
+COL_SCORE_MAX = "score_max"
+RANKING_REPLICATE_COLUMNS = (
+    COL_N_REPLICATES, COL_SCORE_MEAN, COL_SCORE_SD, COL_SCORE_MIN, COL_SCORE_MAX,
+)
diff --git a/autodocker/vspipeline/replicates.py b/autodocker/vspipeline/replicates.py
index 289e539..93fd88d 100644
--- a/autodocker/vspipeline/replicates.py
+++ b/autodocker/vspipeline/replicates.py
@@ -35,27 +35,17 @@ import statistics
 from dataclasses import replace
 from typing import Dict, Iterable, List, Optional, Sequence, Tuple
 
-from vspipeline import outputs, utils
+from vspipeline import defaults, outputs, utils
 
-# Proposed defaults.REPLICATES: one docking per ligand (the 1.x behaviour).
-REPLICATES = 1
+REPLICATES = defaults.REPLICATES
 # Upper bound of --replicates: a typo guard (each replicate is a full docking).
 MAX_REPLICATES = 100
 
-# Proposed outputs.py names.
-REPLICATE_TAG = "_rep"                       # <name>_rep<k>
-REPLICATE_INPUTS_DIR = "replicate_inputs"    # in DOCKED_DIR: staged ligand copies
-REPLICATE_SCORES_CSV = "replicate_scores.csv"
-REPLICATE_SUMMARY_JSON = "replicate_summary.json"
-
-# ranking.csv columns appended after outputs.RANKING_COLUMNS when N > 1.
-COL_N_REPLICATES = "n_replicates"
-COL_SCORE_MEAN = "score_mean"
-COL_SCORE_SD = "score_sd"
-COL_SCORE_MIN = "score_min"
-COL_SCORE_MAX = "score_max"
-RANKING_REPLICATE_COLUMNS = (COL_N_REPLICATES, COL_SCORE_MEAN, COL_SCORE_SD,
-                             COL_SCORE_MIN, COL_SCORE_MAX)
+REPLICATE_TAG = outputs.REPLICATE_TAG
+REPLICATE_INPUTS_DIR = outputs.REPLICATE_INPUTS_DIR
+REPLICATE_SCORES_CSV = outputs.REPLICATE_SCORES_CSV
+REPLICATE_SUMMARY_JSON = outputs.REPLICATE_SUMMARY_JSON
+RANKING_REPLICATE_COLUMNS = outputs.RANKING_REPLICATE_COLUMNS
 REPLICATE_SCORES_COLUMNS = ("Ligand", "Replicate", "Seed", "Status", "Score", "Raw_Score",
                             "Reported", "Docked_File", "Log_File")
 
diff --git a/autodocker/vspipeline/report.py b/autodocker/vspipeline/report.py
index 5bcc1b9..32fcbe9 100644
--- a/autodocker/vspipeline/report.py
+++ b/autodocker/vspipeline/report.py
@@ -78,7 +78,8 @@ _cluster_cache: Dict[Tuple, Dict] = {}
 
 
 def cluster_poses(poses_dir: str, rmsd_threshold: float = defaults.RMSD_THRESHOLD,
-                  ligands: Optional[List[str]] = None) -> Dict[str, List[Dict]]:
+                  ligands: Optional[List[str]] = None,
+                  pose_files: Optional[Dict[str, str]] = None) -> Dict[str, List[Dict]]:
     """Greedy clustering of each ligand's docked modes.
 
     The distance is the heavy-atom RMSD by atom order (``calculate_rmsd``:
@@ -90,9 +91,14 @@ def cluster_poses(poses_dir: str, rmsd_threshold: float = defaults.RMSD_THRESHOL
     Results are cached per (poses_dir, threshold, ligands) so the same
     directory isn't re-scanned when the report and results phases both
     request clustering.
+
+    ``pose_files`` ({ligand: pose file in ``poses_dir``}) names each
+    ligand's reported pose file instead of scanning ``poses_dir``
+    (--replicates: the best replicate's <name>_rep<k>_out.pdbqt).
     """
     keep = None if ligands is None else set(ligands)
-    cache_key = (poses_dir, rmsd_threshold, None if keep is None else frozenset(keep))
+    cache_key = (poses_dir, rmsd_threshold, None if keep is None else frozenset(keep),
+                 None if pose_files is None else frozenset(pose_files.items()))
     if cache_key in _cluster_cache:
         return _cluster_cache[cache_key]
     clusters = defaultdict(list)
@@ -103,7 +109,12 @@ def cluster_poses(poses_dir: str, rmsd_threshold: float = defaults.RMSD_THRESHOL
     # Vina (<ligand>_out.pdbqt) and SMINA (<ligand>_smina.pdbqt) write every
     # mode of a ligand as a MODEL in one file; Vina wins if both exist.
     vina_outputs = {}
-    for file in sorted(os.listdir(poses_dir)):
+    if pose_files is not None:
+        # dock_all writes every pose file to poses_dir: keep the file names.
+        vina_outputs = {lig: os.path.basename(path) for lig, path in sorted(pose_files.items())
+                        if path and os.path.isfile(os.path.join(poses_dir,
+                                                                os.path.basename(path)))}
+    for file in sorted(os.listdir(poses_dir)) if pose_files is None else ():
         for suffix in (outputs.POSE_SUFFIX, outputs.SMINA_POSE_SUFFIX):
             if file.endswith(suffix):
                 lig = file[:-len(suffix)]
@@ -655,8 +666,16 @@ def generate_html_report(results_csv: str, poses_dir: str, output_file: str,
             scatter_data.append({"name": lig.get(outputs.COL_LIGAND, '?'),
                                  "mw": mw, "aff": av})
 
+    # --replicates: cluster each ligand's reported (best-replicate) pose file.
+    pose_files = None
+    if ligands and outputs.COL_N_REPLICATES in ligands[0]:
+        base_dir = os.path.dirname(os.path.abspath(results_csv))
+        pose_files = {l.get(outputs.COL_LIGAND): _pose_file(l, l.get(outputs.COL_LIGAND),
+                                                            poses_dir, base_dir)
+                      for l in ligands}
     clusters = (cluster_poses(poses_dir, rmsd_threshold=rmsd_threshold,
-                              ligands=[l.get(outputs.COL_LIGAND) for l in ligands])
+                              ligands=[l.get(outputs.COL_LIGAND) for l in ligands],
+                              pose_files=pose_files)
                 if os.path.isdir(poses_dir) else {})
 
     css = _REPORT_CSS(bool(ranked_ligands))
diff --git a/autodocker/vspipeline/results.py b/autodocker/vspipeline/results.py
index e877985..f128050 100644
--- a/autodocker/vspipeline/results.py
+++ b/autodocker/vspipeline/results.py
@@ -18,6 +18,7 @@ from typing import Dict, List, Optional, Tuple
 
 from vspipeline import defaults  # noqa: E402
 from vspipeline import outputs  # noqa: E402
+from vspipeline import replicates as replicates_mod  # noqa: E402
 from vspipeline.docking import (SIMSCORE_RMSD_LB_CUTOFF, SIMSCORE_RMSD_UB_CUTOFF,
                                  Status)
 
@@ -61,9 +62,13 @@ def _minimized_text(meta: Dict) -> str:
 class ResultsAnalyzer:
     """Writes the ranking and the text reports of one screen."""
 
-    def __init__(self, outdir: str, top_n: int = defaults.TOP_N):
+    def __init__(self, outdir: str, top_n: int = defaults.TOP_N,
+                 replicates: int = defaults.REPLICATES):
         self.outdir = outdir
         self.top_n = top_n
+        # --replicates N > 1: ranking.csv gets the replicate columns and
+        # Top_hits.txt a replicate line (N = 1 leaves both unchanged).
+        self.replicates = replicates
         self.results_file = os.path.join(outdir, outputs.RESULTS_FULL_TXT)
         self.top_hits_file = os.path.join(outdir, outputs.TOP_HITS_TXT)
         self.csv_file = os.path.join(outdir, outputs.RANKING_CSV)
@@ -145,10 +150,14 @@ class ResultsAnalyzer:
             with open(self.csv_file, "w", newline="", encoding="utf-8") as f:
                 writer = csv.writer(f)
                 if metrics_dict:
-                    writer.writerow(outputs.RANKING_COLUMNS)
+                    extra_columns = (outputs.RANKING_REPLICATE_COLUMNS
+                                     if self.replicates > 1 else ())
+                    writer.writerow(outputs.RANKING_COLUMNS + extra_columns)
                     for lig, status, affinity in rows:
                         metrics = metrics_dict.get(lig, {})
                         meta = self.metadata.get(lig, {})
+                        extra = (replicates_mod.ranking_values(metrics)
+                                 if extra_columns else [])
                         writer.writerow([
                             lig, self._affinity_csv(status, affinity),
                             _blank_if_none(metrics.get("simscore")),
@@ -162,7 +171,7 @@ class ResultsAnalyzer:
                             status.value,
                             metrics.get("engine") or metrics.get("vina_tool", ""),
                             _minimized_text(meta),
-                        ])
+                        ] + extra)
                 else:
                     writer.writerow([outputs.COL_LIGAND, outputs.COL_AFFINITY, outputs.COL_STATUS])
                     for lig, status, affinity in rows:
@@ -253,10 +262,19 @@ class ResultsAnalyzer:
                         if metrics.get("engine"):
                             f.write(f"   Engine: {metrics['engine']}\n")
                         f.write(f"   Docked file: {self._dock_file(name, metrics)}\n")
+                        replicate_line = replicates_mod.top_hits_line(metrics)
+                        if replicate_line:
+                            f.write(f"   {replicate_line}\n")
 
                     f.write(f"   Interpretation: {self._interpretation(status, affinity)}\n")
                     if status is Status.FAILED:
-                        f.write(f"   Log: {outputs.DOCKED_DIR}/{name}{outputs.VINA_LOG_SUFFIX}\n")
+                        metrics = (metrics_dict or {}).get(name) or {}
+                        if replicates_mod.has_replicates(metrics):
+                            f.write(f"   Logs: {outputs.DOCKED_DIR}/{name}"
+                                    f"{outputs.REPLICATE_TAG}<k>{outputs.VINA_LOG_SUFFIX}\n")
+                        else:
+                            f.write(f"   Log: {outputs.DOCKED_DIR}/{name}"
+                                    f"{outputs.VINA_LOG_SUFFIX}\n")
 
                 f.write("\n" + "=" * 70 + "\n")
             logger.info(f"[*] Top hits saved: {self.top_hits_file}")
```
