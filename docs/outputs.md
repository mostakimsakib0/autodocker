# Output files

Everything AutoDocker writes goes to the output directory (`/workspace/output`
in the container, `-o` natively). Text outputs contain no timestamps and no
absolute paths, so two identical runs give byte-identical files; run-level
information (versions, options, checksums, times) is in `run_manifest.json`.
Files that a run does not produce are removed at its start, so a file left by
an earlier run in the same directory is never mistaken for a result.

## Results

| File | Contents |
|---|---|
| `ranking.csv` | One row per ligand, best first (columns below). The file to read. |
| `Top_hits.txt` | The `--top-n` best in-window ligands (default 20) with their scores. |
| `metrics.txt` | Counts (docked, in-window, out of window, failed), a summary of the in-window scores and per-ligand metrics. |
| `Results_full.txt` | Every ligand with all its modes and metrics, in plain text. |
| `results_report.html` | Self-contained report with the ranking, the box and an NGL 3D viewer of the poses; works offline. On by default (`--no-html-report` skips it). |
| `run_manifest.json` | Command line, resolved options, seed, tool paths, versions and SHA-256, the SHA-256 of every input and output, failure counts and status. Checked by `scripts/verify_run.py` (see `reproducibility.md`). |
| `pose_clusters.json` | With `--cluster-poses`: greedy RMSD clustering of each ligand's poses (`--rmsd-threshold`, default 2.0 Å). |
| `consensus_ranking.csv`, `consensus_summary.json` | With `--consensus`: SMINA Vinardo rescoring of Vina's top pose and the rank correlation. |
| `replicate_scores.csv`, `replicate_summary.json` | With `--replicates N`: every replicate's score and the per-ligand summary. |

## Docked poses (`docked/`)

| File | Contents |
|---|---|
| `<ligand>_out.pdbqt` | Vina (or QuickVina) poses, best first, with `REMARK VINA RESULT` lines. |
| `<ligand>_vina.log` | The exact engine command line (`# Command:`), the engine binary and its SHA-256 (`# Engine:`), the engine output and the result line. |
| `<ligand>_smina.pdbqt`, `_smina.log` | With `--smina-only`. |
| `<ligand>_vina_top.pdbqt`, `_smina_rescore.pdbqt` | With `--consensus`. |
| `<ligand>_rep<k>_out.pdbqt` | With `--replicates N`: replicate k, docked with seed `--seed + k`. |

Poses, `receptor.pdbqt` and `grid.conf` are in the docking frame: translated so
that the mean of the receptor heavy atoms is at the origin. `frame_offset.txt`
holds the shift and `scripts/restore_frame.py` maps poses
back onto the input structure.

## Prepared inputs

| File | Contents |
|---|---|
| `protein_clean.pdb` | The receptor after cleaning (chains, first model, one altloc, selenomethionine as methionine). |
| `receptor.pdbqt` | The prepared receptor that was docked. `rigid.pdbqt` and `flex.pdbqt` with flexible residues. |
| `ligands/` | The prepared ligand PDBQT files that were docked (PDBQT inputs are copied unchanged). |
| `ligands_minimized/` | SDF/PubChem ligands after seeded 3D embedding and MMFF94 minimization. |
| `ligand_metadata.json` | Per ligand: source file and format, how its 3D coordinates were made, minimization, and properties (MW, logP, HBA, HBD, TPSA, rotors) when they can be computed (not for PDBQT inputs). |
| `grid.conf` | The docking box (center and size, docking frame). |
| `grid_box.pml` | PyMOL script that draws the box on the receptor: check it before trusting a ranking. |
| `pocket_summary.csv` | fpocket's pockets with scores (the volume column is a Monte Carlo estimate and varies between runs; it is never used to choose a pocket). |
| `receptor_no_reference.pdb`, `reference_ligand.pdb` | With `--box-from-het`: the receptor without the reference group, and the group itself. |

## `ranking.csv` columns

| Column | Meaning |
|---|---|
| `Ligand` | Ligand name (input file name or PubChem name). |
| `Binding_Affinity` | Vina score of the best pose, kcal/mol (with `--replicates`, the best replicate). Not a binding free energy. |
| `SimScore` | Pose consistency, 0 to 1: the mean of the fractions of the other modes within 1.6 Å (Vina's RMSD lower bound) and within 3.2 Å (upper bound) of the best pose. High means the alternative poses agree with the best one. Empty with one mode or SMINA output. |
| `Binding_Modes` | Number of modes written. |
| `Mode_Coverage` | Modes written / modes requested (`--binding-modes`), at most 1. |
| `MW`, `LogP`, `HBA`, `HBD`, `TPSA`, `Rotors` | Ligand properties (Open Babel / RDKit); empty when they cannot be computed (PDBQT inputs). |
| `ZINC_URL` | Historical column name; holds the PubChem URL when the ligand came from PubChem. |
| `Docked_File` | Pose file, relative to the output directory. |
| `Status` | See below. |
| `Engine` | Engine that produced the pose (Vina, QuickVina, SMINA). |
| `Minimized` | `yes` or `no`: whether the ligand was minimized with MMFF94. |
| `n_replicates`, `score_mean`, `score_sd`, `score_min`, `score_max` | Only with `--replicates N > 1`: statistics over the in-window replicates. |

## Status

| Status | Meaning |
|---|---|
| `OK` | Score in the accepted window (−20 ≤ score < `--min-valid-affinity`, default −1 kcal/mol) and below −4 kcal/mol. |
| `WEAK` | Score in the window but −4 kcal/mol or higher: docked, but a weak binder. |
| `OUT_OF_WINDOW` | The engine produced a valid result whose score is outside the window (for example 0.0 when the box misses the receptor). The score is reported but not ranked with the hits, and the run logs a warning. |
| `FAILED` | No score: bad input, engine crash, timeout or unreadable output. More failures than `--max-failure-fraction` (default 0.5) make the run exit with status 3, still writing every output. |

Rows are ordered OK/WEAK by score, then OUT_OF_WINDOW, then FAILED;
`verify_run.py` checks this order.
