# Benchmarks

Validation of the AutoDocker pipeline: re-docking accuracy (J1), screening
enrichment (J2), parity with the raw Vina binary (J3) and throughput (J4).
All benchmarks are run by committed scripts in `scripts/`; the input
datasets (CASF-2016 core-set entries from the RCSB PDB, DUD-E) are
downloaded on demand. Each script writes its tables and a `summary.json`
(protocol, tool versions, host) to its `--outdir`; the defaults below
(`benchmarks/results/...`) are created by the run and are not part of the
repository.

> **Status.** The protocols below were corrected after a review of the 1.1.x
> benchmark scripts (EF definition, bootstrap resampling, failed-ligand
> handling, top-1 RMSD alignment, success-rate denominators, timing
> statistics). Numbers obtained with the 1.1.x scripts are withdrawn. Every
> number on this page is a placeholder until the suite is rerun with
> AutoDocker 2.0.0.

Run the suite:

```bash
# J1 input: CASF-2016 core-set entries rebuilt from the RCSB PDB
python3 scripts/fetch_casf_rcsb.py --outdir data/casf2016_rcsb \
    --sample 40 --sample-seed 42
python3 scripts/benchmark_redock.py --data-dir data/casf2016_rcsb \
    --sample 40 --sample-seed 42 --outdir benchmarks/results/redock \
    --processes 8 --threads 1 --seed 42
# J3 re-runs the Vina commands logged by J1
python3 scripts/benchmark_vina_parity.py \
    --workdir benchmarks/results/redock/work --outdir benchmarks/results/parity
python3 scripts/benchmark_enrichment.py \
    --targets akt1,braf,cdk2,egfr,fa10,hivpr,andr,aldr,hmdh,aces \
    --workdir data/dude --outdir benchmarks/results/enrichment \
    --max-actives 50 --max-decoys 500 --exhaustiveness 4 --processes 8 --seed 42
python3 scripts/benchmark_performance.py \
    --receptor data/dude/akt1/receptor.pdb --ligands data/dude/akt1/ligands \
    --sizes 10,25,50,100 --processes 1,2,4,8 --repeats 3 --threads 1 \
    --exhaustiveness 2 --outdir benchmarks/results/performance
```

## J1 — Re-docking (docking power)

**Input.** `fetch_casf_rcsb.py` rebuilds core-set complexes from the
deposited RCSB entries of the CASF-2016 core-set IDs
(`benchmarks/data/casf2016_core_set.csv`). These are not the
PDBbind-processed files, so results are not directly comparable with
published CASF numbers. Per entry: first model; one alternate conformer per
residue (blank/`A`, else highest occupancy); receptor = protein atoms, the
modified amino acids the pipeline treats as protein
(`vspipeline.protein.MODIFIED_RESIDUES`: MSE, TPO, SEP, ...) and metal ions
(`--hetatm metals`): HETATM residues made of a single metal atom, the rule of
the runner's `--detect-metals`, so the iron of a heme is dropped together
with the heme instead of being kept alone.

Ligand: the crystal instance of the HET group whose heavy-atom formula
matches the core-set SMILES (coordinates from the RCSB ModelServer). One
RDKit molecule is built from it: crystal coordinates and heavy-atom graph;
bond orders, formal charges and hydrogen counts (the protonation state)
from the core-set SMILES (`AllChem.AssignBondOrdersFromTemplate`, which
requires the crystal graph to be isomorphic to the SMILES graph); stereo
from the crystal coordinates; explicit hydrogens added once. In 9 of its 285
SMILES the core-set CSV writes deprotonated phosphates charge-separated
(`[P+](=O)(=O)[O-]`, a valence-6 phosphorus); these are rewritten as the
`P(=O)([O-])[O-]` dianion (same atoms and hydrogens) first. The molecule is
written as `{pdbid}_ligand.sdf` and `{pdbid}_ligand.mol2` (crystal start
and RMSD reference), and the RDKit start conformer is embedded from the
same molecule, so both start modes dock the same topology and protonation.
A complex is skipped with the reason in `manifest.csv` unless the MOL2's
InChIKey equals the SDF's, also after `obabel -h` (the pipeline's MOL2 ->
PDBQT conversion; with every hydrogen explicit it adds none), and matches
the core-set InChIKey in connectivity block and protonation flag; the
start conformer is written only if its InChIKey equals the crystal
MOL2's. `manifest.csv` records the ligand, core-set and CCD InChIKeys, the
full-key (stereo) match, and `identity_check`: whether the CCD bond orders
already gave the core-set connectivity (`ccd-inchikey-connectivity`) or the
SMILES bond orders were needed (`heavy-atom-graph+smiles-bond-orders`, e.g.
1y6r, whose ModelServer SDF lacks a ring double bond). Peptide and
multi-residue ligands cannot be built and are listed with the reason.
RDKit is required. (Data directories written by earlier versions of the
script lost hydrogens in an Open Babel `sdf -> mol2` conversion of the
H-less ModelServer SDF in 12 of 39 complexes; `benchmark_redock.py` refuses
an RDKit start conformer whose InChIKey differs from the crystal MOL2.)

**Protocol.**

| Item | Setting |
|---|---|
| Complexes | seeded random sample of the core set: `--sample 40 --sample-seed 42`, or `--first N` (first N IDs of the CSV in sorted order); both scripts select through the same function, so the same options select the same complexes; `--pdbids` is case-insensitive; no receptor-size filter |
| Start conformer | `--start-conformer rdkit` (default): one ETKDGv3 conformer of the same protonated molecule as the crystal-start MOL2 (InChIKeys must match; randomSeed 42, single thread, then MMFF94, at most 2000 iterations; deterministic for a given RDKit version), so ring conformations, bond lengths and angles do not come from the crystal pose. `--start-conformer crystal` reproduces the 1.1.x protocol for comparison |
| Receptor | all chains; metal ions kept (runner `--detect-metals`) |
| Box | cube centred on the crystal-ligand heavy-atom centroid, edge = max(22.5 Å, ligand extent + 2 × 6 Å), given in input-PDB coordinates |
| Docking | Vina, seed 42, exhaustiveness 8, 9 modes, energy range 3 kcal/mol |
| RMSD | `obrms`: heavy atoms, symmetry-corrected, no superposition, against the crystal pose moved into the docking frame with `frame_offset.txt` |
| Top-1 | mode 1 (first MODEL). If its RMSD cannot be computed the complex is `rmsd_failed`; it is never replaced by mode 2 |
| Success | RMSD ≤ 2.0 Å. Top-1 and best-of-9 successes are both counted among the evaluated complexes (docked, mode-1 RMSD computed) and reported over the evaluated and over the attempted complexes (docking and RMSD failures count as misses); medians are over the evaluated complexes |

A complex is reused on resume only if its stored `params.json` matches the
current docking parameters, the AutoDocker and RDKit versions and the
SHA-256 of the complex's input files (receptor, ligand MOL2/SDF,
precomputed start conformer).

**Results** (AutoDocker 2.0.0, Vina version as recorded in `summary.json`):

| Metric | Value |
|---|---|
| Complexes sampled / built / attempted / evaluated | TBD (rerun with 2.0.0) |
| Excluded complexes (with reasons) | TBD (rerun with 2.0.0) |
| Top-1 success, over evaluated complexes | TBD (rerun with 2.0.0) |
| Top-1 success, over attempted complexes (failures count as misses) | TBD (rerun with 2.0.0) |
| Best-of-9 success, over evaluated / attempted | TBD (rerun with 2.0.0) |
| Median top-1 RMSD | TBD (rerun with 2.0.0) |
| Median best-of-9 RMSD | TBD (rerun with 2.0.0) |
| Top-1 success with `--start-conformer crystal` (comparison) | TBD (rerun with 2.0.0) |

Files: `results.csv` (per complex: status, ranking.csv affinity, mode-1
affinity, top-1/best RMSD, RMSD of every mode, reason for any exclusion) and
`summary.json` (counts with explicit denominators, medians, excluded
complexes, protocol).

## J2 — DUD-E enrichment (screening power)

**Protocol.**

| Item | Setting |
|---|---|
| Targets | akt1, braf, cdk2, egfr (kinases), fa10 (serine protease), hivpr (aspartic protease), andr (nuclear receptor), aldr, hmdh (reductases), aces (esterase) |
| Ligands | per target a seeded sample of 50 unique active and 500 unique decoy compounds (1:10; DUD-E's own ratio is about 1:50). DUD-E lists a compound once per protonation state/tautomer under the same ID (e.g. akt1: 423 active entries for 293 compounds, hivpr 1395 for 536; decoy IDs repeat too), so entries are grouped by compound ID, the sorted IDs are shuffled with the seed, and each sampled compound is docked as ONE entry: the first of its entries in DUD-E's `*_final.mol2.gz`. DUD-E documents no preferred state among them; the first entry is a fixed, reproducible choice, not a claim about the dominant species. Entry and compound counts per file are in `summary.json`, the entry used per compound in `<target>/ligands.csv`. Files are named per compound (`A_<CHEMBL id>.pdbqt`, `D_<ZINC id>.pdbqt`); converted to PDBQT with Open Babel; MOL2 intermediates kept out of the docked folder; conversion failures recorded (`convert_failures`) and replaced by the next compound in the seeded order |
| Receptor | DUD-E `receptor.pdb`, all chains (`--chain all`) |
| Box | `--box-source crystal` (default): cube around DUD-E's `crystal_ligand.mol2`, same rule as J1 (max(22.5 Å, extent + 12 Å)), passed as `--grid-center/--grid-size`. `--box-source fpocket` reproduces the 1.1.x protocol (fpocket pocket; frame alignment in `protein.py::_receptor_centroid`) |
| Docking | Vina, seed 42, exhaustiveness 4, 1 mode |
| Failed ligands | every sampled ligand stays in N; a ligand without a usable score (docking failure, missing from `ranking.csv`) is ranked last, tied with the other failures; failures are counted per class |

**Metrics** (N ligands, n actives, n_top = ⌈xN⌉; lower score = better; ties
resolved by their expected value under random tie-breaking, so no metric
depends on input order). With 50 actives and 500 decoys, N = 550.

Primary (the reported result):

- ROC AUC (tie-aware; equals the Mann-Whitney U statistic / n·(N−n)), with
  its bootstrap CI.
- EF at 5 %: EF = (TP/n_top)/(n/N), n_top = 28, reported with the maximum
  achievable EF_max = (min(n_top, n)/n_top)/(n/N) = 11.

Secondary (computed and written, but only with these caveats; `summary.json`
carries them under `secondary_metrics_caveats`, and per target in
`warnings`):

- EF at 1 %: n_top = 6 ligands (< 10); one active more or less in the top
  changes EF1% by 1.8, so it is unstable at this library size.
- BEDROC, α = 80.5 (Truchon & Bayly, *J. Chem. Inf. Model.* 2007, 47, 488).
  It assumes α·Rₐ ≪ 1; here α·Rₐ = 80.5 × 50/550 = 7.3, so BEDROC saturates
  and is not comparable with values from full DUD-E sets (Rₐ ≈ 0.02).

Also reported:

- ROC enrichment at 1 % and 5 % false-positive rate (TPR/FPR).
- 95 % percentile-bootstrap CIs (`--ci`, `--bootstrap 1000`), stratified:
  each replicate redraws the actives and the decoys separately, with
  replacement, as (score, label) pairs.
- One-sided Mann-Whitney U test (normal approximation, tie-corrected),
  H1: actives score better than decoys. p-values are stored at full double
  precision.

**Results:**

Primary columns: AUC and EF5%. EF1% and BEDROC (†) are secondary, see the
caveats above.

| Target | Family | N (failed act./dec.) | AUC (95 % CI) | EF5% (max) | EF1% (max) † | BEDROC † | p (M-W, one-sided) |
|---|---|---|---|---|---|---|---|
| akt1 | kinase | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| braf | kinase | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| cdk2 | kinase | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| egfr | kinase | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| fa10 | serine protease | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| hivpr | aspartic protease | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| andr | nuclear receptor | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| aldr | reductase | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| hmdh | reductase | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| aces | esterase | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |

Files: `enrichment.csv` (one row per target, all metrics and CIs, DUD-E entry
and compound counts), `<target>/metrics.csv` (ranked ligands incl. failures),
`<target>/ligands.csv` (docked compounds and the DUD-E entry used),
`<target>/ranking.csv` (the pipeline's ranking) and `summary.json` (protocol
incl. the sampling rule and primary/secondary metrics, and per-target records
including warnings and conversion failures).

`--ingest DIR` scores an external `ranking.csv` per target against the same
prepared ligand set, so other docking pipelines can be compared on identical
inputs.

## J3 — Parity with raw Vina

**Protocol.** For every complex docked in J1, the Vina command recorded in
`docked/*_vina.log` (the last `# Command:` line, i.e. the attempt that wrote
the output) is re-executed with the raw binary, writing to a separate file.
A complex passes when

1. every mode's `REMARK VINA RESULT` affinity equals the pipeline output's,
   and
2. the affinity reported in `ranking.csv` equals the raw run's mode-1
   `REMARK VINA RESULT` value, numerically and as a string at the precision
   Vina printed.

A complex that cannot be checked (missing log, binary failure, ligand absent
from `ranking.csv`) counts as a failure, and a run that checks nothing fails.
Whether the two output files are byte-identical is recorded as well. The
binary's `--version` is stored in `summary.json`.

**Results:**

| Metric | Value |
|---|---|
| Complexes checked / passed | TBD (rerun with 2.0.0) |
| Modes compared / matched | TBD (rerun with 2.0.0) |
| ranking.csv equals raw mode 1 | TBD (rerun with 2.0.0) |
| Output files byte-identical | TBD (rerun with 2.0.0) |
| Vina binary and version | TBD (rerun with 2.0.0) |

Files: `parity.csv`, `summary.json`.

## J4 — Throughput and parallel scaling

**Protocol.** The pipeline is run on the akt1 DUD-E receptor (fpocket box)
with seeded subsamples of 10/25/50/100 prepared ligands, exhaustiveness 2,
one mode, at 1/2/4/8 worker processes with one Vina thread per process
(`--threads 1`, the default); each cell is run 3 times (`--repeats 3`),
interleaved across cells. Without `--threads` the runner gives `-p 1` all
cores and each of `-p N` processes cores/N threads, so every cell would use
the same CPU budget and speedup/p would not measure parallel efficiency.
Recorded per run: wall time, split into preparation (start to the runner's
docking step), docking and post-processing from the time-stamped runner
log; peak memory of the whole process tree (summed RSS and, where readable,
summed PSS, sampled every 0.5 s); the largest single process (`ru_maxrss`);
ligands docked successfully; Vina threads per process (requested, and
observed as the `--cpu` value in the logged command) and total threads
(processes × threads). Per cell: median and range of the wall time,
throughput (median over repeats of successfully docked ligands per minute,
not library size per minute), speedup vs. `-p 1` at the same library size
and parallel efficiency against total threads = speedup / (total threads /
total threads at p = 1), i.e. speedup / p with one thread per process
(None when a cell's thread count is unknown). `summary.json` records the
logical and physical core counts, CPU model, Vina and AutoDocker versions.

**Results:**

| Size | p=1 | p=2 | p=4 | p=8 |
|---|---|---|---|---|
| 10 | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| 25 | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| 50 | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |
| 100 | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) | TBD (rerun with 2.0.0) |

(each cell: median wall time [min–max], throughput, speedup vs. p=1)

| Metric | Value |
|---|---|
| Fixed preparation time (median) | TBD (rerun with 2.0.0) |
| Peak process-tree memory (RSS / PSS) | TBD (rerun with 2.0.0) |
| Host (logical / physical cores, CPU model) | TBD (rerun with 2.0.0) |

Files: `runs.csv` (every run), `performance.csv` (per cell), `summary.json`.

## Environment

Record the environment from the `summary.json` files of the rerun: AutoDocker
version (`runner.py --version`), Vina version (`vina --version` of the binary
used; the Docker image builds Vina 1.2.7), Open Babel version, RDKit version
(J1 start conformers), host CPU and core counts. Benchmarks are deterministic
for fixed seeds on a given software stack (see `docs/reproducibility.md`),
except wall times.
