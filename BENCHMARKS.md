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
    --max-actives 50 --max-decoys 100 --exhaustiveness 4 --processes 8 --seed 42
python3 scripts/benchmark_performance.py \
    --receptor data/dude/akt1/receptor.pdb --ligands data/dude/akt1/ligands \
    --sizes 10,25,50,100 --processes 1,2,4,8 --repeats 3 \
    --exhaustiveness 2 --outdir benchmarks/results/performance
```

## J1 — Re-docking (docking power)

**Input.** `fetch_casf_rcsb.py` rebuilds core-set complexes from the
deposited RCSB entries of the CASF-2016 core-set IDs
(`benchmarks/data/casf2016_core_set.csv`). These are not the
PDBbind-processed files, so results are not directly comparable with
published CASF numbers. Per entry: first model; one alternate conformer per
residue (blank/`A`, else highest occupancy); receptor = protein atoms plus
metal ions (`--hetatm metals`); ligand = the crystal instance of the HET group
whose heavy-atom formula matches the core-set SMILES and whose identity is
confirmed by InChIKey connectivity block or heavy-atom graph isomorphism,
with bond orders from the Chemical Component Dictionary (RCSB ModelServer).
Peptide and multi-residue ligands cannot be built and are listed with the
reason in `manifest.csv`.

**Protocol.**

| Item | Setting |
|---|---|
| Complexes | seeded random sample of the core set: `--sample 40 --sample-seed 42` (the same selection in both scripts); no receptor-size filter |
| Start conformer | `--start-conformer rdkit` (default): one ETKDGv3 conformer from the CCD molecule (randomSeed 42, single thread, then MMFF94, at most 2000 iterations; deterministic for a given RDKit version), so ring conformations, bond lengths and angles do not come from the crystal pose. `--start-conformer crystal` reproduces the 1.1.x protocol for comparison |
| Receptor | all chains; metal ions kept (runner `--detect-metals`) |
| Box | cube centred on the crystal-ligand heavy-atom centroid, edge = max(22.5 Å, ligand extent + 2 × 6 Å), given in input-PDB coordinates |
| Docking | Vina, seed 42, exhaustiveness 8, 9 modes, energy range 3 kcal/mol |
| RMSD | `obrms`: heavy atoms, symmetry-corrected, no superposition, against the crystal pose moved into the docking frame with `frame_offset.txt` |
| Top-1 | mode 1 (first MODEL). If its RMSD cannot be computed the complex is `rmsd_failed`; it is never replaced by mode 2 |
| Success | RMSD ≤ 2.0 Å |

A complex is reused on resume only if its stored `params.json` matches the
current docking parameters.

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
| Ligands | per target a seeded sample of 50 actives and 100 decoys from DUD-E, converted to PDBQT with Open Babel; MOL2 intermediates kept out of the docked folder; conversion failures recorded (`convert_failures`) and replaced by the next ligand in the seeded order |
| Receptor | DUD-E `receptor.pdb`, all chains (`--chain all`) |
| Box | `--box-source crystal` (default): cube around DUD-E's `crystal_ligand.mol2`, same rule as J1 (max(22.5 Å, extent + 12 Å)), passed as `--grid-center/--grid-size`. `--box-source fpocket` reproduces the 1.1.x protocol (fpocket pocket; frame alignment in `protein.py::_receptor_centroid`) |
| Docking | Vina, seed 42, exhaustiveness 4, 1 mode |
| Failed ligands | every sampled ligand stays in N; a ligand without a usable score (docking failure, missing from `ranking.csv`) is ranked last, tied with the other failures; failures are counted per class |

**Metrics** (N ligands, n actives, n_top = ⌈xN⌉; lower score = better; ties
resolved by their expected value under random tie-breaking, so no metric
depends on input order):

- ROC AUC (tie-aware; equals the Mann-Whitney U statistic / n·(N−n)).
- EF at 1 % and 5 %: EF = (TP/n_top)/(n/N), reported with the maximum
  achievable EF_max = (min(n_top, n)/n_top)/(n/N). With 150 ligands, n_top is
  2 (1 %) and 8 (5 %); EF on fewer than 10 ligands is unstable and the script
  says so.
- BEDROC, α = 80.5 (Truchon & Bayly, *J. Chem. Inf. Model.* 2007, 47, 488).
  With 50 actives in 150 ligands α·Rₐ ≈ 27 ≫ 1, so BEDROC saturates and is not
  comparable with values from full DUD-E sets; the script warns.
- ROC enrichment at 1 % and 5 % false-positive rate (TPR/FPR).
- 95 % percentile-bootstrap CIs (`--ci`, `--bootstrap 1000`), stratified:
  each replicate redraws the actives and the decoys separately, with
  replacement, as (score, label) pairs.
- One-sided Mann-Whitney U test (normal approximation, tie-corrected),
  H1: actives score better than decoys. p-values are stored at full double
  precision.

**Results:**

| Target | Family | N (failed act./dec.) | AUC (95 % CI) | EF1% (max) | EF5% (max) | BEDROC | p (M-W, one-sided) |
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

Files: `enrichment.csv` (one row per target, all metrics and CIs),
`<target>/metrics.csv` (ranked ligands incl. failures), `<target>/ranking.csv`
(the pipeline's ranking) and `summary.json` (protocol and per-target records
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
one mode, at 1/2/4/8 worker processes; each cell is run 3 times
(`--repeats 3`), interleaved across cells. Recorded per run: wall time,
split into preparation (start to the runner's docking step), docking and
post-processing from the time-stamped runner log; peak memory of the whole
process tree (summed RSS and, where readable, summed PSS, sampled every
0.5 s); the largest single process (`ru_maxrss`); Vina threads per process
(the `--cpu` value in the logged command; by default the runner gives each
of `-p N` processes cores/N threads). Per cell: median and range of the wall
time, throughput from the median, speedup vs. `-p 1` at the same library
size and parallel efficiency (speedup / p). `summary.json` records the
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
