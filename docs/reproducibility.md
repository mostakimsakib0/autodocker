# Reproducibility

AutoDocker aims for the **same input + seed + software versions to give
identical docking results**, on Docker or natively. This page lists what is
fixed, what is not, and how it is tested.

## What is fixed

- Docking seed: `--seed 42` by default (`VS_SEED`), passed to every Vina /
  SMINA call. Seeds are 1 to 2147483647: the engines take a 32-bit seed,
  Vina treats 0 as "pick a random seed" and RDKit rejects negative seeds, so
  other values are refused (exit 2).
- Pinned software in the Docker image: Vina 1.2.7, QuickVina 2, SMINA,
  fpocket 4.2.3 and Open Babel 3.2.1 are compiled from pinned submodules
  (with `-ffp-contract=off`, verified to leave amd64 binaries and results
  unchanged); the other system packages come from a pinned Debian snapshot
  (`DEBIAN_SNAPSHOT`); base images are pinned by digest; Python
  dependencies, including RDKit, are pinned in `requirements.txt`. A native
  install uses whatever versions are on PATH; `run_manifest.json` records
  the path, version and SHA-256 of every tool used, `openbabel_version`, and
  the RDKit version. Open Babel older than 3.2.0 runs with a warning (see
  [`openbabel.md`](openbabel.md)).
- Every docking command is recorded at the top of `docked/*_vina.log` (also
  for failed attempts), so the exact inputs passed to Vina can be recovered
  and re-run (`scripts/benchmark_vina_parity.py` does this).
- Open Babel's random generator. Every Open Babel call runs with
  `OB_RANDOM_SEED` set to the run's `--seed` (Open Babel >= 3.2 reads it;
  3.1 ignores it and seeds from the clock). `run_manifest.json` records the
  value and the policy under `openbabel_random_seed`.
- Receptor hydrogens. Hydrogens in the input receptor are removed and the
  receptor is protonated at `--ph` by Open Babel. Open Babel places some
  hydrogens in a random direction when the neighbouring atoms do not fix
  one (isolated atoms such as waters, two-atom fragments, distorted
  residues). A fixed seed makes that draw repeatable, but it is still a
  random draw: AutoDocker converts once with the seed S and once with S+1
  (and S+2 when those agree), finds every hydrogen that moves with the seed,
  and re-places it deterministically from the heavy atoms, checking that
  Vina still bonds it only to its parent atom. This matters because Vina
  decides which atoms are hydrogen-bond donors from hydrogen positions. The
  prepared receptor is therefore the same for every seed.
- Ligand hydrogens: the same check is applied to every ligand conversion
  that adds or keeps hydrogens.
- 3D coordinates for 2D ligands come from RDKit ETKDGv3 seeded with `--seed`
  (single thread), followed by MMFF94. Without RDKit the pipeline falls back
  to Open Babel `--gen3d` (seeded with `OB_RANDOM_SEED` on Open Babel >= 3.2,
  not reproducible on 3.1) and logs a warning; `ligand_metadata.json`
  records the method per ligand.
- Open Babel always runs single-threaded (`OMP_NUM_THREADS=1`; the image's
  build has OpenMP disabled).
- Pocket choice: pockets are ranked by fpocket score, then druggability, then
  pocket number; fpocket's pocket atoms and scores are deterministic.
- Text outputs and the HTML report (`results_report.html`) contain no
  timestamps or absolute paths; timestamps are only in `run_manifest.json`.
- Benchmark scripts seed every sampling step (ligand subsets, complex
  subsets, bootstrap, RDKit start conformers) and record their parameters
  and tool versions in `summary.json`.

## What is not fixed

- **fpocket volumes.** fpocket estimates pocket volumes by Monte Carlo with a
  clock-based seed, so the `Volume` column of `pocket_summary.csv` varies
  between runs. Pocket selection never uses it.
- **Software versions.** A different Vina, Open Babel or RDKit version can
  change results; the seed fixes results only for a fixed software stack.
  The Open Babel 3.1.1 -> 3.2.1 change of AutoDocker 2.0 is quantified in
  [`openbabel.md`](openbabel.md).
- **CPU architecture.** Results are byte-reproducible on a given
  architecture. Whether linux/amd64 and linux/arm64 give identical results is
  measured in CI (the cross-architecture job) and not guaranteed: glibc's
  math library (used by Vina and by Open Babel's MMFF94 minimization and
  Gasteiger charges; it chooses FMA or non-FMA routines at run time
  depending on the CPU) and SMINA's vectorised code may differ between
  architectures, and in principle between x86 CPUs with and without FMA.

## Thread and process count

Vina seeds each Monte Carlo run from the job seed, independent of how runs
are assigned to threads, so seeded Vina output is not expected to depend on
`--cpu`. This is tested, not assumed:

- `tests/test_e2e.py::test_vina_cpu_count_does_not_change_output` docks one
  ligand with `--cpu 1`, `2` and `4` and requires byte-identical output
  files; a manual check with 1, 2 and 8 threads on one ligand also gave
  byte-identical output.
- `tests/test_e2e.py::test_process_count_does_not_change_results` runs the
  docking step with `-p 1` and `-p 2` (which also halves the threads per Vina
  process) and requires byte-identical `ranking.csv` and poses.

These tests cover the committed fixtures and the Vina build on the test
machine; they are evidence, not a proof for every input.

The threads given to each Vina/QuickVina/SMINA process (`--threads`, else
the available CPUs divided by `-p` when `-p > 1`, else the engine's default
of all cores) are recorded in `run_manifest.json` as `threads_per_process`.

## Resuming

Docking results are checkpointed per ligand under a fingerprint of the
receptor and flexible-residue files, the grid file, the docking parameters
and the primary engine binary (path, SHA-256, version); a change of any of
them discards the checkpoint. The QuickVina fallback (`--qvina-fallback`) is
not part of the fingerprint: a result of the fallback is tagged with the
fallback binary's SHA-256 and docked again on its own when the fallback
changes or is turned off, while the primary engine's results are kept.
Paths in the checkpoint are relative to the output directory, so a resumed
run gives the same `ranking.csv` from any working directory and after the
output directory was moved. Results are collected in library order, never in
the order the workers finish; a ligand whose worker process dies is docked
again in a fresh worker and reported `FAILED` if that worker dies too.

## End-to-end determinism check

```bash
scripts/reproducibility_check.sh protein.pdb ligs/ /tmp/repro \
    -- --grid-center 12.3,4.5,6.7 --grid-size 22.5,22.5,22.5
```

The script runs the full pipeline three times with the same seed (`-p 1`,
`-p 1`, `-p $PROCESSES`) and compares `ranking.csv` (with the run directory
normalised away), every docked pose file and `results_report.html` byte for
byte (the report embeds the NGL viewer from `/ngl/ngl.js` in the image, else
from `tools/ngl` or a checksum-verified download, so all runs must find the
same source). Arguments after
`--` go to every run; an explicit box keeps fpocket out of the comparison.

| Exit code | Meaning |
|---|---|
| 0 | all three runs identical (`DETERMINISTIC`) |
| 1 | results differ (`NON-DETERMINISTIC`, differences listed) |
| 2 | usage error or missing tool |
| 3 | a pipeline run failed (its log is shown) |

`tests/test_e2e.py` asserts the same for the docking stage on committed
fixtures: a binding-site fragment of PDB 2FXS and three drug-like ligands
with rotatable bonds, docked twice with `-p 1` and once with `-p 2`.
`tests/test_runner_e2e.py::test_runs_are_byte_identical_including_the_html_report`
runs the whole pipeline (receptor and ligand preparation included) twice
with `-p 1` and once with `-p 2` and compares every output except
`run_manifest.json`, logs and the docking checkpoint, `results_report.html`
included.

## Versioned results

Pin the pipeline version in your own records:

```bash
python3 autodocker/runner.py --version   # e.g. AutoDocker 2.0.0
```

Container images are tagged by git tag and commit SHA
(`ghcr.io/mostakimsakib0/autodocker:<tag>`).

## Parity with raw Vina

`scripts/benchmark_vina_parity.py` re-executes the exact recorded Vina
command for each re-docked complex and checks that every mode's affinity and
the affinity reported in `ranking.csv` equal raw Vina's output (results:
`BENCHMARKS.md`, J3).
