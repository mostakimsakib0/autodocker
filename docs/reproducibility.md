# Reproducibility

AutoDocker aims for the **same input + seed + software versions to give
identical docking results**, on Docker or natively. This page lists what is
fixed, what is not, and how it is tested.

## What is fixed

- Docking seed: `--seed 42` by default (`VS_SEED`), passed to every Vina /
  SMINA call.
- Pinned software in the Docker image: the Dockerfile compiles Vina 1.2.7
  (tagged release; images up to v1.1.4 used a later commit with identical
  C++ source), QuickVina 2, and fpocket 4.2.3 from pinned submodules, and
  installs Open Babel 3.1.1 from the Debian bookworm repository (the
  `tools/openbabel` submodule is not used by the image); Python and pip
  dependencies are pinned (`python:3.12-slim-bookworm`, exact versions in
  `requirements.txt`). A native install uses whatever versions are on PATH;
  record them (`vina --version`, `obabel -V`).
- Every docking command is recorded in `docked/*_vina.log`, so the exact
  inputs passed to Vina can be recovered and re-run
  (`scripts/benchmark_vina_parity.py` does this).
- Benchmark scripts seed every sampling step (ligand subsets, complex
  subsets, bootstrap, RDKit start conformers) and record their parameters
  and tool versions in `summary.json`.

## What is not fixed by the seed

- **3D generation of ligands without coordinates.** SMILES or 2D inputs are
  embedded with Open Babel `--gen3d`, which is not seeded. Library
  preparation from such inputs is not guaranteed to be identical between
  runs; supply 3D ligands (SDF/MOL2/PDBQT) when bit-identical results
  matter.
- **fpocket.** fpocket seeds its random number generator from the clock, so
  the detected pockets, and therefore the automatically chosen box, can
  differ between runs. For reproducible boxes give the box explicitly
  (`--grid-center X,Y,Z --grid-size SX,SY,SZ`, optionally `--no-fpocket`).
- **Receptor protonation.** Hydrogens are added to the receptor with Open
  Babel at `--ph`; making this step deterministic is being handled in a
  separate change. Until that change is released, compare `receptor.pdbqt`
  between runs (e.g. `sha256sum`) before comparing docking scores.
- **Software versions.** A different Vina, Open Babel or RDKit version can
  change results; the seed fixes results only for a fixed software stack.

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

## End-to-end determinism check

```bash
scripts/reproducibility_check.sh protein.pdb ligs/ /tmp/repro \
    -- --grid-center 12.3,4.5,6.7 --grid-size 22.5,22.5,22.5
```

The script runs the full pipeline three times with the same seed (`-p 1`,
`-p 1`, `-p $PROCESSES`) and compares `ranking.csv` (with the run directory
normalised away) and every docked pose file byte for byte. Arguments after
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
