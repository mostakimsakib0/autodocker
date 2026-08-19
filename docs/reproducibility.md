# Reproducibility

AutoDocker is designed so the **same input + seed reproduces identical
results**, on Docker or natively.

## What is fixed

- Docking seed: `--seed 42` by default (`VS_SEED`).
- Ligand-library sampling, receptor prep, grid placement and docking all
  derive from the seed.
- Pinned software: the Dockerfile compiles Vina 1.2.5, QuickVina 2,
  Open Babel 3.x and fpocket 4.2.3 from pinned submodules; Python and pip
  dependencies are pinned (`python:3.12-slim-bookworm`, exact versions in
  `requirements.txt`).
- Every docking command is recorded in `docked/*_vina.log`, so the exact
  inputs passed to Vina are recoverable after the fact.
- Benchmark scripts are deterministic for a fixed seed (`BENCHMARKS.md`).

## End-to-end determinism check

Run the pipeline twice on identical inputs and confirm `ranking.csv` is
byte-identical:

```bash
scripts/reproducibility_check.sh protein.pdb ligs/ /tmp/repro
```

The script reports `DETERMINISTIC` (exit 0) or `NON-DETERMINISTIC` (exit 1).

Strictly deterministic *values* are also asserted in `tests/test_e2e.py`
(the golden determinism test).

## Versioned results

Pin the pipeline version in your own records:

```bash
python3 autodocker/runner.py --version   # AutoDocker 1.0.0
```

Container images are tagged by git tag and commit SHA
(`ghcr.io/mostakimsakib0/autodocker:<tag>`).

## Bit-for-bit parity with raw Vina

Parity is not just "seeded": `scripts/benchmark_vina_parity.py` re-executes
the exact recorded Vina command and shows every affinity reported by the
pipeline equals raw Vina's score (222/222 modes, 0.0 kcal/mol max diff).