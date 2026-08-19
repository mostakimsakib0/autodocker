# Contributing to AutoDocker

Thanks for contributing! AutoDocker is a scientific software project, so
changes that affect results must not silently break reproducibility.

## Getting started

1. Fork the repository and clone with submodules:
   `git clone --recurse-submodules https://github.com/mostakimsakib0/autodocker.git`
2. Install dev dependencies: `pip install -r requirements.txt pytest`
   You also need the `vina` and `obabel` binaries on `PATH` (the test suite
   imports `autodocker/runner.py`, which refuses to start without them).
3. Run the tests: `python -m pytest tests/ -q`

## Test and benchmark expectations

- Every change must keep `python -m pytest tests/ -q` green.
- Changes that alter docking inputs, grids, or scoring must not regress the
  scientific benchmarks in `BENCHMARKS.md` (re-docking success, DUD-E
  enrichment, raw-Vina parity).
- Reproducibility is a hard requirement: new code paths must be seeded and
  deterministic, and any change in default behavior must be documented.

## Pull requests

- Keep pull requests focused on a single concern and rebase on `main`.
- Add or update tests for new behavior.
- Update `README.md`, `BENCHMARKS.md`, or `paper/paper.md` when user-facing
  behavior or benchmark results change.

## Bug reports

Open an issue with a minimal failing example (input files + the exact
command). Include your environment: OS, Python version, and output of
`python3 autodocker/runner.py --version`.

## License

By contributing, you agree that your contributions are licensed under the
MIT license (see `LICENSE`).