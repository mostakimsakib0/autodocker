# Installation

Two ways to run AutoDocker: **Docker container** (recommended, reproducible)
or **native Python** (lightweight, needs scientific binaries).

## 1. Docker (recommended)

Requires a working Docker installation on Linux, macOS, or Windows (WSL2).

```bash
git clone --recurse-submodules https://github.com/mostakimsakib0/autodocker.git
cd autodocker
docker build -t autodocker .
```

The build compiles Vina, QuickVina 2, Open Babel, and fpocket from pinned
submodules; allow 15-45 minutes on first build. To use a prebuilt image:

```bash
docker pull ghcr.io/mostakimsakib0/autodocker
```

## 2. Native Python

Requires Python 3.9-3.12 and the `vina` and `obabel` executables on `PATH`
(fpocket is optional — the pipeline falls back to a centroid grid).

```bash
# Debian/Ubuntu
sudo apt-get install -y autodock-vina openbabel

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The test suite verifies the installation end-to-end:

```bash
pip install pytest
python -m pytest tests/ -q
```

## First run

```bash
mkdir -p works/ligs
cp protein.pdb works/protein.pdb
cp ligands/*.sdf works/ligs/

# Docker
docker run --rm -v "$PWD/works:/workspace" autodocker

# Native
python3 autodocker/runner.py -r works/protein.pdb -l works/ligs -o works/output
```

Results appear in `works/output/` (`ranking.csv`, `Top_hits.txt`,
`metrics.txt`, `pocket_summary.csv`, `docked/`).

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `Critical tool 'obabel' not found` | `obabel` (Open Babel) not on PATH — install it or use the Docker image |
| `Critical tool 'vina' not found` | `vina` not on PATH — install `autodock-vina` or use Docker |
| `Protein file not found: /workspace/protein.pdb` | container mount wrong — mount the directory containing `protein.pdb` at `/workspace` |
| `No ligand files (.sdf/.pdbqt/...) found` | put ligands in `works/ligs/` (or set `LIGS`) |
| `Receptor has no charges` | receptor PDB was stripped of HETATM/charge info — pass a raw PDB and let the pipeline prepare it |
| Python under 3.9 | upgrade Python or use Docker |

For a reproducible two-run determinism check see `docs/reproducibility.md`.