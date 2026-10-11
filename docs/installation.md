# Installation

Two ways to run AutoDocker: **Docker container** (recommended, reproducible)
or **native Python** (lightweight, needs scientific binaries).

## 1. Docker (recommended)

Requires a working Docker installation on Linux, macOS, or Windows (WSL2).
Images run natively on x86-64 and ARM64 (Apple Silicon via Docker Desktop,
AWS Graviton); `docker pull` selects the right architecture. Results are
byte-reproducible on a given architecture; see `docs/reproducibility.md` for
amd64 versus arm64.

```bash
git clone --recurse-submodules https://github.com/mostakimsakib0/autodocker.git
cd autodocker
docker build -t autodocker \
  --build-arg AUTODOCKER_VERSION="$(cat autodocker/VERSION)" \
  --build-arg VCS_REF="$(git rev-parse HEAD)" .
```

The build arguments set the image's version label and source revision (the
build fails if `AUTODOCKER_VERSION` disagrees with `autodocker/VERSION`).
The build compiles Vina, QuickVina 2, SMINA, fpocket and Open Babel 3.2.1
from pinned submodules (clone with `--recurse-submodules`, or run
`git submodule update --init` first); allow 20-60 minutes on first build. To use a prebuilt image:

```bash
docker pull ghcr.io/mostakimsakib0/autodocker
```

## 2. Native Python

Requires Python 3.9-3.12 and the `vina`, `obabel` and `fpocket` executables
on `PATH` (fpocket is needed unless you pass `--no-fpocket` or an explicit
box). RDKit (in `requirements.txt`) is needed for reproducible 3D coordinates
of 2D ligands.

```bash
# Debian/Ubuntu: Vina (the image builds 1.2.7)
sudo apt-get install -y autodock-vina

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Open Babel

Use **Open Babel 3.2.0 or later**; the image builds 3.2.1, the release
AutoDocker 2.0 is validated with. Open Babel 3.2 reads the `OB_RANDOM_SEED`
variable that AutoDocker sets on every call and no longer adds a hydrogen to
pre-protonated Lys/Arg (see [`openbabel.md`](openbabel.md)). Linux
distributions still ship 3.1.1 (Debian 12/13, Ubuntu 24.04): AutoDocker runs
with it, with a warning, and its results are reproducible but not comparable
with the image's. `run_manifest.json` records `openbabel_version`.

Any of these gives 3.2.1 (check with `obabel -V`):

```bash
# conda-forge
conda install -c conda-forge openbabel=3.2.1

# pip: the upstream wheel (CPython 3.10-3.14) ships the obabel binary
pip install openbabel==3.2.1
# its "obabel" command starts a Python wrapper per call; to use the native
# binary directly (what CI does):
OB="$(python -c 'import os, openbabel; print(os.path.dirname(openbabel.__file__))')"
export PATH="$OB/bin:$PATH" \
       BABEL_LIBDIR="$OB/lib/openbabel/3.2.1" \
       BABEL_DATADIR="$OB/share/openbabel/3.2.1"

# from source: the submodule pinned by this repository, built as the image does
git submodule update --init tools/openbabel/repo
cmake -S tools/openbabel/repo -B build/openbabel -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$HOME/.local/openbabel-3.2.1" \
  -DCMAKE_C_FLAGS=-ffp-contract=off -DCMAKE_CXX_FLAGS=-ffp-contract=off \
  -DBUILD_GUI=OFF -DPYTHON_BINDINGS=OFF -DENABLE_TESTS=OFF -DUSE_RANDOMMT=ON \
  -DWITH_MAEPARSER=OFF -DWITH_COORDGEN=OFF -DWITH_JSON=OFF
cmake --build build/openbabel -j"$(nproc)" && cmake --install build/openbabel
export PATH="$HOME/.local/openbabel-3.2.1/bin:$PATH"
```

The pip wheel and the source build give byte-identical output for the
conversions AutoDocker runs (checked on receptor protonation, ligand `-h`
conversion, MMFF94 minimization and seeded `--gen3d`).

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
| `Open Babel 3.1.1 found; AutoDocker ... needs Open Babel >= 3.2.0` | a distribution's Open Babel; install 3.2.1 (above) to match the image's results |
| `Unable to find OpenBabel plugins` | set `BABEL_LIBDIR` (and `BABEL_DATADIR`) to the install's `lib/openbabel/<version>` and `share/openbabel/<version>` |
| `Critical tool 'vina' not found` | `vina` not on PATH — install `autodock-vina` or use Docker |
| `Protein file not found: /workspace/protein.pdb` | container mount wrong — mount the directory containing `protein.pdb` at `/workspace` |
| `No ligand files (.sdf/.pdbqt/...) found` | put ligands in `works/ligs/` (or set `LIGS`) |
| `Receptor has no charges` | receptor PDB was stripped of HETATM/charge info — pass a raw PDB and let the pipeline prepare it |
| Python under 3.9 | upgrade Python or use Docker |

For a reproducible two-run determinism check see `docs/reproducibility.md`.