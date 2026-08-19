---
title: 'AutoDocker: reproducible, containerized virtual screening with AutoDock Vina'
tags:
  - molecular docking
  - virtual screening
  - computational drug discovery
  - AutoDock Vina
  - OpenBabel
  - fpocket
  - reproducibility
authors:
  - name: Mostakim Sakib
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Independent researcher
    index: 1
date: 20 August 2026
bibliography: paper.bib
---

# Summary

AutoDocker is a containerized virtual-screening pipeline built around
AutoDock Vina (Eberhardt et al., 2021, @eberhardt2021vina; Trott & Olson,
2010, @trott2010vina). It automates the complete screening workflow from raw
input to a ranked, reproducible report: protein preparation with chain
selection and fpocket pocket detection (Le Guilloux et al., 2009,
@leguilloux2009fpocket), ligand preparation from local
SDF/MOL2/PDB/PDBQT files or PubChem-downloaded libraries (Kim et al., 2023,
@kim2023pubchem), multi-engine docking (Vina with an optional QuickVina 2
fallback (Alhossary et al., 2015, @alhossary2015qvina2) and SMINA consensus
scoring (Koes et al., 2013, @koes2013smina)), and output generation ranging
from plain-text rankings to a self-contained HTML report with an embedded
NGL (Rose et al., 2018, @rose2018ngl) 3D viewer. Every step is seeded and
checkpointed, so a given input and seed yields identical results on every
run and interrupted screenings resume in place.

# Statement of need

Modern virtual-screening pipelines are frequently assembled from a loose
collection of scripts, tied to GUI tools, or hidden inside opaque web
services, which makes reproducibility — a core requirement for early-stage
drug discovery — difficult to achieve and verify. Fragmented pipelines also
make it hard to produce auditable, comparable results: reviewers and
collaborators cannot easily reconstruct exactly how a ranked hit list was
produced, which receptor was prepared, which pocket was chosen, or which
docking parameters were used. Existing automation efforts tend to be
GUI-bound (e.g., Raccoon2, MzDOCK) or assume particular library sources and
missing quality controls such as charge validation, pocket-based grid
placement, or checkpointing.

AutoDocker addresses these gaps with a fully local, scriptable,
containerized pipeline that is explicit about every decision it makes. The
container image builds all scientific dependencies (Vina, QuickVina 2,
Open Babel (O'Boyle et al., 2011, @oboyle2011openbabel), fpocket) from
source as pinned submodules, removing environment drift. The pipeline
strictly validates receptor and ligand PDBQT files (nonzero charges, atom
records, plateau-size guards), refuses to silently produce garbage, and
ships open-source scientific benchmarks (PDBbind/CASF-2016 re-docking,
DUD-E enrichment, raw-Vina parity) so that changes are continuously checked
against known-good reference results. AutoDocker is designed to be driven
from a single command or a set of environment variables, which makes it
suitable both for interactive use and for large scripted or parallel
screening campaigns.

# Usage

The quickest path is the container image:

```bash
git clone --recurse-submodules https://github.com/mostakimsakib0/autodocker.git
cd autodocker
docker build -t autodocker .
mkdir -p works/ligs
cp protein.pdb works/protein.pdb
cp ligands/*.sdf works/ligs/
docker run --rm -v "$PWD/works:/workspace" autodocker
```

Run parameters are exposed both on the CLI (e.g., `--chains`, `--pockets`,
`--exhaustiveness`, `--seed`, `--consensus`, `--html-report`) and via
matching environment variables (`CHAIN`, `POCKETS`, `EXHAUSTIVENESS`,
`VS_SEED`, and so on). Results are written to an `output/` directory:
`ranking.csv` with per-ligand binding affinities and descriptors,
`Top_hits.txt`, `pocket_summary.csv`, Docker log files, and (optionally) an
interactive HTML report with the receptor, poses, and affinity statistics.

For reproducibility the pipeline uses a fixed default seed (`--seed 42`),
records every docking command in per-ligand log files, and supports
checkpointed resume, so the exact inputs passed to Vina are always
recoverable. Installing the pipeline natively (instead of via Docker)
requires only Python 3.9+, numpy, requests, and the `vina` and `obabel`
binaries on `PATH`; the test suite (`python -m pytest tests/`) verifies the
pipeline with a real end-to-end docking run.

# Validation

AutoDocker is validated with three reproducible benchmark scripts
(`scripts/benchmark_*.py`, see `BENCHMARKS.md`):

1. **Re-docking power.** 40 PDBbind v2016 core-set complexes are re-docked;
   26 dockable complexes give a 76.9% top-1 success rate within 2 Å RMSD of
   the crystal pose and 96.2% when the best of the 9 predicted modes is
   considered, matching or exceeding the Vina baseline.
2. **Screening power.** DUD-E enrichment (Mysinger et al., 2012,
   @mysinger2012dude) across 10 targets spanning protein kinases,
   proteases, a nuclear receptor, reductases and an esterase, with
   bootstrap confidence intervals and Mann-Whitney U significance tests.
3. **No hidden corruption.** All 222 affinity values reported by the
   pipeline match the raw Vina scores bit-for-bit.

# Availability

- Source code: https://github.com/mostakimsakib0/autodocker
- MIT License
- Archived version: https://doi.org/10.5281/zenodo.0000000 (on first release)