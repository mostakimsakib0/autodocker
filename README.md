# autodocker — Virtual Screening Pipeline

![License](https://img.shields.io/github/license/mostakim_sakib0/autodocker)
![CI](https://github.com/mostakim_sakib0/autodocker/actions/workflows/ci.yml/badge.svg)
![Version](https://img.shields.io/badge/version-1.0.0-blue)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.0000000.svg)](https://doi.org/10.5281/zenodo.0000000)
[![JOSS paper](paper/paper.md)](paper/paper.md)
[![Contributing](CONTRIBUTING.md)](CONTRIBUTING.md)

AutoDocker is a containerized virtual-screening pipeline. It prepares proteins
(chain selection, fpocket pocket detection, receptor PDBQT), prepares ligands
(local SDF/MOL2/PDB/PDBQT or PubChem downloads), docks with AutoDock Vina
(QuickVina fallback), and produces a reproducible ranked report
(`ranking.csv`, `Top_hits.txt`, `metrics.txt`, optional HTML).

---

## Quick start

Build the image (submodule-based; builds Vina, QVina2, OpenBabel, fpocket):

```bash
git clone --recurse-submodules https://github.com/mostakimsakib0/autodocker.git
cd autodocker
docker build -t autodocker .
```

Run a screening on a protein + ligand directory:

```bash
mkdir -p works/ligs
cp protein.pdb works/protein.pdb
cp ligands/*.sdf works/ligs/

docker run --rm \
  -v "$PWD/works:/workspace" \
  autodocker
```

The single-ligand mode downloads a ligand via `INPUT`:

```bash
docker run --rm \
  -v "$PWD/works:/workspace" \
  -e INPUT=/workspace/ligs/single.sdf \
  autodocker
```

Results are written under `/workspace/output/`:
- `ranking.csv` — sorted docking scores with descriptors
- `Top_hits.txt` — the top-N hits with interpretation
- `metrics.txt` — pose statistics and SimScore
- `pocket_summary.csv` — fpocket pocket scores
- `grid_box.py` — PyMOL script visualizing the docking box
- `docked/*_out.pdbqt` — docked poses; `*_vina.log` — per-ligand logs

## CLI reference

Direct use (inside the image or with tools on `PATH`):

| Option | Description | Default |
|---|---|---|
| `-r, --receptor` | Protein PDB file (**required**) | |
| `-l, --ligands` | SDF file, PDBQT/MOL2/PDB file, or directory (**required**) | |
| `-o, --output` | Output directory (**required**) | |
| `--library` | `local` · `fda` · `custom` (PubChem source) | `local` |
| `--no-admet` | Skip Lipinski ADMET filtering | off |
| `--no-minimize` | Skip MMFF94 ligand minimization | off |
| `--chain` | Chains to keep, e.g. `A`, `A,B`, or `all` | auto |
| `--no-fpocket` | Skip fpocket; use centroid grid | off |
| `--pockets` | Pocket numbers, e.g. `1` or `2,3` | best |
| `--padding` | Grid padding around pocket (Å) | `6.0` |
| `-p, --processes` | Parallel workers (`-1` = all cores) | `1` |
| `--exhaustiveness` | Vina exhaustiveness | `8` |
| `--binding-modes` | Vina binding modes | `9` |
| `--seed` | Random seed (reproducibility) | `42` |
| `--no-resume` | Ignore checkpoint and start fresh | off |
| `--top-n` | Hits in `Top_hits.txt` | `20` |
| `--html-report` | Also generate an HTML report | off |
| `--vina-bin` | Custom docking binary | auto |
| `--vina-extra` | Extra args passed to the docking binary | none |
| `--keep-waters` | Keep waters near binding site | off |
| `--detect-metals` | Detect & keep metal ions | off |
| `--detect-cofactors` | Detect & keep cofactors | off |
| `--flexible-residues` | Flexible residues, e.g. `A45,A78` | none |
| `--auto-flexible N` | Auto-select N flexible residues | none |
| `--consensus` | Vina + SMINA consensus scoring (needs SMINA on PATH) | off |
| `--smina-only` | Rank with SMINA only | off |
| `--version` | Print version and exit | — |

## Environment variables

| Variable | Effect | Default |
|---|---|---|
| `INPUT` | Single-file ligand input (entry.sh) | unset |
| `LIGS` | Ligand directory (entry.sh) | `/workspace/ligs` |
| `PDB` | Protein path (entry.sh) | `/workspace/protein.pdb` |
| `CHAIN` | Chain selection | `A` |
| `PROCESSES` / `MAX_WORKERS` | Parallel workers | `-1` |
| `EXHAUSTIVENESS` | Vina exhaustiveness | `8` |
| `BINDING_MODES` | Vina binding modes | `9` |
| `ENERGY_RANGE` | Vina energy range | `3.0` |
| `TOP_N` | Top hits count | `20` |
| `TIMEOUT` | Per-command timeout (s) | `900` |
| `ADMET` | Set `0` to disable ADMET filtering | `1` |
| `RESUME` | Set `0` to disable checkpoint resume | `1` |
| `VS_SEED` | RNG seed | `42` |
| `PUBCHEM_REST` | PubChem REST base (mirrors) | official |

## Reproducibility

Docking is seeded (`--seed 42`); the same input + seed produces identical
rankings (verified by the golden determinism test in `tests/test_e2e.py`
and the standalone `scripts/reproducibility_check.sh` checker).

## Validation

Scientific results are in `BENCHMARKS.md`: a CASF-2016 re-docking benchmark
(76.9% top-1 pose success ≤2 Å), a DUD-E enrichment test (AUC 0.59–0.85
across akt1/braf/aldr), and a raw-Vina parity check (222/222 modes score
identically). All three are reproducible from `scripts/benchmark_*.py`.

## More docs

- [`docs/installation.md`](docs/installation.md) — Docker and native setup, troubleshooting
- [`docs/reproducibility.md`](docs/reproducibility.md) — determinism guarantees and how to verify
- [`docs/user_testing.md`](docs/user_testing.md) — fresh-machine testing template and bug report form
- [`BENCHMARKS.md`](BENCHMARKS.md) — validation, performance, and enrichment results

## Development

```bash
pip install -r requirements.txt pytest
python -m pytest tests/ -q
```

`tests/` covers validation, grid fallbacks, ADMET gating, checkpoint resume,
report generation/HMTL escaping, and an end-to-end smoke dock.

## Citation

If you use AutoDocker, please cite it. A JOSS-style manuscript is in
`paper/paper.md` (with `paper/paper.bib`); machine-readable metadata is in
`CITATION.cff` and `.zenodo.json`.

## License

MIT — see `LICENSE`.