# Autodocker — Virtual Screening Pipeline

![License](https://img.shields.io/github/license/mostakimsakib0/autodocker)
![CI](https://github.com/mostakimsakib0/autodocker/actions/workflows/ci.yml/badge.svg)
![Version](https://img.shields.io/badge/version-2.0.0-blue)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22035568.svg)](https://doi.org/10.5281/zenodo.22035568)

AutoDocker is a containerized virtual-screening pipeline. It prepares proteins
(chain selection, fpocket pocket detection, receptor PDBQT), prepares ligands
(local SDF/MOL2/PDB/PDBQT or PubChem downloads), docks with AutoDock Vina
(QuickVina fallback), and produces a reproducible ranked report
(`ranking.csv`, `Top_hits.txt`, `metrics.txt`, optional HTML).

---

## Quick start

Pull the published image (linux/amd64 and linux/arm64; Docker picks the
right one):

```bash
docker pull ghcr.io/mostakimsakib0/autodocker:v2.0.0
docker tag ghcr.io/mostakimsakib0/autodocker:v2.0.0 autodocker
```

or build it from source (Vina, QuickVina 2, SMINA and fpocket are compiled
from pinned submodules; Open Babel comes from a pinned Debian snapshot):

```bash
git clone --recurse-submodules https://github.com/mostakimsakib0/autodocker.git
cd autodocker
docker build -t autodocker \
  --build-arg AUTODOCKER_VERSION="$(cat autodocker/VERSION)" \
  --build-arg VCS_REF="$(git rev-parse HEAD)" .
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

To screen a PubChem library instead of local files, leave `works/ligs/`
empty and pass the library option. `fda` is a fixed list of 74 well-known
approved drugs (not the full FDA list); `custom` is the same list filtered by
`--mw-min/--mw-max/--logp-max`; `--max-compounds` (default 20) limits how many
are prepared. Each compound's PubChem CID, download URL and retrieval time are
recorded in `ligand_metadata.json`.

```bash
docker run --rm -v "$PWD/works:/workspace" autodocker --library fda
```

The single-ligand mode docks one file given via `INPUT`:

```bash
docker run --rm \
  -v "$PWD/works:/workspace" \
  -e INPUT=/workspace/ligs/single.sdf \
  autodocker
```

Results are written under `/workspace/output/`:
- `ranking.csv` — all ligands ranked by score: in-window results first, then
  `OUT_OF_WINDOW` results (a valid Vina score outside the accepted window,
  shown with its number), then failures. Columns include the descriptors,
  `Status`, `Engine` (which docking binary produced the score) and
  `Minimized`; `Docked_File` is relative to the output directory.
- `Top_hits.txt`, `metrics.txt`, `Results_full.txt` — summaries (no timestamps,
  so repeated runs give identical files)
- `run_manifest.json` — everything needed to audit or repeat the run: version,
  command line and resolved options, seed, tool paths, versions and SHA-256,
  the docking fingerprint, the threads given to each engine process
  (`threads_per_process`), platform, timestamps and failure counts
  (`failures`: ligands that could not be prepared and ligands whose docking
  failed). It is replaced atomically, so it is never seen half-written
- `ligand_metadata.json` — per-ligand source, descriptors, preparation details
- `pocket_summary.csv` — fpocket pocket scores (the `Volume` column is
  fpocket's clock-seeded Monte Carlo estimate and varies between runs; it is
  never used to choose a pocket)
- `grid.conf`, `grid_box.pml` — the docking box and a PyMOL script that draws it
- `docked/*_out.pdbqt` — docked poses; `*_vina.log` — per-ligand logs whose
  first line is the exact command
- `consensus_ranking.csv`, `consensus_summary.json` — with `--consensus`.
  Ligands without both scores follow the ranked ones with `Consensus` `NA`
  (a missing score is `FAILED`); a ligand whose Vina score is
  `OUT_OF_WINDOW` is not rescored and keeps its Vina score, with
  `SMINA_Affinity` `NA` and `Consensus` `OUT_OF_WINDOW`
- `pose_clusters.json` — with `--cluster-poses`
- `frame_offset.txt` — translation back to the input PDB's coordinates (see below)

Rerunning into the same output directory first removes the derived outputs
of the earlier run (the files above, `grid.conf`, `flex.pdbqt`/`rigid.pdbqt`
and fpocket's output; listed in `autodocker/vspipeline/outputs.py`), so
every file present belongs to the latest run; for example
`results_report.html` is gone after a `--no-html-report` rerun, and
`pocket_summary.csv` after a run with an explicit box. Docked poses and the
checkpoint (`docked/`, `.docking_checkpoint*.json`) are kept for resuming;
the checkpoint stores file paths relative to the output directory, so a run
can be resumed from any working directory and after the output directory
was moved.

Only one run at a time can use an output directory: a run takes the lock
`.docking.lock` before it writes anything there (in every mode, including
`--smina-only`) and holds it to the end. A second run into the same
directory stops at once with exit code 2 and changes nothing in it.

Hydrogens present in the input receptor are removed and the receptor is
protonated at pH 7.4 (`--ph`); it is centered at the origin
during preparation, so `receptor.pdbqt`, `grid.conf` and the docked poses
share that centered frame (`--grid-center` is given in the input PDB's
coordinates and converted automatically).
The bundled PyMOL script and HTML report use it consistently, but poses will
not overlay on your original `protein.pdb` (e.g. for MD set-up or comparison
with a crystal ligand) until shifted back:

```bash
python3 scripts/restore_frame.py works/output works/output/docked/*_out.pdbqt
```

## CLI reference

Direct use (inside the image or with tools on `PATH`):

| Option | Description | Default |
|---|---|---|
| `-r, --receptor` | Protein PDB file (**required**) | |
| `-l, --ligands` | Ligand file or directory of SDF/MOL2/PDB/PDBQT files (searched recursively). Same-name files resolve SDF > PDBQT > MOL2 > PDB; the others are reported as superseded | required for `local` |
| `-o, --output` | Output directory (**required**); skipped if it lies inside `-l` | |
| `--library` | `local` · `fda` · `custom` (PubChem) | `local` |
| `--mw-min/--mw-max/--logp-max` | `custom` only: replace Lipinski's MW and logP limits (Open Babel logP) | `200/500/5` |
| `--max-compounds` | PubChem compounds to prepare | `20` |
| `--no-lipinski` | Skip the Lipinski filter (SDF/PubChem inputs; zero violations allowed). `--no-admet` is a deprecated alias | filter on |
| `--no-minimize` | Skip MMFF94 minimization (SDF/PubChem ligands only) | off |
| `--chain` | Chains to keep, e.g. `A`, `A,B`, or `all` (case-sensitive; checked against the receptor before any ligand is prepared) | first chain, with a warning |
| `--ph` | pH for receptor hydrogens (Open Babel `-p`), 0 to 14 | `7.4` |
| `--pockets` | fpocket pocket numbers, e.g. `1` or `2,3` (not with `--no-fpocket`, or with `--grid-center` and `--grid-size` together) | best score |
| `--padding` | Box padding around the pocket (Å) | `6.0` |
| `--no-fpocket` | Skip fpocket; box around the whole receptor | off |
| `--grid-center X,Y,Z` | Explicit box center, input-PDB coordinates; with `--grid-size` too, fpocket is not run (nor required) | auto |
| `--grid-center-frame` | `input` or `docking` (≤1.1.4 used `docking`); needs `--grid-center` | `input` |
| `--grid-size SX,SY,SZ` | Explicit box size (Å) | auto |
| `--keep-waters`, `--detect-metals`, `--detect-cofactors` | Keep waters near the site, metal ions, known cofactors | off |
| `--flexible-residues` | Flexible side chains, e.g. `A45,A78`; each must exist in the docked chains (else exit 2) and be buildable (else exit 1) | none |
| `--auto-flexible N` | The N residues closest to the box center become flexible; if none can be built, docking is rigid (warning). Not with `--flexible-residues` | none |
| `-p, --processes` | Parallel docking processes (`0`/`-1` = all CPUs, `-N` = CPUs ÷ N) | `1` |
| `--threads` | CPU threads per docking process (`--cpu`) | CPUs ÷ processes |
| `--exhaustiveness` | Vina exhaustiveness (≥ 1) | `8` |
| `--binding-modes`, `--energy-range` | Vina output modes and energy range | `9`, `3.0` |
| `--seed` | Random seed, 1 to 2147483647 | `$VS_SEED` or `42` |
| `--min-valid-affinity` | Scores ≥ this (or < −20) are reported as `OUT_OF_WINDOW`; must be > −20 | `-1.0` |
| `--vina-bin`, `--vina-extra` | Custom docking binary; extra engine arguments | auto, none |
| `--qvina-fallback` | Re-dock a ligand with QuickVina 2 when Vina fails on it (the `Engine` column shows which engine scored each ligand). Turning it on or off keeps the checkpointed Vina results; only results of the fallback are docked again | off |
| `--timeout` | Per-command timeout (s) | `$VS_COMMAND_TIMEOUT` or `900` |
| `--max-failure-fraction` | Exit 3 if more than this fraction of ligands fail: (ligands that could not be prepared + ligands whose docking failed) / all ligands (superseded files are not counted; `OUT_OF_WINDOW` is not a failure). Also when all fail; `1.0` never fails | `0.5` |
| `--no-resume` | Ignore the checkpoint (it is also discarded when receptor, box, parameters or engine binaries change) | off |
| `--consensus` | Rescore Vina's top pose with SMINA's Vinardo function; writes `consensus_ranking.csv` and `consensus_summary.json` (needs SMINA) | off |
| `--smina-only` | Dock with SMINA instead of Vina (1 mode, no resume) | off |
| `--top-n` | Top-ranked compounds highlighted in `Top_hits.txt` and the HTML report (`ranking.csv` lists all) | `20` |
| `--no-html-report` | Skip the self-contained HTML report | report on |
| `-v, --verbose` | Debug logging for all modules | off |
| `--version` | Print version and exit | — |

`python3 autodocker/runner.py --help` lists every option with its default.
Long options must be written in full (`--ligands`, not `--lig`); numeric
options must be finite numbers. Settings are logged in `run_manifest.json`
(strict JSON) with the flexible residues actually docked.

## Environment variables

| Variable | Effect | Default |
|---|---|---|
| `INPUT` | Single-file ligand input (entry.sh); ignored, with a warning, by `--library fda/custom` | unset |
| `LIGS` | Ligand directory (entry.sh) | `/workspace/ligs` |
| `PDB` | Protein path (entry.sh) | `/workspace/protein.pdb` |
| `CHAIN` | Chain selection (`A`, `A,B`, `all`) | first chain |
| `PROCESSES` / `MAX_WORKERS` | Parallel workers | `-1` |
| `EXHAUSTIVENESS` | Vina exhaustiveness | `8` |
| `BINDING_MODES` | Vina binding modes | `9` |
| `ENERGY_RANGE` | Vina energy range | `3.0` |
| `TOP_N` | Compounds highlighted in `Top_hits.txt` and the HTML report | `20` |
| `TIMEOUT` | Per-command timeout (s) | `VS_COMMAND_TIMEOUT` or `900` |
| `LIPINSKI` | Set `0` to disable the Lipinski filter (`ADMET=0` is a deprecated alias) | `1` |
| `NO_MINIMIZE` | Set `1` to skip MMFF94 minimization | `0` |
| `PADDING` | Box padding (Å) | `6` |
| `POCKETS` | fpocket pocket number(s) | best |
| `RESUME` | Set `0` to disable checkpoint resume | `1` |
| `OMP_NUM_THREADS` | Recorded in `run_manifest.json` (Open Babel always runs single-threaded) | unset |
| `VS_SEED` | RNG seed (an explicit `--seed` takes precedence) | `42` |
| `PUBCHEM_REST` | PubChem REST base (mirrors) | official |

## Exit status

| Code | Meaning |
|---|---|
| 0 | Screen completed |
| 1 | Pipeline error |
| 2 | Invalid options, input files or environment, a required tool is missing (also from the container's `entry.sh`, e.g. a missing protein or ligand folder), or the output directory is in use by another run (nothing is written then, not even `run_manifest.json`) |
| 3 | More ligands failed than `--max-failure-fraction` allows, including all of them (all outputs are still written) |
| 130 | Interrupted (Ctrl-C) |
| 143 | Terminated by SIGTERM (e.g. `docker stop`) |

On Ctrl-C or SIGTERM the run stops its docking workers and their engine
processes, keeps the checkpoint (rerun the same command to resume), writes
`run_manifest.json` with status `interrupted` and releases the output
directory. A docking worker process that dies (e.g. killed by the
out-of-memory killer) does not stop the screen: its unfinished ligands are
docked again, each in a fresh worker, and a ligand whose worker dies again is
reported as `FAILED`.

For a debug shell as root, run the image with `--root-shell` as the first
argument. `--help` and `--version` work without a workspace.

## Reproducibility

Every random step is seeded or made deterministic: Vina runs with `--seed`;
hydrogens that Open Babel would place in random directions (receptor and
ligand) are re-placed deterministically; 2D ligands get 3D coordinates from
RDKit ETKDG seeded with `--seed`; Open Babel runs single-threaded. Repeating a
run with the same inputs, options and image gives byte-identical
`ranking.csv`, summaries and poses on the same CPU architecture, for any `-p`
or `--threads`. Whether amd64 and arm64 results are identical is measured in
CI and not guaranteed. `run_manifest.json` records what is needed to repeat a
run. See [`docs/reproducibility.md`](docs/reproducibility.md) and
`scripts/reproducibility_check.sh`.

## Validation

`BENCHMARKS.md` describes the bundled benchmarks: CASF-2016 re-docking from
seeded RDKit start conformers, DUD-E enrichment (ROC AUC, EF, BEDROC with
bootstrap confidence intervals) with boxes from the crystal ligands, a
raw-Vina parity check, and throughput scaling. Results for 2.0.0 are pending;
every benchmark is reproducible from `scripts/benchmark_*.py`.

## More docs

- [`docs/installation.md`](docs/installation.md) — Docker and native setup, troubleshooting
- [`docs/reproducibility.md`](docs/reproducibility.md) — determinism guarantees and how to verify
- [`docs/user_testing.md`](docs/user_testing.md) — fresh-machine testing template and bug report form
- [JOSS-style paper](paper/paper.md) — manuscript and references
- [Contributing guide](CONTRIBUTING.md)
- [`BENCHMARKS.md`](BENCHMARKS.md) — validation, performance, and enrichment results

## Development

```bash
pip install -r requirements.txt pytest
python -m pytest tests/ -q
```

`tests/` covers input validation, receptor and ligand preparation
(including repeated-preparation determinism), pocket and box selection, the
Lipinski filter, checkpoint resume, report generation and HTML escaping, the
benchmark statistics, and end-to-end docking runs that must be byte-identical
across repeats and process counts. Tests marked `slow` use the real
Open Babel/Vina binaries for longer; deselect them with `-m "not slow"`.

## Citation

If you use AutoDocker, please cite it. A JOSS-style manuscript is in
`paper/paper.md` (with `paper/paper.bib`); machine-readable metadata is in
`CITATION.cff` and `.zenodo.json`.

## License

MIT — see `LICENSE`.
