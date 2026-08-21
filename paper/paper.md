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
  - name: Md. Mostakim Ahmed Sakib
    orcid: 0009-0001-7176-2016
    affiliation: 1
  - name: Ahmad Hasan Mubashshir
    affiliation: 1
  - name: Md Zahurul Haque
    orcid: 0009-0007-3455-7535
    affiliation: 2
affiliations:
  - name: Independent researcher
    index: 1
  - name: Dept. of CSE, Manarat International University, Dhaka, Bangladesh
    index: 2
date: 20 August 2026
bibliography: paper.bib

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
scoring (Koes et al., 2013, @koes2013smina)), and output ranging from
plain-text rankings to a self-contained HTML report with an embedded NGL
(Rose et al., 2018, @rose2018ngl) 3D viewer. Every step is seeded and
checkpointed, so a given input and seed yields identical results on every
run, and interrupted screenings resume in place.

# Statement of need

Virtual screening is a core step in early-stage drug discovery and
repurposing, but running it *reproducibly* is still unnecessarily hard.
Typical practice assembles a loose chain of heterogeneous
components — a GUI for receptor preparation, a separate tool for pocket
detection, a docking binary invoked by hand, ad-hoc scripts to convert
formats, and more scripts to rank the outputs. Each junction is a place
where parameters (pocket choice, grid box, exhaustiveness, seeds, charge
handling) can be altered, lost, or applied inconsistently, so two otherwise
identical efforts routinely produce different hit lists. The literature
often reports docking scores without the parameters and inputs that
produced them, which makes results difficult to audit or compare across
groups.

The tools that automate screening workflows address this only partially.
GUI wrappers around Vina (e.g., PyRx (Dallakyan & Olson, 2015,
@dallakyan2015pyrx), Raccoon2, MzDOCK) make interactive preparation easier
but are not designed for scripted, repeatable, high-throughput runs and
typically capture little of the decision trail. Complementary engines such
as smina and the QuickVina2 fork (Alhossary et al., 2015, @alhossary2015qvina2)
speed up scoring but leave the surrounding workflow — library curating,
preparation, validation, reporting, resumption — to the user. Script-based
suites such as jamdock-suite (Manso, 2025, @manso2025jamdock) provide a
missing workflow layer but do not include the quality controls, statistics,
or bundled validation that make results citable.

AutoDocker is designed for the gap between a bare docking binary and a GUI:
a fully scriptable, containerized pipeline that makes every preparation and
scoring decision explicit, validates that inputs are chemically usable
before docking (nonzero charges, atom records, size guards), fails loudly
instead of silently producing noise, and ships scientific benchmarks so
that any future change is continuously checked against known-good reference
results. It is targeted at researchers who need auditable, reproducible
screening results at scale — without a GUI, without a proprietary platform,
and with the option to reproduce every result bit-for-bit.

# State of the field

The table below compares AutoDocker with the Vina CLI and representative
existing workflows: PyRx (a widely used GUI wrapper (Dallakyan & Olson,
2015, @dallakyan2015pyrx)), Raccoon2/MzDOCK (GUI automation), the
jamdock-suite script suite (Manso, 2025, @manso2025jamdock), and smina (an
alternative scoring engine (Koes et al., 2013, @koes2013smina)). Rows are
capabilities that affect reproducibility or screening quality.

| Capability | AutoDocker | Vina CLI | PyRx | Raccoon2 / MzDOCK | jamdock-suite | smina |
|---|---|---|---|---|---|---|
| Fully scriptable / headless | ✔ | ✔ (binary only) | ✖ (GUI) | ✖ (GUI) | ✔ | ✔ |
| Pocket detection → grid box (fpocket) | ✔ auto (+centroid fallback) | ✖ manual | ✖ manual box | ✔ | ✔ | ✖ |
| Library source: local files + PubChem/FDA download | ✔ | ✖ | local only | local only | ZINC | ✖ |
| ADMET (Lipinski) filter on library | ✔ | ✖ | ✖ | ✖ | ✖ | ✖ |
| Ligand minimization (MMFF94) | ✔ | ✖ | ✖ | ✖ | ✖ | ✖ |
| Input validation (nonzero charges, atom records, fail-loud) | ✔ | ✖ | partial | partial | ✖ | ✖ |
| Flexible-receptor docking | ✔ (manual + auto residues) | ✔ (manual) | ✖ | ✖ | ✖ | ✖ |
| Consensus scoring (Vina + SMINA) | ✔ optional | ✖ | ✖ | ✖ | ✖ | ✖ (engine only) |
| Checkpoint / resume | ✔ | ✖ | ✖ | ✖ | ✔ (jamresume) | ✖ |
| Seeded determinism + versioned container | ✔ | ✖ | ✖ | ✖ | ✖ | ✖ |
| Interactive HTML report with 3D viewer | ✔ optional | ✖ | ✖ | ✖ | ✖ | ✖ |
| Bundled open benchmarks (re-docking, enrichment, parity) | ✔ | ✖ | ✖ | ✖ | ✖ | ✖ |

AutoDocker's differentiating value is not a new scoring function — it is a
*complete, auditable workflow* on top of Vina: automatic pocket-based grid
placement, chemical validity guards before docking, optional pharmacophore
meta-features (ADMET, minimization, flexible residues, consensus scoring),
single-argument scriptability for HPC, and bundled validation so that
claims are checkable.

# Software design

AutoDocker is two thin layers over well-known scientific binaries. A small
bash entry point (`entry.sh`) maps environment variables to command-line
arguments (both workflows are scored: `INPUT`/`LIGS`/`PDB` and
`-r`/`-l`/`-o`), then a single Python module (`runner.py`) orchestrates the
pipeline in five phases:

1. **Library preparation** — local SDF/MOL2/PDB/PDBQT input or PubChem
   download (FDA-approved and custom libraries), optional Lipinski ADMET
   filtering (Lipinski et al., 2001, @lipinski2001rule), optional MMFF94
   minimization (Halgren, 1996, @halgren1996mmff94), conversion to PDBQT
   with Open Babel (O'Boyle et al., 2011, @oboyle2011openbabel).
2. **Protein preparation** — validate/reject undersized inputs, select
   chains, prepare the receptor PDBQT, detect pockets with fpocket
   (Le Guilloux et al., 2009, @leguilloux2009fpocket; Schmidtke et al.,
   2010, @schmidtke2010fpocket) or fall back to a centroid grid.
3. **Input validation** — before any docking, the receptor and each ligand
   must pass structural guards (atom records present, nonzero partial
   charges, minimum size). Truncated or charge-less inputs are rejected
   instead of being docked silently.
4. **Docking** — Vina (primary) with optional QuickVina 2 fallback, seeded
   for reproducibility, parallelized over ligands with
   `multiprocessing`; optional flexible residues and consensus scoring via
   SMINA.
5. **Result analysis** — ranked CSV, top-hit and metrics text files, pocket
   summary, checkpoint file for resume, and an optional self-contained HTML
   report with affinity statistics, pose clustering, and an embedded NGL
   viewer.

The Docker image is multi-stage: all scientific binaries (Vina, QuickVina
2, Open Babel, fpocket) are compiled from pinned git submodules, which
removes environment drift, and the runtime image is minimal and
non-root. The design goal is that the entire pipeline — from raw receptor
and ligand files to a ranked, citable hit list — is reproducible with a
single command and a fixed seed.

# Validation

AutoDocker ships reproducible benchmark scripts (`scripts/benchmark_*.py`;
results in `BENCHMARKS.md`):

1. **Re-docking power.** 40 PDBbind v2016 core-set complexes are re-docked;
   26 dockable complexes give a **76.9%** top-1 success rate within 2 Å
   RMSD of the crystal pose and **96.2%** when the best of the 9 predicted
   modes is considered, matching or exceeding the Vina baseline.
2. **Screening power.** DUD-E enrichment (Mysinger et al., 2012,
   @mysinger2012dude) on a curated, family-diverse target set (kinases
   akt1/braf, reductase aldr) with seeded sampling (50 actives + 100
   decoys per target), bootstrap 95% confidence intervals and one-sided
   Mann-Whitney U tests. Reported mean AUC 0.74 (akt1 0.76, braf 0.74,
   aldr 0.71), all three actives-vs-decoys comparisons significant
   (p ≤ 1e-5); full per-target results in `BENCHMARKS.md`.
   Benchmark scripts default to the full 10-target curated set; results
   above are from the subset that completed on the legacy DUD-E
   download server ([results][benchmark-results]).
3. **No hidden corruption.** All affinity values reported by the pipeline
   match the raw Vina scores bit-for-bit (mode-by-mode parity).
4. **Performance & scalability.** Wall-clock time, peak RSS, and throughput
   are measured as a function of ligand-library size and process count
   (`scripts/benchmark_performance.py`, see `BENCHMARKS.md`).

# Availability

- Source code: https://github.com/mostakimsakib0/autodocker
- MIT License
- Archived version: https://doi.org/10.5281/zenodo.22035568