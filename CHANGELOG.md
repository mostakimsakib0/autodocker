# Changelog

## 2.0.0

A major release: docking scores, docking boxes and several output formats
change compared with 1.1.x (entries marked **(results)**). Rerun any screen
whose numbers you need to compare. Version 1.1.5 was never released; its
changes (below) are part of 2.0.0.

### Fixed

- **(results)** Receptor protonation is deterministic. Open Babel places
  some hydrogens in a random, clock-seeded direction when the neighbouring
  atoms do not fix one; Vina decides which atoms are hydrogen-bond donors
  from hydrogen positions, so re-preparing the same receptor changed scores.
  Hydrogens in the input receptor are now removed before protonation at
  `--ph` (pre-protonated inputs such as DUD-E receptors previously got an
  extra, randomly placed hydrogen on every Lys NZ and Arg NH1), and any
  hydrogen that still differs between two Open Babel runs is re-placed
  deterministically from the heavy atoms. The 1HSG case-study receptor is
  unchanged.
- **(results)** Ligand hydrogens are deterministic: the same two-run check
  is applied to every ligand conversion that adds or keeps hydrogens (some
  MOL2 inputs got randomly placed amine hydrogens).
- **(results)** 2D ligands get 3D coordinates from RDKit ETKDGv3 seeded with
  `--seed` (Open Babel `--gen3d` gave different coordinates on every run; it
  remains a fallback without RDKit, with a warning).
- **(results)** A 2D SDF with `--no-minimize`, or after an MMFF94 failure,
  was docked flat (all z = 0) because the charge repair rebuilt it from the
  2D input; it is now converted from 3D coordinates and a flat result is
  rejected.
- **(results)** The docking box is centered on the midpoint of the pocket's
  bounding box, so the padding is equal on both sides (it was centered on
  the mean of the pocket atoms). `--no-fpocket` boxes and automatically
  chosen flexible residues ignore hydrogens.
- **(results)** Affinities are read from the output PDBQT at full precision;
  the screen output had rounded them (Vina to 4 significant digits, QuickVina
  and SMINA to 0.1 kcal/mol).
- **(results)** Polar hydrogens of PDB ligands get Gasteiger charges (charges
  were computed before hydrogens were added). Vina does not use charges.
- The run fingerprint includes receptor and flexible-residue hydrogens and
  the SHA-256 and version of each docking binary, so a changed receptor or an
  upgraded engine invalidates the checkpoint. Existing checkpoints are
  discarded once.
- QuickVina receives `--flex` in flexible docking; an engine that cannot
  take it fails the ligand instead of docking against a receptor with the
  side chains removed.
- A valid score outside the accepted window is reported as `OUT_OF_WINDOW`
  with its number, ranked after in-window results and counted separately, not
  as a failed docking. A screen in which every result is out of the window no
  longer stops with "all docking attempts failed".
- All-zero partial charges (e.g. benzene after nonpolar hydrogen merging)
  are a warning, not a failure.
- Pocket ties are broken by pocket number instead of fpocket's random Monte
  Carlo volume.
- No silent fallbacks: a missing or failing fpocket, a receptor without
  pockets, PubChem downloads that yield nothing, and unreadable chain lists
  stop the run with a clear message instead of docking in a whole-receptor
  box, docking the local ligand folder, or assuming chain A.
- Receptor cleaning keeps one alternate location per residue (highest
  occupancy) and the first model, also for kept HETATM residues; kept waters
  are no longer added twice with `--keep-hetero` and `--keep-waters`. Kept
  HETATM residues keep P, F, Cl, Br and I atoms; `--detect-metals` keeps
  metal ions only, not a whole heme.
- V3000 SDF files are read correctly (3D coordinates were discarded).
- A Lipinski-rejected `x.sdf` no longer lets `x.mol2` of the same molecule be
  docked unfiltered; same-name files are resolved by precedence and the
  others are reported as superseded.
- An output directory inside the ligand folder is no longer docked as
  ligands on the next run.
- `--mw-min/--mw-max/--logp-max` are rejected outside `--library custom`
  instead of being ignored; in custom mode they replace Lipinski's limits.
- `-v` enables debug logging in every module.
- `--smina-only` rejects options it cannot honour, and SMINA commands are
  logged like Vina's.
- Checkpoints are written atomically, and a lock prevents two runs from
  sharing an output directory.
- The engine settings reach worker processes explicitly, so `--vina-bin` and
  `--timeout` also work with the spawn start method (macOS, Python 3.14).
- Every docking command is logged before it runs, and failed runs keep
  their stderr and a real traceback.
- The HTML report: the table header no longer sorts with the rows, the
  interpretation guide matches the score colours, and the NGL viewer
  download is pinned and hash-checked.
- The PyMOL box script draws the box in the intended colour.

### Changed

- Consensus scoring rescores Vina's top pose with SMINA's Vinardo function
  (it re-docked with SMINA's default function, which is Vina's, so the
  agreement mostly measured search noise). `consensus_ranking.csv` loses the
  "Agreement (heuristic)" column; the Spearman correlation is written to
  `consensus_summary.json`.
- The "ADMET" filter is named for what it is: the Lipinski rule-of-five
  filter (zero violations allowed; Open Babel descriptors). `--no-lipinski`
  replaces `--no-admet` (kept as a deprecated alias), and `LIPINSKI=0`
  replaces `ADMET=0` in the container.
- Custom PubChem libraries use Open Babel's logP for every filter (the
  PubChem XLogP3 pre-filter is removed). Ranitidine (withdrawn in 2020) is
  removed from the built-in drug list; `--max-compounds` sets the number of
  compounds (default 20).
- `ranking.csv`: new `Engine` and `Minimized` columns; `Docked_File` is
  relative to the output directory; SimScore and modes are blank when
  undefined instead of 0.
- `Top_hits.txt`, `metrics.txt` and `Results_full.txt` contain no timestamps;
  `receptor.pdbqt` contains no output path.
- Chain IDs are case-sensitive.
- Open Babel always runs single-threaded.
- Exit codes: 2 for invalid options, environment or missing tools; 3 when
  more ligands fail than `--max-failure-fraction` (default 0.5) allows;
  130 when interrupted.
- The druggability label uses fpocket's documented threshold of 0.5.
- Default settings are defined once (`autodocker/vspipeline/defaults.py`);
  the container passes an option only when its environment variable is set,
  and the benchmark scripts default to the available CPUs instead of 8.

### Added

- `run_manifest.json` in every output directory: version, command line and
  resolved options, seed, platform, the path, version and SHA-256 of every
  tool, the run fingerprint, RDKit version, timestamps and failure counts.
- Container images for linux/amd64 and linux/arm64 (native builds in one
  manifest). Docking engines and fpocket are compiled with
  `-ffp-contract=off`, which was verified to leave amd64 binaries and results
  unchanged. On arm64, fpocket reads PDB input only, which is all AutoDocker
  uses.
- Reproducible image builds: base images pinned by digest, Debian packages
  from a dated snapshot (Open Babel 3.1.1 asserted), NGL from a lockfile,
  exact Python pins (now including RDKit), and `autodocker/VERSION` as the
  single version source.
- `--max-compounds`, `--max-failure-fraction`, `--no-lipinski`.
- Benchmarks: enrichment boxes from the DUD-E crystal ligand with all chains;
  BEDROC and ROC enrichment; re-docking from seeded RDKit start conformers on
  a seeded sample of the core set; per-mode RMSDs; ligand identity checks in
  `fetch_casf_rcsb.py`; a parity check that `ranking.csv` equals raw Vina;
  repeated, per-phase throughput measurements.
- Tests for determinism of receptor and ligand preparation, and end-to-end
  runs that must be byte-identical across repeats, `-p` and Vina `--cpu`.

### Fixed in the benchmark scripts

- Enrichment factors were computed against the nominal fraction instead of
  the screened fraction and could exceed their maximum; the bootstrap
  collapsed duplicate draws; failed ligands were dropped from the count
  instead of ranked last; `--ci` was ignored; p-values were rounded to zero.
- Re-docking took mode 2 as "top-1" when the mode-1 RMSD failed, reported
  success over docked instead of attempted complexes, and resumed with
  changed parameters.
- The parity check passed when every complex errored.

### CI

- One workflow: Python 3.9 to 3.12 with the official Vina 1.2.7, ruff,
  shellcheck and hadolint; on pull requests the image is built for both
  architectures and a smoke screen must reproduce byte for byte.

## 1.1.5 (never released; superseded by 2.0.0)

Changes marked **(results)** alter docking scores compared with 1.1.4; rerun
screens that need comparable numbers.

### Fixed

- **(results)** The receptor is now protonated (Open Babel `-p`, pH 7.4 by
  default, `--ph`). Earlier versions wrote a receptor without hydrogens, so
  every nitrogen was typed as an acceptor and the protein had no hydrogen-bond
  donors. Preparation fails if a receptor ends up without polar hydrogens.
- **(results)** PDB and MOL2 ligands, and SDF ligands with `--no-minimize`,
  get hydrogens (`-h`); the charge-repair step no longer rebuilds a ligand
  without them (a heavy-atom PDB ligand such as a crystal ligand was docked
  without donors). The repair step also reads MOL2 input as MOL2.
- **(results)** A 2D SDF prepared with `--no-minimize` gets 3D coordinates
  instead of being docked flat.
- Checkpoints record a fingerprint of the receptor, box, docking parameters
  and engine, and a hash of each ligand file; changed inputs are docked again
  and ligands no longer in the library are not reported. Results are saved as
  each ligand finishes, so an interrupted screen keeps its progress.
- `--grid-center` is read in the input PDB's coordinates and converted to the
  centered docking frame (`--grid-center-frame docking` keeps the old
  behaviour).
- Flexible residues: only rotatable side-chain bonds become branches (rings,
  amides, carboxylates and guanidinium groups stay rigid), the backbone stays
  in the rigid receptor, GLY/ALA/PRO are never selected, only residues that
  were actually built are removed from the rigid receptor, auto-selection
  ranks residues by their closest atom, and insertion codes are respected.
  `--flexible-residues`/`--auto-flexible` enable flexible docking on their own.
- Receptor cleaning keeps the first model of multi-model files and the first
  alternate location of each atom, and keeps common modified residues (MSE,
  SEP, TPO, PTR, ...). Kept HETATM records with 4-digit residue numbers are
  written in the correct columns; metals are no longer guessed from atom
  names such as C-alpha.
- A score outside the accepted window is reported as `REJECTED (<score>)`,
  no longer triggers a re-dock with the fallback engine, and is cached so a
  resumed screen does not dock the ligand again.
- `--threads` is passed to Vina, QuickVina and SMINA as `--cpu`; with `-p N`
  each process gets CPUs ÷ N by default, counting the container's CPU quota.
- `--smina-only` docks with SMINA alone and reports SMINA's poses and
  failures; consensus mode lists ligands one engine could not score.
- Ties in `ranking.csv` are ordered by ligand name, so parallel runs give
  identical files.
- HTML report: the docking box is drawn, the viewer shows only the best pose
  (from the file recorded in `ranking.csv`, so SMINA poses appear in
  `--smina-only` mode), and pose clustering works on Vina and SMINA output,
  covers only the ligands of the current screen, ignores flexible side-chain
  atoms and honours `--rmsd-threshold`.
- Multi-molecule SDF/MOL2 files and multi-model PDBQT files are rejected with
  a clear reason; a failed MMFF94 setup falls back to the input structure.
- PubChem modes stop after repeated network failures with a clear message,
  do not retry 404s, and use the 2D record when no 3D conformer exists.
- Container: the exit status reflects the run; `-r` is passed to the
  pipeline (the root debug shell is now `--root-shell`); proxy variables are
  forwarded; `/autodocker` is no longer writable by the unprivileged user;
  under rootless podman the pipeline runs as the (unprivileged) container
  root with all capabilities dropped instead of `chown -R`-ing `/workspace`
  to `dockuser` and back, so a killed run no longer leaves files the user
  cannot modify;
  `entry.sh` honours `-r/-l/-o`, detects extensions case-insensitively and
  stops when `INPUT` points to a missing file.
- The Vina submodule is pinned to the v1.2.7 release tag (same C++ source as
  before).

- SDF files whose molblock contains text such as "404" or "Status:" are no
  longer rejected as PubChem error pages; when Open Babel cannot compute
  descriptors, ADMET properties are read from SDF data items (including
  PubChem tags).
- Checkpoints written by earlier versions (no run fingerprint) are discarded,
  so the first 1.1.5 run on an old output directory docks every ligand again.
- `scripts/patch.sh` stops the image build when a patch no longer applies.

- The run fingerprint ignores receptor hydrogen coordinates. *(Withdrawn in
  2.0.0: the premise that Vina ignores hydrogen positions is wrong. Vina
  decides which atoms are hydrogen-bond donors from hydrogen positions, so
  the random placements changed scores. 2.0.0 makes the placement
  deterministic and fingerprints the full receptor file.)*

### Added

- `--ph` (receptor protonation pH, default 7.4) and `--grid-center-frame`.
- `scripts/fetch_casf_rcsb.py` rebuilds the CASF-2016 core-set complexes
  from public RCSB files (no PDBbind registration needed; two-letter element
  symbols such as CL/BR are normalised); `scripts/benchmark_redock.py` passes
  the box center in input coordinates.

### Removed

- `autodocker/uidmapper.sh` (unused).
