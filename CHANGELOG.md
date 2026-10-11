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
- The container forwards `LIPINSKI`, `OMP_NUM_THREADS` and
  `AUTODOCKER_VERSION` to the pipeline (`LIPINSKI=0` was dropped by
  `setpriv --reset-env` and had no effect); a test checks the forwarded list
  against `entry.sh` and the README.
- `--flexible-residues` naming a residue that is not in the docked chains
  stops the run before ligand preparation (exit 2), and one that cannot be
  built as a flexible side chain stops it with exit 1; the run no longer
  docks rigidly while reporting flexible residues. `--auto-flexible` that
  finds no buildable residue docks rigidly with a warning.
  `run_manifest.json` records the flexible residues actually docked
  (`flexible_residues`, `[]` for rigid docking).
- fpocket is neither run nor required when `--grid-center` and `--grid-size`
  give the whole box; when pockets are detected, a missing fpocket is
  reported before any work (exit 2).
- Invalid values are rejected before any work (exit 2): non-finite numbers
  (`nan`, `inf`) for every numeric option, `--seed`/`VS_SEED` outside 1 to
  2147483647 (the engines take a 32-bit seed, Vina treats 0 as random and
  RDKit rejects negative seeds), `--ph` outside 0-14, `--min-valid-affinity`
  at or below the -20 kcal/mol floor (empty window), malformed `--pockets`,
  `--chain` IDs absent from the receptor (the message lists the available
  chains), and an `-o` that is a file, empty or cannot be created.
  `run_manifest.json` is strict JSON.
- Options that were accepted but had no effect are rejected (exit 2):
  `--auto-flexible` with `--flexible-residues`, `--pockets` with
  `--no-fpocket` or with both `--grid-center` and `--grid-size`, an empty
  `--flexible-residues` list, `--grid-center-frame` without `--grid-center`.
- Outputs of an earlier run in the same output directory (HTML report, pose
  clusters, consensus files, pocket summary, flexible-receptor files, ...;
  `outputs.DERIVED_OUTPUTS`) are removed when a run starts, so the summary
  and the container's closing list show only this run's files. Docked poses
  and the checkpoint are kept. The HTML report also removes an old report
  when it has nothing to write.
- A screen in which every ligand fails writes all outputs (ranking.csv with
  `FAILED` rows, summaries, report, manifest with failure counts) and exits
  with 3 like any other excess of failures (0 with
  `--max-failure-fraction 1.0`), instead of stopping with exit 1.
- `entry.sh`: invalid input (missing `INPUT` file, protein, ligand folder or
  ligand files, an output directory that cannot be created) exits with 2;
  `--help` and `--version` work without a workspace; `-lPATH`, `-l=PATH`,
  `-oPATH`, `-rPATH` and `--opt=value` are recognised; `INPUT` is ignored
  with a warning by `--library fda/custom` (it made the container fail).
- Messages: an interrupt during `--smina-only` docking is reported as such;
  a conflict caused by `VS_ENABLE_QVINA_FALLBACK=1` names the variable; the
  HTML report shows "n/a" instead of 0.00 for the best and average affinity
  when no ligand has an in-window score, drops the dead 3D-viewer link when
  there is no viewer, and clicking the scatter point of a ligand outside the
  top hits says so instead of reloading the current hit.
- One run per output directory: the lock is taken before anything is
  written (it was taken only when docking started, after ligand and receptor
  preparation and the manifest, and never with `--smina-only`, so a second
  run could rewrite the grid, receptor and ligands of a running screen).
  A run that finds the directory in use exits with 2 and writes nothing.
- The checkpoint stores per-ligand file paths relative to the output
  directory: a run resumed from another working directory, or after the
  output directory was moved, wrote a wrong `Docked_File` and `--consensus`
  crashed. A missing Vina pose now fails only that ligand's rescoring.
- A docking or rescoring worker process that dies (e.g. out-of-memory
  killer) no longer hangs the screen: its unfinished ligands are run again,
  each in a fresh worker, and a ligand whose worker dies again is `FAILED`.
- SIGTERM (`docker stop`) is handled like Ctrl-C: workers and their engine
  processes are stopped, the checkpoint is kept, `run_manifest.json` records
  status `interrupted`, and the exit code is 143. Before, the run died
  without a manifest while its workers kept docking into the unlocked
  directory.
- `--max-failure-fraction` counts ligands that could not be prepared, as
  documented ("if more than this fraction of ligands fail"); superseded
  files are not counted. `run_manifest.json` records both counts
  (`failures`).
- `consensus_ranking.csv` labelled `OUT_OF_WINDOW` ligands `FAILED`; they
  keep their Vina score, with `SMINA_Affinity` `NA` and `Consensus`
  `OUT_OF_WINDOW`.
- A non-finite SMINA score (`nan`) fails the ligand in `--smina-only` and
  `--consensus` instead of being reported as a score.
- `run_manifest.json` is written atomically and records the threads given
  to each engine process (`threads_per_process`).
- Turning `--qvina-fallback` on or off no longer discards the checkpoint:
  the fingerprint covers the primary engine only, and results produced by
  the fallback (tagged with its SHA-256) are docked again individually when
  the fallback changes. Checkpoints written by earlier 2.0.0 development
  builds are discarded once.
- Docking results are collected in library order, independent of the order
  in which the workers finish.

### Changed

- **CLI:** abbreviated long options are no longer accepted (`--lig`,
  `--out`, `--mw-ma`, ...; write `--ligands`, `--output`, `--mw-max`). The
  container entry point has to recognise `-r/-l/-o/--library` exactly as the
  runner does, and abbreviations made the two disagree.

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
- QuickVina fallback is the visible option `--qvina-fallback` (the hidden
  `VS_ENABLE_QVINA_FALLBACK=1` still works, with a deprecation warning).
  The undocumented `VS_EXHAUSTIVENESS`, `VS_BINDING_MODES`,
  `VS_ENERGY_RANGE`, `VS_MIN_VALID_AFFINITY` and `VS_THREADS` fallbacks are
  removed; `VS_SEED` and `VS_COMMAND_TIMEOUT` remain.
- The HTML report highlights `--top-n` ligands (default 20) in its top-hits
  table and 3D viewer, like `Top_hits.txt`; it always showed 10.
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
