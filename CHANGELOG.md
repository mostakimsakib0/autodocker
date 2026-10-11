# Changelog

## 1.1.5 (unreleased)

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

### Added

- `--ph` (receptor protonation pH, default 7.4) and `--grid-center-frame`.
- `scripts/fetch_casf_rcsb.py` rebuilds the CASF-2016 core-set complexes
  from public RCSB files (no PDBbind registration needed; two-letter element
  symbols such as CL/BR are normalised); `scripts/benchmark_redock.py` passes
  the box center in input coordinates.

### Removed

- `autodocker/uidmapper.sh` (unused).
