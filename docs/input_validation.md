# Input validation and preparation rules

What AutoDocker accepts as receptor and ligand input, what it changes, and what it
refuses before docking. All of it is deterministic: the same input files give
byte-identical prepared files.

## Receptor (`-r`)

**Chains.**
- Chain IDs are case-sensitive (`A` and `a` are different chains).
- ATOM records with a blank chain column (DUD-E style) are read as chain `A`.
  `get_chains`, `--chain`, the log (`Single chain detected: A (blank chain ID ...)`)
  and residue keys all use that name.
- A file with blank-chain ATOM records next to a real chain `A` is refused, because
  the two cannot be told apart. Give the blank-chain records another chain ID.
- `--chain all` selects every chain that has ATOM records. A chain made only of
  HETATM records is never selected. With `--keep-hetero` its records are listed in a
  warning instead of being dropped silently.

**Selenomethionine.** MSE residues are written as methionine: residue name `MET`,
atom `SE` becomes ` SD `, and the element becomes `S`. The coordinates are kept as
deposited. Two things force this:
- AutoDock Vina has no selenium atom type.
- Open Babel has no Gasteiger parameters for Se, so before this change a SeMet
  receptor produced an empty PDBQT.

The conversion is logged as a warning.

**Metal ions with `--keep-hetero`.** Open Babel has no Gasteiger parameters for Zn,
Fe, Mg, Ca, Mn, Cu, Ni, Co, Hg, Cd, Na or K. AutoDocker handles a single-atom metal
ion like this:
- It stays in `protein_clean.pdb`, so it counts towards the frame offset.
- It is left out of the Open Babel input.
- It is appended to `receptor.pdbqt` with its Vina metal type.

A metal inside a larger residue (the Fe of HEM) or any other element without
Gasteiger parameters still makes Open Babel convert nothing. That is an error that
names the element and the residue.

**Open Babel converting nothing is an error.** Open Babel exits with status 0 when
it skips every molecule (`0 molecules converted`). Every Open Babel call therefore
fails, with Open Babel's messages, in either of these cases:
- it reports 0 molecules converted;
- it leaves an expected output file missing or empty.

**Flexible residues (`--flexible-residues`).** A residue key has one of three forms:

| Form | Examples | Use |
|---|---|---|
| `<letter chain><resnum>[icode]` | `A45`, `a45`, `A100B` | Chain is a letter |
| `<chain>:<resnum>[icode]` | `1:45`, `A:45` | Any chain, required for digit chains |
| `<resnum>[icode]` | `45` | The blank chain |

- Chains match exactly, so `A45` never picks up atoms of chain `a`.
- On a blank-chain receptor, `A45` and `45` name the same residue.
- A residue named twice is built once, and both keys count as built.
- Keys written for auto-detected residues use the same syntax. A digit chain gives
  `1:45`.
- A residue with no atoms in the docked chain(s) stops the run before any ligand is
  prepared (exit 2). A residue that cannot be built as a flexible side chain (missing
  atoms, or none rotatable, e.g. GLY, ALA, PRO) stops it at receptor setup (exit 1).

**Kept HETATM residues** (`--keep-waters`, `--detect-metals`, `--detect-cofactors`)
are identified by chain, residue number with insertion code, and residue name.
`HOH 301` and `HOH 301A` are two residues.

## Ligands (`-l`)

**One molecule per file.** The following are rejected with the reason:
- SDF or MOL2 files with several records;
- PDB files with several `MODEL` records;
- PDBQT files with several `ROOT` blocks, `MODEL` records, or no torsion tree.

Every PDBQT that AutoDocker prepares is checked the same way before docking.

**Salts and other disconnected fragments.** Before conversion, AutoDocker keeps the
fragment with the most heavy atoms, as Open Babel separates them (bonds as written,
or perceived for PDB).
- The removed fragments are logged and recorded in `ligand_metadata.json` as
  `fragments_removed`, as canonical SMILES.
- The kept fragment is written to `ligands_minimized/<name>_fragment.<ext>` with the
  input title.
- For SDF input the Lipinski descriptors are computed on the kept fragment.
- If several fragments tie for the largest size, the file is rejected, because there
  is no way to tell which one is the ligand. Examples are `CCO.CCO` and a 1:1
  co-crystal.

**2D input.** An SDF with every z at 0, or a MOL2 or PDB with more than three atoms
and every z at 0, gets 3D coordinates from RDKit ETKDGv3 (seeded with `--seed`) plus
MMFF94. This is the same recipe for all three formats, and the ligand path never
uses Open Babel `--gen3d`.
- The coordinate generation is recorded as `coordinates` and `coordinates_method`.
- A PDB has no bond orders. When the orders Open Babel perceives from a flat drawing
  cannot be read by RDKit, the ligand is rejected and the reason says so.
- A molecule of up to three atoms (water, CO2, an ion) is planar whatever its
  geometry, so planarity is not held against it.

**Earlier outputs inside `-l`.** AutoDocker treats a directory as an output
directory if it holds `ligand_metadata.json` or `run_manifest.json`, or both the
`ligands/` and `ligands_minimized/` subdirectories.
- An output directory found under `-l` is skipped, and the skip is logged.
- An `-l` that is itself an output directory is refused. Pass its `ligands/`
  subdirectory to re-dock the prepared PDBQTs.

**Names and paths.** The name record of a PDBQT prepared from a PDB is the file
name, not its path, so identical files in different folders give identical PDBQTs.
`source_sdf` of PubChem compounds is relative to the output directory.

**Charges.** Hydrogens that Open Babel adds with `-h` get Gasteiger charges. They
previously kept 0.000, for example for a heavy-atom-only SDF with `--no-minimize` or
for a MOL2 without hydrogens.
- The molecule is protonated first and charged in a separate conversion. Only the
  charge column of the PDBQT changes.
- For a MOL2 the charges in the file, which were computed without those hydrogens,
  are replaced by Gasteiger charges of the protonated molecule.

**Small files and unsupported elements.**
- There is no minimum SDF size; a methane molblock is about 90 bytes. Error pages
  and truncated downloads are caught by the parse checks instead (no `M  END`, no
  atom block).
- A ligand with an element that Open Babel cannot charge (B, Pt, Si, metals, ...)
  fails with the element named.

## Hydrogen placement

Open Babel places some hydrogens in a random direction (waters, isolated atoms,
two-atom fragments, some MOL2 amines). Open Babel >= 3.2 draws that direction from
a generator seeded with `OB_RANDOM_SEED`, which AutoDocker sets to `--seed`; 3.1
seeds it from the clock (`usec % 70092`). Either way the direction is a random
draw, not chemistry, so:
- Each conversion runs with the run seed S and is compared with a probe
  conversion with seed S+1; hydrogens that move between the two are re-placed
  from heavy-atom geometry (the result does not depend on S).
- When the two agree, a third conversion (seed S+2) is compared. With 3.2 two
  different seeds expose every random hydrogen already; the third run covers
  Open Babel 3.1, where two runs can draw the same clock seed.

With Open Babel 3.1 a random hydrogen can therefore survive only if all three
runs drew the same clock seed, about 1 in 70092² (about 2e-10) per conversion.
