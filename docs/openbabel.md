# Open Babel

AutoDocker prepares receptors and ligands with Open Babel (`obabel`): protonation
at `--ph` (`-p`), hydrogens (`-h`), Gasteiger charges, PDBQT writing, MMFF94
minimization and descriptors. Version 2.0.0 uses **Open Babel 3.2.1**; 1.x used
the Debian package of 3.1.1. This page explains why, what AutoDocker relies on,
and how much the results change.

## Where Open Babel comes from

| Setting | Open Babel | How |
|---|---|---|
| Docker image | 3.2.1 | built from the submodule `tools/openbabel/repo` (tag `openbabel-3-2-1`, commit `0e94434fa75c9f61095023e3c12e0d5f2ac035ff`) in the `openbabel_builder` stage, with `-ffp-contract=off` like the other tools, on amd64 and arm64; no GUI, bindings, tests, OpenMP or downloaded components; the build fails unless `obabel -V` reports 3.2.1. SMINA links the same library. |
| CI (pytest, Python 3.9-3.12) | 3.2.1 | the upstream PyPI wheel (CPython 3.12, manylinux x86_64), pinned by SHA-256; its native `obabel` binary is used for every matrix Python |
| Native install | whatever `obabel` is on PATH | conda-forge, pip or source; see [`installation.md`](installation.md) |

The pip wheel and the source build gave byte-identical output for receptor
protonation, ligand `-h` conversion, MMFF94 minimization and seeded `--gen3d`.

## What Open Babel 3.2 changes for AutoDocker

- **Seedable random generator** (upstream PR
  [#2241](https://github.com/openbabel/openbabel/pull/2241); issues
  [#2420](https://github.com/openbabel/openbabel/issues/2420) and
  [#1934](https://github.com/openbabel/openbabel/issues/1934)). Open Babel adds
  some hydrogens in a random direction (waters, isolated atoms, two-atom
  fragments, some amines read from MOL2) and `--gen3d` uses random numbers. 3.1
  seeds that generator from the clock; 3.2 reads `OB_RANDOM_SEED` (unset: a
  random seed).
- **No extra hydrogen on pre-protonated Lys/Arg** (upstream PR
  [#2380](https://github.com/openbabel/openbabel/pull/2380)). 3.1's pH model
  added a hydrogen to every atom it charged, also when the hydrogens were
  already there: a pre-protonated receptor got a fourth hydrogen on Lys NZ and
  a third on Arg NH1/NH2, and the Arg nitrogen was written as acceptor `NA`,
  which Vina types N_DA.
- **PDB header records in PDBQT output.** The 3.2 PDBQT writer copies the
  input's PDB header records (COMPND, AUTHOR, REMARK, ...) after
  `REMARK  Name =`. Vina rejects all of them except REMARK ("Unknown or
  inappropriate tag"), and COMPND held the input's path. AutoDocker removes
  them from every PDBQT Open Babel writes for it (`utils.run_obabel`).
- The rigid (`-xr`) receptor output has no `REMARK x y z ...` column-header
  lines, and every PDBQT atom record carries occupancy 1.00 instead of 0.00.
  Vina reads neither.

## Seeds

- Every Open Babel call of the pipeline (`utils.run_obabel`, and so every
  receptor, ligand, charge, minimization and descriptor call) runs with
  `OMP_NUM_THREADS=1` and `OB_RANDOM_SEED` = the run's `--seed`
  (`defaults.SEED`, 42, for library callers that pass none). The benchmark
  scripts set it too (`_benchlib.obabel_env`).
- A fixed seed makes Open Babel repeatable, but a hydrogen it placed by a
  random draw is still not placed by chemistry, and it would move with the
  seed. Conversions that add hydrogens therefore run with the seed S and are
  compared with a probe conversion with S+1 (then S+2 when the two agree;
  seeds wrap within 1..2147483647). A hydrogen that moves between seeds is
  re-placed from heavy-atom geometry (`hydrogens.py`), as are the structurally
  undefined groups (isolated atoms, two-atom fragments). The prepared receptor
  and ligands are the same for every seed; only the Open Babel `--gen3d`
  fallback (RDKit missing) depends on it.
- With Open Babel 3.1, which ignores the variable, the probes differ by clock
  seed, and the third conversion covers the case of two runs drawing the same
  clock seed.
- `run_manifest.json` records `openbabel_version` (e.g. `"3.2.1"`) and
  `openbabel_random_seed`: the variable, the value, the probe values, whether
  this Open Babel honours it, and the policy in words.

## Older Open Babel

Open Babel < 3.2.0 (e.g. the 3.1.1 that Debian 12/13 and Ubuntu 24.04 ship) is
**accepted with a loud warning**, not refused: the input-hydrogen removal and the
hydrogen re-placement keep its output deterministic run to run (the sweeps below
were also run on 3.1.1: identical across runs), but its MMFF94 minimization
differs from 3.2.1, so results are not comparable with the image's. The version
is in `run_manifest.json`. The test suite skips the 3.2-only tests on such an
installation (`AUTODOCKER_REQUIRE_OPENBABEL_3_2=1`, set in CI, makes them fail
instead).

## Input hydrogens are still removed

The receptor's input hydrogens are removed before protonation at `--ph`, as in
1.x. With 3.1 that was required (the pH model over-protonated Lys/Arg). 3.2 fixes
that, so keeping them was tested on the 39 CASF-2016 and 10 DUD-E receptors
(`--chain all`), 14 of which contain hydrogens:

| Receptor preparation (49 receptors) | N_DA (Vina) | Lys NZ with 4 H | Arg NH typed `NA` | N_D | N_A | HD bonded to != 1 heavy atom |
|---|---:|---:|---:|---:|---:|---:|
| input H removed, Open Babel 3.1.1 (1.x) | 6 | 1 | 4 | 21073 | 919 | 1 |
| input H removed, Open Babel 3.2.1 (2.0) | 6 | 1 | 4 | 21073 | 919 | 1 |
| input H kept, Open Babel 3.1.1 | 426 | 233 | 183 | 20794 | 778 | 4 |
| input H kept, Open Babel 3.2.1 | 6 | 1 | 4 | 21213 | 779 | 2 |

Kept hydrogens on 3.2.1 no longer over-protonate, but they carry the
protonation choices of whoever prepared the file: 145 nitrogens change from
acceptor to donor and 5 from donor to acceptor, all in the 14 pre-protonated
receptors (72 Trp NE1, 69 + 5 His ring N, 4 Gly N; measured before the
Trp/His residue templates below, which now protonate every Trp NE1 and take
the His tautomer from the input's hydrogens), and one more hydrogen ends up
within Vina's bond distance of a second heavy atom (DUD-E fa10). Removal keeps the documented meaning of `--ph` (one protonation
model for every input), makes the receptor independent of how, or whether,
the input was protonated, and leaves the 2.0 receptors identical to 1.x in
every field Vina reads. Keeping the input's hydrogens may become an option
later.

(The remaining N_DA, the Lys NZ with four hydrogens and the HD bonded to two
heavy atoms come from input geometry, identically in every row: in CASF 1e66
the NZ of Lys A325 lies 2.3 Å from its CE and is protonated as an isolated
ammonium; in DUD-E aces the NE of Arg 46 lies 1.95 Å from OE1 of Glu 163, so
Open Babel's NE hydrogen is 1.18 Å from that oxygen.)

## Trp and His ring hydrogens: residue templates

Open Babel's bond-order perception of protein PDB files (3.1.1 and 3.2.1
alike) often fails to kekulize the indole and imidazole rings ("Failed to
kekulize aromatic bonds" in the log), and its pH model then adds no ring
hydrogen. On the 49 receptors below every Trp NE1 (280/280) was written as
acceptor `NA` without a hydrogen, and 206 of the 416 His had no ring hydrogen
(the other 210 got one on ND1). Vina calls an N a donor only when an `HD` atom
lies within its bond distance (`model::assign_bonds`, 1.1 x covalent radii),
so these nitrogens were acceptors, or neither donor nor acceptor, instead of
donors.

Since 2.0.0 these hydrogens come from residue templates
(`vspipeline/residue_templates.py`), applied to the PDBQT Open Babel wrote
(after the hydrogen stabilisation and the residue-label restoration, before the
receptor is centred):

- Every Open Babel hydrogen on Trp NE1 or His ND1/NE2 is removed; the template
  decides all of them, whether or not Open Babel perceived the ring.
- **Trp**: NE1 always gets HE1.
- **His**: the tautomer comes from the first rule that applies:
  1. residue name: HID/HSD (H on ND1), HIE/HSE (H on NE2), HIP/HSP (both); the
     residue is written as HIS;
  2. hydrogens of the input file on ND1 and/or NE2 (read from the input; the
     receptor itself is still protonated without them): ND1 -> HID, NE2 ->
     HIE, both -> HIP;
  3. metal coordination: if exactly one ring N lies within 2.8 Å of a Zn, Fe,
     Cu, Ni, Co, Mn, Cd or Hg HETATM of the input (kept by `--keep-hetero` or
     not), that N stays unprotonated and the other gets the hydrogen;
  4. default: HIE (H on NE2), the more common neutral tautomer and the default
     of Amber/CHARMM preparation tools.

  The choice does not depend on `--ph`: a His is charged (HIP) only when the
  input says so by name or hydrogens. There is no hydrogen-bond network
  optimisation (Reduce/PDB2PQR-style tautomer and flip search).
- The hydrogen lies in the ring plane on the external bisector of the N's two
  ring bonds, at 0.9845 Å (Open Babel's own sp2 N-H length, the length of the
  ND1-H it adds), named HE1/HD1/HE2, typed `HD` with charge +0.000, and is
  written right after its N; atom serials are renumbered. Open Babel computes
  the Gasteiger charges before `-p` adds hydrogens, so all its polar hydrogens
  carry +0.000 and their N the charge of the unprotonated atom; the template
  hydrogens follow that convention (Vina does not read charges).
- A protonated ring N is typed `N` (Vina N_D), an unprotonated His N `NA`
  (N_A). No other atom changes.
- `run_manifest.json` records, under `receptor_residue_templates`, the policy,
  the number of Trp HE1 added, the His states and the rule that decided each
  His (per residue), the His residues renamed, the Open Babel ring hydrogens
  replaced, residues skipped for missing ring atoms, and any template hydrogen
  within Vina's bond distance of a second atom (none in the 49 receptors).

Effect on the 49 receptors (Open Babel 3.2.1, `--chain all`; Vina XS types
from the emulation of `assign_bonds`/`assign_types`):

| | N_D | N_A | N_DA | N_P | O_A | O_DA | HD |
|---|---:|---:|---:|---:|---:|---:|---:|
| before (Open Babel ring hydrogens) | 21073 | 919 | 6 | 741 | 22026 | 2693 | 29291 |
| residue templates | 21580 | 412 | 6 | 741 | 22026 | 2693 | 29798 |

| Change (atoms) | Count |
|---|---:|
| Trp NE1 N_A -> N_D (HE1 added) | 280 |
| His NE2 N_A -> N_D | 350 |
| His ND1 N_A -> N_D | 37 |
| His ND1 N_D -> N_A (Open Babel's ND1-H replaced by HIE) | 160 |
| any other atom (type, coordinates, charge) | 0 |

His states: 329 HIE, 66 HID, 21 HIP; decided by residue name 20 (DUD-E akt1,
cdk2, hivpr), input hydrogens 78 (the pre-protonated CASF 1c5z, 1o3f, 1o5b and
DUD-E receptors), metal coordination 0, default 318. 210 Open Babel ring
hydrogens were replaced. Every Trp NE1 is N_D, every His has the donor/acceptor
pattern of its state, every new HD is bonded (Vina's rule) to its N only, three
preparations of each receptor are byte-identical, and Vina 1.2.5 scores a
ligand against each (`--score_only`).

1HSG case study (`protein.pdb`, ligands 1-4 and indinavir; `--chain all
--pockets 2 --exhaustiveness 2 -p 2`): the receptor gains HE1 on the four Trp,
and His A69/B69 become HIE (Open Babel had put the hydrogen on ND1); 332 HD
instead of 328.

| Ligand | Affinity before (kcal/mol) | Affinity with templates (kcal/mol) | Δ |
|---|---:|---:|---:|
| indinavir_ctrl | -11.576 | -11.568 | +0.008 |
| 3 | -5.81 | -5.81 | 0 |
| 4 | -5.249 | -5.249 | 0 |
| 1 | -4.783 | -4.778 | +0.005 |
| 2 | -4.624 | -4.562 | +0.062 |

The ranking is unchanged, and two runs with templates give byte-identical
outputs.

## Results change from 3.1.1 to 3.2.1

Same pipeline code, Open Babel 3.1.1 (AutoDocker 1.x preparation) against 3.2.1,
on valid inputs.

**Receptors** (39 CASF-2016 + 10 DUD-E, `--chain all`, with and without
`--keep-hetero`): every receptor keeps the same atoms, coordinates (hydrogens
included), Gasteiger charges and AutoDock types, hence the same Vina XS types
(emulation of Vina's `assign_bonds`/`assign_types`, table below). The files
differ only in the occupancy column (1.00) and the two `REMARK x y z` lines.

**Ligands** (the 165-file preparation test set: CASF-2016 MOL2/SDF/PDB, DUD-E
MOL2, 2D, 3D and V3000 SDF, salts, odd file names; 160 prepared, the same 5
rejected for the same reasons):

| Change in the ligand PDBQT | Ligands |
|---|---:|
| occupancy column only | 121 |
| coordinates of MMFF94-minimized SDF ligands (36 of the 69 minimized; heavy-atom RMSD median 0.20 Å, max 0.88 Å) | 36 |
| benzene from a 2D SDF: carbons typed `A` (aromatic) instead of `C` (Vina's XS type is the same) | 1 |
| ethane: the torsion REMARK lists the C-C bond (`TORSDOF 0` unchanged) | 1 |

With `--no-minimize`, 157 differ in the occupancy column only, plus benzene and
ethane as above: the coordinate changes all come from Open Babel's MMFF94
minimization. Atom types, charges and torsion trees are otherwise unchanged,
and every HD is bonded (Vina's 1.1 x covalent-radius rule) to exactly one heavy
atom in all 160 ligands, with either version.

**End to end: 1HSG case study** (`protein.pdb`, ligands 1-4 as PDBQT and
indinavir as PDB; `--chain all --pockets 2 --exhaustiveness 2 -p 2`):

| Ligand | Affinity, 3.1.1 (kcal/mol) | Affinity, 3.2.1 (kcal/mol) | Δ |
|---|---:|---:|---:|
| indinavir_ctrl | -11.576 | -11.576 | 0 |
| 3 | -5.81 | -5.81 | 0 |
| 4 | -5.249 | -5.249 | 0 |
| 1 | -4.783 | -4.783 | 0 |
| 2 | -4.624 | -4.624 | 0 |

`ranking.csv`, `Top_hits.txt`, `metrics.txt`, `Results_full.txt`, `grid.conf`
and the poses of ligands 1-4 are byte-identical. `receptor.pdbqt` (1842 atoms,
328 HD) has the same atoms, coordinates, charges and types; it differs in the
occupancy column and the two REMARK lines, as do the indinavir PDBQT and pose
and the structures embedded in `results_report.html`.

### Vina XS types per receptor

Counts with Open Babel 3.2.1 (input hydrogens removed); "3.1.1 = 3.2.1" means
every XS-type count and the HD count are equal. The last two columns show the
input-hydrogen experiment above.

| Receptor | N_D | N_A | N_DA | N_P | O_DA | O_A | HD | 3.1.1 = 3.2.1 | input H kept, 3.1.1: N_DA / Lys NZ with 4 H | input H kept, 3.2.1: N_DA / Lys NZ with 4 H |
|---|---:|---:|---:|---:|---:|---:|---:|:---:|---:|---:|
| CASF 1bcu | 370 | 15 | 0 | 16 | 36 | 386 | 516 | yes | 0 / 0 | 0 / 0 |
| CASF 1bzc | 386 | 20 | 0 | 12 | 46 | 411 | 542 | yes | 0 / 0 | 0 / 0 |
| CASF 1c5z | 318 | 20 | 0 | 12 | 52 | 319 | 449 | yes | 29 / 15 | 0 / 0 |
| CASF 1e66 | 635 | 32 | 0 | 32 | 83 | 672 | 858 | yes | 0 / 1 | 0 / 1 |
| CASF 1eby | 245 | 7 | 0 | 10 | 18 | 252 | 329 | yes | 0 / 0 | 0 / 0 |
| CASF 1g2k | 243 | 7 | 0 | 10 | 20 | 250 | 327 | yes | 0 / 0 | 0 / 0 |
| CASF 1gpk | 635 | 35 | 0 | 32 | 80 | 675 | 852 | yes | 0 / 0 | 0 / 0 |
| CASF 1gpn | 623 | 32 | 3 | 32 | 81 | 672 | 835 | yes | 3 / 0 | 3 / 0 |
| CASF 1h22 | 639 | 33 | 1 | 32 | 81 | 684 | 857 | yes | 1 / 0 | 1 / 0 |
| CASF 1h23 | 638 | 32 | 1 | 32 | 81 | 690 | 859 | yes | 1 / 0 | 1 / 0 |
| CASF 1k1i | 262 | 9 | 0 | 8 | 54 | 270 | 376 | yes | 0 / 0 | 0 / 0 |
| CASF 1lpg | 369 | 15 | 0 | 9 | 47 | 383 | 514 | yes | 0 / 0 | 0 / 0 |
| CASF 1mq6 | 363 | 15 | 0 | 9 | 46 | 378 | 505 | yes | 0 / 0 | 0 / 0 |
| CASF 1nc1 | 553 | 19 | 0 | 12 | 51 | 610 | 720 | yes | 0 / 0 | 0 / 0 |
| CASF 1nc3 | 553 | 18 | 0 | 12 | 51 | 608 | 720 | yes | 0 / 0 | 0 / 0 |
| CASF 1nvq | 340 | 17 | 0 | 15 | 28 | 368 | 468 | yes | 0 / 0 | 0 / 0 |
| CASF 1o0h | 320 | 14 | 0 | 8 | 62 | 324 | 476 | yes | 0 / 0 | 0 / 0 |
| CASF 1o3f | 262 | 9 | 0 | 8 | 54 | 270 | 376 | yes | 17 / 14 | 0 / 0 |
| CASF 1o5b | 317 | 20 | 0 | 12 | 50 | 318 | 446 | yes | 29 / 15 | 0 / 0 |
| CASF 1owh | 305 | 21 | 0 | 12 | 51 | 310 | 430 | yes | 0 / 0 | 0 / 0 |
| CASF 1oyt | 364 | 16 | 0 | 16 | 37 | 378 | 507 | yes | 0 / 0 | 0 / 0 |
| CASF 1p1n | 307 | 7 | 0 | 7 | 46 | 322 | 419 | yes | 0 / 0 | 0 / 0 |
| CASF 1p1q | 937 | 22 | 0 | 21 | 137 | 981 | 1318 | yes | 0 / 0 | 0 / 0 |
| CASF 1ps3 | 1289 | 87 | 0 | 58 | 167 | 1339 | 1766 | yes | 0 / 0 | 0 / 0 |
| CASF 1pxn | 366 | 18 | 1 | 19 | 40 | 374 | 504 | yes | 1 / 0 | 1 / 0 |
| CASF 1q8t | 467 | 20 | 0 | 14 | 46 | 485 | 654 | yes | 0 / 0 | 0 / 0 |
| CASF 1q8u | 473 | 19 | 0 | 14 | 46 | 488 | 666 | yes | 0 / 0 | 0 / 0 |
| CASF 1qf1 | 384 | 16 | 0 | 8 | 79 | 415 | 539 | yes | 0 / 0 | 0 / 0 |
| CASF 1qkt | 307 | 24 | 0 | 8 | 35 | 323 | 407 | yes | 0 / 0 | 0 / 0 |
| CASF 1r5y | 463 | 15 | 0 | 17 | 49 | 467 | 622 | yes | 0 / 0 | 0 / 0 |
| CASF 1s38 | 482 | 17 | 0 | 17 | 50 | 475 | 656 | yes | 0 / 0 | 0 / 0 |
| CASF 1sqa | 306 | 20 | 0 | 12 | 51 | 310 | 429 | yes | 0 / 0 | 0 / 0 |
| CASF 1syi | 649 | 13 | 0 | 14 | 90 | 676 | 925 | yes | 0 / 0 | 0 / 0 |
| CASF 1u1b | 324 | 10 | 0 | 8 | 62 | 324 | 480 | yes | 0 / 0 | 0 / 0 |
| CASF 1uto | 263 | 8 | 0 | 8 | 54 | 270 | 377 | yes | 17 / 14 | 0 / 0 |
| CASF 1vso | 300 | 8 | 0 | 8 | 48 | 317 | 432 | yes | 0 / 0 | 0 / 0 |
| CASF 1w4o | 320 | 14 | 0 | 8 | 62 | 324 | 476 | yes | 0 / 0 | 0 / 0 |
| CASF 1y6r | 552 | 22 | 0 | 12 | 51 | 610 | 719 | yes | 0 / 0 | 0 / 0 |
| CASF 1yc1 | 262 | 6 | 0 | 4 | 39 | 288 | 365 | yes | 0 / 0 | 0 / 0 |
| DUD-E aces | 530 | 27 | 0 | 25 | 65 | 555 | 759 | yes | 36 / 14 | 0 / 0 |
| DUD-E akt1 | 377 | 17 | 0 | 10 | 43 | 389 | 548 | yes | 41 / 22 | 0 / 0 |
| DUD-E aldr | 357 | 19 | 0 | 21 | 42 | 399 | 511 | yes | 32 / 20 | 0 / 0 |
| DUD-E andr | 290 | 15 | 0 | 10 | 27 | 295 | 399 | yes | 21 / 10 | 0 / 0 |
| DUD-E braf | 405 | 26 | 0 | 11 | 42 | 419 | 584 | yes | 32 / 18 | 0 / 0 |
| DUD-E cdk2 | 344 | 18 | 0 | 18 | 35 | 359 | 477 | yes | 35 / 20 | 0 / 0 |
| DUD-E egfr | 312 | 14 | 0 | 12 | 30 | 320 | 460 | yes | 31 / 19 | 0 / 0 |
| DUD-E fa10 | 302 | 10 | 0 | 5 | 37 | 306 | 431 | yes | 29 / 16 | 0 / 0 |
| DUD-E hivpr | 243 | 7 | 0 | 10 | 20 | 250 | 327 | yes | 20 / 12 | 0 / 0 |
| DUD-E hmdh | 722 | 14 | 0 | 21 | 91 | 718 | 1079 | yes | 51 / 23 | 0 / 0 |

## Verification

- Receptor preparation of the 49 receptors, three times each with Open Babel
  3.2.1 (with and without `--keep-hetero`, which keeps the waters whose
  hydrogens are placed at random) and twice with 3.1.1: byte-identical within
  each version. Every HD is bonded (Vina's rule) to exactly one heavy atom
  except the DUD-E aces Arg 46 hydrogen described above. A water-containing
  test receptor prepares identically with seeds 42, 43 and 2147483647
  (`tests/test_openbabel.py`).
- Ligand preparation of the 165-file test set, three times with 3.2.1 and twice
  with 3.1.1: byte-identical within each version.
- 1HSG end to end with 3.2.1, twice with `-p 2` and once with `-p 1`: every
  output byte-identical, `results_report.html` included (not compared:
  `run_manifest.json`, logs, `pocket_summary.csv`, `protein_clean_out/`, the
  docking checkpoint).
- `tests/test_openbabel.py`: the Open Babel version (CI requires >= 3.2.0),
  `OB_RANDOM_SEED` reproducibility of a water-containing receptor (same seed:
  identical; other seed: different), the version policy and the PDBQT
  header-record removal.

Not covered here and to be confirmed by an image build (CI `docker` job): the
from-source build itself on amd64 and arm64, SMINA linked against Open Babel
3.2.1 inside the image (a local build linked against it scores the 1HSG
indinavir pose identically to the 1.x image's SMINA), and the arm64 results.
