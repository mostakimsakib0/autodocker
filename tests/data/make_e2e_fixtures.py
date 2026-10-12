#!/usr/bin/env python3
"""Regenerate the committed end-to-end docking fixtures in tests/data/e2e/.

Not run by the test suite; kept so the fixtures' provenance is reproducible.

  * receptor: binding-site residues (any atom within 8 A of the ligand) of
    RCSB entry 2FXS (HSP90 N-terminal domain, ligand RDA A1001), first model,
    altloc A, heavy atoms; written as pocket_2fxs.pdb and converted with
    ``obabel -xr -h --partialcharge gasteiger`` to receptor.pdbqt
  * box.json : centroid of the crystal ligand's heavy atoms, 16 A cube
  * ligands  : three drug-like molecules with rotatable bonds, one seeded
    RDKit ETKDGv3 conformer each (randomSeed 42, MMFF94), converted with
    ``obabel --partialcharge gasteiger`` to PDBQT

Usage: python3 tests/data/make_e2e_fixtures.py [--pdb 2fxs.pdb]
"""

import argparse
import json
import math
import os
import subprocess
import urllib.request

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "e2e")
# Small on purpose: the tests dock each ligand three times.
LIGANDS = {
    "lig_ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "lig_paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "lig_sulfonamide": "Cc1ccc(cc1)S(=O)(=O)NCCc1ccccc1",
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", help="local copy of 2FXS (else downloaded)")
    args = ap.parse_args()
    if args.pdb:
        text = open(args.pdb, encoding="utf-8").read()
    else:
        text = urllib.request.urlopen(
            "https://files.rcsb.org/download/2FXS.pdb", timeout=120).read().decode()
    os.makedirs(OUT, exist_ok=True)

    records = []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")) and line[16] in " A":
            if (line[76:78].strip() or line[12:14].strip()).upper() in ("H", "D"):
                continue
            records.append(line[:16] + " " + line[17:])
    lig = [ln for ln in records if ln.startswith("HETATM") and ln[17:20] == "RDA"
           and ln[21] == "A"]
    lig_xyz = [(float(ln[30:38]), float(ln[38:46]), float(ln[46:54])) for ln in lig]
    center = [sum(c[i] for c in lig_xyz) / len(lig_xyz) for i in range(3)]

    near = set()
    for ln in records:
        if not ln.startswith("ATOM"):
            continue
        xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
        if any(math.dist(xyz, c) <= 8.0 for c in lig_xyz):
            near.add((ln[21], ln[22:27]))
    pocket = [ln for ln in records if ln.startswith("ATOM") and (ln[21], ln[22:27]) in near]
    pdb = os.path.join(OUT, "pocket_2fxs.pdb")
    with open(pdb, "w", encoding="utf-8") as f:
        f.write("REMARK binding-site residues of RCSB 2FXS (8 A around RDA A1001)\n")
        f.write("\n".join(pocket) + "\nEND\n")
    # Relative names and cwd=OUT keep local paths out of the REMARK lines.
    subprocess.run(["obabel", "-ipdb", "pocket_2fxs.pdb", "-xr", "-h",
                    "--partialcharge", "gasteiger", "-opdbqt", "-O", "receptor.pdbqt"],
                   check=True, capture_output=True, cwd=OUT)
    with open(os.path.join(OUT, "box.json"), "w", encoding="utf-8") as f:
        json.dump({"center": [round(c, 3) for c in center], "size": 16.0}, f, indent=2)

    from rdkit import Chem
    from rdkit.Chem import AllChem
    for name, smi in LIGANDS.items():
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 1
        assert AllChem.EmbedMolecule(mol, params) == 0
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        sdf = os.path.join(OUT, name + ".tmp.sdf")
        with open(sdf, "w", encoding="utf-8") as f:
            f.write(Chem.MolToMolBlock(mol) + "$$$$\n")
        subprocess.run(["obabel", "-isdf", name + ".tmp.sdf", "-opdbqt", "-O",
                        name + ".pdbqt", "--partialcharge", "gasteiger"],
                       check=True, capture_output=True, cwd=OUT)
        os.remove(sdf)
    print(f"{len(pocket)} pocket atoms, {len(near)} residues, box center {center}")


if __name__ == "__main__":
    main()
