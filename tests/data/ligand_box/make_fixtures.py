#!/usr/bin/env python3
"""Regenerate the committed --box-ligand / --box-from-het fixtures.

Not run by the test suite; kept so the fixtures' provenance is reproducible.

  * 1hsg_site.pdb : PDB entry 1HSG (HIV-1 protease with indinavir, MK1 B 902),
    first model: every residue of chains A and B with an atom within 8 A of
    MK1, the MK1 HETATM records, the waters within 6 A of MK1, then END.
    Records are copied byte for byte (input frame of 1HSG).
  * mk1.pdb       : the MK1 HETATM records alone
  * mk1.sdf, mk1.mol2, mk1.pdbqt : mk1.pdb converted with Open Babel
    (``obabel mk1.pdb -O mk1.<ext>``; heavy atoms only, as in the PDB entry)

Usage: python3 tests/data/ligand_box/make_fixtures.py path/to/1HSG.pdb
"""
import math
import os
import subprocess
import sys

OUT = os.path.dirname(os.path.abspath(__file__))
SITE_RADIUS = 8.0
WATER_RADIUS = 6.0


def _xyz(line):
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def main():
    src = sys.argv[1]
    with open(src, encoding="utf-8") as f:
        lines = [line for line in f]
    end = next((i for i, line in enumerate(lines) if line.startswith("ENDMDL")), len(lines))
    lines = lines[:end]
    mk1 = [line for line in lines if line.startswith("HETATM") and line[17:20] == "MK1"]
    mk1_xyz = [_xyz(line) for line in mk1]

    def near(line, radius):
        p = _xyz(line)
        return any(math.dist(p, q) <= radius for q in mk1_xyz)

    keep_res = {(line[21], line[22:27]) for line in lines
                if line.startswith("ATOM") and near(line, SITE_RADIUS)}
    site = []
    for line in lines:
        if line.startswith("ATOM") and (line[21], line[22:27]) in keep_res:
            site.append(line)
        elif line.startswith("TER"):
            site.append(line)
    waters = [line for line in lines if line.startswith("HETATM")
              and line[17:20] == "HOH" and near(line, WATER_RADIUS)]
    with open(os.path.join(OUT, "1hsg_site.pdb"), "w", encoding="utf-8") as f:
        f.write("REMARK   1HSG binding-site fragment, see make_fixtures.py\n")
        f.writelines(site + mk1 + waters)
        f.write("END\n")
    with open(os.path.join(OUT, "mk1.pdb"), "w", encoding="utf-8") as f:
        f.writelines(mk1)
        f.write("END\n")
    for ext in ("sdf", "mol2", "pdbqt"):
        # relative names: Open Babel writes the input name into the files
        subprocess.run(["obabel", "mk1.pdb", "-O", f"mk1.{ext}"], check=True, cwd=OUT)


if __name__ == "__main__":
    main()
