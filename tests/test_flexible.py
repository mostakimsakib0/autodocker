import os

import pytest

import vspipeline.flexible as fl


def _atom(serial, name, resname, chain, resnum, x, y, z, atom_type):
    line = "ATOM  "
    line += f"{serial:>5d}"
    line += " "
    line += f"{name:<4s}"
    line += " "
    line += f"{resname:<3s}"
    line += " "
    line += chain
    line += f"{resnum:>4s}"
    line = line.ljust(30)
    line += f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
    line = line.ljust(76)
    line += f" {atom_type}"
    return line + "\n"


def _write_receptor(path, atoms_lines):
    path.write_text("".join(atoms_lines))


def test_pdbqt_atom_element():
    assert fl._pdbqt_atom_element("C") == "C"
    assert fl._pdbqt_atom_element("OA") == "O"
    assert fl._pdbqt_atom_element("SA") == "S"
    assert fl._pdbqt_atom_element("HD") == "H"
    assert fl._pdbqt_atom_element("Zn") == "ZN"
    assert fl._pdbqt_atom_element("UNKNOWN") == "C"


def test_parse_residue_key():
    assert fl._parse_residue_key("A45") == ("A", "45")
    assert fl._parse_residue_key("B102") == ("B", "102")
    assert fl._parse_residue_key("bad") is None
    assert fl._parse_residue_key("AB45") is None
    # insertion code and blank chain
    assert fl._parse_residue_key("A100B") == ("A", "100B")
    assert fl._parse_residue_key("45") == ("", "45")


def test_resnum_matches():
    assert fl._resnum_matches("45", "45") is True
    assert fl._resnum_matches("45A", "45") is False
    assert fl._resnum_matches("100B", "100b") is True


def test_detect_flexible_residues(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _write_receptor(rec, [
        _atom(1, "CA", "SER", "A", "45", 0.5, 0.5, 0.5, "C"),
        _atom(2, "CB", "SER", "A", "45", 100.0, 100.0, 100.0, "C"),
        # GLY has no side-chain torsions and is never selected
        _atom(3, "CA", "GLY", "A", "46", 0.1, 0.1, 0.1, "C"),
        # ranked by its closest atom (CB at 1.0), not its first atom (CA at 7.0)
        _atom(4, "CA", "LEU", "A", "47", 7.0, 0.0, 0.0, "C"),
        _atom(5, "CB", "LEU", "A", "47", 1.0, 0.0, 0.0, "C"),
    ])
    flex = fl.detect_flexible_residues(str(rec), (0.0, 0.0, 0.0), radius=8.0)
    assert flex == ["A45", "A47"]


def test_residue_atoms_from_pdbqt(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _write_receptor(rec, [
        _atom(1, "CA", "ALA", "A", "45", 0, 0, 0, "C"),
        _atom(2, "CB", "ALA", "A", "45", 0, 1.5, 0, "C"),
        _atom(3, "CA", "ALA", "B", "10", 0, 0, 0, "C"),
    ])
    atoms = fl._residue_atoms_from_pdbqt(str(rec), "A45")
    assert len(atoms) == 2
    assert {a["name"] for a in atoms} == {"CA", "CB"}


def test_build_flex_residue_block():
    atoms = [dict(a, resname="SER") for a in [
        {"serial": 1, "name": "CA", "atom_type": "C", "element": "C",
         "coords": (0.0, 0.0, 0.0), "line": "CA line\n"},
        {"serial": 2, "name": "N", "atom_type": "N", "element": "N",
         "coords": (-1.0, 0.5, 0.5), "line": "N line\n"},
        {"serial": 3, "name": "C", "atom_type": "C", "element": "C",
         "coords": (1.0, -0.5, 0.5), "line": "C line\n"},
        {"serial": 4, "name": "O", "atom_type": "OA", "element": "O",
         "coords": (1.5, -1.0, 0.5), "line": "O line\n"},
        {"serial": 5, "name": "CB", "atom_type": "C", "element": "C",
         "coords": (0.0, 1.5, 0.0), "line": "CB line\n"},
        {"serial": 6, "name": "OG", "atom_type": "OA", "element": "O",
         "coords": (0.0, 2.5, 0.0), "line": "OG line\n"},
    ]]
    block = fl._build_flex_residue_block(atoms)
    assert block
    text = "".join(block)
    assert "BEGIN_RES" in text and "ROOT" in text and "ENDROOT" in text
    assert "BRANCH" in text and "END_RES" in text
    # Backbone N, C, O stay in the rigid receptor; ROOT is CA alone.
    assert block[1:4] == ["ROOT\n", "CA line\n", "ENDROOT\n"]
    assert "N line" not in text and "O line" not in text
    assert [l.split()[1:] for l in block if l.startswith("BRANCH")] == [["1", "5"], ["5", "6"]]


def test_build_flex_residue_block_no_ca():
    atoms = [{"serial": 1, "name": "CB", "atom_type": "C", "element": "C",
              "coords": (0, 0, 0), "line": "x\n"}]
    assert fl._build_flex_residue_block(atoms) == []


def test_build_flexible_residue_pdbqt(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _write_receptor(rec, [
        _atom(1, "CA", "SER", "A", "45", 0, 0, 0, "C"),
        _atom(2, "CB", "SER", "A", "45", 0, 1.5, 0, "C"),
        _atom(3, "N", "SER", "A", "45", -1, 0.5, 0.5, "N"),
        _atom(4, "OG", "SER", "A", "45", 0, 2.9, 0, "OA"),
        _atom(5, "CA", "GLY", "A", "46", 5, 0, 0, "C"),
    ])
    built = []
    out = fl.build_flexible_residue_pdbqt(str(rec), ["A45", "A46"], str(tmp_path),
                                          built=built)
    assert out is not None
    assert built == ["A45"]  # GLY cannot be flexible and stays rigid
    assert os.path.exists(out)
    assert "BEGIN_RES" in open(out, encoding="utf-8").read()


def test_build_flexible_residue_pdbqt_empty(tmp_path):
    assert fl.build_flexible_residue_pdbqt("x.pdbqt", [], str(tmp_path)) is None


def test_remove_residues_from_pdbqt(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _write_receptor(rec, [
        _atom(1, "CA", "ALA", "A", "45", 0, 0, 0, "C"),
        _atom(2, "CB", "ALA", "A", "45", 0, 1.5, 0, "C"),
        _atom(3, "CA", "ALA", "B", "10", 0, 0, 0, "C"),
    ])
    out = tmp_path / "rigid.pdbqt"
    fl._remove_residues_from_pdbqt(str(rec), ["A45"], str(out))
    content = out.read_text()
    assert "B" in content and "A45" not in content.split("\n")[0]
    # A45 lines removed
    assert "A  45" not in content


def _phe_atoms():
    # Phenylalanine side chain: CA-CB-CG and a planar ring (1.39 A bonds).
    import math
    ring = []
    for k, name in enumerate(["CG", "CD1", "CE1", "CZ", "CE2", "CD2"]):
        # CG (k=0) points back at CB: ring centre 1.51 + 1.39 A beyond CB.
        ang = math.radians(180 + 60 * k)
        ring.append((name, 4.43 + 1.39 * math.cos(ang), 1.39 * math.sin(ang)))
    atoms = [("N", -1.0, 1.0, 0.0, "N"), ("CA", 0.0, 0.0, 0.0, "C"),
             ("C", -0.5, -1.4, 0.0, "C"), ("O", -1.5, -1.9, 0.0, "OA"),
             ("CB", 1.53, 0.0, 0.0, "C")]
    atoms += [(n, x + 0.0, y, 0.0, "A") for n, x, y in ring]
    return atoms


def test_aromatic_ring_is_one_rigid_fragment(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _write_receptor(rec, [_atom(i, n, "PHE", "A", "53", x, y, z, t)
                          for i, (n, x, y, z, t) in enumerate(_phe_atoms(), 1)])
    block = fl._build_flex_residue_block(fl._residue_atoms_from_pdbqt(str(rec), "A53"))
    branches = [l for l in block if l.startswith("BRANCH")]
    assert len(branches) == 2  # chi1 (CA-CB) and chi2 (CB-CG) only
    atom_lines = [l for l in block if l.startswith("ATOM")]
    assert len(atom_lines) == 8  # CA, CB and the six ring carbons


def test_proline_is_not_flexible(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _write_receptor(rec, [
        _atom(1, "CA", "PRO", "A", "81", 0, 0, 0, "C"),
        _atom(2, "CB", "PRO", "A", "81", 1.5, 0, 0, "C"),
    ])
    assert fl._build_flex_residue_block(fl._residue_atoms_from_pdbqt(str(rec), "A81")) == []


def test_remove_residues_keeps_backbone(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _write_receptor(rec, [
        _atom(1, "N", "SER", "A", "45", -1, 0, 0, "N"),
        _atom(2, "CA", "SER", "A", "45", 0, 0, 0, "C"),
        _atom(3, "CB", "SER", "A", "45", 0, 1.5, 0, "C"),
        _atom(4, "C", "SER", "A", "45", 1, 0, 0, "C"),
        _atom(5, "O", "SER", "A", "45", 1, 1, 0, "OA"),
    ])
    out = tmp_path / "rigid.pdbqt"
    fl._remove_residues_from_pdbqt(str(rec), ["A45"], str(out))
    names = [l[12:16].strip() for l in out.read_text().splitlines()]
    assert names == ["N", "C", "O"]


def test_flex_files_are_accepted_by_vina(tmp_path):
    """Real Vina must parse the rigid/flex pair built from a protonated receptor."""
    import shutil
    import subprocess
    vina = shutil.which("vina")
    obabel = shutil.which("obabel")
    if not (vina and obabel):
        pytest.skip("vina/obabel not installed")
    # SER A37 and PHE A53 of HIV-1 protease (PDB 1HSG), protonated like a receptor.
    (tmp_path / "rec.pdb").write_text(
        "ATOM    280  N   SER A  37      12.249  30.179  26.852  1.00 34.23           N  \n"
        "ATOM    281  CA  SER A  37      11.163  31.036  27.460  1.00 33.81           C  \n"
        "ATOM    282  C   SER A  37      11.347  32.442  26.960  1.00 32.96           C  \n"
        "ATOM    283  O   SER A  37      12.475  32.929  27.007  1.00 32.62           O  \n"
        "ATOM    284  CB  SER A  37      11.300  31.071  28.976  1.00 34.02           C  \n"
        "ATOM    285  OG  SER A  37       9.948  31.148  29.559  1.00 38.20           O  \n"
        "ATOM    395  N   PHE A  53       3.806  18.788   9.925  1.00 31.50           N  \n"
        "ATOM    396  CA  PHE A  53       3.593  19.996  10.731  1.00 30.76           C  \n"
        "ATOM    397  C   PHE A  53       4.015  19.881  12.194  1.00 31.86           C  \n"
        "ATOM    398  O   PHE A  53       3.930  18.810  12.781  1.00 36.37           O  \n"
        "ATOM    399  CB  PHE A  53       2.121  20.351  10.670  1.00 28.77           C  \n"
        "ATOM    400  CG  PHE A  53       1.760  21.152   9.484  1.00 30.51           C  \n"
        "ATOM    401  CD1 PHE A  53       1.725  20.567   8.216  1.00 27.94           C  \n"
        "ATOM    402  CD2 PHE A  53       1.556  22.518   9.607  1.00 30.74           C  \n"
        "ATOM    403  CE1 PHE A  53       1.500  21.332   7.075  1.00 29.27           C  \n"
        "ATOM    404  CE2 PHE A  53       1.327  23.302   8.496  1.00 32.36           C  \n"
        "ATOM    405  CZ  PHE A  53       1.290  22.718   7.212  1.00 28.01           C  \n"
        "END\n")
    # The real receptor preparation: protonation at pH 7.4, labels restored
    # (two separate residues make Open Babel re-perceive and rename them).
    import vspipeline.protein as protein
    prep = protein.ProteinPreparation(str(tmp_path / "rec.pdb"), str(tmp_path))
    prep.prepare_receptor(chain="all")
    rec = prep.receptor_pdbqt
    built = []
    flex = fl.build_flexible_residue_pdbqt(str(rec), ["A53", "A37"], str(tmp_path), built=built)
    assert built == ["A53", "A37"]
    rigid = fl._remove_residues_from_pdbqt(str(rec), built, str(tmp_path / "rigid.pdbqt"),
                                           flex_file=flex)
    # Every atom is in exactly one of the two files.
    n_atoms = lambda p: sum(1 for l in open(p, encoding="utf-8") if l.startswith(("ATOM", "HETATM")))
    assert n_atoms(rigid) + n_atoms(flex) == n_atoms(rec)
    # The serine hydroxyl hydrogen (named just "H" by Open Babel) moves with OG.
    ser_block = open(flex, encoding="utf-8").read().split("BEGIN_RES")[2]
    assert " SER " in ser_block
    assert any(l[77:79].strip() == "HD" for l in ser_block.splitlines()
               if l.startswith("ATOM"))
    lig = tmp_path / "lig.pdbqt"
    subprocess.run([obabel, "-:CCO", "--gen3d", "-opdbqt", "-O", str(lig)],
                   check=True, capture_output=True)
    proc = subprocess.run([vina, "--receptor", rigid, "--flex", flex, "--ligand", str(lig),
                           "--center_x", "0", "--center_y", "0", "--center_z", "0",
                           "--size_x", "30", "--size_y", "30", "--size_z", "30",
                           "--score_only"], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_detect_flexible_residues_ignores_hydrogens(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _write_receptor(rec, [
        _atom(1, "CA", "SER", "A", "45", 5.0, 0.0, 0.0, "C"),
        _atom(2, "OG", "SER", "A", "45", 5.5, 0.0, 0.0, "OA"),
        _atom(3, "H", "SER", "A", "45", 0.1, 0.0, 0.0, "HD"),  # misplaced H near site
        _atom(4, "CA", "THR", "A", "46", 3.0, 0.0, 0.0, "C"),
    ])
    assert fl.detect_flexible_residues(str(rec), (0.0, 0.0, 0.0), radius=8.0) == ["A46", "A45"]


def test_build_flexible_duplicate_keys_built_once(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _write_receptor(rec, [
        _atom(1, "N", "SER", "A", "45", 0.0, 0.0, 0.0, "N"),
        _atom(2, "CA", "SER", "A", "45", 1.45, 0.0, 0.0, "C"),
        _atom(3, "C", "SER", "A", "45", 2.0, 1.4, 0.0, "C"),
        _atom(4, "CB", "SER", "A", "45", 2.0, -0.8, 1.2, "C"),
        _atom(5, "OG", "SER", "A", "45", 3.4, -0.8, 1.2, "OA"),
    ])
    built = []
    flex = fl.build_flexible_residue_pdbqt(str(rec), ["A45", "A45", "A99"], str(tmp_path),
                                           built=built)
    text = open(flex, encoding="utf-8").read()
    assert built == ["A45"]
    assert text.count("BEGIN_RES") == 1
    assert text.splitlines()[0] == "REMARK flexible residues: A45"
