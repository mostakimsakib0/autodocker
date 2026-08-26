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
    assert fl._parse_residue_key("A45X") is None


def test_resnum_matches():
    assert fl._resnum_matches("45", "45") is True
    assert fl._resnum_matches("45A", "45") is False


def test_detect_flexible_residues(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _write_receptor(rec, [
        _atom(1, "CA", "ALA", "A", "45", 0.5, 0.5, 0.5, "C"),
        _atom(2, "CB", "ALA", "A", "45", 100.0, 100.0, 100.0, "C"),
    ])
    flex = fl.detect_flexible_residues(str(rec), (0.0, 0.0, 0.0), radius=8.0)
    assert flex == ["A45"]


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
    atoms = [
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
    ]
    block = fl._build_flex_residue_block(atoms)
    assert block
    text = "".join(block)
    assert "BEGIN_RES" in text and "ROOT" in text and "ENDROOT" in text
    assert "BRANCH" in text and "END_RES" in text


def test_build_flex_residue_block_no_ca():
    atoms = [{"serial": 1, "name": "CB", "atom_type": "C", "element": "C",
              "coords": (0, 0, 0), "line": "x\n"}]
    assert fl._build_flex_residue_block(atoms) == []


def test_build_flexible_residue_pdbqt(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _write_receptor(rec, [
        _atom(1, "CA", "ALA", "A", "45", 0, 0, 0, "C"),
        _atom(2, "CB", "ALA", "A", "45", 0, 1.5, 0, "C"),
        _atom(3, "N", "ALA", "A", "45", -1, 0.5, 0.5, "N"),
    ])
    out = fl.build_flexible_residue_pdbqt(str(rec), ["A45"], str(tmp_path))
    assert out is not None
    assert os.path.exists(out)
    assert "BEGIN_RES" in open(out).read()


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
