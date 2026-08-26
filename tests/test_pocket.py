import os

import pytest

import vspipeline.pocket as pk


def _het(atom, resname, chain, resnum, x, y, z, element="O", serial=1):
    line = "HETATM"
    line += f"{serial:>5d}"
    line += " "
    line += f"{atom:<4s}"
    line += " "
    line += f"{resname:<3s}"
    line += " "
    line += chain
    line += f"{resnum:>4s}"
    line = line.ljust(30)
    line += f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
    line = line.ljust(76)
    line += f"{element:>2s}"
    return line


def _pdb(path, lines):
    path.write_text("\n".join(lines) + "\n")


def test_parse_hetatm():
    d = _het("O", "HOH", "A", "1", 1.0, 2.0, 3.0, element="O")
    p = pk._parse_hetatm(d)
    assert p["resname"] == "HOH"
    assert p["element"] == "O"
    assert p["coords"] == (1.0, 2.0, 3.0)
    assert p["chain"] == "A"


def test_parse_hetatm_non_hetatm():
    assert pk._parse_hetatm("ATOM      1  C   ALA A   1") is None


def test_parse_hetatm_malformed():
    assert pk._parse_hetatm("HETATM short") is None


def test_detect_water_molecules(tmp_path):
    p = tmp_path / "prot.pdb"
    _pdb(p, [
        _het("O", "HOH", "A", "1", 0.0, 0.0, 0.0, element="O"),
        _het("O", "HOH", "A", "2", 100.0, 100.0, 100.0, element="O"),
    ])
    waters = pk.detect_water_molecules(str(p), (0.0, 0.0, 0.0), distance_threshold=4.0)
    assert len(waters) == 1
    assert waters[0]["resnum"] == "1"


def test_detect_metal_ions(tmp_path):
    p = tmp_path / "prot.pdb"
    _pdb(p, [
        _het("ZN", "ZN", "A", "1", 1.0, 1.0, 1.0, element="ZN"),
        _het("O", "HOH", "A", "2", 2.0, 2.0, 2.0, element="O"),
    ])
    metals = pk.detect_metal_ions(str(p))
    assert len(metals) == 1
    assert metals[0]["name"] == "Zinc"


def test_detect_cofactors(tmp_path):
    p = tmp_path / "prot.pdb"
    _pdb(p, [
        _het("NAD", "NAD", "A", "1", 1.0, 1.0, 1.0, element="C"),
        _het("O", "HOH", "A", "2", 2.0, 2.0, 2.0, element="O"),
    ])
    cofactors = pk.detect_cofactors(str(p))
    assert len(cofactors) == 1
    assert cofactors[0]["name"] == "NAD"


def test_max_pdbqt_serial(tmp_path):
    p = tmp_path / "r.pdbqt"
    _pdb(p, [
        "ATOM      1  C   ALA A   1      0.0 0.0 0.0  0.0 0.0    -0.1 C",
        "ATOM     42  N   ALA A   2      0.0 0.0 0.0  0.0 0.0    -0.3 N",
        "HETATM   10  O   HOH A   1      0.0 0.0 0.0  0.0 0.0    -0.2 O",
    ])
    assert pk._max_pdbqt_serial(str(p)) == 42


def test_hetatm_to_pdbqt_line_types():
    base = pk._parse_hetatm(_het("O", "HOH", "A", "1", 1.0, 2.0, 3.0, element="O"))
    line = pk._hetatm_to_pdbqt_line(base, 1)
    assert "OA" in line and line.startswith("HETATM")

    zn = pk._parse_hetatm(_het("ZN", "ZN", "A", "1", 0, 0, 0, element="ZN"))
    assert "Zn" in pk._hetatm_to_pdbqt_line(zn, 2)

    unsup = pk._parse_hetatm(_het("CL", "CL", "A", "1", 0, 0, 0, element="CL"))
    assert pk._hetatm_to_pdbqt_line(unsup, 3) is None


def test_append_hetatm_to_receptor(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _pdb(rec, ["ATOM      1  C   ALA A   1      0.0 0.0 0.0  0.0 0.0    -0.1 C"])
    pdb = tmp_path / "prot.pdb"
    _pdb(pdb, [
        _het("O", "HOH", "A", "1", 0.0, 0.0, 0.0, element="O"),
        _het("O", "HOH", "B", "1", 0.0, 0.0, 0.0, element="O"),
    ])
    keep = {("A", "1", "HOH")}
    appended = pk._append_hetatm_to_receptor(str(rec), str(pdb), keep, {"A"})
    assert appended == 1
    content = rec.read_text()
    assert "HETATM" in content
    assert "OA" in content
    assert "B" not in content  # chain B excluded
