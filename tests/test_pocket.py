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

    cl = pk._parse_hetatm(_het("CL", "CL", "A", "1", 0, 0, 0, element="CL"))
    assert pk._hetatm_to_pdbqt_line(cl, 3)[77:79] == "Cl"

    unsup = pk._parse_hetatm(_het("SE", "SEX", "A", "1", 0, 0, 0, element="SE"))
    assert pk._hetatm_to_pdbqt_line(unsup, 4) is None


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


def test_append_hetatm_applies_frame_offset(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _pdb(rec, ["ATOM      1  C   ALA A   1      0.0 0.0 0.0  0.0 0.0    -0.1 C"])
    pdb = tmp_path / "prot.pdb"
    _pdb(pdb, [_het("O", "HOH", "A", "1", 11.0, 22.0, 33.0, element="O")])
    pk._append_hetatm_to_receptor(str(rec), str(pdb), {("A", "1", "HOH")},
                                  {"A"}, offset=(10.0, 20.0, 30.0))
    het = [l for l in rec.read_text().splitlines() if l.startswith("HETATM")][0]
    assert (float(het[30:38]), float(het[38:46]), float(het[46:54])) == (1.0, 2.0, 3.0)


def test_kept_hetatm_columns_with_four_digit_resnum():
    parsed = pk._parse_hetatm(
        "HETATM 1601  O   HOH A1001      10.000  20.000  30.000  1.00 20.00           O\n")
    line = pk._hetatm_to_pdbqt_line(parsed, 123456)
    assert line[22:26] == "1001"
    assert (float(line[30:38]), float(line[38:46]), float(line[46:54])) == (10.0, 20.0, 30.0)
    assert line[6:11] == "23456"
    assert line[77:79].strip() == "OA"


def test_element_from_blank_element_column():
    assert pk._element_from_atom_name(" CA ") == "C"   # C-alpha
    assert pk._element_from_atom_name("CA  ") == "CA"  # calcium
    assert pk._element_from_atom_name("HG21") == "H"   # methyl hydrogen
    assert pk._element_from_atom_name("ZN  ") == "ZN"
    assert pk._element_from_atom_name(" NA ") == "N"   # heme nitrogen


@pytest.mark.parametrize("element,atom_type", [
    ("P", "P"), ("F", "F"), ("CL", "Cl"), ("BR", "Br"), ("I", "I"), ("S", "SA"),
    ("MG", "Mg"), ("NA", "Na"),
])
def test_hetatm_supported_elements(element, atom_type):
    parsed = pk._parse_hetatm(_het(element, "LIG", "A", "1", 0, 0, 0, element=element))
    assert pk._hetatm_to_pdbqt_line(parsed, 1)[77:79].strip() == atom_type


def test_resname_tables_fit_the_pdb_resname_field():
    for names in (pk.WATER_RESNAMES, pk.COFACTOR_RESNAMES):
        assert all(len(n) <= 3 for n in names)
    assert {"NAP", "NDP", "NAI"} <= pk.COFACTOR_RESNAMES
    # Every detected metal can be written to the receptor PDBQT.
    assert set(pk.METAL_ELEMENTS) <= set(pk.VINA_METAL_TYPE_MAP)


def test_metal_in_multi_atom_residue_is_not_an_ion(tmp_path, caplog):
    p = tmp_path / "prot.pdb"
    _pdb(p, [
        _het("FE", "HEM", "A", "500", 0.0, 0.0, 0.0, element="FE", serial=1),
        _het("NA", "HEM", "A", "500", 2.0, 0.0, 0.0, element="N", serial=2),
        _het("ZN", "ZN", "A", "600", 9.0, 9.0, 9.0, element="ZN", serial=3),
    ])
    import logging
    with caplog.at_level(logging.WARNING):
        metals = pk.detect_metal_ions(str(p))
    assert [(m["resname"], m["atom"]) for m in metals] == [("ZN", "ZN")]
    assert "HEM" in caplog.text


def test_append_uses_first_model_and_one_altloc(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _pdb(rec, ["ATOM      1  C   ALA A   1      0.0 0.0 0.0  0.0 0.0    -0.1 C"])
    pdb = tmp_path / "prot.pdb"
    a = _het("O", "HOH", "A", "5", 1.0, 1.0, 1.0, element="O")
    b = _het("O", "HOH", "A", "5", 1.5, 1.0, 1.0, element="O")
    a = a[:16] + "A" + a[17:54] + "  0.60" + a[60:]
    b = b[:16] + "B" + b[17:54] + "  0.40" + b[60:]
    later_model = _het("O", "HOH", "A", "5", 5.0, 5.0, 5.0, element="O")
    _pdb(pdb, ["MODEL        1", a, b, "ENDMDL", "MODEL        2", later_model, "ENDMDL"])
    assert pk._append_hetatm_to_receptor(str(rec), str(pdb), {("A", "5", "HOH")}, {"A"}) == 1
    het = [l for l in rec.read_text().splitlines() if l.startswith("HETATM")]
    assert len(het) == 1 and float(het[0][30:38]) == 1.0


def test_append_skips_residues_already_in_receptor(tmp_path):
    # --keep-hetero already put the water in receptor.pdbqt; --keep-waters
    # must not add it a second time.
    rec = tmp_path / "r.pdbqt"
    _pdb(rec, ["ATOM      1  C   ALA A   1       0.000   0.000   0.000  0.00  0.00    -0.100 C",
               "HETATM    2  O   HOH A   9       1.000   2.000   3.000  0.00  0.00    -0.400 OA"])
    pdb = tmp_path / "prot.pdb"
    _pdb(pdb, [_het("O", "HOH", "A", "9", 11.0, 22.0, 33.0, element="O")])
    n = pk._append_hetatm_to_receptor(str(rec), str(pdb), {("A", "9", "HOH")}, {"A"},
                                      offset=(10.0, 20.0, 30.0))
    assert n == 0
    assert rec.read_text().count("HOH") == 1


def test_append_rejects_residue_with_untypeable_atom(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _pdb(rec, ["ATOM      1  C   ALA A   1      0.0 0.0 0.0  0.0 0.0    -0.1 C"])
    pdb = tmp_path / "prot.pdb"
    _pdb(pdb, [_het("C1", "LIG", "A", "1", 0.0, 0.0, 0.0, element="C", serial=1),
               _het("SE1", "LIG", "A", "1", 1.9, 0.0, 0.0, element="SE", serial=2)])
    with pytest.raises(ValueError, match="SE"):
        pk._append_hetatm_to_receptor(str(rec), str(pdb), {("A", "1", "LIG")}, {"A"})
    assert "LIG" not in rec.read_text()


def test_append_types_cofactor_atoms(tmp_path):
    rec = tmp_path / "r.pdbqt"
    _pdb(rec, ["ATOM      1  C   ALA A   1      0.0 0.0 0.0  0.0 0.0    -0.1 C"])
    pdb = tmp_path / "prot.pdb"
    _pdb(pdb, [
        # N bonded to three carbons (no lone pair) and a pyridine-like N
        _het("N1", "LIG", "A", "1", 0.0, 0.0, 0.0, element="N", serial=1),
        _het("C1", "LIG", "A", "1", 1.45, 0.0, 0.0, element="C", serial=2),
        _het("C2", "LIG", "A", "1", -0.7, 1.25, 0.0, element="C", serial=3),
        _het("C3", "LIG", "A", "1", -0.7, -1.25, 0.0, element="C", serial=4),
        _het("N2", "LIG", "A", "1", 2.2, 1.2, 0.0, element="N", serial=5),
        # hydroxyl O with its H, and a methyl H on C3
        _het("O1", "LIG", "A", "1", -1.4, -2.4, 0.0, element="O", serial=6),
        _het("HO1", "LIG", "A", "1", -2.3, -2.2, 0.0, element="H", serial=7),
        _het("H31", "LIG", "A", "1", -0.2, -1.7, 0.9, element="H", serial=8),
    ])
    n = pk._append_hetatm_to_receptor(str(rec), str(pdb), {("A", "1", "LIG")}, {"A"})
    types = {l[12:16].strip(): l[77:79].strip()
             for l in rec.read_text().splitlines() if l.startswith("HETATM")}
    assert n == 7
    assert types == {"N1": "N", "C1": "C", "C2": "C", "C3": "C", "N2": "NA",
                     "O1": "OA", "HO1": "HD"}
