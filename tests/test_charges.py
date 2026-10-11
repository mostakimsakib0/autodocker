import os

import pytest

import vspipeline.charges as ch
from vspipeline import tools, utils


def _write(path, lines):
    path.write_text("\n".join(lines) + "\n")


ATOM_NONZERO = "ATOM      1  C   ALA A   1      1.0   2.0   3.0  0.00  0.00    -0.120 C"
ATOM_ZERO = "ATOM      1  C   ALA A   1      1.0   2.0   3.0  0.00  0.00    0.000 C"


def test_extract_pdbqt_charge():
    assert ch._extract_pdbqt_charge(ATOM_NONZERO) == -0.12
    assert ch._extract_pdbqt_charge(ATOM_ZERO) == 0.0
    assert ch._extract_pdbqt_charge("ATOM foo bar") is None
    assert ch._extract_pdbqt_charge("HETATM    1 O   X   1 1 2 3 0.5") == 0.5


def test_ensure_pdbqt_has_charges(tmp_path):
    p = tmp_path / "a.pdbqt"
    _write(p, [ATOM_NONZERO])
    assert ch._ensure_pdbqt_has_charges(str(p)) is True


def test_ensure_pdbqt_has_charges_zero_only(tmp_path):
    p = tmp_path / "z.pdbqt"
    _write(p, [ATOM_ZERO, ATOM_ZERO])
    assert ch._ensure_pdbqt_has_charges(str(p)) is False


def test_ensure_pdbqt_has_charges_no_atoms(tmp_path):
    p = tmp_path / "e.pdbqt"
    _write(p, ["REMARK nothing"])
    assert ch._ensure_pdbqt_has_charges(str(p)) is False


def test_pdbqt_has_atoms(tmp_path):
    p = tmp_path / "a.pdbqt"
    _write(p, [ATOM_NONZERO])
    assert ch._pdbqt_has_atoms(str(p)) is True
    q = tmp_path / "b.pdbqt"
    _write(q, ["REMARK x"])
    assert ch._pdbqt_has_atoms(str(q)) is False


def test_sanitize_receptor_pdbqt(tmp_path):
    p = tmp_path / "r.pdbqt"
    _write(p, [
        "ROOT",
        "ATOM      1  C   ALA A   1      0.0 0.0 0.0  0.0 0.0    -0.1 C",
        "BRANCH   1 2",
        "HETATM    2  O   HOH A   2      0.0 0.0 0.0  0.0 0.0    -0.2 O",
        "ENDBRANCH   1 2",
        "TER",
        "TORSDOF 3",
        "WEIRDLINE foo",
    ])
    ch._sanitize_receptor_pdbqt(str(p))
    out = p.read_text()
    assert "ROOT" not in out
    assert "BRANCH" not in out
    assert "ENDBRANCH" not in out
    assert "TORSDOF" not in out
    assert "ATOM" in out
    assert "HETATM" in out
    assert "TER" in out
    assert "REMARK WEIRDLINE" in out


def test_fix_pdbqt_charges_success(tmp_path, monkeypatch):
    out = tmp_path / "out.pdbqt"
    src = tmp_path / "src.pdb"
    src.write_text("REMARK source")

    def fake_run(cmd, capture=False, **kwargs):
        # Simulate OpenBabel producing a charged PDBQT
        out.write_text(
            "ATOM      1  C   ALA A   1      0.0 0.0 0.0  0.0 0.0    -0.250 C\n")
        return ("", "", 0)

    monkeypatch.setattr(utils, "run", fake_run)
    tools.configure(obabel="obabel")
    ch._fix_pdbqt_charges(str(out), str(src))
    assert ch._ensure_pdbqt_has_charges(str(out)) is True


def test_fix_pdbqt_charges_failure_raises(tmp_path, monkeypatch):
    out = tmp_path / "out.pdbqt"
    src = tmp_path / "src.pdb"
    src.write_text("REMARK source")

    def fake_run(cmd, capture=False, **kwargs):
        raise RuntimeError("obabel missing")

    monkeypatch.setattr(utils, "run", fake_run)
    tools.configure(obabel="obabel")
    with pytest.raises(RuntimeError):
        ch._fix_pdbqt_charges(str(out), str(src))


def test_new_name_and_backwards_compatible_alias():
    assert ch._ensure_pdbqt_has_charges is ch._pdbqt_has_nonzero_charges
    assert ch._ensure_pdbqt_has_charges is ch._pdbqt_has_nonzero_charges


ATOM_3D = "ATOM      1  C   LIG A   1       1.000   2.000   3.000  0.00  0.00    +0.000 C \n"
ATOM_2D = "ATOM      1  C   LIG A   1       1.000   2.000   0.000  0.00  0.00    +0.000 C \n"


def _capture_run(monkeypatch, write=None):
    calls = []

    def fake_run(cmd, capture=False, **kwargs):
        calls.append(cmd)
        target = cmd[cmd.index("-O") + 1]
        with open(target, "w", encoding="utf-8") as f:
            f.write(write or "ATOM      1  C   LIG A   1       1.000   2.000   3.000  0.00  0.00    -0.250 C \n")
        return ("", "", 0)

    monkeypatch.setattr(utils, "run", fake_run)
    tools.configure(obabel="obabel")
    return calls


def test_fix_charges_rebuilds_3d_pdbqt_in_place(tmp_path, monkeypatch):
    pdbqt = tmp_path / "l.pdbqt"
    pdbqt.write_text(ATOM_3D + ATOM_3D.replace("1.000", "2.000", 1))
    src = tmp_path / "l.sdf"
    src.write_text("2D source that must not be used\n")
    calls = _capture_run(monkeypatch)
    ch._fix_pdbqt_charges(str(pdbqt), str(src))
    cmd = calls[0]
    assert cmd[1:3] == ["-ipdbqt", str(pdbqt)]   # geometry from the 3D PDBQT
    assert "--gen3d" not in cmd
    assert "-0.250" in pdbqt.read_text()
    assert not [p for p in os.listdir(tmp_path) if ".tmp" in p]


def test_fix_charges_generates_3d_for_flat_source(tmp_path, monkeypatch):
    pdbqt = tmp_path / "l.pdbqt"
    pdbqt.write_text(ATOM_2D)                      # flat first conversion
    src = tmp_path / "l.sdf"
    src.write_text("x\n  RDKit          2D\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
                   "    0.0000    0.0000    0.0000 C   0  0\n"
                   "    1.5000    0.0000    0.0000 O   0  0\n"
                   "  1  2  1  0\nM  END\n$$$$\n")
    calls = _capture_run(monkeypatch)
    ch._fix_pdbqt_charges(str(pdbqt), str(src))
    assert calls[0][1:3] == ["-isdf", str(src)]
    assert "--gen3d" in calls[0]


def test_fix_charges_3d_source_without_gen3d(tmp_path, monkeypatch):
    src = tmp_path / "l.mol2"
    src.write_text("@<TRIPOS>MOLECULE\nx\n@<TRIPOS>ATOM\n"
                   "      1 C1   0.0 0.0 1.2 C.3 1 MOL 0.0\n"
                   "      2 O1   1.4 0.0 0.0 O.3 1 MOL 0.0\n@<TRIPOS>BOND\n")
    calls = _capture_run(monkeypatch)
    ch._fix_pdbqt_charges(str(tmp_path / "missing_yet.pdbqt"), str(src))
    assert calls[0][1] == "-imol2" and "--gen3d" not in calls[0]


def test_fix_charges_accepts_all_zero_hydrocarbon(tmp_path, monkeypatch, caplog):
    pdbqt = tmp_path / "benzene.pdbqt"
    pdbqt.write_text(ATOM_3D)
    _capture_run(monkeypatch, write=ATOM_3D)       # Gasteiger gives 0.000 after H merge
    with caplog.at_level("WARNING"):
        ch._fix_pdbqt_charges(str(pdbqt), str(pdbqt))
    assert "All Gasteiger charges are zero" in caplog.text


@pytest.mark.skipif(tools.settings().obabel is None, reason="obabel not available")
def test_fix_charges_real_benzene_keeps_coordinates(tmp_path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    p = AllChem.ETKDGv3()
    p.randomSeed = 42
    AllChem.EmbedMolecule(m, p)
    src = str(tmp_path / "benzene.mol")
    Chem.MolToMolFile(m, src)
    out = str(tmp_path / "benzene.pdbqt")
    utils.run([tools.settings().obabel, src, "-opdbqt", "-O", out])
    before = ch._read_coordinates(out)
    ch._fix_pdbqt_charges(out, src)               # must not raise for a hydrocarbon
    after = ch._read_coordinates(out)
    assert ch._has_3d_coordinates(out)
    assert sorted(before) == sorted(after)        # heavy-atom geometry kept
