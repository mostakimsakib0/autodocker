import pytest

import runner
import vspipeline.charges as ch


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

    def fake_run(cmd, capture=False):
        # Simulate OpenBabel producing a charged PDBQT
        out.write_text(
            "ATOM      1  C   ALA A   1      0.0 0.0 0.0  0.0 0.0    -0.250 C\n")
        return ("", "", 0)

    monkeypatch.setattr(runner, "run", fake_run)
    monkeypatch.setattr(runner, "OBABEL", "obabel")
    ch._fix_pdbqt_charges(str(out), str(src))
    assert ch._ensure_pdbqt_has_charges(str(out)) is True


def test_fix_pdbqt_charges_failure_raises(tmp_path, monkeypatch):
    out = tmp_path / "out.pdbqt"
    src = tmp_path / "src.pdb"
    src.write_text("REMARK source")

    def fake_run(cmd, capture=False):
        raise RuntimeError("obabel missing")

    monkeypatch.setattr(runner, "run", fake_run)
    monkeypatch.setattr(runner, "OBABEL", "obabel")
    with pytest.raises(RuntimeError):
        ch._fix_pdbqt_charges(str(out), str(src))


def test_assign_simple_charges(tmp_path):
    p = tmp_path / "s.pdbqt"
    line = "ATOM" + "A" * 66 + " 0.000" + " " + " O"  # charge 0.0, atom type O at 77-78
    _write(p, [line])
    ch._assign_simple_charges(str(p))
    content = p.read_text()
    # Oxygen maps to -0.5; charge field should no longer be 0.000
    assert "0.000" not in content
    assert "-0.500" in content or "+0.500" in content or "0.500" in content
