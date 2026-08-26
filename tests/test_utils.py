import os

import pytest

import runner
import vspipeline.utils as u


def test_run_capture_success():
    out, err, code = u.run(["echo", "hello"], capture=True)
    assert code == 0
    assert out.strip() == "hello"
    assert err == ""


def test_run_no_capture(tmp_path):
    marker = tmp_path / "ok.txt"
    u.run(["sh", "-c", f"echo done > '{marker}'"], capture=False)
    assert marker.read_text().strip() == "done"


def test_run_failure_raises():
    with pytest.raises(Exception):
        u.run(["false"], capture=True)


def test_run_missing_binary_raises():
    with pytest.raises(FileNotFoundError):
        u.run(["this_binary_does_not_exist_xyz"], capture=True)


def _pdb(path, lines):
    path.write_text("\n".join(lines) + "\n")


def _atom(x, y, z, rec="ATOM"):
    # Pad record name so coordinates land in PDB columns 31-54 (0-indexed 30-53).
    prefix = rec + " " * (30 - len(rec))
    return f"{prefix}{x:>8.3f}{y:>8.3f}{z:>8.3f} X"


def test_parse_pdb_coords():
    import tempfile, pathlib
    d = tempfile.TemporaryDirectory()
    f = os.path.join(d.name, "x.pdb")
    _pdb(pathlib.Path(f), [
        _atom(1.0, 2.0, 3.0),
        _atom(4.0, 5.0, 6.0, rec="HETATM"),
        _atom(7.0, 8.0, 9.0),
    ])
    xs, ys, zs = u.parse_pdb_coords(f)
    assert xs == [1.0, 4.0, 7.0]
    assert ys == [2.0, 5.0, 8.0]
    assert zs == [3.0, 6.0, 9.0]


def test_parse_pdb_coords_empty_raises():
    import tempfile, pathlib
    d = tempfile.TemporaryDirectory()
    f = os.path.join(d.name, "empty.pdb")
    _pdb(pathlib.Path(f), ["REMARK nothing here"])
    with pytest.raises(ValueError):
        u.parse_pdb_coords(f)


def test_safe_ligand_id():
    assert u._safe_ligand_id("Aspirin 123!!") == "aspirin_123__"
    assert u._safe_ligand_id("  ") == "compound"
    assert u._safe_ligand_id("Caffeine-1.0") == "caffeine-1.0"


def test_get_chains(tmp_path):
    p = tmp_path / "c.pdb"
    _pdb(p, [
        "ATOM      1  N   ALA A   1      0.0 0.0 0.0",
        "ATOM      2  N   ALA B   1      0.0 0.0 0.0",
        "ATOM      3  N   ALA A   2      0.0 0.0 0.0",
    ])
    assert u.get_chains(str(p)) == ["A", "B"]


def test_parse_chain_selection_all():
    assert u._parse_chain_selection("all", ["A", "B", "C"]) == ["A", "B", "C"]
    assert u._parse_chain_selection("*", ["A", "B"]) == ["A", "B"]
    assert u._parse_chain_selection("", ["A", "B"]) == ["A", "B"]


def test_parse_chain_selection_valid():
    assert u._parse_chain_selection("A, C", ["A", "B", "C"]) == ["A", "C"]


def test_parse_chain_selection_invalid():
    with pytest.raises(ValueError):
        u._parse_chain_selection("Z", ["A", "B"])
