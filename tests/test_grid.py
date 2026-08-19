import pytest

import runner


def _write_grid(path, text):
    with open(path, "w") as f:
        f.write(text)
    return str(path)


def test_parse_grid_config_ok(tmp_path):
    p = _write_grid(tmp_path / "grid.conf",
                    "center_x = 1.5\ncenter_y = 2\ncenter_z = 3\n"
                    "size_x = 24\nsize_y = 24\nsize_z = 24\n")
    cfg = runner._parse_grid_config(p)
    assert cfg["center_x"] == 1.5
    assert cfg["size_z"] == 24.0


def test_parse_grid_config_missing_params(tmp_path):
    p = _write_grid(tmp_path / "grid.conf", "center_x = 1\ncenter_y = 2\n")
    with pytest.raises(ValueError, match="missing required"):
        runner._parse_grid_config(p)


def test_parse_grid_config_bad_size(tmp_path):
    p = _write_grid(tmp_path / "grid.conf",
                    "center_x=0\ncenter_y=0\ncenter_z=0\n"
                    "size_x=-5\nsize_y=24\nsize_z=24\n")
    with pytest.raises(ValueError, match="size"):
        runner._parse_grid_config(p)


def test_parse_grid_config_too_big(tmp_path):
    p = _write_grid(tmp_path / "grid.conf",
                    "center_x=0\ncenter_y=0\ncenter_z=0\n"
                    "size_x=200\nsize_y=24\nsize_z=24\n")
    with pytest.raises(ValueError, match="80"):
        runner._parse_grid_config(p)


def _prep(protein_pdb, outdir):
    return runner.ProteinPreparation(protein_pdb, outdir)


def test_parse_pocket_selection_default(protein_pdb, outdir):
    pockets = [{"number": 3}, {"number": 7}]
    assert _prep(protein_pdb, outdir)._parse_pocket_selection(
        None, pockets) == [3]


def test_parse_pocket_selection_explicit(protein_pdb, outdir):
    pockets = [{"number": 3}, {"number": 7}]
    assert _prep(protein_pdb, outdir)._parse_pocket_selection(
        "7,3", pockets) == [7, 3]


def test_parse_pocket_selection_invalid(protein_pdb, outdir):
    pockets = [{"number": 3}]
    with pytest.raises(ValueError, match="not found"):
        _prep(protein_pdb, outdir)._parse_pocket_selection("9", pockets)


def test_parse_pocket_selection_non_numeric(protein_pdb, outdir):
    pockets = [{"number": 3}]
    with pytest.raises(ValueError, match="Invalid pocket"):
        _prep(protein_pdb, outdir)._parse_pocket_selection("abc", pockets)


def test_centroid_fallback_grid(protein_pdb, outdir, monkeypatch):
    prep = runner.ProteinPreparation(protein_pdb, outdir)
    monkeypatch.setattr(runner.ProteinPreparation, "detect_pocket", None)
    cx, cy, cz, sx, sy, sz = prep._protein_centroid_grid(padding=6.0)
    # protein_pdb has two chains, x in [0,9], y in [0,6], z fixed at 5.
    assert (4.5, 3.0, 5.0) == pytest.approx((cx, cy, cz))
    assert sx >= 24.0 and sy >= 24.0 and sz >= 24.0


def test_write_grid_produces_parsable_conf(protein_pdb, outdir):
    prep = runner.ProteinPreparation(protein_pdb, outdir)
    prep.write_grid(1.0, 2.0, 3.0, 24.0, 25.0, 26.0)
    cfg = runner._parse_grid_config(prep.grid_conf)
    assert cfg["center_x"] == 1.0
    assert cfg["size_z"] == 26.0


def test_write_grid_box_script(protein_pdb, outdir):
    prep = runner.ProteinPreparation(protein_pdb, outdir)
    prep.write_grid_box_script(0, 0, 0, 24, 24, 24)
    with open(prep.grid_box_script) as f:
        assert "VERTEX" in f.read()


def test_parse_grid_triplet_ok():
    assert runner._parse_grid_triplet("12.5,4,6.75", "grid-center") == (
        12.5, 4.0, 6.75)


def test_parse_grid_triplet_bad_arity():
    with pytest.raises(ValueError, match="three comma-separated"):
        runner._parse_grid_triplet("1,2", "grid-center")


def test_parse_grid_triplet_non_numeric():
    with pytest.raises(ValueError, match="must be floats"):
        runner._parse_grid_triplet("1,abc,3", "grid-size")


def test_parse_grid_triplet_nan():
    with pytest.raises(ValueError, match="finite"):
        runner._parse_grid_triplet("1,nan,3", "grid-center")