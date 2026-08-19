import runner


def test_is_valid_sdf_ok(small_sdf):
    ok, err = runner._is_valid_sdf(small_sdf)
    assert ok is True, err


def test_is_valid_sdf_missing_file(tmp_path):
    ok, err = runner._is_valid_sdf(str(tmp_path / "nope.sdf"))
    assert ok is False
    assert "not found" in err


def test_is_valid_sdf_html_like(tmp_path):
    p = tmp_path / "bad.sdf"
    p.write_text("Status: 404 Not Found\n")
    ok, err = runner._is_valid_sdf(str(p))
    assert ok is False


def test_is_valid_sdf_missing_end_marker(tmp_path):
    p = tmp_path / "bad.sdf"
    p.write_text(
        "X" * 150 + "\n" + " " * 50 + "0.0   1.0   2.0\n" + "M  ENDish\n")
    ok, err = runner._is_valid_sdf(str(p))
    assert ok is False


def test_pdbqt_has_atoms(charged_pdbqt, empty_pdbqt):
    assert runner._pdbqt_has_atoms(charged_pdbqt) is True
    assert runner._pdbqt_has_atoms(empty_pdbqt) is False


def test_ensure_pdbqt_has_charges(charged_pdbqt, zero_charge_pdbqt):
    assert runner._ensure_pdbqt_has_charges(charged_pdbqt) is True
    assert runner._ensure_pdbqt_has_charges(zero_charge_pdbqt) is False


def test_ensure_pdbqt_has_charges_missing_file(tmp_path):
    assert runner._ensure_pdbqt_has_charges(
        str(tmp_path / "missing.pdbqt")) is False


def test_extract_pdbqt_charge():
    line = ("ATOM      1  C   UNK A   1       0.000   0.000   0.000 "
            "1.00 20.00      C   0.120")
    assert runner._extract_pdbqt_charge(line) == 0.120
    assert runner._extract_pdbqt_charge("END") is None