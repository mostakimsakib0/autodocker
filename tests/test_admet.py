import os

import vspipeline.admet as admet

VALID_SDF = (
    "ethanol\n"
    "comment\n"
    "  3  2  0  0  0\n"
    "   1.00000   2.00000   3.00000\n"
    "   4.00000   5.00000   6.00000\n"
    "M  END\n"
)


def test_is_valid_sdf_valid(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text(VALID_SDF)
    ok, msg = admet._is_valid_sdf(str(p))
    assert ok and msg == ""


def test_is_valid_sdf_missing(tmp_path):
    ok, msg = admet._is_valid_sdf(str(tmp_path / "nope.sdf"))
    assert not ok


def test_is_valid_sdf_too_small(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text("x")
    ok, msg = admet._is_valid_sdf(str(p))
    assert not ok


def test_is_valid_sdf_no_mend(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text("mol\n\n  3  2  0  0  0\n   1.0   2.0   3.0\nno end here\n")
    ok, msg = admet._is_valid_sdf(str(p))
    assert not ok


def test_is_valid_sdf_error_content(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text("404 NotFound short")
    ok, msg = admet._is_valid_sdf(str(p))
    assert not ok


def test_check_lipinski_pass():
    ok, v = admet.ADMETFilter.check_lipinski(
        {"mw": 300, "logp": 2, "hba": 2, "hbd": 1})
    assert ok and v == []


def test_check_lipinski_fail():
    ok, v = admet.ADMETFilter.check_lipinski(
        {"mw": 600, "logp": 6, "hba": 12, "hbd": 6})
    assert not ok and len(v) == 4


def test_check_lipinski_missing():
    ok, v = admet.ADMETFilter.check_lipinski(None)
    assert not ok and "Invalid SDF" in v[0]


def test_parse_sdf_properties(monkeypatch, tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text(VALID_SDF)
    monkeypatch.setattr(
        admet, "_obabel_descriptors",
        lambda f: {"mw": 300.0, "logp": 2.0, "tpsa": 60.0,
                   "rotors": 3.0, "hbd": 1.0, "hba": 2.0})
    props = admet.ADMETFilter.parse_sdf_properties(str(p))
    assert props["mw"] == 300.0 and props["hba"] == 2.0


def test_obabel_descriptors_real(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text(VALID_SDF)
    result = admet._obabel_descriptors(str(p))
    assert result is None or isinstance(result, dict)
