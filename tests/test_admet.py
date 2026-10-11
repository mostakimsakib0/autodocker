import os

import pytest

import runner
import vspipeline.admet as admet

DATA = os.path.join(os.path.dirname(__file__), "data", "runner")

VALID_SDF = (
    "ethanol\n"
    "  test\n"
    "\n"
    "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
    "   -0.8883    0.1670   -0.0273 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    0.4658   -0.5116   -0.0368 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "    1.4311    0.3229    0.5867 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
    "  1  2  1  0\n"
    "  2  3  1  0\n"
    "M  END\n"
)


def test_is_valid_sdf_valid(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text(VALID_SDF)
    ok, msg = admet._is_valid_sdf(str(p))
    assert ok and msg == ""


def test_is_valid_sdf_v3000():
    ok, msg = admet._is_valid_sdf(os.path.join(DATA, "acetaminophen_v3000_3d.sdf"))
    assert ok, msg


def test_is_valid_sdf_missing(tmp_path):
    ok, msg = admet._is_valid_sdf(str(tmp_path / "nope.sdf"))
    assert not ok


def test_is_valid_sdf_garbage_and_empty(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text("x")
    ok, msg = admet._is_valid_sdf(str(p))
    assert not ok and "M  END" in msg
    p.write_text("")
    ok, msg = admet._is_valid_sdf(str(p))
    assert not ok and "empty" in msg


def test_is_valid_sdf_accepts_tiny_molecule(tmp_path):
    """L6: a complete molblock of under 100 bytes (methane) is valid; there
    is no byte threshold, only the parse checks."""
    p = tmp_path / "methane.sdf"
    p.write_text("CH4\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                 "    0.0000    0.0000    0.0000 C   0  0\nM  END\n$$$$\n")
    assert p.stat().st_size < 100
    ok, msg = admet._is_valid_sdf(str(p))
    assert ok, msg
    assert not hasattr(admet, "MIN_SDF_BYTES")


def test_is_valid_sdf_no_mend(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text(VALID_SDF.replace("M  END", "no end here") + "x" * 50)
    ok, msg = admet._is_valid_sdf(str(p))
    assert not ok and "M  END" in msg


def test_is_valid_sdf_error_page(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text("Status: 404\nMessage: NotFound\n" + "x" * 120)
    ok, msg = admet._is_valid_sdf(str(p))
    assert not ok


def test_is_valid_sdf_without_atom_block(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text("mol\n  test\n\n  garbage counts line\n" + "x" * 120 + "\nM  END\n")
    ok, msg = admet._is_valid_sdf(str(p))
    assert not ok and "coordinates" in msg


def test_sdf_coordinates_v2000_and_v3000(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text(VALID_SDF)
    coords = admet.sdf_coordinates(str(p))
    assert coords[2] == (1.4311, 0.3229, 0.5867)
    v3000 = admet.sdf_coordinates(os.path.join(DATA, "acetaminophen_v3000_3d.sdf"))
    assert len(v3000) == 20
    assert any(abs(z) > 0.1 for _, _, z in v3000)
    assert admet.sdf_coordinates(str(tmp_path / "missing.sdf")) is None


def test_check_lipinski_pass():
    ok, v = admet.LipinskiFilter.check_lipinski(
        {"mw": 300, "logp": 2, "hba": 2, "hbd": 1})
    assert ok and v == []


def test_check_lipinski_fail_messages_unchanged():
    ok, v = admet.LipinskiFilter.check_lipinski(
        {"mw": 600, "logp": 6, "hba": 12, "hbd": 6})
    assert not ok
    assert v == ["MW=600.0 (>500)", "LogP=6.00 (>5)", "HBD=6 (>5)", "HBA=12 (>10)"]


def test_check_lipinski_zero_violations_allowed():
    ok, v = admet.LipinskiFilter.check_lipinski({"mw": 501, "logp": 2, "hba": 2, "hbd": 1})
    assert not ok and len(v) == 1


def test_check_lipinski_custom_window():
    props = {"mw": 600, "logp": 6, "hba": 2, "hbd": 1}
    ok, v = admet.LipinskiFilter.check_lipinski(props, mw_max=700, logp_max=7)
    assert ok
    ok, v = admet.LipinskiFilter.check_lipinski(props, mw_min=650, mw_max=700, logp_max=7)
    assert not ok and v == ["MW=600.0 (<650)"]
    ok, v = admet.LipinskiFilter.check_lipinski(
        {"mw": 300, "logp": 2, "hba": 20, "hbd": 9}, hbd_max=None, hba_max=None)
    assert ok


def test_check_lipinski_missing():
    ok, v = admet.LipinskiFilter.check_lipinski(None)
    assert not ok and "Invalid SDF" in v[0]


def test_admet_filter_is_a_deprecated_alias():
    assert admet.ADMETFilter is admet.LipinskiFilter
    assert runner.ADMETFilter is runner.LipinskiFilter


def test_parse_sdf_properties(monkeypatch, tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text(VALID_SDF)
    monkeypatch.setattr(
        admet, "_obabel_descriptors",
        lambda f: {"mw": 300.0, "logp": 2.0, "tpsa": 60.0,
                   "rotors": 3.0, "hbd": 1.0, "hba": 2.0})
    props = admet.LipinskiFilter.parse_sdf_properties(str(p))
    assert props["mw"] == 300.0 and props["hba"] == 2.0
    assert props["properties_source"] == "obabel_descriptors"
    assert props["logp_model"] == admet.LOGP_MODEL_OBABEL


def test_obabel_descriptors_real(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text(VALID_SDF)
    result = admet._obabel_descriptors(str(p))
    if runner.OBABEL is None:
        assert result is None
    else:
        assert result["mw"] == pytest.approx(46.07, abs=0.05)
        assert result["hbd"] == 1.0 and result["hba"] == 1.0


def test_obabel_descriptors_single_threaded(monkeypatch, tmp_path):
    calls = []

    def fake_run(cmd, capture=False, extra_env=None, **kw):
        calls.append((cmd, extra_env))
        return "CCO\tethanol\t46.07\t-0.0014\t20.23\t0\t1\t1\n", "", 0

    monkeypatch.setattr(runner, "run", fake_run)
    p = tmp_path / "a.sdf"
    p.write_text(VALID_SDF)
    props = admet._obabel_descriptors(str(p))
    assert props == {"mw": 46.07, "logp": -0.0014, "tpsa": 20.23, "rotors": 0.0,
                     "hbd": 1.0, "hba": 1.0}
    assert calls[0][1] == {"OMP_NUM_THREADS": "1"}
    assert "HBA2" in calls[0][0][-1]


def test_small_valid_sdf_containing_404_is_accepted(tmp_path):
    p = tmp_path / "a.sdf"
    p.write_text("mol\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                 "    0.4040    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                 "M  END\n$$$$\n")
    ok, msg = admet._is_valid_sdf(str(p))
    assert ok, msg


def test_sdf_data_items_used_when_obabel_fails(tmp_path, monkeypatch):
    monkeypatch.setattr(admet, "_obabel_descriptors", lambda p: None)
    p = tmp_path / "a.sdf"
    p.write_text("mol\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n"
                 "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
                 "M  END\n"
                 "> <PUBCHEM_MOLECULAR_WEIGHT>\n180.16\n\n"
                 "> <PUBCHEM_XLOGP3>\n1.2\n\n"
                 "> <PUBCHEM_CACTVS_HBOND_DONOR>\n1\n\n"
                 "> <PUBCHEM_CACTVS_HBOND_ACCEPTOR>\n4\n\n$$$$\n")
    props = admet.LipinskiFilter.parse_sdf_properties(str(p))
    assert (props["mw"], props["logp"], props["hbd"], props["hba"]) == (180.16, 1.2, 1, 4)
    assert props["properties_source"] == "sdf_data_items"
    assert props["logp_model"] == admet.LOGP_MODEL_SDF
