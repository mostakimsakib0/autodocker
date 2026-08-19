import pytest

import runner


VINA_STDOUT = """mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1        -7.5       0.000      0.000
   2        -6.8       1.213      2.147
   3        -5.9       2.334      3.981
"""


def test_extract_score_mode_table():
    assert runner._extract_score(VINA_STDOUT) == pytest.approx(-7.5)


def test_extract_score_vina_result_line():
    out = "REMARK VINA RESULT:    -9.2   0.000   0.000\n"
    assert runner._extract_score(out) == pytest.approx(-9.2)


def test_extract_score_empty():
    with pytest.raises(ValueError):
        runner._extract_score("")


def test_extract_score_no_valid():
    with pytest.raises(ValueError, match="No valid affinity"):
        runner._extract_score("garbage output without scores\n")


def test_parse_vina_modes():
    modes = runner._parse_vina_modes(VINA_STDOUT)
    assert len(modes) == 3
    assert modes[0]["mode"] == 1
    assert modes[0]["affinity"] == pytest.approx(-7.5)
    assert modes[1]["rmsd_lb"] == pytest.approx(1.213)


def test_parse_vina_modes_ignores_garbage():
    assert runner._parse_vina_modes("not a mode row\n") == []


def test_score_bounds_accepted():
    for s in (-1.5, -9.2, -19.5):
        assert -20 < s < -1


def test_score_bounds_rejected():
    for s in (-0.5, 0.0, 5.0, -25.0, -99.0):
        assert not (-20 < s < -1)


def test_metrics_modes_and_simscore():
    metrics = runner._calculate_metrics("lig.pdbqt", "out.pdbqt",
                                        VINA_STDOUT, expected_modes=9)
    assert metrics["binding_modes"] == 3
    assert 0.0 <= metrics["simscore"] <= 1.0
    assert 0 < metrics["mode_coverage"] <= 1.0