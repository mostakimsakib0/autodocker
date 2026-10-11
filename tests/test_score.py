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
    from vspipeline import docking
    for s in (-1.5, -9.2, -19.5, -20.0):
        assert docking.affinity_in_window(s, -1.0)


def test_score_bounds_rejected():
    from vspipeline import docking
    for s in (-1.0, -0.5, 0.0, 5.0, -20.5, -25.0, -99.0):
        assert not docking.affinity_in_window(s, -1.0)


VINA_PDBQT = """MODEL 1
REMARK VINA RESULT:    -7.512      0.000      0.000
ENDMDL
MODEL 2
REMARK VINA RESULT:    -6.803      1.213      2.147
ENDMDL
MODEL 3
REMARK VINA RESULT:    -5.900      2.334      3.981
ENDMDL
"""


def test_metrics_modes_and_simscore(tmp_path):
    out = tmp_path / "out.pdbqt"
    out.write_text(VINA_PDBQT)
    metrics = runner._calculate_metrics(str(out), expected_modes=9)
    assert metrics["binding_modes"] == 3
    # mode 2 passes both cut-offs, mode 3 neither -> (1/2 + 1/2) / 2
    assert metrics["simscore"] == pytest.approx(0.5)
    assert metrics["mode_coverage"] == pytest.approx(3 / 9)


def test_simscore_single_mode_is_undefined(tmp_path):
    out = tmp_path / "out.pdbqt"
    out.write_text(VINA_PDBQT.split("MODEL 2")[0])
    assert runner._calculate_metrics(str(out))["simscore"] is None


def test_simscore_smina_output_is_undefined(tmp_path):
    out = tmp_path / "out.pdbqt"
    out.write_text("MODEL 1\nREMARK minimizedAffinity -7.1\nENDMDL\n"
                   "MODEL 2\nREMARK minimizedAffinity -6.1\nENDMDL\n")
    assert runner._calculate_metrics(str(out))["simscore"] is None


def test_simscore_thresholds_are_strict():
    from vspipeline import docking
    modes = [{"mode": 1, "rmsd_lb": 0.0, "rmsd_ub": 0.0},
             {"mode": 2, "rmsd_lb": docking.SIMSCORE_RMSD_LB_CUTOFF,
              "rmsd_ub": docking.SIMSCORE_RMSD_UB_CUTOFF}]
    assert docking._calculate_protocol_simscore(modes) == (0.0, 0.0, 0.0)
