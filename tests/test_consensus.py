import os

import pytest

import runner
import vspipeline.consensus as cons


def test_rankdata_no_ties():
    assert cons._rankdata([3.0, 1.0, 2.0]) == [3.0, 1.0, 2.0]


def test_rankdata_ties():
    assert cons._rankdata([1.0, 2.0, 2.0, 3.0]) == [1.0, 2.5, 2.5, 4.0]


def test_pearson_perfect():
    assert abs(cons._pearson([1, 2, 3, 4], [10, 20, 30, 40]) - 1.0) < 1e-9


def test_pearson_anti():
    assert abs(cons._pearson([1, 2, 3, 4], [40, 30, 20, 10]) + 1.0) < 1e-9


def test_pearson_constant():
    assert cons._pearson([1, 1, 1], [1, 2, 3]) == 0.0


def test_spearman_perfect():
    rho, n = cons.scorer_agreement_spearman([1, 2, 3, 4], [10, 20, 30, 40])
    assert n == 4 and abs(rho - 1.0) < 1e-9


def test_spearman_anti():
    rho, n = cons.scorer_agreement_spearman([1, 2, 3, 4], [40, 30, 20, 10])
    assert n == 4 and abs(rho + 1.0) < 1e-9


def test_spearman_ties():
    # monotonic-ish with a tie -> partial negative correlation
    rho, n = cons.scorer_agreement_spearman([1, 2, 2, 3], [5, 9, 9, 1])
    assert n == 4 and rho is not None


def test_spearman_undefined_const():
    assert cons.scorer_agreement_spearman([1, 1, 1, 1], [1, 2, 3, 4]) == (None, 4)


def test_spearman_undefined_single():
    assert cons.scorer_agreement_spearman([1], [2]) == (None, 1)


def test_consensus_rank_heuristic():
    row = cons.consensus_rank(-7.0, -7.5)
    assert row["consensus"] == -7.25
    assert abs(row["agreement"] - 95.0) < 1e-9
    assert row["smina"] == -7.5


def test_consensus_rank_no_smina():
    row = cons.consensus_rank(-7.0, None)
    assert row["smina"] is None
    assert row["agreement"] == 0
    assert row["consensus"] == -7.0


def test_parse_smina_affinity_table():
    stdout = (
        "Some header\n"
        "mode |   affinity | dist from best mode\n"
        "  0 |  -7.500  |  0.000\n"
        "  1 |  -7.200  |  1.500\n"
    )
    assert cons._parse_smina_affinity(stdout) == -7.5


def test_parse_smina_affinity_fallback():
    stdout = "1 -5.3 something\n2 -4.1 other\n"
    assert cons._parse_smina_affinity(stdout) == -5.3


def test_parse_smina_affinity_none():
    assert cons._parse_smina_affinity("no numbers here") is None


SMINA_TABLE = (
    "mode |   affinity | dist from best mode\n"
    "  0 |  -8.100  |  0.000\n"
)


@pytest.fixture
def fake_smina(monkeypatch):
    monkeypatch.setattr(runner, "SMINA", "smina")
    calls = {}

    class FakeResult:
        returncode = 0
        stdout = SMINA_TABLE
        stderr = ""

    def fake_run(cmd, capture_output=False, text=False, timeout=None):
        calls["cmd"] = cmd
        # SMINA would write the output PDBQT; create it so the check passes.
        if "-o" in cmd:
            out_idx = cmd.index("-o") + 1
            with open(cmd[out_idx], "w") as fh:
                fh.write("MODEL        1\n"
                         "ATOM      1  C   ALA A   1      0.0 0.0 0.0\n"
                         "ENDMDL\n")
        return FakeResult()

    monkeypatch.setattr(cons.subprocess, "run", fake_run)
    return calls


def test_run_smina_docking(fake_smina, tmp_path):
    rec = tmp_path / "r.pdbqt"
    lig = tmp_path / "l.pdbqt"
    out = tmp_path / "o.pdbqt"
    rec.write_text("ATOM\n")
    lig.write_text("ATOM\n")
    score = cons.run_smina_docking(
        str(rec), str(lig), str(out), 1, 2, 3, 20, 20, 20,
        exhaustiveness=16, seed=42)
    assert score == -8.1
    assert "--exhaustiveness" in fake_smina["cmd"]
    assert "16" in fake_smina["cmd"]
    assert "--seed" in fake_smina["cmd"]
    assert "42" in fake_smina["cmd"]


def test_run_smina_docking_no_binary(monkeypatch):
    monkeypatch.setattr(runner, "SMINA", None)
    assert cons.run_smina_docking("r", "l", "o", 1, 2, 3, 4, 5, 6) is None


def test_run_smina_scoring(fake_smina, tmp_path):
    rec = tmp_path / "r.pdbqt"
    rec.write_text("ATOM\n")
    ligs = []
    for i in range(3):
        p = tmp_path / f"lig{i}.pdbqt"
        p.write_text("ATOM\n")
        ligs.append(str(p))
    scores = cons._run_smina_scoring(
        str(rec), ligs, str(tmp_path / "out"), (1, 2, 3, 20, 20, 20),
        {"exhaustiveness": 8, "seed": 1}, num_processes=1)
    assert len(scores) == 3
    assert all(abs(v + 8.1) < 1e-9 for v in scores.values())


def test_run_smina_scoring_parallel(fake_smina, tmp_path):
    rec = tmp_path / "r.pdbqt"
    rec.write_text("ATOM\n")
    ligs = []
    for i in range(4):
        p = tmp_path / f"lig{i}.pdbqt"
        p.write_text("ATOM\n")
        ligs.append(str(p))
    scores = cons._run_smina_scoring(
        str(rec), ligs, str(tmp_path / "out"), (1, 2, 3, 20, 20, 20),
        {"exhaustiveness": 8}, num_processes=2)
    assert len(scores) == 4
