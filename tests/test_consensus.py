import os
from dataclasses import dataclass

import pytest

import runner
import vspipeline.consensus as cons
from vspipeline import docking


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


def test_parse_smina_affinity_table():
    stdout = (
        "Some header\n"
        "mode |   affinity | dist from best mode\n"
        "  0 |  -7.500  |  0.000\n"
        "  1 |  -7.200  |  1.500\n"
    )
    assert cons._parse_smina_affinity(stdout) == -7.5


def test_parse_smina_affinity_real_output():
    # Verbatim SMINA 2020 output: units row, separator, whitespace columns.
    stdout = (
        "mode |   affinity | dist from best mode\n"
        "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n"
        "-----+------------+----------+----------\n"
        "1       -9.1       0.000      0.000    \n"
        "2       -8.4       2.136      2.961    \n"
        "Refine time 1.978\n"
        "Loop time 3.186\n"
    )
    assert cons._parse_smina_affinity(stdout) == -9.1


def test_smina_score_outside_window_is_rejected(fake_smina, tmp_path):
    rec = tmp_path / "r.pdbqt"
    rec.write_text("ATOM\n")
    lig = tmp_path / "l.pdbqt"
    lig.write_text("ATOM\n")
    scores = cons._run_smina_scoring(
        str(rec), [str(lig)], str(tmp_path / "out"), (1, 2, 3, 20, 20, 20),
        {"exhaustiveness": 8, "min_valid_affinity": -9.0}, num_processes=1)
    assert scores == {"l": None}  # -8.1 >= -9.0, so outside the window


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
            with open(cmd[out_idx], "w", encoding="utf-8") as fh:
                fh.write("MODEL        1\n"
                         "REMARK minimizedAffinity -8.12345\n"
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
    assert score == -8.12345   # PDBQT value, not the rounded stdout -8.1
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
    assert all(abs(v + 8.12345) < 1e-9 for v in scores.values())


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


def test_smina_docking_writes_log(fake_smina, tmp_path):
    rec = tmp_path / "r.pdbqt"
    rec.write_text("ATOM\n")
    lig = tmp_path / "l.pdbqt"
    lig.write_text("ATOM\n")
    cons._run_smina_scoring(str(rec), [str(lig)], str(tmp_path / "out"),
                            (1, 2, 3, 20, 20, 20), {"exhaustiveness": 8}, num_processes=1)
    log = (tmp_path / "out" / "docked" / "l_smina.log").read_text()
    assert log.startswith("# Command: smina -r ")


VINA_OUT = """MODEL 1
REMARK VINA RESULT:    -7.512      0.000      0.000
ROOT
ATOM      1  C   LIG A   1       0.000   0.000   0.000  0.00  0.00    +0.000 C
ENDROOT
TORSDOF 0
ENDMDL
MODEL 2
REMARK VINA RESULT:    -6.803      1.213      2.147
ROOT
ATOM      1  C   LIG A   1       1.000   0.000   0.000  0.00  0.00    +0.000 C
ENDROOT
TORSDOF 0
ENDMDL
"""


@pytest.fixture
def fake_rescore(monkeypatch):
    calls = []

    class R:
        def __init__(self, rc=0, out="", err=""):
            self.returncode, self.stdout, self.stderr = rc, out, err

    def fake_run(cmd, capture_output=False, text=False, timeout=None):
        calls.append(cmd)
        if "__autodocker_list__" in cmd:
            return R(1, "", "Usage error: Invalid builtin scoring function: x. Options are:\n"
                             "ad4_scoring\ndefault\ndkoes_scoring\nvina\nvinardo\n")
        lig = cmd[cmd.index("-l") + 1]
        if "fail" in os.path.basename(lig):
            return R(1, "", "boom")
        with open(cmd[cmd.index("-o") + 1], "w", encoding="utf-8") as fh:
            fh.write("MODEL 1\nREMARK minimizedAffinity -5.43219\n"
                     "REMARK VINA RESULT:    -7.512      0.000      0.000\nENDMDL\n")
        return R(0, "Affinity: -5.43219 (kcal/mol)\n")

    monkeypatch.setattr(cons.subprocess, "run", fake_run)
    return calls


def _dock_outputs(tmp_path, names):
    d = tmp_path / "out" / "docked"
    d.mkdir(parents=True, exist_ok=True)
    for n in names:
        (d / f"{n}_out.pdbqt").write_text(VINA_OUT)
    return str(tmp_path / "out")


def test_consensus_rescores_vina_top_pose_with_vinardo(fake_rescore, tmp_path):
    outdir = _dock_outputs(tmp_path, ["a", "b"])
    results = [("a", -7.512), ("b", -6.0), ("c", None)]
    scores, scoring = cons.rescore_with_smina(
        "rec.pdbqt", results, outdir, smina="smina", timeout=10)
    assert scoring == "vinardo"
    assert scores == {"a": -5.43219, "b": -5.43219, "c": None}
    cmd = [c for c in fake_rescore if "__autodocker_list__" not in c][0]
    assert cmd[:2] == ["smina", "-r"]
    assert "--score_only" in cmd
    assert cmd[cmd.index("--scoring") + 1] == "vinardo"
    assert "--exhaustiveness" not in cmd     # no re-docking
    pose = open(cmd[cmd.index("-l") + 1], encoding="utf-8").read()
    assert "1.000" not in pose and "MODEL" not in pose   # first model only
    log = open(os.path.join(outdir, "docked", "a_smina.log"), encoding="utf-8").read()
    assert log.startswith("# Command: smina -r rec.pdbqt")


def test_rescoring_falls_back_to_dkoes_when_vinardo_missing(monkeypatch):
    monkeypatch.setitem(cons._SMINA_SCORING_CACHE, "old_smina", ("default", "dkoes_scoring"))
    assert cons.select_rescoring_function("old_smina", "r", "l", 5) == "dkoes_scoring"


def test_rescoring_flex_residues_passed(fake_rescore, tmp_path):
    outdir = _dock_outputs(tmp_path, [])
    flex_out = VINA_OUT.replace("TORSDOF 0\nENDMDL", "TORSDOF 0\nBEGIN_RES LYS A  10\n"
                                "ATOM      2  NZ  LYS A  10       2.000   0.000   0.000  0.00  0.00    -0.300 N\n"
                                "END_RES LYS A  10\nENDMDL", 1)
    (tmp_path / "out" / "docked" / "f_out.pdbqt").write_text(flex_out)
    scores, _ = cons.rescore_with_smina("rec.pdbqt", [("f", -7.0)], outdir,
                                        flex_file="flex.pdbqt", scoring="vinardo",
                                        smina="smina", timeout=10)
    cmd = fake_rescore[-1]
    flex = open(cmd[cmd.index("--flex") + 1], encoding="utf-8").read()
    assert flex.startswith("BEGIN_RES") and "NZ  LYS" in flex
    assert "NZ  LYS" not in open(cmd[cmd.index("-l") + 1], encoding="utf-8").read()
    assert scores["f"] == -5.43219


def test_consensus_outputs_summary_json_and_no_agreement_column(tmp_path):
    import csv
    import json
    results = [("a", -9.0), ("b", -8.0), ("c", -7.0), ("d", None)]
    smina = {"a": -6.0, "b": -5.0, "c": -4.0, "d": None}
    agreement = cons.write_consensus_outputs(str(tmp_path), results, smina, "vinardo")
    assert agreement["n"] == 3 and abs(agreement["rho"] - 1.0) < 1e-9
    with open(tmp_path / "consensus_ranking.csv", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[0] == ["Ligand", "Vina_Affinity", "SMINA_Affinity", "Consensus"]
    assert not any(r and r[0].startswith("#") for r in rows)
    assert rows[-1] == ["d", "FAILED", "FAILED", "NA"]
    summary = json.loads((tmp_path / "consensus_summary.json").read_text())
    assert summary["n"] == 3 and summary["metric"] == "spearman_rho"
    assert "vinardo" in summary["scorer_b"]


@dataclass(frozen=True)
class _Job:
    name: str
    marks: str


def _dying_worker(job):
    """Module level (picklable). Records each attempt; "boom" kills its
    worker process, like the out-of-memory killer would."""
    with open(os.path.join(job.marks, job.name), "a", encoding="utf-8") as f:
        f.write("x")
    if job.name == "boom":
        os._exit(1)
    return job.name, -1.0, None


def test_dead_worker_fails_only_its_ligand(tmp_path):
    """A worker process that dies breaks the pool instead of hanging it: the
    unfinished ligands run again in fresh processes, the one that kills its
    worker again is FAILED, every other ligand gets its result."""
    marks = str(tmp_path)
    jobs = [_Job(n, marks) for n in ("a", "b", "boom", "c", "d")]
    out = cons._collect_with_retry(_dying_worker, jobs, 2, "t")
    assert out["boom"] == (None, docking.WORKER_DIED_ERROR)
    assert {n: out[n] for n in "abcd"} == {n: (-1.0, None) for n in "abcd"}
    assert (tmp_path / "boom").read_text() == "xx"     # retried once, not forever


def test_consensus_probe_missing_pose_fails_only_that_ligand(fake_rescore, tmp_path):
    """The scoring-function probe used the first ligand's pose outside the
    per-ligand error handling, so one missing pose crashed --consensus."""
    outdir = _dock_outputs(tmp_path, ["b"])
    results = [("a", -8.0), ("b", -7.512)]   # a's pose file is missing
    scores, scoring = cons.rescore_with_smina(
        "rec.pdbqt", results, outdir, smina="smina", timeout=10)
    assert scoring == "vinardo"
    assert scores["a"] is None and scores["b"] is not None


def test_smina_rescoring_nan_is_failed(tmp_path, monkeypatch):
    """A non-finite SMINA score is a failure (as in the Vina path), never a
    number in consensus_ranking.csv."""
    smina = tmp_path / "smina"
    smina.write_text("#!/bin/bash\nwhile [ $# -gt 0 ]; do [ \"$1\" = \"-o\" ] && out=$2; "
                     "shift; done\nprintf 'MODEL 1\\nREMARK minimizedAffinity nan\\n"
                     "ATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  0.00"
                     "    +0.000 C\\nENDMDL\\n' > \"$out\"\n")
    smina.chmod(0o755)
    outdir = _dock_outputs(tmp_path, ["a"])
    scores, _ = cons.rescore_with_smina("rec.pdbqt", [("a", -7.512)], outdir,
                                        scoring="vinardo", smina=str(smina), timeout=10)
    assert scores == {"a": None}
    log = (tmp_path / "out" / "docked" / "a_smina.log").read_text()
    assert "Non-finite score nan" in log
    with pytest.raises(RuntimeError, match="Non-finite"):
        cons.run_smina_docking("rec.pdbqt", "lig.pdbqt", str(tmp_path / "o.pdbqt"),
                               0, 0, 0, 20, 20, 20, smina=str(smina), timeout=10)


def test_consensus_ranking_out_of_window_keeps_its_score(tmp_path):
    """OUT_OF_WINDOW Vina results were written as FAILED in
    consensus_ranking.csv; they keep their score, as in ranking.csv."""
    import csv
    results = [("a", -9.0), ("b", -8.0), ("oow", None), ("d", None)]
    smina = {"a": -6.0, "b": -5.0}
    metrics = {"a": {"status": "OK"}, "b": {"status": "OK"}, "d": {"status": "FAILED"},
               "oow": {"status": "OUT_OF_WINDOW", "raw_score": -0.5}}
    cons.write_consensus_outputs(str(tmp_path), results, smina, "vinardo",
                                 metrics_dict=metrics)
    with open(tmp_path / "consensus_ranking.csv", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[1:] == [["a", "-9.0", "-6.0", "-7.5"], ["b", "-8.0", "-5.0", "-6.5"],
                        ["d", "FAILED", "FAILED", "NA"],
                        ["oow", "-0.5", "NA", "OUT_OF_WINDOW"]]
