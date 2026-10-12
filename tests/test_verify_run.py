"""Tests of scripts/verify_run.py (run verification, re-docking, comparison).

tests/data/verify/run_small is a real AutoDocker 2.0.0 run (3 ligands,
exhaustiveness 1, 3 modes, seed 3; Vina 1.2.7, Open Babel 3.2.1) with the
machine paths replaced by /work/... and /opt/...; the receptor files are left
out. The re-docking test makes its own tiny run with the tools on PATH.
"""
import importlib.util
import json
import os
import shutil
import subprocess
import sys

import pytest

from vspipeline import defaults, docking, outputs, results

ROOT = os.path.join(os.path.dirname(__file__), "..")
FIXTURE = os.path.join(os.path.dirname(__file__), "data", "verify", "run_small")


def _load(name):
    spec = importlib.util.spec_from_file_location(
        f"bench_{name}", os.path.join(ROOT, "scripts", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


vr = _load("verify_run")


def _run(argv, capsys):
    code = vr.main(list(argv) + ["--json"])
    doc = json.loads(capsys.readouterr().out)
    return code, doc


def _failed(doc):
    return {c["check"] for c in doc["checks"] if c["status"] == "fail"}


@pytest.fixture
def run_copy(tmp_path):
    dst = tmp_path / "run"
    shutil.copytree(FIXTURE, dst)
    return dst


def _edit(path, old, new, count=1):
    text = path.read_text(encoding="utf-8")
    assert old in text, f"{old!r} not in {path}"
    path.write_text(text.replace(old, new, count), encoding="utf-8")


# ---------------------------------------------------------------------------
# Names and rules mirrored from the pipeline
# ---------------------------------------------------------------------------

def test_mirrored_names_match_the_pipeline():
    assert vr.RANKING_COLUMNS == outputs.RANKING_COLUMNS
    assert vr.RUN_MANIFEST_JSON == outputs.RUN_MANIFEST_JSON
    assert vr.RECEPTOR_PDBQT == outputs.RECEPTOR_PDBQT
    assert vr.POCKET_SUMMARY_CSV == outputs.POCKET_SUMMARY_CSV
    assert vr.SMINA_POSE_SUFFIX == outputs.SMINA_POSE_SUFFIX
    assert vr.SMINA_LOG_SUFFIX == outputs.SMINA_LOG_SUFFIX
    assert set(vr.STATUSES) == {s.value for s in docking.Status}
    assert vr.STATUS_GROUP == {s.value: g for s, g in results._GROUP.items()}
    assert vr.AFFINITY_FLOOR == docking.AFFINITY_FLOOR
    assert vr.WEAK_AFFINITY_THRESHOLD == docking.WEAK_AFFINITY_THRESHOLD
    assert vr.DEFAULT_MIN_VALID_AFFINITY == defaults.MIN_VALID_AFFINITY


def test_expected_status_follows_the_pipeline():
    for score in (-25.0, -20.0, -7.0, -4.0, -3.9, -1.0, -0.5, 2.0):
        if docking.affinity_in_window(score, -1.0):
            want = docking.status_for_score(score).value
        else:
            want = "OUT_OF_WINDOW"
        assert vr.expected_status(score, -1.0) == want


# ---------------------------------------------------------------------------
# (a) consistency and tamper detection
# ---------------------------------------------------------------------------

def test_fixture_run_verifies(capsys):
    code, doc = _run([FIXTURE], capsys)
    assert code == 0 and doc["result"] == "PASS", _failed(doc)
    assert doc["autodocker_version"] == "2.0.0"
    checks = {c["check"]: c for c in doc["checks"]}
    assert checks["ligand[acetaminophen].score"]["detail"] == "-5.124"
    assert checks["ligand[butanoic_acid].status"]["status"] == "pass"  # WEAK row
    # The logged binary is not on this machine: skipped, not failed.
    assert checks["engine[/opt/vina/bin/vina].sha256"]["status"] == "skip"
    assert any("SHA-256 of the output files" in g for g in doc["manifest_gaps"])


@pytest.mark.parametrize("tamper, check", [
    # ranking.csv score no longer equals the pose file's REMARK VINA RESULT
    (lambda d: _edit(d / "ranking.csv", "-5.124", "-5.125"), "ligand[acetaminophen].score"),
    # pose file score edited
    (lambda d: _edit(d / "docked" / "phenylethanol_out.pdbqt",
                     "REMARK VINA RESULT:    -4.437", "REMARK VINA RESULT:    -4.937"),
     "ligand[phenylethanol].score"),
    # status inconsistent with the score
    (lambda d: _edit(d / "ranking.csv", "docked/butanoic_acid_out.pdbqt,WEAK",
                     "docked/butanoic_acid_out.pdbqt,OK"), "ligand[butanoic_acid].status"),
    # pose file removed
    (lambda d: os.remove(d / "docked" / "acetaminophen_out.pdbqt"),
     "ligand[acetaminophen].pose"),
    # engine log removed
    (lambda d: os.remove(d / "docked" / "acetaminophen_vina.log"),
     "ligand[acetaminophen].log"),
    # logged command used another seed than the manifest
    (lambda d: _edit(d / "docked" / "phenylethanol_vina.log", "--seed 3", "--seed 4"),
     "ligand[phenylethanol].log_command"),
    # manifest is not strict JSON
    (lambda d: _edit(d / "run_manifest.json", '"exit_code": 0', '"exit_code": NaN'),
     "manifest.strict_json"),
    # run did not finish
    (lambda d: _edit(d / "run_manifest.json", '"status": "completed"',
                     '"status": "interrupted"'), "manifest.status"),
    # manifest counts disagree with the ranking
    (lambda d: _edit(d / "run_manifest.json", '"docked": 3', '"docked": 4'),
     "manifest.counts"),
    # rows out of ranking order
    (lambda d: _swap_rows(d / "ranking.csv"), "ranking.order"),
])
def test_tampering_is_detected(run_copy, capsys, tamper, check):
    tamper(run_copy)
    code, doc = _run([str(run_copy)], capsys)
    assert code == 1 and doc["result"] == "FAIL"
    assert check in _failed(doc)


def _swap_rows(path):
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    lines[1], lines[2] = lines[2], lines[1]
    path.write_text("".join(lines), encoding="utf-8")


def test_missing_manifest_and_ranking_fail(tmp_path, capsys):
    code, doc = _run([str(tmp_path)], capsys)
    assert code == 1
    assert {"manifest.present", "ranking.present"} <= _failed(doc)


def test_orphan_pose_is_a_warning(run_copy, capsys):
    shutil.copy(run_copy / "docked" / "acetaminophen_out.pdbqt",
                run_copy / "docked" / "stale_out.pdbqt")
    code, doc = _run([str(run_copy)], capsys)
    assert code == 0
    warn = [c for c in doc["checks"] if c["check"] == "docked.orphans"]
    assert warn and warn[0]["status"] == "warn" and "stale_out.pdbqt" in warn[0]["detail"]


def test_engine_sha256_checked_against_manifest(tmp_path):
    binary = tmp_path / "vina"
    binary.write_bytes(b"not really vina")
    manifest = {"docking_engines": [{"path": str(binary), "sha256": "0" * 64}]}
    report = vr.Report()
    assert vr.check_binary(str(binary), manifest, report) == vr.FAIL
    manifest["docking_engines"][0]["sha256"] = vr._benchlib.file_sha256(str(binary))
    assert vr.check_binary(str(binary), manifest, report) == vr.PASS
    # a binary the manifest does not know: a warning in (a), a failure for --redock
    assert vr.check_binary("/elsewhere/vina", manifest, report) == vr.WARN
    assert vr.check_binary("/elsewhere/vina", manifest, report, required=True) == vr.FAIL


def test_relocate_command_maps_the_run_directory(tmp_path):
    cmd = ["/opt/vina", "--receptor", "/old/run/receptor.pdbqt",
           "--ligand", "/old/run/ligands/x.pdbqt", "--out", "/old/run/docked/x_out.pdbqt",
           "--seed", "3", "--center_x", "1.0"]
    (tmp_path / "ligands").mkdir()
    (tmp_path / "receptor.pdbqt").write_text("x", encoding="utf-8")
    new, missing = vr.relocate_command(cmd, str(tmp_path), "/tmp/new_out.pdbqt")
    assert new[2] == str(tmp_path / "receptor.pdbqt")
    assert new[4] == str(tmp_path / "ligands" / "x.pdbqt")
    assert new[6] == "/tmp/new_out.pdbqt"
    assert new[7:] == ["--seed", "3", "--center_x", "1.0"]
    assert missing == [f"--ligand {tmp_path / 'ligands' / 'x.pdbqt'}"]


def test_parse_log_takes_the_last_attempt(tmp_path):
    log = tmp_path / "x_vina.log"
    log.write_text("# Attempt: Vina (primary)\n# Command: /a/vina --out 'p q/x.pdbqt'\n"
                   "# Error: crashed\n\n# Attempt: QuickVina (fallback)\n"
                   "# Command: /b/qvina --out /r/docked/x_out.pdbqt\n"
                   "# Result: WEAK - Affinity: -3.200 kcal/mol\n", encoding="utf-8")
    cmd, result = vr.parse_log(str(log))
    assert cmd == ["/b/qvina", "--out", "/r/docked/x_out.pdbqt"]
    assert vr.parse_result_line(result) == ("WEAK", -3.2)
    assert vr.original_outdir(cmd) == "/r"


# ---------------------------------------------------------------------------
# (c) comparing two runs
# ---------------------------------------------------------------------------

def test_spearman():
    assert vr.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert vr.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    assert vr.spearman([1, 2, 2, 3], [1, 2, 3, 4]) == pytest.approx(0.9486833)
    assert vr.spearman([1], [1]) is None
    assert vr.spearman([1, 1], [1, 2]) is None


def _place(run_dir):
    """Write the run's own absolute path into its files, as a real run does."""
    for path in (run_dir / "run_manifest.json", run_dir / "docked" / "acetaminophen_vina.log"):
        _edit(path, "/work/run_small", str(run_dir), count=-1)


def test_compare_identical_runs(tmp_path, capsys):
    a, b = tmp_path / "a", tmp_path / "b"
    shutil.copytree(FIXTURE, a)
    shutil.copytree(FIXTURE, b)
    _place(a)
    _place(b)
    _edit(b / "run_manifest.json", '"finished_at": "', '"finished_at": "1999')
    _edit(b / "docked" / "acetaminophen_vina.log", "AutoDock", "autodock")  # logs skipped
    (b / ".docking.lock").write_text("123", encoding="utf-8")
    for d, vol in ((a, "101.5"), (b, "99.25")):  # fpocket's Monte-Carlo volume
        (d / "pocket_summary.csv").write_text(
            "Pocket,Score,Druggability_Score,Volume,Box_Volume,Could_Be_Druggable,File\n"
            f"1,0.5,0.4,{vol},4096.000,no,p1.pdb\n", encoding="utf-8")
    code, doc = _run([str(a), str(b)], capsys)
    assert code == 0, doc["compare"]["different"]
    cmp_ = doc["compare"]
    assert "run_manifest.json" in cmp_["identical"] and "pocket_summary.csv" in cmp_["identical"]
    assert "docked/acetaminophen_vina.log" not in cmp_["identical"] + cmp_["different"]
    assert cmp_["ranking"]["spearman_rho"] == pytest.approx(1.0)
    assert cmp_["ranking"]["max_abs_diff"] == 0


def test_compare_different_runs(tmp_path, capsys):
    a, b = tmp_path / "a", tmp_path / "b"
    shutil.copytree(FIXTURE, a)
    shutil.copytree(FIXTURE, b)
    # B: phenylethanol scores better than acetaminophen; pose and ranking agree.
    _edit(b / "docked" / "phenylethanol_out.pdbqt",
          "REMARK VINA RESULT:    -4.437", "REMARK VINA RESULT:    -5.437")
    _edit(b / "ranking.csv", "-4.437", "-5.437")
    os.remove(b / "Top_hits.txt")
    code, doc = _run([str(a), str(b)], capsys)
    assert code == 1
    cmp_ = doc["compare"]
    assert set(cmp_["different"]) == {"docked/phenylethanol_out.pdbqt", "ranking.csv"}
    assert cmp_["only_in_a"] == ["Top_hits.txt"]
    rk = cmp_["ranking"]
    assert rk["spearman_rho"] == pytest.approx(0.5)
    assert rk["max_abs_diff"] == pytest.approx(1.0)
    by = {d["ligand"]: d for d in rk["per_ligand"]}
    assert by["phenylethanol"]["diff"] == pytest.approx(-1.0)
    assert by["acetaminophen"]["diff"] == 0


def test_cli_usage_errors(tmp_path):
    with pytest.raises(SystemExit) as e:
        vr.main([str(tmp_path), str(tmp_path), "--redock"])
    assert e.value.code == 2
    with pytest.raises(SystemExit) as e:
        vr.main([str(tmp_path / "missing")])
    assert e.value.code == 2


# ---------------------------------------------------------------------------
# (b) re-docking parity on a real tiny run
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def real_run(tmp_path_factory):
    if not (shutil.which("vina") and shutil.which("obabel")):
        pytest.skip("vina/obabel not available on PATH")
    root = tmp_path_factory.mktemp("verify_real")
    ligs = root / "ligs"
    ligs.mkdir()
    shutil.copy(os.path.join(ROOT, "tests", "data", "runner", "acetaminophen_v3000_3d.sdf"),
                ligs / "acetaminophen.sdf")
    out = root / "run"
    box = json.loads(open(os.path.join(ROOT, "tests", "data", "e2e", "box.json"),
                          encoding="utf-8").read())
    cx, cy, cz = box["center"]
    s = box["size"]
    proc = subprocess.run(
        [sys.executable, os.path.join(ROOT, "autodocker", "runner.py"),
         "-r", os.path.join(ROOT, "tests", "data", "e2e", "pocket_2fxs.pdb"),
         "-l", str(ligs), "-o", str(out), "--no-fpocket",
         f"--grid-center={cx},{cy},{cz}", f"--grid-size={s},{s},{s}",
         "--exhaustiveness", "1", "--binding-modes", "2", "--seed", "3",
         "--threads", "1", "--minimize-steps", "20", "--no-html-report"],
        capture_output=True, text=True, timeout=600,
        env={**os.environ, "OMP_NUM_THREADS": "1"})
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    return out


def test_redock_reproduces_the_run_byte_for_byte(real_run, capsys):
    code, doc = _run([str(real_run), "--redock", "--sample", "0"], capsys)
    assert code == 0, _failed(doc)
    assert doc["mode"] == "redock"
    (rec,) = doc["redock"]
    assert rec["identical"] is True and rec["ranking_match"] is True
    assert rec["modes"] == rec["modes_matched"] >= 1
    checks = {c["check"]: c["status"] for c in doc["checks"]}
    assert checks["redock[acetaminophen]"] == "pass"
    assert any(k.startswith("engine[") and v == "pass" for k, v in checks.items())


def test_redock_of_a_moved_run_and_tampered_pose(real_run, tmp_path, capsys):
    moved = tmp_path / "moved"
    shutil.copytree(real_run, moved)
    code, doc = _run([str(moved), "--redock"], capsys)
    assert code == 0, _failed(doc)
    assert str(moved) in doc["redock"][0]["command"]
    assert str(real_run) not in doc["redock"][0]["command"]
    # A coordinate edit keeps the run self-consistent (scores untouched) ...
    pose = moved / "docked" / "acetaminophen_out.pdbqt"
    lines = pose.read_text(encoding="utf-8").splitlines(keepends=True)
    i = next(n for n, line in enumerate(lines) if line.startswith(("ATOM", "HETATM")))
    x = float(lines[i][30:38])
    lines[i] = lines[i][:30] + f"{x + 0.5:8.3f}" + lines[i][38:]
    pose.write_text("".join(lines), encoding="utf-8")
    code, doc = _run([str(moved)], capsys)
    assert code == 0
    # ... but re-docking exposes it.
    code, doc = _run([str(moved), "--redock"], capsys)
    assert code == 1
    assert "redock[acetaminophen]" in _failed(doc)
    assert doc["redock"][0]["identical"] is False
