import json
import os

import pytest

import runner
from vspipeline import docking


def test_duplicate_names_rejected(tmp_path, charged_pdbqt):
    lig = charged_pdbqt
    outdir = str(tmp_path / "out")
    os.makedirs(outdir, exist_ok=True)
    with pytest.raises(ValueError, match="Duplicate ligand names"):
        runner.dock_all(
            receptor=lig, ligands=[lig, lig], grid_file=lig,
            outdir=outdir, num_processes=1, resume=False)


def test_dock_all_empty_ligands(tmp_path, charged_pdbqt):
    outdir = str(tmp_path / "out")
    os.makedirs(outdir, exist_ok=True)
    results, checkpoint, metrics = runner.dock_all(
        receptor=charged_pdbqt, ligands=[], grid_file=charged_pdbqt,
        outdir=outdir, num_processes=1, resume=False)
    assert results == []
    assert checkpoint.completed == {}


def test_all_fail_reports_failed_rows_with_pointer(tmp_path, empty_pdbqt, caplog):
    """When every ligand fails, dock_all logs an error that points to the
    logs and returns FAILED rows (the runner writes all outputs and exits 3,
    see test_runner.py::test_all_ligands_failed_writes_outputs_and_exits_3)."""
    outdir = str(tmp_path / "out")
    os.makedirs(outdir, exist_ok=True)
    receptor = str(tmp_path / "rec.pdbqt")
    with open(receptor, "w", encoding="utf-8") as f:
        f.write("END\n")
    with caplog.at_level("ERROR"):
        results, _checkpoint, metrics = runner.dock_all(
            receptor=receptor, ligands=[empty_pdbqt], grid_file=receptor,
            outdir=outdir, num_processes=1, resume=False,
            vina_params={"seed": 42})
    assert results == [("empty", None)]
    assert metrics["empty"]["status"] == "FAILED"
    assert "All docking attempts failed" in caplog.text and "_vina.log" in caplog.text


def test_dock_ligand_missing_ligand_file(tmp_path, charged_pdbqt):
    name, score, metrics = runner.dock_ligand((
        charged_pdbqt, str(tmp_path / "missing.pdbqt"),
        charged_pdbqt, str(tmp_path), {}, None))
    assert score is None
    assert metrics["status"] == "FAILED"
    assert "missing" in metrics["error"].lower()


def test_dock_ligand_missing_receptor(tmp_path, charged_pdbqt):
    name, score, metrics = runner.dock_ligand((
        str(tmp_path / "missing_rec.pdbqt"), charged_pdbqt,
        charged_pdbqt, str(tmp_path), {}, None))
    assert score is None
    assert metrics["status"] == "FAILED"


def test_zero_charge_receptor_and_ligand_only_warn(tmp_path, zero_charge_pdbqt, caplog):
    """Vina ignores partial charges: all-zero charges are a warning, not FAILED."""
    eng = _fake_engine(tmp_path)
    lig = _ligand(tmp_path, charge="0.000")
    with caplog.at_level("WARNING"):
        name, score, metrics = docking.dock_ligand(_job(tmp_path, eng, receptor=zero_charge_pdbqt,
                                                        ligand=lig))
    assert metrics["status"] == "OK", metrics
    assert "partial charges are zero" in caplog.text


def test_build_vina_command_standard():
    cmd = runner._build_vina_command(
        "/usr/bin/vina", "rec.pdbqt", "lig.pdbqt", "out.pdbqt",
        {"center_x": 1, "center_y": 2, "center_z": 3,
         "size_x": 24, "size_y": 24, "size_z": 24},
        {"exhaustiveness": 8, "binding_modes": 9, "energy_range": 3.0,
         "seed": 42},
        grid_file="grid.conf", flex_file="flex.pdbqt")
    assert "--center_x" in cmd
    assert "--flex" in cmd
    assert "--seed" in cmd
    assert "--config" not in cmd


def test_build_vina_command_quickvina_uses_config():
    cmd = runner._build_vina_command(
        "qvina2", "rec.pdbqt", "lig.pdbqt", "out.pdbqt",
        {"center_x": 1}, {"exhaustiveness": 8, "binding_modes": 9,
                          "energy_range": 3.0, "seed": 42},
        grid_file="grid.conf", flex_file="flex.pdbqt")
    assert "--config" in cmd
    assert "--center_x" not in cmd
    # QuickVina 2 accepts --flex; omitting it would dock against a receptor
    # whose flexible side chains were removed.
    assert cmd[cmd.index("--flex") + 1] == "flex.pdbqt"


def test_build_vina_command_threads(monkeypatch):
    monkeypatch.setattr(runner, "_vina_cpu_flag", lambda p: "--cpu")
    cmd = runner._build_vina_command(
        "/usr/bin/vina", "rec.pdbqt", "lig.pdbqt", "out.pdbqt",
        {"center_x": 1, "size_x": 24},
        {"exhaustiveness": 8, "binding_modes": 9, "energy_range": 3.0,
         "seed": 42, "threads": 2})
    assert "--cpu" in cmd
    assert cmd[cmd.index("--cpu") + 1] == "2"


def test_build_vina_command_threads_unsupported(monkeypatch):
    monkeypatch.setattr(runner, "_vina_cpu_flag", lambda p: None)
    cmd = runner._build_vina_command(
        "/usr/bin/vina", "rec.pdbqt", "lig.pdbqt", "out.pdbqt",
        {"center_x": 1, "size_x": 24},
        {"exhaustiveness": 8, "binding_modes": 9, "energy_range": 3.0,
         "seed": 42, "threads": 2})
    assert "--cpu" not in cmd and "--threads" not in cmd


@pytest.mark.parametrize("help_text,expected", [
    ("  --cpu arg (=0)   the number of CPUs to use", "--cpu"),
    ("  --threads arg    number of threads", "--threads"),
    ("  --seed arg       explicit random seed", None),
])
def test_vina_cpu_flag_detection(monkeypatch, help_text, expected):
    class Result:
        stdout = help_text
        stderr = ""
    monkeypatch.setattr(docking.subprocess, "run", lambda *a, **k: Result())
    docking._VINA_CPU_FLAG_CACHE.clear()
    assert docking._vina_cpu_flag("/fake/vina") == expected
    docking._VINA_CPU_FLAG_CACHE.clear()


def test_build_vina_command_no_threads_by_default():
    cmd = runner._build_vina_command(
        "/usr/bin/vina", "rec.pdbqt", "lig.pdbqt", "out.pdbqt",
        {"center_x": 1, "size_x": 24},
        {"exhaustiveness": 8, "binding_modes": 9, "energy_range": 3.0,
         "seed": 42})
    assert "--threads" not in cmd and "--cpu" not in cmd


# Serial vs parallel docking of real ligands (byte-identical poses and
# scores) is tested in test_e2e.py (test_process_count_does_not_change_results).


# ---------------------------------------------------------------------------
# Fake engine helpers (no real docking)
# ---------------------------------------------------------------------------

FAKE_ENGINE = r"""#!/usr/bin/env python3
import os, sys, time
args = sys.argv[1:]
if "--help" in args or "--help_advanced" in args:
    print(os.environ.get("FAKE_HELP", "AutoDock Vina fake\n  --flex arg\n  --cpu arg"))
    sys.exit(0)
if "--version" in args:
    print("AutoDock Vina vFAKE"); sys.exit(0)
with open(os.path.join(os.path.dirname(sys.argv[0]), "calls.log"), "a") as f:
    f.write(" ".join(args) + "\n")
mode = os.environ.get("FAKE_MODE_" + os.path.basename(sys.argv[0]), "ok")
if mode == "crash":
    sys.stderr.write("simulated crash\n"); sys.exit(3)
if mode == "sleep":
    time.sleep(30)
killer = os.environ.get("FAKE_KILL_PARENT_FOR")
if killer and killer in os.path.basename(args[args.index("--ligand") + 1]):
    os.kill(os.getppid(), 9)   # the docking worker dies (as by the OOM killer)
    time.sleep(30)
out = args[args.index("--out") + 1]
aff = os.environ.get("FAKE_AFF", "-11.521")
table = os.environ.get("FAKE_TABLE", "-11.52")
with open(out, "w") as f:
    f.write("MODEL 1\nREMARK VINA RESULT:   %s      0.000      0.000\n"
            "ATOM      1  C   LIG A   1       0.000   0.000   0.000  0.00  0.00    +0.000 C\n"
            "ENDMDL\nMODEL 2\nREMARK VINA RESULT:   -9.000      1.000      2.000\n"
            "ATOM      1  C   LIG A   1       1.000   0.000   0.000  0.00  0.00    +0.000 C\n"
            "ENDMDL\n" % aff)
print("mode |   affinity | dist from best mode")
print("     | (kcal/mol) | rmsd l.b.| rmsd u.b.")
print("-----+------------+----------+----------")
print("   1     %s          0          0" % table)
print("   2       -9          1          2")
sys.stderr.write("WARNING: fake engine warning\n")
"""


def _fake_engine(tmp_path, name="vina"):
    d = tmp_path / "bin"
    d.mkdir(exist_ok=True)
    exe = d / name
    exe.write_text(FAKE_ENGINE)
    exe.chmod(0o755)
    docking._VINA_TYPE_CACHE.pop(str(exe), None)
    docking._FLEX_SUPPORT_CACHE.pop(str(exe), None)
    return str(exe)


def _ligand(tmp_path, charge="+0.120", name="lig"):
    p = tmp_path / f"{name}.pdbqt"
    p.write_text(
        "REMARK  test ligand padded to pass the minimum size check.............\n"
        "ROOT\n"
        f"ATOM      1  C   LIG A   1       0.000   0.000   0.000  0.00  0.00    {charge} C \n"
        f"ATOM      2  O   LIG A   1       1.400   0.000   0.000  0.00  0.00    {charge} OA\n"
        "ENDROOT\nTORSDOF 0\n")
    return str(p)


def _job(tmp_path, primary, fallback=None, receptor=None, ligand=None, flex=None,
         params=None, timeout=60):
    grid = tmp_path / "grid.conf"
    grid.write_text("center_x=0\ncenter_y=0\ncenter_z=0\nsize_x=20\nsize_y=20\nsize_z=20\n")
    if receptor is None:
        receptor = str(tmp_path / "rec.pdbqt")
        with open(receptor, "w", encoding="utf-8") as f:
            f.write("ATOM      1  N   LYS A   1       5.000   0.000   0.000  0.00  0.00    -0.300 N \n")
    dock_dir = tmp_path / "docked"
    dock_dir.mkdir(exist_ok=True)
    return docking.DockJob(
        receptor, ligand or _ligand(tmp_path), str(grid), str(dock_dir),
        params or {"seed": 42, "binding_modes": 2, "min_valid_affinity": -1.0},
        flex, docking.EngineConfig(primary, fallback or primary, None, timeout))


def test_score_read_from_pdbqt_not_rounded_stdout(tmp_path):
    eng = _fake_engine(tmp_path)
    name, score, metrics = docking.dock_ligand(_job(tmp_path, eng))
    assert score == -11.521          # stdout table says -11.52
    assert metrics["status"] == "OK"
    assert metrics["engine"] == "Vina"
    assert metrics["binding_modes"] == 2
    assert metrics["simscore"] == pytest.approx(1.0)


def test_stdout_disagreement_is_logged(tmp_path, monkeypatch, caplog):
    eng = _fake_engine(tmp_path)
    monkeypatch.setenv("FAKE_TABLE", "-10.9")
    with caplog.at_level("WARNING"):
        _, score, _ = docking.dock_ligand(_job(tmp_path, eng))
    assert score == -11.521
    assert "disagrees with the stdout" in caplog.text


def test_parse_pdbqt_results(tmp_path):
    p = tmp_path / "o.pdbqt"
    p.write_text("MODEL 1\nREMARK VINA RESULT:   -11.521      0.000      0.000\nENDMDL\n"
                 "MODEL 2\nREMARK VINA RESULT:   -10.004      1.234      2.345\nENDMDL\n")
    assert docking.parse_pdbqt_results(str(p)) == [(-11.521, 0.0, 0.0),
                                                   (-10.004, 1.234, 2.345)]
    s = tmp_path / "s.pdbqt"   # smina --score_only keeps the input's Vina remark
    s.write_text("MODEL 1\nREMARK minimizedAffinity -3.53028774\nREMARK minimizedRMSD -1\n"
                 "REMARK VINA RESULT:    -2.717      0.000      0.000\nENDMDL\n")
    assert docking.parse_pdbqt_results(str(s)) == [(-3.53028774, None, None)]
    q = tmp_path / "q.pdbqt"   # QuickVina prints one decimal
    q.write_text("MODEL 1\nREMARK VINA RESULT:      -2.7      0.000      0.000\nENDMDL\n")
    assert docking.parse_pdbqt_results(str(q)) == [(-2.7, 0.0, 0.0)]


def test_scores_agree_uses_print_precision():
    assert docking._scores_agree("-11.521", "-11.52")
    assert docking._scores_agree("-2.7", "-2.7")
    assert not docking._scores_agree("-11.521", "-11.4")


def test_out_of_window_result_is_reported_not_failed(tmp_path, monkeypatch):
    eng = _fake_engine(tmp_path)
    monkeypatch.setenv("FAKE_AFF", "-0.512")
    monkeypatch.setenv("FAKE_TABLE", "-0.512")
    name, score, metrics = docking.dock_ligand(_job(tmp_path, eng))
    assert score is None
    assert metrics["status"] == "OUT_OF_WINDOW"
    assert metrics["raw_score"] == -0.512
    assert "outside the accepted window" in metrics["error"]
    log = open(metrics["log_file"], encoding="utf-8").read()
    assert "# Result: OUT_OF_WINDOW" in log


@pytest.mark.parametrize("score,inside", [
    (-1.0, False), (-1.001, True), (-20.0, True), (-20.001, False), (0.0, False)])
def test_affinity_window_bounds(score, inside):
    assert docking.affinity_in_window(score, docking.DEFAULT_MIN_VALID_AFFINITY) is inside


def test_command_logged_before_running_and_stderr_kept(tmp_path, monkeypatch):
    eng = _fake_engine(tmp_path)
    _, _, metrics = docking.dock_ligand(_job(tmp_path, eng))
    log = open(metrics["log_file"], encoding="utf-8").read().splitlines()
    assert log[0] == "# Attempt: Vina (primary)"
    assert log[1].startswith("# Command: " + eng + " --receptor ")
    # plain paths: identical to the old ' '.join() text the parity script splits
    assert log[1].split(":", 1)[1].strip().split()[0] == eng
    assert "WARNING: fake engine warning" in "\n".join(log)

    monkeypatch.setenv("FAKE_MODE_vina", "crash")
    _, score, metrics = docking.dock_ligand(_job(tmp_path, eng))
    assert score is None and metrics["status"] == "FAILED"
    text = open(metrics["log_file"], encoding="utf-8").read()
    assert "# Command: " in text and "simulated crash" in text


def test_traceback_is_real(tmp_path, monkeypatch):
    eng = _fake_engine(tmp_path)

    def boom(*a, **k):
        raise KeyError("unexpected")

    monkeypatch.setattr(docking, "_build_metrics", boom)
    _, score, metrics = docking.dock_ligand(_job(tmp_path, eng))
    assert metrics["status"] == "FAILED"
    text = open(metrics["log_file"], encoding="utf-8").read()
    assert "# Traceback:\nTraceback (most recent call last):" in text
    assert "KeyError" in text


def test_fallback_engine_is_recorded_and_warned(tmp_path, monkeypatch, caplog):
    primary = _fake_engine(tmp_path, "vina")
    fallback = _fake_engine(tmp_path, "qvina2")
    monkeypatch.setenv("FAKE_MODE_vina", "crash")
    with caplog.at_level("WARNING"):
        _, score, metrics = docking.dock_ligand(_job(tmp_path, primary, fallback))
    assert score == -11.521
    assert metrics["engine"] == "QuickVina"
    assert metrics["fallback_used"] is True
    assert "fallback engine" in caplog.text


def test_timeout_comes_from_job_config(tmp_path, monkeypatch):
    eng = _fake_engine(tmp_path)
    monkeypatch.setenv("FAKE_MODE_vina", "sleep")
    monkeypatch.setattr(runner, "COMMAND_TIMEOUT", 999)
    t0 = __import__("time").time()
    _, score, metrics = docking.dock_ligand(_job(tmp_path, eng, timeout=1))
    assert __import__("time").time() - t0 < 20
    assert "timeout after 1s" in metrics["error"]


def test_qvina_flexible_docking_passes_flex(tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_HELP", "QuickVina 2.1\n  --flex arg\n  --config arg")
    eng = _fake_engine(tmp_path, "qvina2")
    assert docking._detect_vina_type(eng) == "quickvina"
    flex = tmp_path / "flex.pdbqt"
    flex.write_text("BEGIN_RES LYS A   1\nEND_RES LYS A   1\n")
    _, score, metrics = docking.dock_ligand(_job(tmp_path, eng, flex=str(flex)))
    assert score is not None, metrics
    calls = open(os.path.join(os.path.dirname(eng), "calls.log"), encoding="utf-8").read()
    assert "--flex " + str(flex) in calls


def test_engine_without_flex_support_fails_clearly(tmp_path, monkeypatch):
    eng = _fake_engine(tmp_path, "vina")
    monkeypatch.setenv("FAKE_HELP", "AutoDock Vina fake\n  --cpu arg\n  --seed arg")
    flex = tmp_path / "flex.pdbqt"
    flex.write_text("BEGIN_RES LYS A   1\nEND_RES LYS A   1\n")
    _, score, metrics = docking.dock_ligand(_job(tmp_path, eng, flex=str(flex)))
    assert metrics["status"] == "FAILED"
    assert "do not accept --flex" in metrics["error"]


def test_simscore_undefined_is_none(tmp_path):
    p = tmp_path / "o.pdbqt"
    p.write_text("MODEL 1\nREMARK VINA RESULT:   -7.000      0.000      0.000\nENDMDL\n")
    m = docking._calculate_metrics(str(p), expected_modes=9)
    assert m["simscore"] is None and m["binding_modes"] == 1
    assert docking._calculate_protocol_simscore([]) == (None, None, None)


def test_receptor_cache_keyed_by_content_stamp(tmp_path):
    rec = tmp_path / "r.pdbqt"
    rec.write_text("ATOM      1  N   LYS A   1       5.000   0.000   0.000  0.00  0.00    +0.000 N \n")
    assert docking._cached_receptor_check(str(rec)) == (True, False)
    rec.write_text("ATOM      1  N   LYS A   1       5.000   0.000   0.000  0.00  0.00    -0.300 N \n"
                   "REMARK rewritten\n")
    assert docking._cached_receptor_check(str(rec)) == (True, True)


def test_explicit_params_beat_environment(tmp_path, monkeypatch):
    # Only $VS_SEED is read; other undocumented VS_* variables are ignored.
    monkeypatch.setenv("VS_EXHAUSTIVENESS", "32")
    monkeypatch.setenv("VS_SEED", "7")
    seen = {}

    def fake(job):
        seen.update(job.vina_params)
        return runner.ligand_name(job.ligand), -7.0, {"status": "OK"}

    monkeypatch.setattr(docking, "dock_ligand", fake)
    lig = _ligand(tmp_path)
    grid = tmp_path / "g.conf"
    grid.write_text("center_x=0\ncenter_y=0\ncenter_z=0\nsize_x=20\nsize_y=20\nsize_z=20\n")
    docking.dock_all(lig, [lig], str(grid), str(tmp_path / "o"), resume=False,
                     vina_params={"exhaustiveness": 4})
    assert seen["exhaustiveness"] == 4      # explicit wins
    assert seen["seed"] == 7                # missing key filled from $VS_SEED


def test_engine_config_reaches_worker_explicitly(tmp_path, monkeypatch):
    """dock_all passes its EngineConfig inside each job (spawn-safe)."""
    seen = []

    def fake(job):
        seen.append(job.engine)
        return runner.ligand_name(job.ligand), -7.0, {"status": "OK"}

    monkeypatch.setattr(docking, "dock_ligand", fake)
    lig = _ligand(tmp_path)
    grid = tmp_path / "g.conf"
    grid.write_text("center_x=0\ncenter_y=0\ncenter_z=0\nsize_x=20\nsize_y=20\nsize_z=20\n")
    cfg = docking.EngineConfig("/opt/x/vina", "/opt/x/vina", None, 33)
    docking.dock_all(lig, [lig], str(grid), str(tmp_path / "o"), resume=False,
                     vina_params={}, engine=cfg)
    assert seen == [cfg]


SPAWN_SCRIPT = r"""
import json, multiprocessing, sys
sys.path.insert(0, sys.argv[1])
multiprocessing.set_start_method("spawn", force=True)
from vspipeline import docking          # imported before runner on purpose
if __name__ == "__main__":
    engine, rec, lig1, lig2, grid, out = sys.argv[2:8]
    cfg = docking.EngineConfig(engine, engine, None, 60)
    results, cp, metrics = docking.dock_all(
        rec, [lig1, lig2], grid, out, num_processes=2, resume=False,
        vina_params={"seed": 42, "binding_modes": 2}, engine=cfg)
    print(json.dumps({"results": sorted(results),
                      "engines": sorted(m["engine_path"] for m in metrics.values())}))
"""


def test_dock_all_under_spawn_uses_explicit_engine(tmp_path):
    """Workers started with 'spawn' re-import runner (whose globals would
    point at the system vina); the EngineConfig passed in each job must win."""
    import json
    import subprocess
    import sys
    eng = _fake_engine(tmp_path)
    job = _job(tmp_path, eng)
    lig2 = _ligand(tmp_path, name="lig2")
    script = tmp_path / "spawn_run.py"
    script.write_text(SPAWN_SCRIPT)
    pkg = os.path.join(os.path.dirname(__file__), "..", "autodocker")
    proc = subprocess.run(
        [sys.executable, str(script), pkg, eng, job.receptor, job.ligand, lig2,
         job.grid_file, str(tmp_path / "spawn_out")],
        capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]
    data = json.loads(proc.stdout.strip().splitlines()[-1])
    assert data["results"] == [["lig", -11.521], ["lig2", -11.521]]
    assert data["engines"] == [eng, eng]


# ---------------------------------------------------------------------------
# Worker processes that die; fallback engine and the checkpoint
# ---------------------------------------------------------------------------

def _grid(tmp_path):
    grid = tmp_path / "g.conf"
    grid.write_text("center_x=0\ncenter_y=0\ncenter_z=0\nsize_x=20\nsize_y=20\nsize_z=20\n")
    return str(grid)


def test_dead_worker_fails_its_ligand_and_the_run_finishes(tmp_path, monkeypatch):
    """A worker killed while docking (kill -9, OOM killer) used to hang the
    screen forever ("Progress 8/9"). The ligand whose worker dies twice is
    FAILED; every other ligand is docked; results are in library order."""
    eng = _fake_engine(tmp_path)
    monkeypatch.setenv("FAKE_KILL_PARENT_FOR", "killer")
    job = _job(tmp_path, eng)
    ligs = [_ligand(tmp_path, name=n) for n in ("l1", "killer", "l2", "l3", "l4")]
    results, cp, metrics = docking.dock_all(
        job.receptor, ligs, job.grid_file, str(tmp_path / "o"), num_processes=2,
        resume=False, vina_params={"seed": 42, "binding_modes": 2},
        engine=docking.EngineConfig(eng, eng, None, 60))
    assert results == [("l1", -11.521), ("killer", None), ("l2", -11.521),
                       ("l3", -11.521), ("l4", -11.521)]
    assert metrics["killer"]["status"] == "FAILED"
    assert metrics["killer"]["error"] == docking.WORKER_DIED_ERROR
    assert "killer" not in cp.completed
    attempts = [c for c in (tmp_path / "bin" / "calls.log").read_text().splitlines()
                if "killer" in c]
    assert len(attempts) == 2      # first pool, then once in a fresh worker


def _fallback_engines(tmp_path):
    primary = _fake_engine(tmp_path, "vina")
    fallback = _fake_engine(tmp_path, "qvina2")
    return primary, fallback


def _calls(tmp_path):
    path = tmp_path / "bin" / "calls.log"
    return path.read_text().count("\n") if path.exists() else 0


def test_toggling_fallback_keeps_primary_results(tmp_path):
    """The fallback engine is not part of the run fingerprint: turning
    --qvina-fallback on or off keeps every result the primary produced."""
    primary, fallback = _fallback_engines(tmp_path)
    job = _job(tmp_path, primary)
    ligs = [_ligand(tmp_path, name=n) for n in ("a", "b")]
    out = str(tmp_path / "o")
    params = {"seed": 42, "binding_modes": 2}
    docking.dock_all(job.receptor, ligs, job.grid_file, out, vina_params=params,
                     engine=docking.EngineConfig(primary, primary, None, 60))
    assert _calls(tmp_path) == 2
    results, _, _ = docking.dock_all(
        job.receptor, ligs, job.grid_file, out, vina_params=params,
        engine=docking.EngineConfig(primary, fallback, None, 60))
    assert _calls(tmp_path) == 2 and results == [("a", -11.521), ("b", -11.521)]


def test_fallback_result_redocked_when_fallback_turned_off(tmp_path, monkeypatch):
    """A result of the fallback engine is tagged with its SHA-256 and docked
    again (alone) once that fallback is no longer configured."""
    primary, fallback = _fallback_engines(tmp_path)
    job = _job(tmp_path, primary)
    ligs = [_ligand(tmp_path, name=n) for n in ("a",)]
    out = str(tmp_path / "o")
    params = {"seed": 42, "binding_modes": 2}
    with_fallback = docking.EngineConfig(primary, fallback, None, 60)
    monkeypatch.setenv("FAKE_MODE_vina", "crash")
    _, _, metrics = docking.dock_all(job.receptor, ligs, job.grid_file, out,
                                     vina_params=params, engine=with_fallback)
    assert metrics["a"]["fallback_used"] is True
    assert metrics["a"]["fallback_sha256"] == docking._file_sha256(fallback)
    n = _calls(tmp_path)
    docking.dock_all(job.receptor, ligs, job.grid_file, out, vina_params=params,
                     engine=with_fallback)
    assert _calls(tmp_path) == n          # same fallback: cached
    monkeypatch.delenv("FAKE_MODE_vina")
    results, _, metrics = docking.dock_all(
        job.receptor, ligs, job.grid_file, out, vina_params=params,
        engine=docking.EngineConfig(primary, primary, None, 60))
    assert _calls(tmp_path) == n + 1      # docked again with the primary
    assert results == [("a", -11.521)] and metrics["a"]["fallback_used"] is False


def test_checkpointed_paths_survive_cwd_change_and_move(tmp_path, monkeypatch):
    """dock_file/log_file were stored relative to the working directory of
    the run that docked; a resume from elsewhere (or of a moved directory)
    pointed at files that do not exist."""
    eng = _fake_engine(tmp_path)
    job = _job(tmp_path, eng)
    lig = _ligand(tmp_path, name="a")
    params = {"seed": 42, "binding_modes": 2}
    cfg = docking.EngineConfig(eng, eng, None, 60)
    monkeypatch.chdir(tmp_path)
    docking.dock_all(job.receptor, [lig], job.grid_file, "o", vina_params=params,
                     engine=cfg)
    with open(tmp_path / "o" / ".docking_checkpoint.json", encoding="utf-8") as f:
        stored = json.load(f)["a"]["metrics"]
    assert stored["dock_file"] == os.path.join("docked", "a_out.pdbqt")
    assert stored["log_file"] == os.path.join("docked", "a_vina.log")
    from vspipeline import checkpoint as ckpt
    ckpt.release_lock("o")
    (tmp_path / "moved").mkdir()
    os.rename(tmp_path / "o", tmp_path / "moved" / "o")
    monkeypatch.chdir(tmp_path / "bin")
    new_out = str(tmp_path / "moved" / "o")
    _, _, metrics = docking.dock_all(job.receptor, [lig], job.grid_file, new_out,
                                     vina_params=params, engine=cfg)
    assert metrics["a"]["dock_file"] == os.path.join(new_out, "docked", "a_out.pdbqt")
    assert os.path.exists(metrics["a"]["dock_file"])
    ckpt.release_lock(new_out)
