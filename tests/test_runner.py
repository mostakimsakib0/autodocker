import csv
import hashlib
import json
import logging
import os
import sys

import pytest

import runner


@pytest.fixture
def cli(tmp_path, protein_pdb):
    """Valid base arguments: an existing receptor and ligand directory."""
    ligs = tmp_path / "ligs"
    ligs.mkdir()
    (ligs / "x.pdbqt").write_text("ROOT\nENDROOT\nTORSDOF 0\n")
    out = tmp_path / "out"
    return ["-r", protein_pdb, "-l", str(ligs), "-o", str(out)], out


@pytest.fixture
def fake_tools(monkeypatch):
    """Use the Python interpreter as a stand-in for every external tool."""
    for name in ("OBABEL", "VINA_PRIMARY", "VINA_FALLBACK"):
        monkeypatch.setattr(runner, name, sys.executable)
    for name in ("SMINA", "FPOCKET", "REDUCE"):
        monkeypatch.setattr(runner, name, None)


@pytest.fixture(autouse=True)
def _restore_state(monkeypatch):
    monkeypatch.setattr(runner, "_ENV_ERRORS", [])
    monkeypatch.setattr(runner, "COMMAND_TIMEOUT", runner.COMMAND_TIMEOUT)
    yield
    runner._configure_logging(False)


def _exit_code(argv):
    with pytest.raises(SystemExit) as exc:
        runner.main(argv)
    return exc.value.code


def test_find_tool_present():
    assert runner.find_tool("python3") is not None


def test_find_tool_missing():
    assert runner.find_tool("definitely-not-a-real-binary-xyz") is None


def test_find_tool_prefers_first():
    path = runner.find_tool("python3", "sh")
    assert path is not None and path.endswith("python3")


def test_default_seed_is_42_without_env(monkeypatch):
    monkeypatch.delenv("VS_SEED", raising=False)
    assert runner.default_seed() == 42


def test_default_seed_reads_vs_seed(monkeypatch):
    monkeypatch.setenv("VS_SEED", "7")
    assert runner.default_seed() == 7


def test_default_seed_rejects_non_integer(monkeypatch):
    monkeypatch.setenv("VS_SEED", "abc")
    with pytest.raises(runner.ConfigError, match="VS_SEED"):
        runner.default_seed()


# ---------------------------------------------------------------------------
# Version
# ---------------------------------------------------------------------------

def test_version_sources(tmp_path, monkeypatch):
    monkeypatch.delenv("AUTODOCKER_VERSION", raising=False)
    assert runner._read_version(str(tmp_path)) == runner._FALLBACK_VERSION
    monkeypatch.setenv("AUTODOCKER_VERSION", "9.8.7")
    assert runner._read_version(str(tmp_path)) == "9.8.7"
    (tmp_path / "VERSION").write_text("1.2.3\n")
    assert runner._read_version(str(tmp_path)) == "1.2.3"


def test_version_flag_survives_bad_environment(monkeypatch, capsys):
    monkeypatch.setenv("VS_SEED", "abc")
    monkeypatch.setattr(runner, "_ENV_ERRORS", ["VS_COMMAND_TIMEOUT must be ..."])
    assert _exit_code(["--version"]) == 0
    assert "AutoDocker" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# Fail-loud validation (exit code 2, before any ligand is prepared)
# ---------------------------------------------------------------------------

def test_bad_vs_seed_is_a_clean_usage_error(monkeypatch, capsys, cli, fake_tools):
    argv, out = cli
    monkeypatch.setenv("VS_SEED", "abc")
    assert _exit_code(argv) == runner.EXIT_USAGE
    err = capsys.readouterr().err
    assert "VS_SEED must be an integer" in err and "Traceback" not in err
    assert not out.exists()


def test_bad_command_timeout_env_is_a_clean_usage_error(monkeypatch, capsys, cli, fake_tools):
    monkeypatch.setenv("VS_COMMAND_TIMEOUT", "ten")
    assert runner._env_command_timeout() == runner.DEFAULT_COMMAND_TIMEOUT
    assert "VS_COMMAND_TIMEOUT" in runner._ENV_ERRORS[0]
    argv, _ = cli
    assert _exit_code(argv) == runner.EXIT_USAGE
    assert "VS_COMMAND_TIMEOUT must be a positive integer" in capsys.readouterr().err


@pytest.mark.parametrize("extra, message", [
    (["--exhaustiveness", "0"], "--exhaustiveness must be >= 1"),
    (["--timeout", "0"], "--timeout must be >= 1"),
    (["--auto-flexible", "0"], "--auto-flexible must be >= 1"),
    (["--threads", "0"], "--threads must be >= 1"),
    (["--binding-modes", "0"], "--binding-modes must be >= 1"),
    (["--energy-range", "0"], "--energy-range must be > 0"),
    (["--vina-bin", "/nonexistent/vina"], "does not exist or is not executable"),
    (["--mw-max", "600"], "only apply to --library custom"),
    (["--mw-ma", "600"], "only apply to --library custom"),  # abbreviation
    (["--library", "fda", "--logp-max", "4"], "only apply to --library custom"),
    (["--max-compounds", "5"], "--max-compounds only applies"),
    (["--library", "custom", "--mw-min", "600"], "--mw-min (600) must be"),
    (["--max-failure-fraction", "1.5"], "--max-failure-fraction must be between"),
    (["--grid-size", "1,2"], "three comma-separated values"),
    (["--grid-size", "0,20,20"], "--grid-size must be positive"),
    (["--smina-only", "--consensus"], "--smina-only cannot be combined with --consensus"),
    (["--smina-only", "--binding-modes", "5"], "--binding-modes"),
    (["--smina-only", "--energy-range", "4"], "--energy-range"),
    (["--smina-only", "--vina-extra=--local_only"], "--vina-extra"),
    (["--smina-only", "--vina-bin", sys.executable], "--vina-bin"),
    (["--qvina-fallback", "--vina-bin", sys.executable], "--qvina-fallback cannot be combined with --vina-bin"),
    (["--qvina-fallback", "--smina-only"], "--qvina-fallback cannot be combined with --smina-only"),
])
def test_invalid_options_exit_2(cli, fake_tools, capsys, extra, message):
    argv, out = cli
    assert _exit_code(argv + extra) == runner.EXIT_USAGE
    assert message in capsys.readouterr().err
    assert not out.exists()  # rejected before any output/ligand preparation


def test_ligands_required_for_local_only(tmp_path, protein_pdb, fake_tools, capsys, monkeypatch):
    out = tmp_path / "out"
    assert _exit_code(["-r", protein_pdb, "-o", str(out)]) == runner.EXIT_USAGE
    assert "-l/--ligands is required" in capsys.readouterr().err
    monkeypatch.setattr(runner, "_run_pipeline", lambda args, state: runner.EXIT_OK)
    runner.main(["-r", protein_pdb, "-o", str(out), "--library", "fda",
                 "--max-compounds", "3"])


def test_missing_receptor_exit_2(tmp_path, fake_tools, capsys):
    assert _exit_code(["-r", str(tmp_path / "none.pdb"), "-l", str(tmp_path),
                       "-o", str(tmp_path / "o")]) == runner.EXIT_USAGE
    assert "no such file" in capsys.readouterr().err


def test_missing_tool_exit_2(cli, fake_tools, monkeypatch, caplog):
    monkeypatch.setattr(runner, "OBABEL", None)
    argv, out = cli
    with caplog.at_level(logging.ERROR):
        assert _exit_code(argv) == runner.EXIT_USAGE
    assert "Required tool(s) not found in PATH: obabel" in caplog.text
    assert not out.exists()


def test_consensus_requires_smina(cli, fake_tools, caplog):
    argv, _ = cli
    with caplog.at_level(logging.ERROR):
        assert _exit_code(argv + ["--consensus"]) == runner.EXIT_USAGE
    assert "smina" in caplog.text


def test_processes_resolution(cli, fake_tools, monkeypatch):
    seen = {}

    def pipeline(args, state):
        seen["p"] = args.processes
        return runner.EXIT_OK

    monkeypatch.setattr(runner, "_run_pipeline", pipeline)
    monkeypatch.setattr(runner, "available_cpus", lambda: 8)
    argv, _ = cli
    for value, expected in (("3", 3), ("0", 8), ("-1", 8), ("-4", 2)):
        runner.main(argv + ["-p", value])
        assert seen["p"] == expected


def test_no_admet_is_a_deprecated_alias(cli, fake_tools, monkeypatch, caplog):
    seen = {}
    monkeypatch.setattr(runner, "_run_pipeline",
                        lambda args, state: seen.setdefault("v", args.no_lipinski) and 0)
    argv, _ = cli
    with caplog.at_level(logging.WARNING):
        runner.main(argv + ["--no-admet"])
    assert seen["v"] is True
    assert "--no-admet is deprecated" in caplog.text


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def test_verbose_reaches_vspipeline_loggers():
    runner._configure_logging(True)
    assert logging.getLogger("vspipeline.library").getEffectiveLevel() == logging.DEBUG
    assert logging.getLogger("vspipeline.docking").getEffectiveLevel() == logging.DEBUG
    runner._configure_logging(False)
    assert logging.getLogger("vspipeline.library").getEffectiveLevel() == logging.INFO


def test_verbose_flag_in_main(cli, fake_tools, monkeypatch):
    levels = {}

    def pipeline(args, state):
        levels["lib"] = logging.getLogger("vspipeline.library").getEffectiveLevel()
        return runner.EXIT_OK

    monkeypatch.setattr(runner, "_run_pipeline", pipeline)
    argv, _ = cli
    runner.main(argv + ["-v"])
    assert levels["lib"] == logging.DEBUG


# ---------------------------------------------------------------------------
# Exit codes and run manifest
# ---------------------------------------------------------------------------

def _manifest(out):
    with open(out / "run_manifest.json", encoding="utf-8") as f:
        return json.load(f)


def test_run_manifest_on_success(cli, fake_tools, monkeypatch):
    monkeypatch.setenv("VS_SEED", "7")
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    monkeypatch.setattr(runner, "_run_pipeline", lambda args, state: runner.EXIT_OK)
    argv, out = cli
    runner.main(argv + ["-p", "0"])
    text = (out / "run_manifest.json").read_text()
    m = json.loads(text)
    assert list(m) == sorted(m)  # deterministic key order
    assert m["status"] == "completed" and m["exit_code"] == 0 and m["error"] is None
    assert m["autodocker_version"] == runner.__version__
    assert m["argv"][1:] == argv + ["-p", "0"]
    assert m["seed"] == 7 and m["args"]["seed"] == 7
    assert m["args"]["processes"] == runner.available_cpus()  # resolved value
    assert m["environment"]["VS_SEED"] == "7"
    assert m["environment"]["OMP_NUM_THREADS"] == "1"
    assert m["python"]["version"] == ".".join(map(str, sys.version_info[:3]))
    assert m["cpu_count"]["available"] >= 1
    assert m["started_at"] <= m["finished_at"]
    obabel = m["tools"]["obabel"]
    assert obabel["path"] == sys.executable
    assert obabel["version"].startswith("Python")
    with open(os.path.realpath(sys.executable), "rb") as f:
        assert obabel["sha256"] == hashlib.sha256(f.read()).hexdigest()
    assert "vina" in m["tools"] and m["tools"]["fpocket"]["path"] is None


def test_run_manifest_and_exit_code_on_failure(cli, fake_tools, monkeypatch, caplog):
    def boom(args, state):
        raise ValueError("boom")

    monkeypatch.setattr(runner, "_run_pipeline", boom)
    argv, out = cli
    with caplog.at_level(logging.ERROR):
        assert _exit_code(argv) == runner.EXIT_FAILURE
    assert "Pipeline failed: ValueError: boom" in caplog.text
    m = _manifest(out)
    assert m["status"] == "failed" and m["error"] == "ValueError: boom"


@pytest.mark.parametrize("started, expected", [
    (False, "before docking started"), (True, "checkpoint")])
def test_ctrl_c_exits_130(cli, fake_tools, monkeypatch, caplog, started, expected):
    def interrupted(args, state):
        state["docking_started"] = started
        raise KeyboardInterrupt

    monkeypatch.setattr(runner, "_run_pipeline", interrupted)
    argv, out = cli
    with caplog.at_level(logging.WARNING):
        assert _exit_code(argv) == runner.EXIT_INTERRUPTED
    assert expected in caplog.text
    if not started:
        assert "checkpoint" not in caplog.text
    assert _manifest(out)["status"] == "interrupted"


def test_too_many_failures_exit_3(cli, fake_tools, monkeypatch):
    monkeypatch.setattr(runner, "_run_pipeline",
                        lambda args, state: runner.EXIT_TOO_MANY_FAILURES)
    argv, out = cli
    assert _exit_code(argv) == runner.EXIT_TOO_MANY_FAILURES
    assert _manifest(out)["status"] == "too_many_failures"


def test_failure_fraction_rules():
    results = [("a", -7.0), ("b", None), ("c", None), ("d", None)]
    metrics = {"b": {"status": "FAILED"}, "c": {"status": "REJECTED (score)"},
               "d": {"status": "OUT_OF_WINDOW"}}
    assert runner._failed_docking(results, metrics) == ["b"]
    assert runner._failure_exit_code(1, 4, 0.5) == runner.EXIT_OK
    assert runner._failure_exit_code(3, 4, 0.5) == runner.EXIT_TOO_MANY_FAILURES
    assert runner._failure_exit_code(4, 4, 1.0) == runner.EXIT_OK
    assert runner._failure_exit_code(0, 0, 0.0) == runner.EXIT_OK


# ---------------------------------------------------------------------------
# Outputs and helpers
# ---------------------------------------------------------------------------

def test_consensus_csv_has_only_data_rows(tmp_path):
    path = tmp_path / "consensus_ranking.csv"
    agreement = runner.write_consensus_outputs(
        str(tmp_path), [("a", -8.0), ("b", -6.0), ("c", -7.0), ("d", None)],
        {"a": -7.5, "b": -6.5, "c": -7.2, "d": -5.0}, "vinardo")
    rows = list(csv.reader(open(path, encoding="utf-8")))
    assert rows[0] == ["Ligand", "Vina_Affinity", "SMINA_Affinity", "Consensus"]
    assert [r[0] for r in rows[1:]] == ["a", "c", "b", "d"]
    assert rows[-1] == ["d", "FAILED", "-5.0", "NA"]
    assert not any(line.startswith("#") for line in open(path, encoding="utf-8"))
    assert agreement["n"] == 3 and agreement["rho"] == pytest.approx(1.0)


def test_keep_residue_key_uses_resname_or_name():
    assert runner._keep_residue_key({"chain": "A", "resnum": 5, "resname": "HOH"}) == ("A", 5, "HOH")
    assert runner._keep_residue_key({"chain": "A", "resnum": 9, "name": "HEM"}) == ("A", 9, "HEM")


def test_help_and_docstring_are_neutral():
    text = runner.build_parser().format_help() + (runner.__doc__ or "")
    for word in ("PHASE", "ENHANCED", "professional", "Smart preprocessing",
                 "Original Features", "ADMET filtering"):
        assert word not in text
    assert "[1/5]" in text and "SDF" in text and "MOL2" in text and "PDBQT" in text
    assert "(default: 0.5)" in text  # ArgumentDefaultsHelpFormatter


@pytest.mark.parametrize("module", ["runner.py", "vspipeline/library.py",
                                    "vspipeline/admet.py", "vspipeline/net.py",
                                    "vspipeline/utils.py"])
def test_sources_are_ascii(module):
    """No emoji or check-mark log markers."""
    path = os.path.join(os.path.dirname(runner.__file__), module)
    text = open(path, encoding="utf-8").read()
    assert all(ord(c) < 128 for c in text), module


def test_cli_defaults_come_from_the_defaults_module():
    """Every scientific default is defined once, in vspipeline/defaults.py."""
    from vspipeline import defaults
    p = runner.build_parser()
    got = {a.dest: a.default for a in p._actions}
    expected = {
        "exhaustiveness": defaults.EXHAUSTIVENESS,
        "binding_modes": defaults.BINDING_MODES,
        "energy_range": defaults.ENERGY_RANGE,
        "processes": defaults.PROCESSES,
        "min_valid_affinity": defaults.MIN_VALID_AFFINITY,
        "max_failure_fraction": defaults.MAX_FAILURE_FRACTION,
        "ph": defaults.PH,
        "padding": defaults.PADDING,
        "water_distance": defaults.WATER_DISTANCE,
        "minimize_steps": defaults.MINIMIZE_STEPS,
        "max_compounds": defaults.MAX_COMPOUNDS,
        "mw_min": defaults.CUSTOM_MW_MIN,
        "mw_max": defaults.CUSTOM_MW_MAX,
        "logp_max": defaults.CUSTOM_LOGP_MAX,
        "top_n": defaults.TOP_N,
        "rmsd_threshold": defaults.RMSD_THRESHOLD,
    }
    assert {k: got[k] for k in expected} == expected
    assert runner.DEFAULT_SEED == defaults.SEED
    assert runner.AUTO_FLEX_RADIUS == defaults.AUTO_FLEX_RADIUS


def test_qvina_fallback_needs_quickvina(cli, fake_tools, monkeypatch, capsys):
    monkeypatch.setattr(runner, "QVINA", None)
    argv, out = cli
    assert _exit_code(argv + ["--qvina-fallback"]) == runner.EXIT_USAGE
    assert "needs QuickVina 2" in capsys.readouterr().err


def test_qvina_fallback_is_off_without_the_option():
    """No hidden switch: the fallback equals the primary engine unless
    --qvina-fallback (or the deprecated variable) asks for QuickVina."""
    assert runner.VINA_FALLBACK == runner.VINA_PRIMARY


def test_deprecated_qvina_env_enables_the_option(cli, fake_tools, monkeypatch, caplog):
    monkeypatch.setattr(runner, "QVINA", None)
    monkeypatch.setenv("VS_ENABLE_QVINA_FALLBACK", "1")
    argv, out = cli
    # With QuickVina missing, the (deprecated) request is validated like the option.
    assert _exit_code(argv) == runner.EXIT_USAGE
    assert "deprecated; use --qvina-fallback" in caplog.text
