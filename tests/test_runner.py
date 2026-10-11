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
    monkeypatch.setattr(runner, "FPOCKET", sys.executable)
    for name in ("SMINA", "REDUCE"):
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
    assert "vina" in m["tools"] and m["tools"]["fpocket"]["path"] == sys.executable


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


# ---------------------------------------------------------------------------
# Review fixes: validation of values and combinations (exit 2, nothing written)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("extra, message", [
    # finite numbers only (argparse's float() accepts nan and inf)
    (["--ph", "nan"], "--ph must be a finite number"),
    (["--padding", "nan"], "--padding must be a finite number"),
    (["--padding", "inf"], "--padding must be a finite number"),
    (["--energy-range", "inf"], "--energy-range must be a finite number"),
    (["--min-valid-affinity", "nan"], "--min-valid-affinity must be a finite number"),
    (["--rmsd-threshold", "nan"], "--rmsd-threshold must be a finite number"),
    (["--water-distance", "inf"], "--water-distance must be a finite number"),
    (["--max-failure-fraction", "nan"], "--max-failure-fraction must be a finite number"),
    (["--library", "custom", "--mw-max", "inf"], "--mw-max must be a finite number"),
    (["--library", "custom", "--mw-min", "nan"], "--mw-min must be a finite number"),
    (["--library", "custom", "--logp-max", "nan"], "--logp-max must be a finite number"),
    (["--grid-center", "nan,0,0"], "finite"),
    (["--grid-size", "20,inf,20"], "finite"),
    # seed: Vina/QuickVina/SMINA take a 32-bit int, Vina treats 0 as random,
    # RDKit rejects negative seeds
    (["--seed", str(2**31)], "--seed must be between 1 and 2147483647"),
    (["--seed", "99999999999999"], "--seed must be between 1 and 2147483647"),
    (["--seed", "0"], "--seed must be between 1 and 2147483647"),
    (["--seed", "-1"], "--seed must be between 1 and 2147483647"),
    (["--ph", "14.5"], "--ph must be between 0 and 14"),
    (["--ph", "-1"], "--ph must be between 0 and 14"),
    (["--min-valid-affinity", "-20"], "--min-valid-affinity must be > -20"),
    (["--min-valid-affinity", "-100"], "--min-valid-affinity must be > -20"),
    # options that would be accepted and ignored
    (["--auto-flexible", "3", "--flexible-residues", "A5"], "cannot be combined"),
    (["--flexible-residues", ",,"], "--flexible-residues lists no residue"),
    (["--flexible-residues", "A5,45-"], "invalid residue key(s) 45-"),
    (["--pockets", "1", "--no-fpocket"], "--pockets cannot be combined with --no-fpocket"),
    (["--pockets", "1", "--grid-center", "1,2,3", "--grid-size", "20,20,20"],
     "--pockets cannot be combined with --grid-center and --grid-size"),
    (["--grid-center-frame", "docking"], "--grid-center-frame only applies to --grid-center"),
    (["--grid-center", ""], "three comma-separated values"),
    # --pockets syntax and --chain against the receptor (before ligand prep)
    (["--pockets", "abc"], "--pockets must list fpocket pocket numbers"),
    (["--pockets", "0"], "--pockets must list fpocket pocket numbers"),
    (["--pockets", ","], "--pockets must list fpocket pocket numbers"),
    (["--pockets", "1,x"], "--pockets must list fpocket pocket numbers"),
    (["--chain", "Z"], "Invalid chain(s): Z. Available chains: A, B"),
    (["--chain", "a"], "did you mean 'a' -> 'A'"),
    # flexible residues must exist in the docked chain(s)
    (["--flexible-residues", "A9999"], "no such residue in the docked chain(s) A"),
    (["--flexible-residues", "A5,B80"], "docked chain(s) A of"),  # first chain only
    (["--chain", "B", "--flexible-residues", "A5"], "no such residue in the docked chain(s) B"),
])
def test_review_invalid_options_exit_2(cli, fake_tools, capsys, extra, message):
    argv, out = cli
    assert _exit_code(argv + extra) == runner.EXIT_USAGE
    err = capsys.readouterr().err
    assert message in err and "Traceback" not in err
    assert not out.exists()


@pytest.mark.parametrize("extra", [
    ["--seed", "1"], ["--seed", str(2**31 - 1)], ["--ph", "0"], ["--ph", "14"],
    ["--min-valid-affinity", "-19.5"], ["--pockets", "2,3"], ["--pockets", "auto"],
    ["--chain", "all"], ["--chain", "A,B", "--flexible-residues", "A5,B80"],
    ["--flexible-residues", " A5 , A6 "], ["--auto-flexible", "2"],
    ["--grid-center", "1,2,3", "--grid-center-frame", "docking"],
])
def test_review_valid_options_accepted(cli, fake_tools, monkeypatch, extra):
    monkeypatch.setattr(runner, "_run_pipeline", lambda args, state: runner.EXIT_OK)
    argv, out = cli
    runner.main(argv + extra)
    assert _manifest(out)["status"] == "completed"


def test_vs_seed_out_of_range_exit_2(cli, fake_tools, monkeypatch, capsys):
    monkeypatch.setenv("VS_SEED", "0")
    argv, out = cli
    assert _exit_code(argv) == runner.EXIT_USAGE
    assert "$VS_SEED must be between 1 and" in capsys.readouterr().err
    assert not out.exists()


def test_long_option_abbreviations_are_rejected(cli, fake_tools, capsys):
    """entry.sh must see the same -r/-l/-o/--library as the runner, so
    abbreviated long options (--lig, --out, --mw-ma) are not accepted."""
    argv, out = cli
    assert _exit_code(argv + ["--mw-ma", "600"]) == runner.EXIT_USAGE
    assert "unrecognized arguments: --mw-ma" in capsys.readouterr().err
    receptor, ligands = argv[1], argv[3]
    assert _exit_code(["--rec", receptor, "--lig", ligands, "--out", str(out)]) \
        == runner.EXIT_USAGE
    assert "required" in capsys.readouterr().err
    assert not out.exists()


def test_attached_short_option_values(cli, fake_tools, monkeypatch):
    monkeypatch.setattr(runner, "_run_pipeline", lambda args, state: runner.EXIT_OK)
    argv, out = cli
    runner.main(["-r" + argv[1], "-l" + argv[3], "-o" + str(out)])
    assert _manifest(out)["args"]["ligands"] == argv[3]


# ---------------------------------------------------------------------------
# -o validation (exit 2 with a message, no traceback)
# ---------------------------------------------------------------------------

def test_output_that_is_a_file_exit_2(cli, fake_tools, capsys, tmp_path):
    argv, _ = cli
    afile = tmp_path / "afile"
    afile.write_text("x")
    assert _exit_code(argv[:4] + ["-o", str(afile)]) == runner.EXIT_USAGE
    err = capsys.readouterr().err
    assert "exists and is not a directory" in err and "Traceback" not in err
    assert _exit_code(argv[:4] + ["-o", str(afile / "sub")]) == runner.EXIT_USAGE
    assert "cannot be created" in capsys.readouterr().err
    assert afile.read_text() == "x"


def test_empty_output_exit_2(cli, fake_tools, capsys):
    argv, _ = cli
    assert _exit_code(argv[:4] + ["-o", ""]) == runner.EXIT_USAGE
    assert "-o/--output must not be empty" in capsys.readouterr().err


@pytest.mark.skipif(os.geteuid() == 0, reason="root can write anywhere")
def test_unwritable_output_parent_exit_2(cli, fake_tools, capsys, tmp_path):
    argv, _ = cli
    locked = tmp_path / "locked"
    locked.mkdir()
    locked.chmod(0o500)
    try:
        assert _exit_code(argv[:4] + ["-o", str(locked / "out")]) == runner.EXIT_USAGE
        assert "is not writable" in capsys.readouterr().err
    finally:
        locked.chmod(0o700)


def test_output_creation_failure_after_validation_exit_2(cli, fake_tools, monkeypatch, caplog):
    def fail(path, exist_ok=False):
        raise PermissionError(13, "Permission denied", path)

    monkeypatch.setattr(runner.os, "makedirs", fail)
    argv, out = cli
    with caplog.at_level(logging.ERROR):
        assert _exit_code(argv) == runner.EXIT_USAGE
    assert "-o/--output: cannot create" in caplog.text


def test_tool_probe_failure_is_a_clean_pipeline_error(cli, fake_tools, monkeypatch, caplog):
    def boom(args):
        raise RuntimeError("probe exploded")

    monkeypatch.setattr(runner, "_probe_tools", boom)
    argv, out = cli
    with caplog.at_level(logging.ERROR):
        assert _exit_code(argv) == runner.EXIT_FAILURE
    m = _manifest(out)
    assert m["status"] == "failed" and "probe exploded" in m["error"]


# ---------------------------------------------------------------------------
# fpocket: needed only when pockets are detected; checked up front
# ---------------------------------------------------------------------------

def test_missing_fpocket_exit_2_before_any_work(cli, fake_tools, monkeypatch, caplog):
    monkeypatch.setattr(runner, "FPOCKET", None)
    argv, out = cli
    with caplog.at_level(logging.ERROR):
        assert _exit_code(argv) == runner.EXIT_USAGE
    assert "Required tool(s) not found in PATH: fpocket" in caplog.text
    assert "--no-fpocket" in caplog.text
    assert not out.exists()


@pytest.mark.parametrize("extra", [
    ["--no-fpocket"], ["--grid-center", "1,2,3", "--grid-size", "20,20,20"]])
def test_fpocket_not_needed_without_pocket_detection(cli, fake_tools, monkeypatch, extra):
    monkeypatch.setattr(runner, "FPOCKET", None)
    monkeypatch.setattr(runner, "_run_pipeline", lambda args, state: runner.EXIT_OK)
    argv, out = cli
    runner.main(argv + extra)
    m = _manifest(out)
    assert m["status"] == "completed" and "fpocket" not in m["tools"]


def test_fpocket_still_needed_with_only_one_explicit_triplet(cli, fake_tools, monkeypatch):
    monkeypatch.setattr(runner, "FPOCKET", None)
    argv, _ = cli
    assert _exit_code(argv + ["--grid-center", "1,2,3"]) == runner.EXIT_USAGE


# ---------------------------------------------------------------------------
# run_manifest.json is strict JSON
# ---------------------------------------------------------------------------

def _strict_json(text):
    def reject(constant):
        raise ValueError(f"non-standard JSON constant {constant}")
    return json.loads(text, parse_constant=reject)


def test_run_manifest_is_strict_json(tmp_path):
    path = tmp_path / "run_manifest.json"
    runner.write_run_manifest(str(path), {
        "a": float("nan"), "b": [float("inf"), 1.5], "c": {"d": float("-inf")}})
    m = _strict_json(path.read_text())
    assert m == {"a": "nan", "b": ["inf", 1.5], "c": {"d": "-inf"}}


def test_run_manifest_of_a_run_is_strict_json(cli, fake_tools, monkeypatch):
    def pipeline(args, state):
        state["consensus"] = {"rho": float("nan"), "n": 1}
        return runner.EXIT_OK

    monkeypatch.setattr(runner, "_run_pipeline", pipeline)
    argv, out = cli
    runner.main(argv)
    m = _strict_json((out / "run_manifest.json").read_text())
    assert m["consensus"]["rho"] == "nan"


# ---------------------------------------------------------------------------
# Outputs of an earlier run in the same -o
# ---------------------------------------------------------------------------

def test_stale_derived_outputs_removed_but_docking_results_kept(cli, fake_tools, monkeypatch,
                                                                caplog):
    from vspipeline import outputs
    argv, out = cli
    out.mkdir()
    for name in outputs.DERIVED_OUTPUTS:
        if name == outputs.FPOCKET_OUT_DIR:
            (out / name / "pockets").mkdir(parents=True)
        else:
            (out / name).write_text("stale")
    (out / "docked").mkdir()
    (out / "docked" / "x_out.pdbqt").write_text("pose")
    (out / ".docking_checkpoint.json").write_text("{}")
    (out / "ligand_metadata.json").write_text("{}")
    seen = {}

    def pipeline(args, state):
        seen["left"] = sorted(n for n in outputs.DERIVED_OUTPUTS if (out / n).exists())
        return runner.EXIT_OK

    monkeypatch.setattr(runner, "_run_pipeline", pipeline)
    with caplog.at_level(logging.INFO):
        runner.main(argv + ["--no-html-report"])
    assert seen["left"] == []
    assert "Removed outputs of an earlier run" in caplog.text
    assert outputs.HTML_REPORT in caplog.text
    assert (out / "docked" / "x_out.pdbqt").read_text() == "pose"
    assert (out / ".docking_checkpoint.json").exists()
    assert (out / "run_manifest.json").exists()


def test_derived_outputs_are_defined_once():
    from vspipeline import outputs
    for name in ("results_report.html", "pose_clusters.json", "consensus_ranking.csv",
                 "consensus_summary.json", "pocket_summary.csv", "ranking.csv"):
        assert name in outputs.DERIVED_OUTPUTS
    for kept in (outputs.DOCKED_DIR, outputs.RUN_MANIFEST_JSON, outputs.LIGANDS_DIR):
        assert kept not in outputs.DERIVED_OUTPUTS


# ---------------------------------------------------------------------------
# Messages
# ---------------------------------------------------------------------------

def test_ctrl_c_during_smina_only_docking(cli, fake_tools, monkeypatch, caplog):
    monkeypatch.setattr(runner, "SMINA", sys.executable)

    def interrupted(args, state):
        state["docking_started"] = True
        raise KeyboardInterrupt

    monkeypatch.setattr(runner, "_run_pipeline", interrupted)
    argv, out = cli
    with caplog.at_level(logging.WARNING):
        assert _exit_code(argv + ["--smina-only"]) == runner.EXIT_INTERRUPTED
    assert "interrupted during SMINA docking" in caplog.text
    assert "before docking started" not in caplog.text


@pytest.mark.parametrize("extra", [["--vina-bin", sys.executable], ["--smina-only"]])
def test_qvina_env_conflict_names_the_variable(cli, fake_tools, monkeypatch, capsys, extra):
    monkeypatch.setenv("VS_ENABLE_QVINA_FALLBACK", "1")
    monkeypatch.setattr(runner, "QVINA", sys.executable)
    argv, _ = cli
    assert _exit_code(argv + extra) == runner.EXIT_USAGE
    err = capsys.readouterr().err
    assert "$VS_ENABLE_QVINA_FALLBACK=1 in the environment" in err
    assert "cannot be combined" in err


# ---------------------------------------------------------------------------
# Output-directory lock, SIGTERM, failure fraction, manifest (review fixes)
# ---------------------------------------------------------------------------

def _hold_lock(outdir):
    """Another process holding the output directory's lock (as a running
    screen does); returns the Popen (release with .communicate(""))."""
    import subprocess
    code = ("import sys; sys.path.insert(0, sys.argv[2]);"
            "from vspipeline import checkpoint as c;"
            "c.acquire_lock(sys.argv[1]); print('locked', flush=True); sys.stdin.read()")
    pkg = os.path.join(os.path.dirname(__file__), "..", "autodocker")
    holder = subprocess.Popen([sys.executable, "-c", code, str(outdir), pkg],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    assert holder.stdout.readline().strip() == "locked"
    return holder


def _tree_digest(root):
    """{relative path: (sha256, mtime_ns)} of every file under root."""
    digest = {}
    for dirpath, _, files in os.walk(root):
        for name in files:
            path = os.path.join(dirpath, name)
            with open(path, "rb") as f:
                digest[os.path.relpath(path, root)] = (
                    hashlib.sha256(f.read()).hexdigest(), os.stat(path).st_mtime_ns)
    return digest


@pytest.mark.parametrize("extra", [[], ["--smina-only"]])
def test_locked_output_dir_is_left_untouched(cli, fake_tools, monkeypatch, caplog, extra):
    """A second run into the -o of a running screen used to rewrite ligands/,
    receptor.pdbqt, grid.conf and run_manifest.json before it hit the lock
    (and --smina-only never locked). It now stops with exit 2 before writing
    anything, in every mode."""
    monkeypatch.setattr(runner, "SMINA", sys.executable)

    def must_not_run(args, state):
        raise AssertionError("the pipeline must not start")

    monkeypatch.setattr(runner, "_run_pipeline", must_not_run)
    argv, out = cli
    out.mkdir()
    for name in ("run_manifest.json", "grid.conf", "ranking.csv"):
        (out / name).write_text("from the running screen\n")
    (out / "docked").mkdir()
    (out / "docked" / "a_out.pdbqt").write_text("pose\n")
    holder = _hold_lock(out)
    try:
        before = _tree_digest(out)
        with caplog.at_level(logging.ERROR):
            assert _exit_code(argv + extra) == runner.EXIT_USAGE
        assert _tree_digest(out) == before
    finally:
        holder.communicate("")
    assert "in use by another AutoDocker run" in caplog.text


def test_lock_is_held_for_the_whole_run_and_released(cli, fake_tools, monkeypatch):
    from vspipeline import checkpoint as ckpt
    seen = {}

    def pipeline(args, state):
        seen["held"] = ckpt.holds_lock(args.output)
        # The checkpoint re-uses the run's lock and must not release it.
        cp = runner.DockingCheckpoint(args.output)
        cp.close()
        seen["after_close"] = ckpt.holds_lock(args.output)
        return runner.EXIT_OK

    monkeypatch.setattr(runner, "_run_pipeline", pipeline)
    argv, out = cli
    runner.main(argv)
    assert seen == {"held": True, "after_close": True}
    assert not ckpt.holds_lock(str(out))


def test_sigterm_exits_143_with_manifest(cli, fake_tools, monkeypatch, caplog):
    """SIGTERM (docker stop) takes the interrupt path: manifest written with
    status interrupted, exit 143, lock released, previous handler restored."""
    import signal
    from vspipeline import checkpoint as ckpt
    previous = signal.getsignal(signal.SIGTERM)

    def pipeline(args, state):
        state["docking_started"] = True
        os.kill(os.getpid(), signal.SIGTERM)
        raise AssertionError("SIGTERM must interrupt the run")

    monkeypatch.setattr(runner, "_run_pipeline", pipeline)
    argv, out = cli
    with caplog.at_level(logging.WARNING):
        assert _exit_code(argv) == runner.EXIT_TERMINATED == 143
    m = _manifest(out)
    assert m["status"] == "interrupted" and m["error"] == "SIGTERM"
    assert m["exit_code"] == 143
    assert "terminated (SIGTERM)" in caplog.text
    assert signal.getsignal(signal.SIGTERM) == previous
    assert not ckpt.holds_lock(str(out))


def test_failure_counts_include_preparation_failures():
    """README: exit 3 "if more than this fraction of ligands fail": ligands
    that could not be prepared are failures too."""
    counts = runner._failure_counts(n_prep_failed=2, n_dock_failed=1, n_docked=3,
                                    max_fraction=0.5)
    assert counts == {"preparation_failed": 2, "docking_failed": 1, "failed": 3,
                      "ligands": 5, "fraction": 0.6, "max_failure_fraction": 0.5}
    assert runner._failure_exit_code(counts["failed"], counts["ligands"], 0.5) \
        == runner.EXIT_TOO_MANY_FAILURES
    # Docking alone (1/3) would have passed.
    assert runner._failure_exit_code(1, 3, 0.5) == runner.EXIT_OK


def test_run_manifest_write_is_atomic(tmp_path, monkeypatch):
    path = tmp_path / "run_manifest.json"
    runner.write_run_manifest(str(path), {"a": 1})
    old = path.read_text()

    def boom(*a, **k):
        raise OSError("disk full")

    monkeypatch.setattr(runner.json, "dump", boom)
    with pytest.raises(OSError):
        runner.write_run_manifest(str(path), {"a": 2})
    assert path.read_text() == old
    assert os.listdir(tmp_path) == ["run_manifest.json"]


def test_threads_per_process_recorded(cli, fake_tools, monkeypatch):
    monkeypatch.setattr(runner, "available_cpus", lambda: 8)
    seen = {}

    def fake_args(processes, threads):
        return type("A", (), {"processes": processes, "threads": threads})()

    assert runner._threads_per_process(fake_args(4, None)) == 2
    assert runner._threads_per_process(fake_args(1, None)) is None
    assert runner._threads_per_process(fake_args(4, 3)) == 3

    def pipeline(args, state):
        state["threads_per_process"] = {"threads": 2, "source": "x", "processes": 4}
        seen["ok"] = True
        return runner.EXIT_OK

    monkeypatch.setattr(runner, "_run_pipeline", pipeline)
    argv, out = cli
    runner.main(argv)
    assert _manifest(out)["threads_per_process"]["threads"] == 2
