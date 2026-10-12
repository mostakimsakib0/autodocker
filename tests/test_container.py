"""The container entry chain: aclmgr's environment forwarding and entry.sh's
reading of the command line (with a stub in place of runner.py)."""
import os
import re
import shutil
import signal
import subprocess
import time

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
ACLMGR_C = os.path.join(ROOT, "tools", "aclmgr.c")
ENTRY_SH = os.path.join(ROOT, "autodocker", "entry.sh")
README = os.path.join(ROOT, "README.md")


def _read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def _forwarded_env():
    """Names in aclmgr.c's forwarded_env[] list."""
    source = _read(ACLMGR_C)
    body = re.search(r"forwarded_env\[\]\s*=\s*\{(.*?)\};", source, re.S).group(1)
    body = re.sub(r"/\*.*?\*/", "", body, flags=re.S)
    return set(re.findall(r'"([A-Za-z_][A-Za-z0-9_]*)"', body))


def _forwarded(name):
    # aclmgr also forwards every VS_* variable by prefix.
    return name.startswith("VS_") or name in _forwarded_env()


def _entry_sh_env_reads():
    """Variables entry.sh reads from the environment: under ``set -u`` every
    such read has a default, ``${NAME:-...}`` / ``${NAME-...}``."""
    return set(re.findall(r"\$\{([A-Z_][A-Z0-9_]*):?-", _read(ENTRY_SH)))


def _readme_env_table():
    text = _read(README)
    section = text.split("## Environment variables", 1)[1].split("\n## ", 1)[0]
    names = set()
    for line in section.splitlines():
        if line.startswith("|") and not line.startswith(("| Variable", "|---")):
            names.update(re.findall(r"`([A-Z_][A-Z0-9_]*)`", line.split("|")[1]))
    return names


def test_forwarded_env_parsed():
    names = _forwarded_env()
    assert {"ADMET", "PDB", "HTTP_PROXY"} <= names


def test_aclmgr_forwards_every_variable_entry_sh_reads():
    reads = _entry_sh_env_reads()
    assert {"LIPINSKI", "ADMET", "CHAIN", "PDB"} <= reads  # the parse finds them
    missing = sorted(name for name in reads if not _forwarded(name))
    assert not missing, f"aclmgr.c drops {missing} (setpriv --reset-env)"


def test_aclmgr_forwards_every_documented_variable():
    documented = _readme_env_table()
    assert {"LIPINSKI", "INPUT", "VS_SEED", "PUBCHEM_REST"} <= documented
    missing = sorted(name for name in documented if not _forwarded(name))
    assert not missing, f"aclmgr.c drops documented variable(s) {missing}"


def test_aclmgr_forwards_what_the_manifest_records():
    from vspipeline import manifest
    source = open(manifest.__file__, encoding="utf-8").read()
    snapshot = re.search(r"def _environment_snapshot.*?\n\n", source, re.S).group(0)
    names = set(re.findall(r'"([A-Z_][A-Z0-9_]*)"', snapshot))
    assert "OMP_NUM_THREADS" in names
    assert not sorted(n for n in names if not _forwarded(n))


@pytest.mark.skipif(shutil.which("gcc") is None, reason="gcc not available")
def test_aclmgr_compiles(tmp_path):
    proc = subprocess.run(["gcc", ACLMGR_C, "-o", str(tmp_path / "aclmgr")],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


# ---------------------------------------------------------------------------
# entry.sh with a stub runner
# ---------------------------------------------------------------------------

STUB_RUNNER = """#!/bin/bash
printf 'STUB ARGS:'; printf ' [%s]' "$@"; echo
exit "${STUB_RC:-0}"
"""


@pytest.fixture
def entry(tmp_path):
    """Run a copy of entry.sh whose workspace is tmp_path/ws; returns
    run(*args, env=...) -> (exit code, output)."""
    if shutil.which("bash") is None:
        pytest.skip("bash not available")
    bindir = tmp_path / "bin"
    bindir.mkdir()
    script = _read(ENTRY_SH).replace('WORK="/workspace"', f'WORK="{tmp_path}/ws"')
    script = script.replace('LIGS="${LIGS:-/workspace/ligs}"', 'LIGS="${LIGS:-$WORK/ligs}"')
    (bindir / "entry.sh").write_text(script, encoding="utf-8")
    (bindir / "runner.py").write_text(STUB_RUNNER, encoding="utf-8")
    for f in ("entry.sh", "runner.py"):
        (bindir / f).chmod(0o755)
    ws = tmp_path / "ws"
    (ws / "ligs").mkdir(parents=True)
    (ws / "protein.pdb").write_text("ATOM\n")
    (ws / "ligs" / "a.sdf").write_text("x\n")
    (tmp_path / "other").mkdir()
    (tmp_path / "other" / "b.pdbqt").write_text("x\n")

    def environment(env=None):
        base = {k: v for k, v in os.environ.items()
                if k not in ("INPUT", "LIGS", "PDB", "LIPINSKI", "ADMET")}
        return {**base, **(env or {})}

    def run(*args, env=None):
        proc = subprocess.run(["bash", str(bindir / "entry.sh"), *args],
                              env=environment(env), capture_output=True,
                              text=True, timeout=60)
        return proc.returncode, proc.stdout + proc.stderr

    run.ws = ws
    run.tmp = tmp_path
    run.bindir = bindir
    run.environment = environment
    return run


def _stub_args(output):
    line = next(line for line in output.splitlines() if line.startswith("STUB ARGS:"))
    return re.findall(r"\[([^\]]*)\]", line)


@pytest.mark.parametrize("flag", ["--help", "-h", "--version"])
def test_entry_help_and_version_need_no_workspace(entry, flag):
    (entry.ws / "protein.pdb").unlink()
    code, out = entry(flag)
    assert code == 0
    assert _stub_args(out) == [flag]


@pytest.mark.parametrize("form", [
    ["-l{d}"], ["-l={d}"], ["--ligands={d}"], ["--ligands", "{d}"], ["-vl{d}"]])
def test_entry_reads_every_ligand_option_form(entry, form):
    other = str(entry.tmp / "other")
    code, out = entry(*[f.format(d=other) for f in form])
    assert code == 0
    assert f"Batch mode: scanning {other}" in out
    assert "PDBQT: 1" in out


@pytest.mark.parametrize("form", [["-o{d}"], ["--output={d}"], ["-o", "{d}"]])
def test_entry_banner_uses_the_runner_output(entry, form):
    outdir = str(entry.tmp / "results")
    code, out = entry(*[f.format(d=outdir) for f in form])
    assert code == 0
    assert f"Results saved to: {outdir}" in out
    assert os.path.isdir(outdir)


@pytest.mark.parametrize("form", [["-r{p}"], ["--receptor={p}"]])
def test_entry_reads_attached_receptor(entry, form):
    pdb = entry.tmp / "other" / "rec.pdb"
    pdb.write_text("ATOM\n")
    (entry.ws / "protein.pdb").unlink()
    code, out = entry(*[f.format(p=pdb) for f in form])
    assert code == 0
    assert "Protein loaded: rec.pdb" in out


@pytest.mark.parametrize("form", [["--library=fda"], ["--library", "custom"]])
def test_entry_detects_pubchem_mode(entry, form):
    for f in (entry.ws / "ligs").iterdir():
        f.unlink()
    code, out = entry(*form)
    assert code == 0
    assert "PubChem library mode" in out


def test_entry_input_is_ignored_in_pubchem_mode(entry):
    single = entry.ws / "ligs" / "a.sdf"
    code, out = entry("--library", "fda", env={"INPUT": str(single)})
    assert code == 0
    assert f"INPUT={single} is ignored with --library fda" in out
    assert "PubChem library mode" in out
    assert single.is_file()  # never mkdir'ed over or touched


@pytest.mark.parametrize("env, args, message", [
    ({"INPUT": "/nonexistent/x.sdf"}, [], "INPUT is set but the file does not exist"),
    ({"PDB": "/nonexistent/p.pdb"}, [], "Protein file not found"),
    ({"LIGS": "/nonexistent/ligs"}, [], "Ligands directory not found"),
    ({}, ["-l", "/nonexistent/ligs"], "Ligands directory not found"),
    ({"INPUT": "{ws}/protein.txt"}, [], "Unsupported ligand file type"),
    ({"LIGS": "{tmp}/empty"}, [], "No ligand files"),
    ({}, ["-o", "{ws}/protein.pdb/sub"], "Cannot create the output directory"),
])
def test_entry_invalid_input_exits_2(entry, env, args, message):
    (entry.ws / "protein.txt").write_text("x\n")
    (entry.tmp / "empty").mkdir()
    fmt = {"ws": entry.ws, "tmp": entry.tmp}
    env = {k: v.format(**fmt) for k, v in env.items()}
    code, out = entry(*[a.format(**fmt) for a in args], env=env)
    assert code == 2, out
    assert message in out
    assert "STUB ARGS" not in out


def test_entry_passes_runner_exit_code(entry):
    code, out = entry(env={"STUB_RC": "3"})
    assert code == 3
    assert "too many ligands failed" in out
    code, out = entry(env={"STUB_RC": "2"})
    assert code == 2


def test_entry_lipinski_switch(entry):
    code, out = entry(env={"LIPINSKI": "0"})
    assert code == 0 and "--no-lipinski" in _stub_args(out)


# A runner that reports when it is running and records the signal it gets.
SIGNAL_STUB_RUNNER = """#!/bin/bash
trap 'echo TERM > "$STUB_DIR/signal"; exit 143' TERM
trap 'echo INT > "$STUB_DIR/signal"; exit 130' INT
touch "$STUB_DIR/running"
for _ in $(seq 600); do sleep 0.1; done
exit 0
"""


@pytest.mark.parametrize("sig, code", [(signal.SIGTERM, 143), (signal.SIGINT, 130)])
def test_entry_forwards_signals_to_the_runner(entry, sig, code):
    """docker stop sends SIGTERM to entry.sh (through aclmgr); bash would not
    pass it on to a foreground child, so the runner never saw it. entry.sh
    forwards SIGTERM and SIGINT to the runner as SIGTERM (a background child
    ignores SIGINT) and exits with 143 or 130."""
    runner_stub = entry.bindir / "runner.py"
    runner_stub.write_text(SIGNAL_STUB_RUNNER, encoding="utf-8")
    stub_dir = entry.tmp / "stub"
    stub_dir.mkdir()
    proc = subprocess.Popen(["bash", str(entry.bindir / "entry.sh")],
                            env=entry.environment({"STUB_DIR": str(stub_dir)}),
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    try:
        deadline = time.monotonic() + 30
        while not (stub_dir / "running").exists():
            assert proc.poll() is None and time.monotonic() < deadline, "runner did not start"
            time.sleep(0.05)
        proc.send_signal(sig)
        out, _ = proc.communicate(timeout=30)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.communicate()
    assert (stub_dir / "signal").read_text().strip() == "TERM"
    assert proc.returncode == code, out
    assert "Pipeline Completed" not in out
