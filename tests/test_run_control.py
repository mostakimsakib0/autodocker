"""Run control end to end, with the real Open Babel, Vina and SMINA: resuming
from another working directory or a moved output directory, SIGTERM, and
ligand-preparation failures in --max-failure-fraction."""
import json
import os
import shutil
import signal
import subprocess
import sys
import time

import pytest

import runner

pytestmark = pytest.mark.skipif(
    not (shutil.which("vina") and shutil.which("obabel")),
    reason="vina/obabel not available on PATH")

RUNNER = os.path.join(os.path.dirname(__file__), "..", "autodocker", "runner.py")
ENV = {**os.environ, "OMP_NUM_THREADS": "1"}


def _ligands(tmp_path, smiles=(("CCO", "ethanol"), ("CCCO", "propanol"))):
    ligs = tmp_path / "ligs"
    ligs.mkdir(exist_ok=True)
    for smi, name in smiles:
        subprocess.run([shutil.which("obabel"), f"-:{smi} {name}", "--gen3d", "-osdf",
                        "-O", str(ligs / f"{name}.sdf")], check=True, capture_output=True,
                       env=ENV)
    return ligs


def _argv(protein_pdb, ligs, out):
    return ["-r", protein_pdb, "-l", str(ligs), "-o", str(out), "--no-fpocket",
            "--exhaustiveness", "1", "--binding-modes", "2", "--seed", "3",
            "--minimize-steps", "20", "--no-html-report"]


def _run(argv, cwd, timeout=600):
    return subprocess.run([sys.executable, os.path.abspath(RUNNER)] + argv, cwd=cwd,
                          capture_output=True, text=True, env=ENV, timeout=timeout)


def _read(path):
    with open(path, "rb") as f:
        return f.read()


def _manifest(out):
    with open(os.path.join(out, "run_manifest.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.mark.skipif(not shutil.which("smina"), reason="smina not available")
def test_resume_from_another_cwd_and_after_a_move(tmp_path, protein_pdb):
    """dock_file/log_file were stored relative to the working directory: a
    resume from elsewhere wrote a wrong Docked_File and --consensus crashed.
    A resumed run (other cwd, absolute -o; then the directory moved) gives
    the same ranking.csv and consensus_ranking.csv as a fresh run."""
    ligs = _ligands(tmp_path)
    work = tmp_path / "work"
    work.mkdir()
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    extra = ["--consensus", "-p", "2"]

    fresh = tmp_path / "fresh"
    assert _run(_argv(protein_pdb, ligs, fresh) + extra, cwd=tmp_path).returncode == 0

    first = _run(_argv(protein_pdb, ligs, "rel_out") + extra, cwd=work)
    assert first.returncode == 0, first.stderr[-3000:]
    out = work / "rel_out"
    resumed = _run(_argv(protein_pdb, ligs, out) + extra, cwd=elsewhere)
    assert resumed.returncode == 0, resumed.stderr[-3000:]
    assert "All ligands already docked (checkpoint)" in resumed.stderr
    for name in ("ranking.csv", "consensus_ranking.csv", "Top_hits.txt", "metrics.txt"):
        assert _read(out / name) == _read(fresh / name), name

    moved = tmp_path / "moved" / "out"
    moved.parent.mkdir()
    os.rename(out, moved)
    again = _run(_argv(protein_pdb, ligs, moved) + extra, cwd=elsewhere)
    assert again.returncode == 0, again.stderr[-3000:]
    assert "All ligands already docked (checkpoint)" in again.stderr
    for name in ("ranking.csv", "consensus_ranking.csv"):
        assert _read(moved / name) == _read(fresh / name), name


def _descendants(pid):
    """PIDs of every process below ``pid`` (from /proc)."""
    children = {}
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open(f"/proc/{entry}/stat", encoding="utf-8") as f:
                ppid = int(f.read().rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            continue
        children.setdefault(ppid, []).append(int(entry))
    found, todo = [], [pid]
    while todo:
        for child in children.get(todo.pop(), []):
            found.append(child)
            todo.append(child)
    return found


def _alive(pid):
    try:
        with open(f"/proc/{pid}/stat", encoding="utf-8") as f:
            return f.read().rsplit(")", 1)[1].split()[0] != "Z"
    except OSError:
        return False


@pytest.mark.skipif(not os.path.isdir("/proc"), reason="needs /proc")
@pytest.mark.parametrize("processes", ["1", "2"])
def test_sigterm_stops_workers_and_engines(tmp_path, protein_pdb, processes):
    """SIGTERM (docker stop) used to kill the runner without a manifest while
    its workers and their Vina processes kept writing into the (now unlocked)
    output directory. Now: exit 143, status "interrupted", no process left."""
    vina = tmp_path / "slow_vina"
    vina.write_text("#!/bin/sh\ncase \"$1\" in --help*|--version) exec vina \"$@\";; esac\n"
                    "exec sleep 300\n")
    vina.chmod(0o755)
    ligs = _ligands(tmp_path)
    out = tmp_path / "out"
    argv = _argv(protein_pdb, ligs, out) + ["--vina-bin", str(vina), "-p", processes]
    proc = subprocess.Popen([sys.executable, os.path.abspath(RUNNER)] + argv,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                            env=ENV, cwd=tmp_path)
    try:
        deadline = time.time() + 120
        engines = []
        while time.time() < deadline and len(engines) < int(processes):
            time.sleep(0.2)
            engines = [p for p in _descendants(proc.pid) if _cmdline(p).startswith("sleep")]
        assert len(engines) == int(processes), "the engines did not start"
        tree = _descendants(proc.pid)
        proc.send_signal(signal.SIGTERM)
        output, _ = proc.communicate(timeout=60)
    finally:
        if proc.poll() is None:
            proc.kill()
    assert proc.returncode == runner.EXIT_TERMINATED, output[-3000:]
    m = _manifest(out)
    assert m["status"] == "interrupted" and m["error"] == "SIGTERM"
    time.sleep(0.5)
    assert [p for p in tree if _alive(p)] == []
    # The lock is free again.
    from vspipeline import checkpoint as ckpt
    assert ckpt.acquire_lock(str(out))
    ckpt.release_lock(str(out))


def _cmdline(pid):
    try:
        with open(f"/proc/{pid}/cmdline", "rb") as f:
            return f.read().replace(b"\0", b" ").decode(errors="replace")
    except OSError:
        return ""


def test_preparation_failures_count_toward_max_failure_fraction(tmp_path, protein_pdb):
    """A ligand that cannot be prepared is a failed ligand (README: exit 3
    "if more than this fraction of ligands fail")."""
    ligs = _ligands(tmp_path, (("CCO", "ethanol"),))
    (ligs / "broken.sdf").write_text("not a molecule\n")
    out = tmp_path / "out"
    strict = _run(_argv(protein_pdb, ligs, out) + ["--max-failure-fraction", "0.4"],
                  cwd=tmp_path)
    assert strict.returncode == runner.EXIT_TOO_MANY_FAILURES, strict.stderr[-3000:]
    m = _manifest(out)
    assert m["failures"] == {"preparation_failed": 1, "docking_failed": 0, "failed": 1,
                             "ligands": 2, "fraction": 0.5, "max_failure_fraction": 0.4}
    assert m["docking_failures"] == {"failed": 0, "docked": 1}
    assert m["threads_per_process"] == {"threads": None, "processes": 1,
                                        "source": "engine default (all cores)"}
    lenient = _run(_argv(protein_pdb, ligs, out), cwd=tmp_path)   # 0.5 is not > 0.5
    assert lenient.returncode == runner.EXIT_OK, lenient.stderr[-3000:]
