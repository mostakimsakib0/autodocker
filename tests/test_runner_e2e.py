"""End-to-end run of runner.main() with the real tools (tiny and fast)."""
import csv
import json
import logging
import os
import shutil
import subprocess

import pytest

import runner

pytestmark = pytest.mark.skipif(
    not (shutil.which("vina") and shutil.which("obabel")),
    reason="vina/obabel not available on PATH")


def test_main_end_to_end_writes_manifest(tmp_path, protein_pdb, caplog):
    ligs = tmp_path / "ligs"
    ligs.mkdir()
    subprocess.run([shutil.which("obabel"), "-:CCO ethanol", "--gen3d", "-osdf",
                    "-O", str(ligs / "ethanol.sdf")], check=True, capture_output=True,
                   env={**os.environ, "OMP_NUM_THREADS": "1"})
    out = tmp_path / "out"
    argv = ["-r", protein_pdb, "-l", str(ligs), "-o", str(out), "--no-fpocket",
            "--exhaustiveness", "1", "--binding-modes", "1", "--seed", "3",
            "--no-html-report", "--minimize-steps", "20", "--max-failure-fraction", "1"]
    try:
        with caplog.at_level(logging.INFO):
            runner.main(argv)
    finally:
        runner._configure_logging(False)
    assert "obabel:" in caplog.text and "Open Babel" in caplog.text
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["status"] == "completed"
    assert manifest["seed"] == 3
    assert "Open Babel" in manifest["tools"]["obabel"]["version"]
    assert "Vina" in manifest["tools"]["vina"]["version"]
    assert len(manifest["tools"]["vina"]["sha256"]) == 64
    assert manifest["ligand_preparation"] == {"failed": 0, "superseded": 0}
    assert (out / "ranking.csv").exists()
    meta = json.loads((out / "ligand_metadata.json").read_text())
    assert meta["ethanol"]["minimized"] is True


# ---------------------------------------------------------------------------
# Review fixes, end to end with the real Open Babel and Vina
# ---------------------------------------------------------------------------

def _ethanol_dir(tmp_path):
    ligs = tmp_path / "ligs"
    ligs.mkdir()
    subprocess.run([shutil.which("obabel"), "-:CCO ethanol", "--gen3d", "-osdf",
                    "-O", str(ligs / "ethanol.sdf")], check=True, capture_output=True,
                   env={**os.environ, "OMP_NUM_THREADS": "1"})
    return ligs


def _base_argv(tmp_path, protein_pdb, modes=True):
    argv = ["-r", protein_pdb, "-l", str(_ethanol_dir(tmp_path)), "-o",
            str(tmp_path / "out"), "--exhaustiveness", "1", "--seed", "3",
            "--minimize-steps", "20"]
    return argv + (["--binding-modes", "1"] if modes else [])


def _main(argv):
    """runner.main(argv) -> exit code."""
    try:
        runner.main(argv)
        return runner.EXIT_OK
    except SystemExit as e:
        return e.code
    finally:
        runner._configure_logging(False)


def _manifest(out):
    with open(os.path.join(out, "run_manifest.json"), encoding="utf-8") as f:
        return json.load(f)


def test_explicit_box_needs_no_fpocket_and_stale_outputs_go(tmp_path, protein_pdb,
                                                            monkeypatch, caplog):
    """--grid-center + --grid-size: fpocket is neither required nor run; the
    outputs an earlier run left and this one does not write are removed."""
    from vspipeline import outputs

    def no_fpocket(*args, **kwargs):
        raise AssertionError("fpocket must not run with an explicit box")

    monkeypatch.setattr(runner, "FPOCKET", None)
    monkeypatch.setattr(runner.ProteinPreparation, "detect_pocket", no_fpocket)
    argv = _base_argv(tmp_path, protein_pdb)
    out = tmp_path / "out"
    out.mkdir()
    stale = (outputs.HTML_REPORT, outputs.POCKET_SUMMARY_CSV, outputs.POSE_CLUSTERS_JSON,
             outputs.CONSENSUS_RANKING_CSV, outputs.CONSENSUS_SUMMARY_JSON, outputs.FLEX_PDBQT)
    for name in stale:
        (out / name).write_text("stale")
    with caplog.at_level(logging.INFO):
        code = _main(argv + ["--grid-center", "4.5,3,5", "--grid-size", "20,20,20",
                             "--no-html-report"])
    assert code == runner.EXIT_OK
    assert "fpocket not run" in caplog.text
    for name in stale:
        assert not (out / name).exists(), name
    assert (out / outputs.RANKING_CSV).exists() and (out / outputs.GRID_CONF).exists()
    assert "Pocket summary" not in caplog.text and "HTML report:" not in caplog.text
    m = _manifest(out)
    assert "fpocket" not in m["tools"] and m["flexible_residues"] == []


def test_unbuildable_flexible_residue_is_an_error(tmp_path, protein_pdb, caplog):
    """A residue that exists but cannot be built as a flexible side chain
    (the test receptor has CA atoms only) stops the run (exit 1) instead of
    docking rigidly while claiming flexible residues."""
    argv = _base_argv(tmp_path, protein_pdb)
    with caplog.at_level(logging.INFO):
        code = _main(argv + ["--no-fpocket", "--flexible-residues", "A5"])
    assert code == runner.EXIT_FAILURE
    m = _manifest(tmp_path / "out")
    assert m["status"] == "failed"
    assert "--flexible-residues: could not build A5" in m["error"]
    assert m["flexible_residues"] == []
    assert not (tmp_path / "out" / "ranking.csv").exists()


def test_auto_flexible_without_buildable_residue_docks_rigid(tmp_path, protein_pdb, caplog):
    argv = _base_argv(tmp_path, protein_pdb)
    with caplog.at_level(logging.INFO):
        code = _main(argv + ["--no-fpocket", "--auto-flexible", "2"])
    assert code == runner.EXIT_OK
    assert "docking with a RIGID receptor" in caplog.text
    assert "Flexible residues: none (rigid docking)" in caplog.text
    assert _manifest(tmp_path / "out")["flexible_residues"] == []
    assert not (tmp_path / "out" / "flex.pdbqt").exists()


def _failing_binary(tmp_path, name):
    path = tmp_path / name
    path.write_text("#!/bin/sh\necho 'simulated engine crash' >&2\nexit 1\n")
    path.chmod(0o755)
    return str(path)


@pytest.mark.parametrize("fraction, expected", [
    (None, runner.EXIT_TOO_MANY_FAILURES), ("1.0", runner.EXIT_OK)])
def test_all_ligands_failed_writes_outputs_and_exits_3(tmp_path, protein_pdb, caplog,
                                                       fraction, expected):
    """Every ligand failing is the extreme failure fraction: all outputs are
    written (FAILED rows) and the exit code follows --max-failure-fraction
    (3 by default; 1.0 never fails)."""
    argv = _base_argv(tmp_path, protein_pdb) + [
        "--no-fpocket", "--vina-bin", _failing_binary(tmp_path, "vina")]
    if fraction:
        argv += ["--max-failure-fraction", fraction]
    with caplog.at_level(logging.INFO):
        code = _main(argv)
    assert code == expected
    out = tmp_path / "out"
    with open(out / "ranking.csv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert [(r["Ligand"], r["Status"]) for r in rows] == [("ethanol", "FAILED")]
    for name in ("Top_hits.txt", "metrics.txt", "Results_full.txt", "results_report.html"):
        assert (out / name).exists(), name
    report = (out / "results_report.html").read_text(encoding="utf-8")
    assert "n/a" in report and 'href="#viewer"' not in report
    m = _manifest(out)
    assert m["docking_failures"] == {"failed": 1, "docked": 1}
    assert m["status"] == ("completed" if expected == 0 else "too_many_failures")
    assert "All docking attempts failed" in caplog.text


def test_smina_only_all_failed_exits_3(tmp_path, protein_pdb, monkeypatch):
    monkeypatch.setattr(runner, "SMINA", _failing_binary(tmp_path, "smina"))
    argv = _base_argv(tmp_path, protein_pdb, modes=False) + ["--no-fpocket", "--smina-only"]
    assert _main(argv) == runner.EXIT_TOO_MANY_FAILURES
    out = tmp_path / "out"
    assert (out / "ranking.csv").exists()
    assert _manifest(out)["docking_failures"] == {"failed": 1, "docked": 1}
