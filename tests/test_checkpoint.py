import json
import os

import pytest

import runner


def test_checkpoint_roundtrip(outdir):
    cp = runner.DockingCheckpoint(outdir)
    cp.save_result("ligA", -7.5, {"status": "OK"})
    cp.save_result("ligB", -6.2, {"status": "OK"})

    cp2 = runner.DockingCheckpoint(outdir)
    assert cp2.is_completed("ligA")
    assert cp2.is_completed("ligB")
    assert not cp2.is_completed("ligC")

    results = dict(cp2.get_results())
    assert results == {"ligA": -7.5, "ligB": -6.2}


def test_checkpoint_corrupt_json(outdir):
    with open(runner.DockingCheckpoint(outdir).checkpoint_file, "w") as f:
        f.write("{this is not json")
    cp = runner.DockingCheckpoint(outdir)
    assert cp.completed == {}
    assert cp.is_completed("ligA") is False


def test_checkpoint_missing_file(outdir):
    cp = runner.DockingCheckpoint(outdir)
    assert cp.completed == {}


def test_checkpoint_get_metrics(outdir):
    cp = runner.DockingCheckpoint(outdir)
    cp.save_result("ligA", -7.5, {"status": "OK", "simscore": 0.8})
    metrics = cp.get_metrics()
    assert metrics["ligA"]["status"] == "OK"


def test_no_resume_clears_checkpoint(outdir):
    cp = runner.DockingCheckpoint(outdir)
    cp.save_result("ligA", -7.5, {})
    assert os.path.exists(cp.checkpoint_file)
    cp_file = cp.checkpoint_file
    os.remove(cp_file)
    assert not os.path.exists(cp_file)


def test_resume_skips_completed(outdir, charged_pdbqt):
    lig = charged_pdbqt
    cp = runner.DockingCheckpoint(outdir)
    cp.save_result("charged", -7.5, {})
    ligands_todo = [
        lig for lig in [lig]
        if not cp.is_completed("charged")
    ]
    assert ligands_todo == []


def test_resume_with_only_failing_ligands_is_partial_not_fatal(
        tmp_path, charged_pdbqt, zero_charge_pdbqt, monkeypatch, outdir):
    """Regression: a resumed run whose remaining ligands are permanent
    outliers (e.g. unrealistic affinities) must still succeed when valid
    scores are already checkpointed. Before the fix the old guard
    `if ligands_todo and success_count == 0` raised, treating a partial
    resume as a full pipeline failure."""
    grid = tmp_path / "grid.txt"
    grid.write_text("center_x=1.0\ncenter_y=2.0\ncenter_z=3.0\n"
                    "size_x=20.0\nsize_y=20.0\nsize_z=20.0\n")

    # Simulate a docking binary that always fails (no valid score produced).
    def _failing_run(cmd, capture=False):
        raise RuntimeError("docking crashed for this ligand")

    monkeypatch.setattr(runner, "run", _failing_run)

    # First "run": dock the good ligand into the checkpoint so a valid
    # score exists. Use a ligand that passes pre-flight checks.
    prep = runner.DockingCheckpoint(outdir)
    prep.save_result("good", -8.4, {"status": "OK", "vina_tool": "vina q"})

    # Second "run": only the failing ligand remains to be docked.
    result, checkpoint, metrics = runner.dock_all(
        charged_pdbqt, [zero_charge_pdbqt, charged_pdbqt],
        str(grid), outdir,
        num_processes=1, resume=True, vina_params={"seed": 42,
                                                   "extra_args": []})
    scores = [s for _, s in result if s is not None]
    assert scores == [-8.4]  # cached success retained; failure tolerated
    assert checkpoint.completed["good"]["score"] == -8.4