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