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
    from vspipeline import docking
    grid = tmp_path / "grid.txt"
    grid.write_text("center_x=1.0\ncenter_y=2.0\ncenter_z=3.0\n"
                    "size_x=20.0\nsize_y=20.0\nsize_z=20.0\n")
    ligands = [zero_charge_pdbqt, charged_pdbqt]
    good = runner.ligand_name(charged_pdbqt)

    def first_run(args):
        name = runner.ligand_name(args[1])
        if name == good:
            return name, -8.4, {"status": "OK"}
        return name, None, {"status": "FAILED"}

    monkeypatch.setattr(docking, "dock_ligand", first_run)
    runner.dock_all(charged_pdbqt, ligands, str(grid), outdir, num_processes=1,
                    resume=True, vina_params={"seed": 42, "extra_args": []})

    # Second run: the cached ligand is skipped and the other still fails.
    docked = []

    def second_run(args):
        docked.append(runner.ligand_name(args[1]))
        return runner.ligand_name(args[1]), None, {"status": "FAILED"}

    monkeypatch.setattr(docking, "dock_ligand", second_run)
    result, checkpoint, metrics = runner.dock_all(
        charged_pdbqt, ligands, str(grid), outdir, num_processes=1,
        resume=True, vina_params={"seed": 42, "extra_args": []})
    assert docked == [runner.ligand_name(zero_charge_pdbqt)]
    scores = [s for _, s in result if s is not None]
    assert scores == [-8.4]  # cached success retained; failure tolerated
    assert checkpoint.completed[good]["score"] == -8.4


def _ok_run(args):
    return runner.ligand_name(args[1]), -7.0, {"status": "OK"}


def test_checkpoint_invalidated_when_parameters_change(
        tmp_path, charged_pdbqt, monkeypatch, outdir):
    from vspipeline import docking
    grid = tmp_path / "grid.txt"
    grid.write_text("center_x=1.0\ncenter_y=2.0\ncenter_z=3.0\n"
                    "size_x=20.0\nsize_y=20.0\nsize_z=20.0\n")
    monkeypatch.setattr(docking, "dock_ligand", _ok_run)
    runner.dock_all(charged_pdbqt, [charged_pdbqt], str(grid), outdir,
                    vina_params={"seed": 42, "exhaustiveness": 8})
    docked = []

    def tracking(args):
        docked.append(args[1])
        return _ok_run(args)

    monkeypatch.setattr(docking, "dock_ligand", tracking)
    # Same settings: cached. Different exhaustiveness: docked again.
    runner.dock_all(charged_pdbqt, [charged_pdbqt], str(grid), outdir,
                    vina_params={"seed": 42, "exhaustiveness": 8})
    assert docked == []
    runner.dock_all(charged_pdbqt, [charged_pdbqt], str(grid), outdir,
                    vina_params={"seed": 42, "exhaustiveness": 32})
    assert docked == [charged_pdbqt]


def test_checkpoint_rechecks_edited_ligand_and_drops_removed_ones(
        tmp_path, charged_pdbqt, monkeypatch, outdir):
    import shutil
    from vspipeline import docking
    grid = tmp_path / "grid.txt"
    grid.write_text("center_x=1.0\ncenter_y=2.0\ncenter_z=3.0\n"
                    "size_x=20.0\nsize_y=20.0\nsize_z=20.0\n")
    other = str(tmp_path / "other.pdbqt")
    shutil.copy(charged_pdbqt, other)
    monkeypatch.setattr(docking, "dock_ligand", _ok_run)
    runner.dock_all(charged_pdbqt, [charged_pdbqt, other], str(grid), outdir,
                    vina_params={"seed": 42})
    with open(other, "a") as f:
        f.write("REMARK edited\n")
    docked = []

    def tracking(args):
        docked.append(args[1])
        return _ok_run(args)

    monkeypatch.setattr(docking, "dock_ligand", tracking)
    result, _, _ = runner.dock_all(charged_pdbqt, [other], str(grid), outdir,
                                   vina_params={"seed": 42})
    assert docked == [other]  # content changed, so it is docked again
    assert [n for n, _ in result] == ["other"]  # the removed ligand is not reported


def test_results_are_checkpointed_as_they_arrive(
        tmp_path, charged_pdbqt, monkeypatch, outdir):
    import shutil
    import pytest
    from vspipeline import docking
    grid = tmp_path / "grid.txt"
    grid.write_text("center_x=1.0\ncenter_y=2.0\ncenter_z=3.0\n"
                    "size_x=20.0\nsize_y=20.0\nsize_z=20.0\n")
    ligs = []
    for n in ("a", "b", "c"):
        ligs.append(str(tmp_path / f"{n}.pdbqt"))
        shutil.copy(charged_pdbqt, ligs[-1])

    def interrupted(args):
        if args[1].endswith("c.pdbqt"):
            raise KeyboardInterrupt
        return _ok_run(args)

    monkeypatch.setattr(docking, "dock_ligand", interrupted)
    with pytest.raises(KeyboardInterrupt):
        runner.dock_all(charged_pdbqt, ligs, str(grid), outdir, vina_params={"seed": 42})
    assert set(runner.DockingCheckpoint(outdir).completed) == {"a", "b"}

def test_failed_redock_of_edited_ligand_drops_stale_score(
        tmp_path, charged_pdbqt, monkeypatch, outdir):
    import shutil
    from vspipeline import docking
    grid = tmp_path / "grid.txt"
    grid.write_text("center_x=1.0\ncenter_y=2.0\ncenter_z=3.0\n"
                    "size_x=20.0\nsize_y=20.0\nsize_z=20.0\n")
    lig = str(tmp_path / "x.pdbqt")
    keep = str(tmp_path / "y.pdbqt")
    shutil.copy(charged_pdbqt, lig)
    shutil.copy(charged_pdbqt, keep)
    monkeypatch.setattr(docking, "dock_ligand", _ok_run)
    runner.dock_all(charged_pdbqt, [lig, keep], str(grid), outdir, vina_params={"seed": 42})
    with open(lig, "a") as f:
        f.write("REMARK edited\n")
    monkeypatch.setattr(docking, "dock_ligand",
                        lambda a: ("x", None, {"status": "FAILED"}))
    result, _, _ = runner.dock_all(charged_pdbqt, [lig, keep], str(grid), outdir,
                                   vina_params={"seed": 42})
    # y stays cached; x is reported once, as failed, not also with its old score
    assert sorted(result, key=lambda r: r[0]) == [("x", None), ("y", -7.0)]


def test_rejected_score_is_cached(tmp_path, charged_pdbqt, monkeypatch, outdir):
    import shutil
    from vspipeline import docking
    grid = tmp_path / "grid.txt"
    grid.write_text("center_x=1.0\ncenter_y=2.0\ncenter_z=3.0\n"
                    "size_x=20.0\nsize_y=20.0\nsize_z=20.0\n")
    weak = str(tmp_path / "weak.pdbqt")
    good = str(tmp_path / "good.pdbqt")
    shutil.copy(charged_pdbqt, weak)
    shutil.copy(charged_pdbqt, good)
    calls = []

    def dock(args):
        calls.append(runner.ligand_name(args[1]))
        if args[1] == weak:
            return "weak", None, {"status": "REJECTED (-0.50)", "raw_score": -0.5}
        return _ok_run(args)

    monkeypatch.setattr(docking, "dock_ligand", dock)
    for _ in range(2):
        result, _, metrics = runner.dock_all(charged_pdbqt, [weak, good], str(grid),
                                             outdir, vina_params={"seed": 42})
        assert sorted(result, key=lambda r: r[0]) == [("good", -7.0), ("weak", None)]
        assert metrics["weak"]["status"].startswith("REJECTED")
    assert sorted(calls) == ["good", "weak"]  # the second run docked nothing
