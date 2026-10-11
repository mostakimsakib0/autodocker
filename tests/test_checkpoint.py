import json
import os

import pytest

from vspipeline import checkpoint as checkpoint_mod
from vspipeline import utils
def test_checkpoint_roundtrip(outdir):
    cp = checkpoint_mod.DockingCheckpoint(outdir)
    cp.save_result("ligA", -7.5, {"status": "OK"})
    cp.save_result("ligB", -6.2, {"status": "OK"})

    cp2 = checkpoint_mod.DockingCheckpoint(outdir)
    assert cp2.is_completed("ligA")
    assert cp2.is_completed("ligB")
    assert not cp2.is_completed("ligC")

    results = dict(cp2.get_results())
    assert results == {"ligA": -7.5, "ligB": -6.2}


def test_checkpoint_corrupt_json(outdir):
    with open(checkpoint_mod.DockingCheckpoint(outdir).checkpoint_file, "w", encoding="utf-8") as f:
        f.write("{this is not json")
    cp = checkpoint_mod.DockingCheckpoint(outdir)
    assert cp.completed == {}
    assert cp.is_completed("ligA") is False


def test_checkpoint_missing_file(outdir):
    cp = checkpoint_mod.DockingCheckpoint(outdir)
    assert cp.completed == {}


def test_checkpoint_get_metrics(outdir):
    cp = checkpoint_mod.DockingCheckpoint(outdir)
    cp.save_result("ligA", -7.5, {"status": "OK", "simscore": 0.8})
    metrics = cp.get_metrics()
    assert metrics["ligA"]["status"] == "OK"


def test_no_resume_clears_checkpoint(outdir):
    cp = checkpoint_mod.DockingCheckpoint(outdir)
    cp.save_result("ligA", -7.5, {})
    assert os.path.exists(cp.checkpoint_file)
    cp_file = cp.checkpoint_file
    os.remove(cp_file)
    assert not os.path.exists(cp_file)


def test_resume_skips_completed(outdir, charged_pdbqt):
    lig = charged_pdbqt
    cp = checkpoint_mod.DockingCheckpoint(outdir)
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
    good = utils.ligand_name(charged_pdbqt)

    def first_run(args):
        name = utils.ligand_name(args.ligand)
        if name == good:
            return name, -8.4, {"status": "OK"}
        return name, None, {"status": "FAILED"}

    monkeypatch.setattr(docking, "dock_ligand", first_run)
    docking.dock_all(charged_pdbqt, ligands, str(grid), outdir, num_processes=1,
                    resume=True, vina_params={"seed": 42, "extra_args": []})

    # Second run: the cached ligand is skipped and the other still fails.
    docked = []

    def second_run(args):
        docked.append(utils.ligand_name(args.ligand))
        return utils.ligand_name(args.ligand), None, {"status": "FAILED"}

    monkeypatch.setattr(docking, "dock_ligand", second_run)
    result, checkpoint, metrics = docking.dock_all(
        charged_pdbqt, ligands, str(grid), outdir, num_processes=1,
        resume=True, vina_params={"seed": 42, "extra_args": []})
    assert docked == [utils.ligand_name(zero_charge_pdbqt)]
    scores = [s for _, s in result if s is not None]
    assert scores == [-8.4]  # cached success retained; failure tolerated
    assert checkpoint.completed[good]["score"] == -8.4


def _ok_run(args):
    return utils.ligand_name(args.ligand), -7.0, {"status": "OK"}


def test_checkpoint_invalidated_when_parameters_change(
        tmp_path, charged_pdbqt, monkeypatch, outdir):
    from vspipeline import docking
    grid = tmp_path / "grid.txt"
    grid.write_text("center_x=1.0\ncenter_y=2.0\ncenter_z=3.0\n"
                    "size_x=20.0\nsize_y=20.0\nsize_z=20.0\n")
    monkeypatch.setattr(docking, "dock_ligand", _ok_run)
    docking.dock_all(charged_pdbqt, [charged_pdbqt], str(grid), outdir,
                    vina_params={"seed": 42, "exhaustiveness": 8})
    docked = []

    def tracking(args):
        docked.append(args.ligand)
        return _ok_run(args)

    monkeypatch.setattr(docking, "dock_ligand", tracking)
    # Same settings: cached. Different exhaustiveness: docked again.
    docking.dock_all(charged_pdbqt, [charged_pdbqt], str(grid), outdir,
                    vina_params={"seed": 42, "exhaustiveness": 8})
    assert docked == []
    docking.dock_all(charged_pdbqt, [charged_pdbqt], str(grid), outdir,
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
    docking.dock_all(charged_pdbqt, [charged_pdbqt, other], str(grid), outdir,
                    vina_params={"seed": 42})
    with open(other, "a", encoding="utf-8") as f:
        f.write("REMARK edited\n")
    docked = []

    def tracking(args):
        docked.append(args.ligand)
        return _ok_run(args)

    monkeypatch.setattr(docking, "dock_ligand", tracking)
    result, _, _ = docking.dock_all(charged_pdbqt, [other], str(grid), outdir,
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
        if args.ligand.endswith("c.pdbqt"):
            raise KeyboardInterrupt
        return _ok_run(args)

    monkeypatch.setattr(docking, "dock_ligand", interrupted)
    with pytest.raises(KeyboardInterrupt):
        docking.dock_all(charged_pdbqt, ligs, str(grid), outdir, vina_params={"seed": 42})
    assert set(checkpoint_mod.DockingCheckpoint(outdir).completed) == {"a", "b"}

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
    docking.dock_all(charged_pdbqt, [lig, keep], str(grid), outdir, vina_params={"seed": 42})
    with open(lig, "a", encoding="utf-8") as f:
        f.write("REMARK edited\n")
    monkeypatch.setattr(docking, "dock_ligand",
                        lambda a: ("x", None, {"status": "FAILED"}))
    result, _, _ = docking.dock_all(charged_pdbqt, [lig, keep], str(grid), outdir,
                                   vina_params={"seed": 42})
    # y stays cached; x is reported once, as failed, not also with its old score
    assert sorted(result, key=lambda r: r[0]) == [("x", None), ("y", -7.0)]


def test_out_of_window_score_is_cached(tmp_path, charged_pdbqt, monkeypatch, outdir):
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
        calls.append(utils.ligand_name(args.ligand))
        if args.ligand == weak:
            return "weak", None, {"status": "OUT_OF_WINDOW", "raw_score": -0.5}
        return _ok_run(args)

    monkeypatch.setattr(docking, "dock_ligand", dock)
    for _ in range(2):
        result, _, metrics = docking.dock_all(charged_pdbqt, [weak, good], str(grid),
                                             outdir, vina_params={"seed": 42})
        assert sorted(result, key=lambda r: r[0]) == [("good", -7.0), ("weak", None)]
        assert metrics["weak"]["status"] == "OUT_OF_WINDOW"
        assert metrics["weak"]["raw_score"] == -0.5
    assert sorted(calls) == ["good", "weak"]  # the second run docked nothing


def test_all_out_of_window_is_not_fatal(tmp_path, charged_pdbqt, monkeypatch, outdir):
    from vspipeline import docking
    grid = tmp_path / "grid.txt"
    grid.write_text("center_x=1.0\ncenter_y=2.0\ncenter_z=3.0\n"
                    "size_x=20.0\nsize_y=20.0\nsize_z=20.0\n")
    monkeypatch.setattr(docking, "dock_ligand", lambda a: (
        utils.ligand_name(a.ligand), None,
        {"status": "OUT_OF_WINDOW", "raw_score": -0.4}))
    result, _, metrics = docking.dock_all(charged_pdbqt, [charged_pdbqt], str(grid),
                                         outdir, vina_params={"seed": 42})
    assert result == [("charged", None)]
    assert metrics["charged"]["raw_score"] == -0.4


def _fp_files(tmp_path):
    grid = tmp_path / "grid.txt"
    grid.write_text("center_x=0\ncenter_y=0\ncenter_z=0\nsize_x=20\nsize_y=20\nsize_z=20\n")
    heavy = "ATOM      1  N   LYS A   1       1.000   2.000   3.000  0.00  0.00    -0.300 N \n"
    h_at = lambda x: f"ATOM      2  H   LYS A   1    {x:8.3f}   2.000   3.000  0.00  0.00    +0.000 HD\n"
    a, b, c = (tmp_path / n for n in ("a.pdbqt", "b.pdbqt", "c.pdbqt"))
    a.write_text(heavy + h_at(1.5))
    b.write_text(heavy + h_at(0.4))                      # hydrogen moved
    c.write_text(heavy.replace("1.000", "1.200") + h_at(1.5))  # heavy atom moved
    return grid, a, b, c


def test_fingerprint_includes_hydrogen_coordinates(tmp_path):
    """Vina perceives donor hydrogens by bond distance, so moving a polar H can
    change the result: the fingerprint must cover H coordinates (receptor and
    flex file)."""
    from vspipeline import docking
    grid, a, b, c = _fp_files(tmp_path)
    fp = lambda r, flex=None: docking._run_fingerprint(str(r), str(grid), flex, {"seed": 42})
    assert fp(a) != fp(b)
    assert fp(a) != fp(c)
    assert fp(a) == fp(a)
    assert fp(a, str(a)) != fp(a, str(b))  # flex file hashed in full too


def test_fingerprint_covers_engine_binary(tmp_path):
    from vspipeline import docking
    grid, a, _, _ = _fp_files(tmp_path)
    exe1 = tmp_path / "vina1"
    exe2 = tmp_path / "vina2"
    exe1.write_text("#!/bin/sh\necho 'AutoDock Vina v1.2.5'\n")
    exe2.write_text("#!/bin/sh\necho 'AutoDock Vina v1.2.6'\n")
    os.chmod(exe1, 0o755)
    os.chmod(exe2, 0o755)
    ident = docking._engine_identity(str(exe1))
    assert ident["version"] == "AutoDock Vina v1.2.5"
    assert ident["sha256"] == docking._file_sha256(str(exe1))
    fp = lambda exe: docking._run_fingerprint(
        str(a), str(grid), None, {"seed": 42},
        docking.EngineConfig(vina_primary=str(exe), vina_fallback=str(exe)))
    assert fp(exe1) != fp(exe2)       # different binary (same path layout)
    first = fp(exe1)
    exe1.write_text("#!/bin/sh\necho 'AutoDock Vina v1.2.5' # rebuilt\n")
    assert fp(exe1) != first          # same path, binary replaced


def test_checkpoint_flush_is_atomic_and_truncates_journal(outdir, monkeypatch):
    cp = checkpoint_mod.DockingCheckpoint(outdir)
    cp.save_result("a", -7.0, {"status": "OK"})
    cp.save_result("b", -6.0, {"status": "OK"})
    assert os.path.getsize(cp.journal_file) > 0
    cp.flush()
    assert os.path.getsize(cp.journal_file) == 0
    with open(cp.checkpoint_file, encoding="utf-8") as f:
        assert set(json.load(f)) == {"a", "b"}
    assert not [p for p in os.listdir(outdir) if ".tmp." in p]
    # A crash during the write leaves the previous checkpoint intact.
    from vspipeline import checkpoint as ckpt
    cp.save_result("c", -5.0, {"status": "OK"})

    def boom(*a, **k):
        raise OSError("disk full")

    monkeypatch.setattr(ckpt.json, "dump", boom)
    cp.flush()
    monkeypatch.undo()
    with open(cp.checkpoint_file, encoding="utf-8") as f:
        assert set(json.load(f)) == {"a", "b"}
    assert os.path.getsize(cp.journal_file) > 0  # not truncated after a failed write
    assert set(checkpoint_mod.DockingCheckpoint(outdir).completed) == {"a", "b", "c"}


def test_checkpoint_lock_rejects_second_process(outdir):
    import subprocess
    import sys
    from vspipeline import checkpoint as ckpt
    code = (
        "import sys; sys.path.insert(0, sys.argv[2]);"
        "from vspipeline import checkpoint as c;"
        "c.acquire_lock(sys.argv[1]); print('locked', flush=True); sys.stdin.read()")
    pkg = os.path.join(os.path.dirname(__file__), "..", "autodocker")
    holder = subprocess.Popen([sys.executable, "-c", code, outdir, pkg],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    try:
        assert holder.stdout.readline().strip() == "locked"
        with pytest.raises(ckpt.OutputDirLockedError, match="in use by another"):
            checkpoint_mod.DockingCheckpoint(outdir)
    finally:
        holder.communicate("")
    checkpoint_mod.DockingCheckpoint(outdir)  # free again once the other run exits
    ckpt.release_lock(outdir)


def test_checkpoint_timestamp_not_in_metrics(outdir):
    cp = checkpoint_mod.DockingCheckpoint(outdir)
    cp.save_result("a", -7.0, {"status": "OK"})
    assert "timestamp" not in cp.get_metrics()["a"]
