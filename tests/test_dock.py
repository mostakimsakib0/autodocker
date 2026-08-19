import os
import shutil

import pytest

import runner


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


def test_all_fail_raises_with_pointer(tmp_path, empty_pdbqt):
    """Regression for C8: when every ligand fails, raise with a pointer to logs."""
    outdir = str(tmp_path / "out")
    os.makedirs(outdir, exist_ok=True)
    receptor = str(tmp_path / "rec.pdbqt")
    with open(receptor, "w") as f:
        f.write("END\n")
    with pytest.raises(RuntimeError, match="vina.log"):
        runner.dock_all(
            receptor=receptor, ligands=[empty_pdbqt], grid_file=receptor,
            outdir=outdir, num_processes=1, resume=False,
            vina_params={"seed": 42})


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


def test_dock_ligand_rejects_zero_charge_receptor(tmp_path, charged_pdbqt, zero_charge_pdbqt):
    name, score, metrics = runner.dock_ligand((
        zero_charge_pdbqt, charged_pdbqt,
        charged_pdbqt, str(tmp_path), {}, None))
    assert score is None
    assert metrics["status"] == "FAILED"
    assert "charges" in metrics["error"].lower()


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
    assert "--flex" not in cmd


def test_build_vina_command_threads(monkeypatch):
    monkeypatch.setattr(runner, "_vina_supports_threads", lambda p: True)
    cmd = runner._build_vina_command(
        "/usr/bin/vina", "rec.pdbqt", "lig.pdbqt", "out.pdbqt",
        {"center_x": 1, "size_x": 24},
        {"exhaustiveness": 8, "binding_modes": 9, "energy_range": 3.0,
         "seed": 42, "threads": 2})
    assert "--threads" in cmd
    assert cmd[cmd.index("--threads") + 1] == "2"


def test_build_vina_command_no_threads_by_default():
    cmd = runner._build_vina_command(
        "/usr/bin/vina", "rec.pdbqt", "lig.pdbqt", "out.pdbqt",
        {"center_x": 1, "size_x": 24},
        {"exhaustiveness": 8, "binding_modes": 9, "energy_range": 3.0,
         "seed": 42})
    assert "--threads" not in cmd


def test_parallel_equals_serial(tmp_path, protein_pdb):
    """Regression for C4: same inputs + same seed must give same result
    regardless of process count. Uses a full docking when vina exists."""
    vina = shutil.which("vina")
    obabel = shutil.which("obabel")
    if vina is None or obabel is None:
        pytest.skip("vina/obabel not available")

    receptor = str(tmp_path / "rec.pdbqt")
    runner.run([runner.OBABEL, "-ipdb", protein_pdb, "-opdbqt",
                "-O", receptor, "-xr", "-c",
                "--partialcharge", "gasteiger"])

    def run_dock(np_):
        outdir = str(tmp_path / f"out_{np_}")
        os.makedirs(outdir, exist_ok=True)
        ligdir = str(tmp_path / f"ligs_{np_}")
        os.makedirs(ligdir, exist_ok=True)
        lig = os.path.join(ligdir, "lig.pdbqt")
        runner.run([runner.OBABEL, "-:CCO", "-opdbqt", "-O", lig,
                    "--partialcharge", "gasteiger"])
        grid = str(tmp_path / f"grid_{np_}.conf")
        with open(grid, "w") as f:
            f.write("center_x = 4.5\ncenter_y = 6.5\ncenter_z = 5.0\n"
                    "size_x = 24\nsize_y = 24\nsize_z = 24\n")
        vina_params = {
            "exhaustiveness": 2, "binding_modes": 1, "energy_range": 1.0,
            "seed": 42, "min_valid_affinity": -1.0}
        results, _, metrics = runner.dock_all(
            receptor, [lig], grid, outdir, num_processes=np_,
            resume=False, vina_params=vina_params)
        assert len(results) == 1
        assert results[0][1] is not None
        return results[0][1]

    score_serial = run_dock(1)
    score_parallel = run_dock(2)
    assert score_serial == score_parallel