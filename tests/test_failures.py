"""Failure-mode and edge-case tests for the autodocker pipeline."""

import argparse
import os
import sys

import pytest

import runner

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "scripts"))


# ---------------------------------------------------------------------------
# Versioning
# ---------------------------------------------------------------------------

def test_version_attribute_and_flag(capsys, monkeypatch):
    assert isinstance(runner.__version__, str)
    assert runner.__version__.count(".") == 2
    monkeypatch.setattr(sys, "argv", ["runner", "--version"])
    with pytest.raises(SystemExit) as exc:
        runner.main()
    assert exc.value.code == 0
    assert "AutoDocker" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# Input rejection (fail-loud behaviour)
# ---------------------------------------------------------------------------

def test_tiny_receptor_rejected(tiny_protein_pdb, outdir):
    prep = runner.ProteinPreparation(tiny_protein_pdb, outdir)
    with pytest.raises(ValueError, match="need >= 100"):
        prep.validate()


def test_empty_pdbqt_has_no_atoms(empty_pdbqt):
    assert runner._pdbqt_has_atoms(empty_pdbqt) is False


def test_zero_charge_pdbqt_detected(charged_pdbqt, empty_pdbqt):
    assert runner._ensure_pdbqt_has_charges(charged_pdbqt) is True


def test_charge_check_empty_file(empty_pdbqt):
    assert runner._ensure_pdbqt_has_charges(empty_pdbqt) is False


def _grid_config(tmp_path):
    grid = tmp_path / "grid.txt"
    grid.write_text("center_x=1.0\ncenter_y=2.0\ncenter_z=3.0\n"
                    "size_x=20.0\nsize_y=20.0\nsize_z=20.0\n")
    return grid


def test_missing_receptor_fails_dock(tmp_path, charged_pdbqt, monkeypatch, outdir):
    """dock_ligand must return a FAILED record when inputs are missing."""
    grid = _grid_config(tmp_path)
    result = runner.dock_ligand((
        str(tmp_path / "does_not_exist_receptor.pdbqt"),
        charged_pdbqt,
        str(grid), outdir,
        {"seed": 42, "extra_args": []}, None))
    name, score, meta = result
    assert score is None
    assert meta["status"] == "FAILED"
    assert "Receptor file missing" in meta["error"]


def test_all_tools_failed_is_fatal(tmp_path, charged_pdbqt, monkeypatch, outdir):
    """When every docking binary fails, dock_ligand returns FAILED (no silence)."""
    def _boom(*a, **k):
        raise RuntimeError("simulated docking binary crash")

    monkeypatch.setattr(runner, "run", _boom)
    grid = _grid_config(tmp_path)
    result = runner.dock_ligand((
        charged_pdbqt, charged_pdbqt,
        str(grid), outdir,
        {"seed": 42, "extra_args": []}, None))
    name, score, meta = result
    assert score is None
    assert meta["status"] == "FAILED"


def test_duplicate_ligand_names_rejected(tmp_path, outdir, charged_pdbqt):
    lig_dir = tmp_path / "ligs"
    lig_dir.mkdir()
    (lig_dir / "dup.pdbqt").write_text(open(charged_pdbqt).read())
    (lig_dir / "dup.pdbqt").touch()
    ligs = [str(lig_dir / "dup.pdbqt")] * 3
    with pytest.raises(ValueError, match="Duplicate ligand names"):
        runner.dock_all(charged_pdbqt, ligs, tmp_path / "g.txt", outdir,
                        resume=False, vina_params={})


# ---------------------------------------------------------------------------
# Grid parsing edge cases
# ---------------------------------------------------------------------------

def test_grid_missing_axis(tmp_path):
    cfg = tmp_path / "bad_grid.txt"
    cfg.write_text("cx 1.0\n")
    with pytest.raises(Exception):
        runner._parse_grid_config(str(cfg))


# ---------------------------------------------------------------------------
# Enrichment / performance helper edge cases (benchmark scripts)
# ---------------------------------------------------------------------------

def _load_benchmark(name):
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, os.path.join(
        os.path.dirname(__file__), "..", "scripts", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_enrichment_auc_degenerates_to_none():
    be = _load_benchmark("benchmark_enrichment")
    # no actives -> None
    scores = {f"D_{i}": -5.0 + i for i in range(10)}
    labels = {n: n.startswith("A_") for n in scores}
    assert be.auc_and_ef(scores, labels) == (None, None, None)
    # perfect separation still works
    scores = {f"A_{i}": -8.0 - i for i in range(5)}
    scores.update({f"D_{i}": -2.0 for i in range(5)})
    auc, ef1, _ = be.auc_and_ef(scores, labels={n: n.startswith("A_")
                                                for n in scores})
    assert auc == 1.0


def test_enrichment_mannwhitney_ties():
    be = _load_benchmark("benchmark_enrichment")
    u1, p1, _ = be.mannwhitney([-1.0, -1.0], [-1.0, -1.0])
    assert u1 is not None
    assert p1 is not None


def test_performance_build_library_requires_pdbqt(tmp_path):
    bp = _load_benchmark("benchmark_performance")
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(RuntimeError, match="no \\*.pdbqt"):
        bp.build_library(str(empty), 5, 42, str(tmp_path))


# ---------------------------------------------------------------------------
# CLI parse sanity (a wrong option must fail, not silently pass)
# ---------------------------------------------------------------------------

def test_unknown_argument_fails(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["runner", "--definitely-not-an-option"])
    with pytest.raises(SystemExit) as exc:
        runner.main()
    assert exc.value.code != 0