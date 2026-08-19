"""Golden determinism test (ROADMAP H7): same input + seed -> identical ranking.csv.

This test performs a real end-to-end smoke dock through the pipeline entrypoints.
It is skipped when vina or obabel are not on PATH.
"""
import hashlib
import os
import shutil
import subprocess
import sys

import pytest

runner = None
pytestmark = pytest.mark.skipif(
    not (shutil.which("vina") and shutil.which("obabel")),
    reason="vina/obabel not available on PATH")


def _hash_file(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def _hash_ranking(path):
    """Hash only the reproducibility-critical columns (paths differ by outdir)."""
    import csv
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f):
            rows.append((
                row["Ligand"], row["Binding_Affinity"], row["SimScore"],
                row["Status"]))
    return hashlib.sha256(repr(sorted(rows)).encode()).hexdigest()


@pytest.fixture
def receptor(tmp_path, protein_pdb):
    out = str(tmp_path / "rec.pdbqt")
    subprocess.run([shutil.which("obabel"), "-ipdb", protein_pdb, "-opdbqt",
                    "-O", out, "-xr", "-c", "--partialcharge", "gasteiger"],
                   check=True, capture_output=True)
    return out


@pytest.fixture
def ligands(tmp_path):
    d = tmp_path / "ligs"
    d.mkdir()
    for name, smiles in (("eth", "CCO"), ("prop", "CCC")):
        out = str(d / f"{name}.pdbqt")
        subprocess.run([shutil.which("obabel"), f"-:{smiles}",
                        "-opdbqt", "-O", out, "--partialcharge", "gasteiger"],
                       check=True, capture_output=True)
    return str(d)


def _run_pipeline(receptor, ligands, outdir):
    import runner  # imported lazily so the skip marker works
    import logging
    logging.disable(logging.CRITICAL)
    lib = runner.LibraryManager(outdir, ligands)
    prepared = lib._prepare_local_sdf(apply_admet=False)
    assert len(prepared) == 2

    prep = runner.ProteinPreparation(receptor, outdir)
    cx, cy, cz, sx, sy, sz = prep._protein_centroid_grid(padding=6.0)
    prep.write_grid(cx, cy, cz, sx, sy, sz)

    vina_params = {
        "exhaustiveness": 2, "binding_modes": 1, "energy_range": 1.0,
        "seed": 42, "min_valid_affinity": -1.0,
    }
    results, checkpoint, metrics = runner.dock_all(
        receptor, prepared, prep.grid_conf, outdir,
        num_processes=1, resume=False, vina_params=vina_params)

    analyzer = runner.ResultsAnalyzer(outdir)
    analyzer.save_ranking(results, mode="extended", metrics_dict=metrics)
    analyzer.save_metrics_report(results, metrics)
    return os.path.join(outdir, "ranking.csv")


def test_golden_determinism(tmp_path, receptor, ligands):
    csv1 = _run_pipeline(receptor, ligands, str(tmp_path / "run1"))
    csv2 = _run_pipeline(receptor, ligands, str(tmp_path / "run2"))
    assert _hash_ranking(csv1) == _hash_ranking(csv2)


def test_smoke_dock_produces_all_outputs(tmp_path, receptor, ligands):
    outdir = str(tmp_path / "run")
    csv_path = _run_pipeline(receptor, ligands, outdir)
    assert os.path.exists(csv_path)
    assert os.path.exists(os.path.join(outdir, "ranking.csv"))
    assert os.path.exists(os.path.join(outdir, "Top_hits.txt"))
    assert os.path.exists(os.path.join(outdir, "metrics.txt"))
    assert os.path.exists(os.path.join(outdir, "docked"))