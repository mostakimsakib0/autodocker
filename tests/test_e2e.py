"""Golden determinism tests: same input + seed -> byte-identical results.

Real docking on committed fixtures (tests/data/e2e, see make_e2e_fixtures.py):
a binding-site fragment of PDB 2FXS prepared as PDBQT, and three drug-like
ligands with rotatable bonds and 3D coordinates from a seeded RDKit
embedding. The receptor and ligands are committed as prepared PDBQT, so these
tests check docking, ranking and output writing, not receptor protonation.

Skipped when vina is not on PATH.
"""
import glob
import json
import logging
import os
import shutil
import subprocess

import pytest

DATA = os.path.join(os.path.dirname(__file__), "data", "e2e")
LIGANDS = sorted(glob.glob(os.path.join(DATA, "lig_*.pdbqt")))
VINA = shutil.which("vina")

pytestmark = [
    pytest.mark.skipif(VINA is None, reason="vina not available on PATH"),
    pytest.mark.slow,
]


def _box():
    with open(os.path.join(DATA, "box.json"), encoding="utf-8") as f:
        box = json.load(f)
    return box["center"], box["size"]


def _write_grid(path):
    (cx, cy, cz), s = _box()
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"center_x = {cx}\ncenter_y = {cy}\ncenter_z = {cz}\n"
                f"size_x = {s}\nsize_y = {s}\nsize_z = {s}\n")
    return path


def _run_pipeline(outdir, processes):
    # imported lazily so the skip marker works
    from vspipeline import docking, results as results_mod
    logging.disable(logging.CRITICAL)
    try:
        os.makedirs(outdir)
        receptor = os.path.join(DATA, "receptor.pdbqt")
        grid = _write_grid(os.path.join(outdir, "grid.conf"))
        vina_params = {"exhaustiveness": 1, "binding_modes": 3, "energy_range": 3.0,
                       "seed": 42, "min_valid_affinity": -1.0}
        results, _checkpoint, metrics = docking.dock_all(
            receptor, LIGANDS, grid, outdir, num_processes=processes,
            resume=False, vina_params=vina_params)
        analyzer = results_mod.ResultsAnalyzer(outdir)
        analyzer.save_ranking(results, mode="extended", metrics_dict=metrics)
        analyzer.save_metrics_report(results, metrics)
    finally:
        logging.disable(logging.NOTSET)
    return outdir


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("golden")
    return {
        "p1_a": _run_pipeline(str(root / "p1_a"), 1),
        "p1_b": _run_pipeline(str(root / "p1_b"), 1),
        "p2": _run_pipeline(str(root / "p2"), 2),
    }


def _ranking_bytes(outdir):
    """ranking.csv with the run directory normalised away (versions that
    record absolute Docked_File paths)."""
    with open(os.path.join(outdir, "ranking.csv"), "rb") as f:
        return f.read().replace(os.path.abspath(outdir).encode() + b"/", b"<RUN>/")


def _poses(outdir):
    out = {}
    for path in sorted(glob.glob(os.path.join(outdir, "docked", "*_out.pdbqt"))):
        with open(path, "rb") as f:
            out[os.path.basename(path)] = f.read()
    return out


def test_fixtures_are_realistic():
    assert len(LIGANDS) >= 3
    for lig in LIGANDS:
        lines = open(lig, encoding="utf-8").read().splitlines()
        coords = {ln[30:54] for ln in lines if ln.startswith(("ATOM", "HETATM"))}
        assert len(coords) > 5  # real 3D coordinates, not all at the origin
        torsdof = next(int(ln.split()[1]) for ln in lines if ln.startswith("TORSDOF"))
        assert torsdof >= 2
    n_rec = sum(1 for ln in open(os.path.join(DATA, "receptor.pdbqt"), encoding="utf-8") if ln.startswith("ATOM"))
    assert n_rec > 300


def test_all_ligands_docked(runs):
    import csv
    with open(os.path.join(runs["p1_a"], "ranking.csv"), encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == len(LIGANDS)
    for row in rows:
        assert float(row["Binding_Affinity"]) < 0
    assert len(_poses(runs["p1_a"])) == len(LIGANDS)
    for name in ("ranking.csv", "Top_hits.txt", "metrics.txt"):
        assert os.path.exists(os.path.join(runs["p1_a"], name))


def test_repeat_run_is_byte_identical(runs):
    assert _ranking_bytes(runs["p1_a"]) == _ranking_bytes(runs["p1_b"])
    assert _poses(runs["p1_a"]) == _poses(runs["p1_b"])


def test_process_count_does_not_change_results(runs):
    # -p 2 also gives each Vina process half the cores (--cpu), so this
    # covers the thread count as well.
    assert _ranking_bytes(runs["p1_a"]) == _ranking_bytes(runs["p2"])
    assert _poses(runs["p1_a"]) == _poses(runs["p2"])


def test_vina_cpu_count_does_not_change_output(tmp_path):
    (cx, cy, cz), s = _box()
    ligand = next(lig for lig in LIGANDS if "paracetamol" in lig)
    outputs = []
    for cpu in (1, 2, 4):
        out = str(tmp_path / f"cpu{cpu}.pdbqt")
        subprocess.run(
            [VINA, "--receptor", os.path.join(DATA, "receptor.pdbqt"), "--ligand", ligand,
             "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
             "--size_x", str(s), "--size_y", str(s), "--size_z", str(s),
             "--exhaustiveness", "2", "--num_modes", "3", "--seed", "42",
             "--cpu", str(cpu), "--out", out],
            check=True, capture_output=True, timeout=600)
        with open(out, "rb") as f:
            outputs.append(f.read())
    assert outputs[0] == outputs[1] == outputs[2]
