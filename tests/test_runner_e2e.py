"""End-to-end run of runner.main() with the real tools (tiny and fast)."""
import csv
import json
import logging
import os
import shutil
import subprocess

import pytest

import runner
from vspipeline import report as report_mod
from vspipeline import manifest as manifest_mod
from vspipeline import protein, tools, utils
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
    # Open Babel version and OB_RANDOM_SEED policy (docs/openbabel.md).
    version = utils.obabel_version()
    assert manifest["openbabel_version"] == ".".join(map(str, version))
    ob_seed = manifest["openbabel_random_seed"]
    assert ob_seed["variable"] == "OB_RANDOM_SEED" and ob_seed["value"] == 3
    assert ob_seed["hydrogen_probe_values"] == [4, 5]
    assert ob_seed["honoured"] is (version >= utils.OBABEL_SEEDABLE_VERSION)
    assert (out / "ranking.csv").exists()
    meta = json.loads((out / "ligand_metadata.json").read_text())
    assert meta["ethanol"]["minimized"] is True
    # Input and output hashes (scripts/verify_run.py checks them).
    assert manifest["inputs"]["receptor"] == {
        "path": os.path.abspath(protein_pdb), "sha256": manifest_mod._sha256(protein_pdb)}
    lig = manifest["inputs"]["ligands"]["ethanol"]
    assert lig["path"] == str(ligs / "ethanol.sdf") and len(lig["sha256"]) == 64
    assert "box_ligand" not in manifest["inputs"]
    hashed = manifest["outputs"]
    assert list(hashed) == sorted(hashed)
    for rel in ("ranking.csv", "receptor.pdbqt", "protein_clean.pdb", "grid.conf",
                "ligand_metadata.json", "docked/ethanol_out.pdbqt", "docked/ethanol_vina.log",
                "ligands/ethanol.pdbqt"):
        assert hashed[rel] == manifest_mod._sha256(str(out / rel)), rel
    assert not any(rel.startswith(".") or "/." in rel or rel == "run_manifest.json"
                   for rel in hashed)
    log = (out / "docked" / "ethanol_vina.log").read_text()
    assert f"# Engine: {manifest['tools']['vina']['path']} sha256=" in log


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


# Run records excluded from the byte comparison of two runs: the manifest
# (times, paths), logs, the fpocket output and the docking bookkeeping.
_RUN_RECORDS = ("run_manifest.json", "pocket_summary.csv")


def _output_bytes(out):
    """{relative path: bytes} of every output of a run except _RUN_RECORDS,
    *.log, protein_clean_out/ and .docking* files."""
    files = {}
    for root, dirs, names in os.walk(out):
        dirs[:] = [d for d in dirs if d != "protein_clean_out" and not d.startswith(".docking")]
        for name in names:
            if name in _RUN_RECORDS or name.endswith(".log") or name.startswith(".docking"):
                continue
            path = os.path.join(root, name)
            with open(path, "rb") as f:
                files[os.path.relpath(path, out)] = f.read()
    return files


@pytest.mark.slow
def test_runs_are_byte_identical_including_the_html_report(tmp_path, protein_pdb,
                                                          monkeypatch):
    """Two runs with -p 1 and one with -p 2 write the same bytes in every
    output, results_report.html included (it has no timestamp)."""
    monkeypatch.setattr(report_mod, "_fetch_ngljs", lambda: b"/* fake ngl */")
    ligs = _ethanol_dir(tmp_path)
    subprocess.run([shutil.which("obabel"), "-:CC(=O)NCCO acetamide", "--gen3d", "-osdf",
                    "-O", str(ligs / "acetamide.sdf")], check=True, capture_output=True,
                   env={**os.environ, "OMP_NUM_THREADS": "1", "OB_RANDOM_SEED": "1"})
    runs = []
    for name, procs in (("a", "1"), ("b", "1"), ("p2", "2")):
        out = str(tmp_path / name)
        argv = ["-r", protein_pdb, "-l", str(ligs), "-o", out, "--no-fpocket",
                "--exhaustiveness", "1", "--binding-modes", "2", "--seed", "3",
                "--minimize-steps", "20", "-p", procs]
        assert _main(argv) == runner.EXIT_OK
        runs.append(_output_bytes(out))
    assert "results_report.html" in runs[0] and "receptor.pdbqt" in runs[0]
    for other in runs[1:]:
        assert sorted(other) == sorted(runs[0])
        assert [k for k in runs[0] if runs[0][k] != other[k]] == []


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

    tools.configure(fpocket=None)
    monkeypatch.setattr(protein.ProteinPreparation, "detect_pocket", no_fpocket)
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
    tools.configure(smina=_failing_binary(tmp_path, "smina"))
    argv = _base_argv(tmp_path, protein_pdb, modes=False) + ["--no-fpocket", "--smina-only"]
    assert _main(argv) == runner.EXIT_TOO_MANY_FAILURES
    out = tmp_path / "out"
    assert (out / "ranking.csv").exists()
    assert _manifest(out)["docking_failures"] == {"failed": 1, "docked": 1}


# ---------------------------------------------------------------------------
# --box-from-het / --box-ligand and --replicates through runner.main
# ---------------------------------------------------------------------------

_LIGAND_BOX_DATA = os.path.join(os.path.dirname(__file__), "data", "ligand_box")
_SITE_PDB = os.path.join(_LIGAND_BOX_DATA, "1hsg_site.pdb")


def _site_argv(tmp_path, out):
    return ["-r", _SITE_PDB, "-l", str(_ethanol_dir(tmp_path)), "-o", str(out),
            "--chain", "all", "--exhaustiveness", "1", "--binding-modes", "1",
            "--seed", "3", "--minimize-steps", "20", "--no-html-report"]


def _grid_conf(out):
    values = {}
    for line in (out / "grid.conf").read_text().splitlines():
        key, _, value = line.partition("=")
        values[key.strip()] = float(value)
    return values


def test_box_from_het_mk1_end_to_end(tmp_path, caplog):
    """--box-from-het MK1 on the 1HSG site: MK1 is removed from the receptor
    (also with --keep-hetero), the box is the reference box in the docking
    frame, and the manifest records it."""
    out = tmp_path / "out"
    with caplog.at_level(logging.INFO):
        code = _main(_site_argv(tmp_path, out) + ["--box-from-het", "MK1", "--keep-hetero"])
    assert code == runner.EXIT_OK
    assert "MK1" not in (out / "receptor.pdbqt").read_text()
    assert "MK1" not in (out / "receptor_no_reference.pdb").read_text()
    assert (out / "reference_ligand.pdb").read_text().count("HETATM") == 45
    box = _manifest(out)["box"]
    assert box["source"] == "het_group" and box["reference"] == "MK1:B:902"
    assert box["reference_heavy_atoms"] == 45 and box["reference_atoms_in_box"] == 45
    assert box["size"] == box["size_computed"] and box["size_source"] == "reference ligand"
    grid = _grid_conf(out)
    assert [grid["center_x"], grid["center_y"], grid["center_z"]] == pytest.approx(
        box["center_docking_frame"], abs=1e-3)
    assert [grid["size_x"], grid["size_y"], grid["size_z"]] == pytest.approx(box["size"],
                                                                            abs=1e-3)
    offset = [float(v) for v in (out / "frame_offset.txt").read_text().split()
              if v.lstrip("-").replace(".", "", 1).isdigit()][:3]
    assert [c - o for c, o in zip(box["center_input_frame"], offset)] == pytest.approx(
        box["center_docking_frame"], abs=2e-3)
    assert "--box-from-het MK1:B:902: removed its 45 records" in caplog.text
    if tools.settings().fpocket:
        assert box["fpocket"] is not None and "top_ranked_pocket" in box["fpocket"]
        assert (out / "pocket_summary.csv").exists()


def test_box_from_het_without_fpocket_installed(tmp_path, caplog):
    """--box-from-het without fpocket on PATH: the box comes from the
    reference, so the run goes ahead and only the pocket report is skipped."""
    tools.configure(fpocket=None)
    out = tmp_path / "out"
    with caplog.at_level(logging.INFO):
        code = _main(_site_argv(tmp_path, out) + ["--box-from-het", "MK1"])
    assert code == runner.EXIT_OK
    assert "fpocket not found: the pocket report" in caplog.text
    m = _manifest(out)
    assert m["box"]["source"] == "het_group" and m["box"]["fpocket"] is None
    assert "fpocket" not in m["tools"]
    assert not (out / "pocket_summary.csv").exists()


def test_box_ligand_with_grid_size_and_no_fpocket(tmp_path, caplog):
    """--box-ligand mk1.sdf --grid-size: the center from the reference, the
    size given; --no-fpocket skips the fpocket report. With --keep-hetero the
    receptor keeps its own MK1, with a warning naming --box-from-het."""
    out = tmp_path / "out"
    with caplog.at_level(logging.INFO):
        code = _main(_site_argv(tmp_path, out) + [
            "--box-ligand", os.path.join(_LIGAND_BOX_DATA, "mk1.sdf"), "--no-fpocket",
            "--grid-size", "22.5,22.5,22.5", "--keep-hetero"])
    assert code == runner.EXIT_OK
    box = _manifest(out)["box"]
    assert box["source"] == "ligand_file" and box["reference_sha256"]
    assert box["size"] == [22.5, 22.5, 22.5] and box["size_source"] == "--grid-size"
    assert box["fpocket"] is None and not (out / "pocket_summary.csv").exists()
    assert "MK1" in (out / "receptor.pdbqt").read_text()
    assert "use --box-from-het MK1:B:902 to remove it" in caplog.text


def test_replicates_end_to_end_and_resume(tmp_path, caplog):
    """--replicates 2: two seeded dockings per ligand, the best one reported
    with the replicate columns; a rerun docks nothing and writes the same
    bytes; another N re-docks every (ligand, replicate)."""
    out = tmp_path / "out"
    argv = _site_argv(tmp_path, out) + ["--box-from-het", "MK1", "--no-fpocket",
                                        "--replicates", "2"]
    with caplog.at_level(logging.INFO):
        assert _main(argv) == runner.EXIT_OK
    assert "Docking 2 new (ligand, replicate) pairs (0 cached)" in caplog.text
    for k in (0, 1):
        assert (out / "docked" / f"ethanol_rep{k}_out.pdbqt").exists()
        log = (out / "docked" / f"ethanol_rep{k}_vina.log").read_text()
        assert f"--seed {3 + k}" in log
    with open(out / "ranking.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1 and rows[0]["n_replicates"] == "2"
    assert float(rows[0]["score_min"]) == float(rows[0]["Binding_Affinity"])
    with open(out / "replicate_scores.csv", newline="", encoding="utf-8") as f:
        reps = list(csv.DictReader(f))
    assert [(r["Replicate"], r["Seed"]) for r in reps] == [("0", "3"), ("1", "4")]
    assert sum(r["Reported"] == "yes" for r in reps) == 1
    summary = json.loads((out / "replicate_summary.json").read_text())
    assert summary["seeds"] == [3, 4] and _manifest(out)["replicates"] == summary
    first = _output_bytes(out)

    caplog.clear()
    with caplog.at_level(logging.INFO):
        assert _main(argv) == runner.EXIT_OK
    assert "Docking 0 new (ligand, replicate) pairs (2 cached)" in caplog.text
    assert _output_bytes(out) == first

    caplog.clear()
    with caplog.at_level(logging.INFO):
        assert _main(argv[:-1] + ["3"]) == runner.EXIT_OK
    assert "Docking 3 new (ligand, replicate) pairs (0 cached)" in caplog.text
