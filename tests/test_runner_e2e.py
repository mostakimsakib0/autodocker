"""End-to-end run of runner.main() with the real tools (tiny and fast)."""
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
