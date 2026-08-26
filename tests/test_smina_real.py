import os

import pytest

import runner
import vspipeline.consensus as consensus


SMINA_AVAILABLE = runner.SMINA is not None


@pytest.fixture
def tiny_system(tmp_path):
    """Create a tiny receptor + ligand PDBQT with obabel (must be on PATH)."""
    if runner.OBABEL is None:
        pytest.skip("obabel not available")
    rec = str(tmp_path / "receptor.pdbqt")
    lig = str(tmp_path / "ligand.pdbqt")
    out = str(tmp_path / "out.pdbqt")
    runner.run([runner.OBABEL, "-:c1ccccc1", "-opdbqt", "-O", rec])
    runner.run([runner.OBABEL, "-:CO", "-opdbqt", "-O", lig])
    assert os.path.exists(rec) and os.path.exists(lig)
    return rec, lig, out


@pytest.mark.skipif(not SMINA_AVAILABLE, reason="smina binary not installed")
def test_run_smina_docking_real(tiny_system):
    rec, lig, out = tiny_system
    aff = consensus.run_smina_docking(
        rec, lig, out,
        cx=0.0, cy=0.0, cz=0.0,
        sx=30.0, sy=30.0, sz=30.0,
        exhaustiveness=1,
    )
    assert isinstance(aff, float)
    assert os.path.exists(out)
    assert any(line.startswith("MODEL") for line in open(out))


@pytest.mark.skipif(not SMINA_AVAILABLE, reason="smina binary not installed")
def test_parse_smina_affinity_real_stdout():
    # Captured from a real smina run; verifies parser on genuine output.
    stdout = (
        "mode |   affinity | dist from best mode\n"
        "     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n"
        "-----+------------+----------+----------\n"
        "1       -0.4       0.000      0.000    \n"
    )
    assert consensus._parse_smina_affinity(stdout) == -0.4
