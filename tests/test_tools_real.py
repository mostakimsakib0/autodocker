import shutil
import subprocess
import importlib.util

import pytest

REDUCE = shutil.which("reduce")
MKP = shutil.which("mk_prepare_ligand.py")
QVINA = shutil.which("qvina2") or shutil.which("qvina02") or shutil.which("qvina")
AUTODOCK = shutil.which("autodock4")
RDKIT = importlib.util.find_spec("rdkit") is not None
OBABEL = shutil.which("obabel")


def _run_ok(cmd):
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        return result.returncode == 0
    except (FileNotFoundError, subprocess.CalledProcessError, OSError):
        return False


# ---------------------------------------------------------------------------
# reduce (protonation) -- binary works locally
# ---------------------------------------------------------------------------
@pytest.mark.skipif(REDUCE is None, reason="reduce binary not installed")
@pytest.mark.skipif(OBABEL is None, reason="obabel not installed")
def test_reduce_adds_hydrogens(tmp_path):
    pdb = tmp_path / "eth.pdb"
    subprocess.run([OBABEL, "-:CCO", "-opdb", "-O", str(pdb)],
                   check=True, stderr=subprocess.DEVNULL)
    out = subprocess.run([REDUCE, str(pdb)], capture_output=True, text=True)
    assert out.returncode == 0
    # Protonated output must contain at least one hydrogen atom line.
    assert " H" in out.stdout


# ---------------------------------------------------------------------------
# mk_prepare_ligand.py (meeko) -- needs rdkit (not installed here)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(MKP is None or not RDKIT, reason="mk_prepare_ligand.py/rdkit not available")
@pytest.mark.skipif(OBABEL is None, reason="obabel not installed")
def test_mk_prepare_ligand_real(tmp_path):
    mol2 = tmp_path / "lig.mol2"
    subprocess.run([OBABEL, "-:CCO", "-omol2", "-O", str(mol2)],
                   check=True, stderr=subprocess.DEVNULL)
    pdbqt = tmp_path / "lig.pdbqt"
    r = subprocess.run([MKP, "-i", str(mol2), "-o", str(pdbqt)],
                       capture_output=True, text=True)
    assert r.returncode == 0
    assert pdbqt.exists() and pdbqt.stat().st_size > 0


# ---------------------------------------------------------------------------
# qvina -- present but needs libboost_filesystem (not installed, no root)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(QVINA is None, reason="qvina binary not installed")
@pytest.mark.skipif(not _run_ok([QVINA, "--help"]),
                    reason="qvina present but missing libboost (cannot run)")
def test_qvina_docking_real(tmp_path):
    rec = tmp_path / "receptor.pdbqt"
    lig = tmp_path / "ligand.pdbqt"
    out = tmp_path / "out.pdbqt"
    subprocess.run([OBABEL, "-:c1ccccc1", "-opdbqt", "-O", str(rec)],
                   check=True, stderr=subprocess.DEVNULL)
    subprocess.run([OBABEL, "-:CO", "-opdbqt", "-O", str(lig)],
                   check=True, stderr=subprocess.DEVNULL)
    r = subprocess.run([QVINA, "--receptor", str(rec), "--ligand", str(lig),
                       "--out", str(out), "--center_x", "0", "--center_y", "0",
                       "--center_z", "0", "--size_x", "30", "--size_y", "30",
                       "--size_z", "30", "--exhaustiveness", "1"],
                      capture_output=True, text=True)
    assert r.returncode == 0
    assert out.exists()


# ---------------------------------------------------------------------------
# autodock -- precompiled tarball behind a registration wall; not wired in
# ---------------------------------------------------------------------------
@pytest.mark.skipif(AUTODOCK is None, reason="autodock4 binary not installed")
def test_autodock_runs():
    out = subprocess.run([AUTODOCK], capture_output=True, text=True)
    assert out.returncode in (0, 1)
