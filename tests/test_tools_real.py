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


def _has_real_charges(pdbqt):
    for line in open(pdbqt):
        if line.startswith(("ATOM", "HETATM")):
            chunk = line[54:60].strip()
            try:
                if float(chunk) != 0.0:
                    return True
            except ValueError:
                pass
    return False


def _make_pdbqt(smiles, path, charged=False):
    cmd = [OBABEL, "-:" + smiles, "--gen3d", "-h", "-opdbqt", "-O", str(path)]
    if charged:
        cmd += ["--partialcharge", "gasteiger"]
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)
    return _has_real_charges(str(path))


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
    assert " H" in out.stdout


# ---------------------------------------------------------------------------
# mk_prepare_ligand.py (meeko) -- needs rdkit
# ---------------------------------------------------------------------------
@pytest.mark.skipif(MKP is None or not RDKIT,
                    reason="mk_prepare_ligand.py needs rdkit (pip install rdkit "
                           "times out on this link / no root)")
@pytest.mark.skipif(OBABEL is None, reason="obabel not installed")
def test_mk_prepare_ligand_real(tmp_path):
    mol2 = tmp_path / "lig.mol2"
    subprocess.run([OBABEL, "-:CCO", "-omol2", "-O", str(mol2)],
                   check=True, stderr=subprocess.DEVNULL)
    pdbqt = tmp_path / "lig.pdbqt"
    r = subprocess.run([MKP, "-i", str(mol2), "-o", str(pdbqt)],
                       capture_output=True, text=True)
    assert r.returncode == 0 and pdbqt.exists() and pdbqt.stat().st_size > 0


# ---------------------------------------------------------------------------
# qvina -- present but needs libboost_filesystem (not installed, no root)
# ---------------------------------------------------------------------------
@pytest.mark.skipif(QVINA is None, reason="qvina binary not installed")
@pytest.mark.skipif(not _run_ok([QVINA, "--help"]),
                    reason="qvina present but libboost_filesystem.so.1.84.0 missing "
                           "(no root to install Boost; anaconda throttles boost-cpp)")
@pytest.mark.skipif(OBABEL is None, reason="obabel not installed")
def test_qvina_docking_real(tmp_path):
    rec, lig, out = (tmp_path / "receptor.pdbqt", tmp_path / "ligand.pdbqt",
                     tmp_path / "out.pdbqt")
    subprocess.run([OBABEL, "-:c1ccccc1", "-opdbqt", "-O", str(rec)],
                   check=True, stderr=subprocess.DEVNULL)
    subprocess.run([OBABEL, "-:CO", "-opdbqt", "-O", str(lig)],
                   check=True, stderr=subprocess.DEVNULL)
    r = subprocess.run([QVINA, "--receptor", str(rec), "--ligand", str(lig),
                       "--out", str(out), "--center_x", "0", "--center_y", "0",
                       "--center_z", "0", "--size_x", "30", "--size_y", "30",
                       "--size_z", "30", "--exhaustiveness", "1"],
                      capture_output=True, text=True)
    assert r.returncode == 0 and out.exists()


# ---------------------------------------------------------------------------
# autodock -- binary installed, but a real docking needs a CHARGED receptor
# PDBQT. obabel in this env cannot compute partial charges (gasteiger returns
# '??') and MGLTools is absent, so a genuine docking run is not possible here.
# The test attempts a real run and skips cleanly when charges can't be made.
# ---------------------------------------------------------------------------
@pytest.mark.skipif(AUTODOCK is None, reason="autodock4 binary not installed")
@pytest.mark.skipif(OBABEL is None, reason="obabel not installed")
def test_autodock_real_docking(tmp_path):
    rec = tmp_path / "rec.pdbqt"
    lig = tmp_path / "lig.pdbqt"
    if not _make_pdbqt("O=C(O)c1ccccc1", rec, charged=True) or \
       not _make_pdbqt("CO", lig, charged=True):
        pytest.skip("obabel cannot compute partial charges here (MGLTools needed "
                    "to build a charged receptor for AutoDock)")

    (tmp_path / "rec.gpf").write_text(
        "receptor_file rec.pdbqt\n"
        "gridcenter 0.0 0.0 0.0\nnpts 20 20 20\nspacing 1.0\n"
        "gridfld rec.maps.fld\nmap rec.C.map\nmap rec.A.map\nmap rec.HD.map\n"
        "map rec.OA.map\nmap rec.N.map\nmap rec.NA.map\nmap rec.SA.map\n"
        "map rec.e.map\nelecmap rec.e.map\ndesolvmap rec.d.map\n"
        "dielectric -0.1465\n")
    subprocess.run([AUTODOCK.replace("autodock4", "autogrid4"),
                    "-p", str(tmp_path / "rec.gpf"), "-l", str(tmp_path / "rec.glg")],
                   check=True, capture_output=True, text=True)
    assert (tmp_path / "rec.maps.fld").exists()

    (tmp_path / "rec.dpf").write_text(
        "receptor_file rec.pdbqt\nligand_file lig.pdbqt\nmaps rec\n"
        "flexible_ligand 1\nga_run 10\noutlev 1\n")
    r = subprocess.run([AUTODOCK, "-p", str(tmp_path / "rec.dpf"),
                        "-l", str(tmp_path / "rec.dlg")],
                       capture_output=True, text=True)
    assert r.returncode == 0
    assert "Estimated Free Energy" in (tmp_path / "rec.dlg").read_text()
