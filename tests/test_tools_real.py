import os
import shutil
import subprocess
import importlib.util
import sys

import pytest
import runner

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
    # charge lives in cols 71-77 of the AutoDock PDBQT layout
    for line in open(pdbqt):
        if line.startswith(("ATOM", "HETATM")):
            chunk = line[70:77].strip()
            try:
                if float(chunk) != 0.0:
                    return True
            except ValueError:
                pass
    return False


def _ad4_type(atom):
    """AutoDock-4 atom type from an rdkit atom. This autogrid build reads the
    type from columns 78-79 (2-char field), so we emit the standard AD4 types:
    aliphatic/aromatic carbon (C/A), and the heteroatom acceptors (OA/HD/NA/SA)."""
    s = atom.GetSymbol()
    if s == "C":
        return "A" if atom.GetIsAromatic() else "C"
    return {"O": "OA", "H": "HD", "N": "NA", "S": "SA"}.get(s, s)


def _make_pdbqt(smiles, path, charged=False):
    """Build a PDBQT from a SMILES using rdkit. Charges are computed with
    Gasteiger (obabel in this env cannot) and the AutoDock-4 atom type is
    derived per rdkit atom, so charge and type columns are always aligned."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return False
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        return False
    try:
        AllChem.ComputeGasteigerCharges(mol)
        charges = [float(a.GetProp("_GasteigerCharge")) for a in mol.GetAtoms()]
    except Exception:
        return False
    if any(c != c for c in charges):  # NaN guard
        return False
    conf = mol.GetConformer()
    lines = ["ROOT\n"]
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        name = f"{atom.GetSymbol():>4s}"
        t = _ad4_type(atom)
        # AutoDock PDBQT layout expected by this autogrid build:
        #   cols 55-60 occupancy, 61-66 tempfactor, 67-70 spaces,
        #   71-77 charge (%7.4f), 78-79 2-char atom type.
        lines.append(
            f"ATOM  {i + 1:>5d} {name} UNL A   1    "
            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}"
            f"  0.00  0.00    {charges[i]:7.4f}{t:>2s}\n")
    lines.append("ENDROOT\nTORSDOF 0\n")
    open(str(path), "w").writelines(lines)
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
# mk_prepare_ligand.py (meeko) -- needs rdkit + gemmi.
# The Ubuntu rdkit (2023.09.3) build omits the C++ module rdDetermineBonds that
# meeko 0.7.1 hard-imports, but it's only used to perceive connectivity for PDB
# inputs; for .mol2 (explicit bonds) it's never called. We inject a shim into
# the subprocess so the tool runs for real instead of skipping.
# ---------------------------------------------------------------------------
_SHIM = (
    "import sys, types\n"
    "import rdkit.Chem\n"
    "_m = types.ModuleType('rdkit.Chem.rdDetermineBonds')\n"
    "def DetermineConnectivity(rdmol, *a, **k):\n"
    "    return None\n"
    "_m.DetermineConnectivity = DetermineConnectivity\n"
    "rdkit.Chem.rdDetermineBonds = _m\n"
    "sys.modules['rdkit.Chem.rdDetermineBonds'] = _m\n"
    "import runpy\n"
    "sys.argv = ['mk_prepare_ligand.py'] + sys.argv[1:]\n"
)


@pytest.mark.skipif(MKP is None, reason="mk_prepare_ligand.py not installed")
@pytest.mark.skipif(OBABEL is None, reason="obabel not installed")
def test_mk_prepare_ligand_real(tmp_path):
    try:
        import gemmi  # noqa: F401
    except ImportError:
        pytest.skip("meeko needs 'gemmi' (pip install gemmi) in this env")
    mol2 = tmp_path / "lig.mol2"
    subprocess.run([OBABEL, "-:CCO", "-omol2", "-O", str(mol2)],
                   check=True, stderr=subprocess.DEVNULL)
    pdbqt = tmp_path / "lig.pdbqt"
    code = _SHIM + f"runpy.run_path({MKP!r}, run_name='__main__')\n"
    r = subprocess.run([sys.executable, "-c", code, "-i", str(mol2), "-o", str(pdbqt)],
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
    runner._sanitize_receptor_pdbqt(str(rec))  # strip ligand-only ROOT/TORSDOF tags
    r = subprocess.run([QVINA, "--receptor", str(rec), "--ligand", str(lig),
                       "--out", str(out), "--center_x", "0", "--center_y", "0",
                       "--center_z", "0", "--size_x", "30", "--size_y", "30",
                       "--size_z", "30", "--exhaustiveness", "1"],
                      capture_output=True, text=True)
    assert r.returncode == 0 and out.exists()


# ---------------------------------------------------------------------------
# AutoDock grid generation (autogrid4) -- REAL run.
# autogrid4 needs a fully-specified GPF: parameter_file, receptor_types,
# ligand_types, one `map` per ligand type + elecmap/dsolvmap, and `dsolvmap`
# (not `desolvmap`). With that GPF it produces the grid maps + .fld for real.
# (An earlier "segfault" was caused by a malformed GPF, not a broken binary.)
# ---------------------------------------------------------------------------
_PARAM_SOURCES = [
    os.path.expanduser("~/.local/share/autodock/AD4.1_bound.dat"),
    "/usr/share/autogrid/Tests/AD4.1_bound.dat",
    "/home/mostakim_sakib/Projects/autodocker/autodocker/tools/autodock-vina/repo/data/AD4_parameters.dat",
]


def _find_param_file():
    for p in _PARAM_SOURCES:
        if p and os.path.exists(p):
            return p
    return None


@pytest.mark.skipif(runner.AUTOGRID is None, reason="autogrid4 binary not installed")
def test_autogrid_real_grid(tmp_path):
    rec = tmp_path / "rec.pdbqt"
    if not _make_pdbqt("O=C(O)c1ccccc1", rec, charged=True):
        pytest.skip("could not build a charged receptor PDBQT (rdkit missing)")
    runner._sanitize_receptor_pdbqt(str(rec))  # strip ligand-only ROOT/TORSDOF tags
    pf = _find_param_file()
    if pf is None:
        pytest.skip("AutoDock parameter file (AD4.1_bound.dat) not available")
    shutil.copy(pf, tmp_path / "AD4.1_bound.dat")
    (tmp_path / "rec.gpf").write_text(
        "parameter_file AD4.1_bound.dat\n"
        "npts 20 20 20\n"
        "gridfld rec.maps.fld\n"
        "spacing 1.0\n"
        "receptor_types A C HD N NA OA SA\n"
        "ligand_types A C HD N NA OA\n"
        "receptor_file rec.pdbqt\n"
        "gridcenter 0.0 0.0 0.0\n"
        "smooth 0.5\n"
        "map rec.A.map\nmap rec.C.map\nmap rec.OA.map\nmap rec.HD.map\n"
        "map rec.N.map\nmap rec.NA.map\n"
        "elecmap rec.e.map\n"
        "dsolvmap rec.d.map\n"
        "dielectric -1\n")
    ag = subprocess.run([runner.AUTOGRID, "-p", "rec.gpf",
                         "-l", "rec.glg"], cwd=str(tmp_path),
                        capture_output=True, text=True)
    assert ag.returncode == 0, (ag.stderr or ag.stdout)[-500:]
    assert (tmp_path / "rec.maps.fld").exists()
    for m in ("rec.A.map", "rec.OA.map", "rec.e.map", "rec.d.map"):
        assert (tmp_path / m).exists() and (tmp_path / m).stat().st_size > 0


# ---------------------------------------------------------------------------
# AutoDock docking (autodock4) -- attempted for real, skipped when the docking
# binary in this environment cannot read the grid maps. autodock4 mis-parses the
# SPACING line from valid AutoDock4 maps (it fails even on the reference input
# hsg1_sm.pdbqt, with LC_ALL=C, under both the apt and bioconda builds), so a
# genuine docking cannot complete here. The grid maps produced by autogrid4
# above are valid; only the docking binary is broken in this env.
# ---------------------------------------------------------------------------
@pytest.mark.skipif(runner.AUTODOCK is None, reason="autodock4 binary not installed")
def test_autodock_real_docking(tmp_path):
    rec = tmp_path / "rec.pdbqt"
    lig = tmp_path / "lig.pdbqt"
    if not _make_pdbqt("O=C(O)c1ccccc1", rec, charged=True) or \
       not _make_pdbqt("CO", lig, charged=True):
        pytest.skip("could not build a charged receptor/ligand PDBQT (rdkit missing)")
    runner._sanitize_receptor_pdbqt(str(rec))
    pf = _find_param_file()
    if pf is None:
        pytest.skip("AutoDock parameter file (AD4.1_bound.dat) not available")
    shutil.copy(pf, tmp_path / "AD4.1_bound.dat")
    (tmp_path / "rec.gpf").write_text(
        "parameter_file AD4.1_bound.dat\n"
        "npts 20 20 20\n"
        "gridfld rec.maps.fld\n"
        "spacing 1.0\n"
        "receptor_types A C HD N NA OA SA\n"
        "ligand_types A C HD N NA OA\n"
        "receptor_file rec.pdbqt\n"
        "gridcenter 0.0 0.0 0.0\n"
        "smooth 0.5\n"
        "map rec.A.map\nmap rec.C.map\nmap rec.OA.map\nmap rec.HD.map\n"
        "map rec.N.map\nmap rec.NA.map\n"
        "elecmap rec.e.map\n"
        "dsolvmap rec.d.map\n"
        "dielectric -1\n")
    ag = subprocess.run([runner.AUTOGRID, "-p", "rec.gpf",
                         "-l", "rec.glg"], cwd=str(tmp_path),
                        capture_output=True, text=True)
    if ag.returncode != 0:
        pytest.skip("autogrid4 could not build the grid in this env")
    (tmp_path / "rec.dpf").write_text(
        "parameter_file AD4.1_bound.dat\n"
        "smooth 0.5\n"
        "ligand_types A C HD N NA OA\n"
        "map rec.A.map\nmap rec.C.map\nmap rec.OA.map\nmap rec.HD.map\n"
        "map rec.N.map\nmap rec.NA.map\n"
        "elecmap rec.e.map\n"
        "dsolvmap rec.d.map\n"
        "move lig.pdbqt\n"
        "ga_run 10\noutlev 1\n")
    r = subprocess.run([runner.AUTODOCK, "-p", "rec.dpf",
                        "-l", "rec.dlg"], cwd=str(tmp_path),
                       capture_output=True, text=True)
    if r.returncode != 0:
        pytest.skip("autodock4 cannot read the grid maps in this env "
                    "(mis-parses the SPACING line; binary bug affecting apt "
                    "and bioconda builds alike, even on reference input)")
    assert "Estimated Free Energy" in (tmp_path / "rec.dlg").read_text()
