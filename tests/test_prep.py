import os
import shutil

import pytest

import runner

obabel = shutil.which("obabel")
pytestmark = pytest.mark.skipif(
    obabel is None, reason="Open Babel not available on PATH")


def _make_ethanol_sdf(tmp_path):
    """Build an ethanol SDF via obabel (has nonzero partial charges)."""
    path = str(tmp_path / "ethanol.sdf")
    runner.run([runner.OBABEL, "-:CCO", "-osdf", "-O", path, "--gen3d"])
    return path


def test_sdf_to_pdbqt_with_charges(tmp_path):
    sdf = _make_ethanol_sdf(tmp_path)
    out = str(tmp_path / "ligand.pdbqt")
    runner.run([runner.OBABEL, "-isdf", sdf, "-opdbqt", "-O", out,
                "--partialcharge", "gasteiger"])
    assert os.path.getsize(out) > 0
    assert runner._pdbqt_has_atoms(out)
    assert runner._ensure_pdbqt_has_charges(out)


def test_mol2_to_pdbqt_with_charges(tmp_path):
    mol2 = str(tmp_path / "ethanol.mol2")
    runner.run([runner.OBABEL, "-:CCO", "-omol2", "-O", mol2, "--gen3d"])
    out = str(tmp_path / "ligand.pdbqt")
    runner.run([runner.OBABEL, "-imol2", mol2, "-opdbqt", "-O", out,
                "--partialcharge", "gasteiger"])
    assert os.path.getsize(out) > 0
    assert runner._pdbqt_has_atoms(out)
    assert runner._ensure_pdbqt_has_charges(out)


def test_pdb_to_pdbqt_with_charges(protein_pdb, tmp_path):
    out = str(tmp_path / "ligand.pdbqt")
    runner.run([runner.OBABEL, "-ipdb", protein_pdb, "-opdbqt", "-O", out])
    assert os.path.getsize(out) > 0
    assert runner._pdbqt_has_atoms(out)


def test_fix_pdbqt_charges_assigns_charges(zero_charge_pdbqt, tmp_path):
    sdf = _make_ethanol_sdf(tmp_path)
    runner._fix_pdbqt_charges(zero_charge_pdbqt, sdf)
    assert runner._ensure_pdbqt_has_charges(zero_charge_pdbqt) is True


def test_fix_pdbqt_charges_raises_on_bad_source(tmp_path):
    bad = str(tmp_path / "bad.pdbqt")
    with open(bad, "w") as f:
        f.write("ATOM  junk\n")
    with pytest.raises(Exception) as exc:
        runner._fix_pdbqt_charges(bad, str(tmp_path / "missing.sdf"))
    assert "Cannot assign valid charges" in str(exc.value)


def test_library_prepare_local_sdf_end_to_end(tmp_path):
    d = tmp_path / "ligands"
    d.mkdir()
    sdf = _make_ethanol_sdf(tmp_path)
    shutil.copy(sdf, str(d / "ethanol.sdf"))
    lm = runner.LibraryManager(str(tmp_path / "out"), str(d))
    out = lm._prepare_local_sdf(apply_admet=False)
    assert len(out) == 1
    assert runner._pdbqt_has_atoms(out[0])
    assert runner._ensure_pdbqt_has_charges(out[0])


def test_library_prepare_local_sdf_batch_recursion(tmp_path):
    d = tmp_path / "ligands" / "sub" / "nested"
    d.mkdir(parents=True)
    sdf = _make_ethanol_sdf(tmp_path)
    shutil.copy(sdf, str(d / "nested.sdf"))
    lm = runner.LibraryManager(str(tmp_path / "out"), str(tmp_path / "ligands"))
    out = lm._prepare_local_sdf(apply_admet=False)
    assert len(out) == 1
    assert runner._pdbqt_has_atoms(out[0])


def test_blank_chain_pdb_kept_for_default_chain(protein_pdb, tmp_path):
    """PDBs with blank chain IDs (common in DUD-E-style receptors) must keep
    their atoms when the default chain 'A' is selected."""
    lines = []
    with open(protein_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                lines.append(line[:21] + " " + line[22:])
            else:
                lines.append(line.rstrip("\n"))
    pdbfile = tmp_path / "blankchain.pdb"
    pdbfile.write_text("\n".join(lines) + "\n")
    prep = runner.ProteinPreparation(str(pdbfile), str(tmp_path))
    prep.prepare_receptor("A")
    n = sum(1 for l in open(prep.receptor_pdbqt)
            if l.startswith(("ATOM", "HETATM")))
    assert n > 0


def test_receptor_centroid_offsets_frame_shift(protein_pdb, tmp_path):
    """The fpocket grid is computed on pdb_clean (original frame) while the
    receptor PDBQT is centered at the origin; _receptor_centroid must return
    the original-frame centroid so detect_pocket can shift the grid center."""
    lines = []
    with open(protein_pdb) as f:
        for line in f:
            if line.startswith("ATOM"):
                x = float(line[30:38]) + 100.0
                y = float(line[38:46]) + 200.0
                z = float(line[46:54]) + 300.0
                lines.append(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
            else:
                lines.append(line.rstrip("\n"))
    pdbfile = tmp_path / "offset.pdb"
    pdbfile.write_text("\n".join(lines) + "\n")
    prep = runner.ProteinPreparation(str(pdbfile), str(tmp_path))
    prep.prepare_receptor("A")

    def centroid(path):
        xs, ys, zs = [], [], []
        for l in open(path):
            if l.startswith("ATOM"):
                xs.append(float(l[30:38]))
                ys.append(float(l[38:46]))
                zs.append(float(l[46:54]))
        return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)

    exp = centroid(prep.pdb_clean)
    got = prep._receptor_centroid()
    assert abs(got[0] - exp[0]) < 1e-3
    assert abs(got[1] - exp[1]) < 1e-3
    assert abs(got[2] - exp[2]) < 1e-3
    assert abs(exp[0] - 104.5) < 1e-3

    rctr = centroid(prep.receptor_pdbqt)
    for c in rctr:
        assert abs(c) < 1e-3