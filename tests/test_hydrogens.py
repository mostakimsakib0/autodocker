"""Deterministic hydrogen placement for Open Babel output (ligands and
generic files) and deterministic 3D coordinates for 2D ligands."""
import logging
import math
import os
import shutil

import pytest

import vspipeline.conformers as conf
import vspipeline.hydrogens as hyd
from vspipeline import tools, utils

DATA = os.path.join(os.path.dirname(__file__), "data", "hydrogens")
# CASF-2016 1nc3 ligand: Open Babel (-h, MOL2 input) puts the two hydrogens
# of its amino N (atoms 23-24 in the PDBQT) at random: clock-seeded in 3.1.1,
# from OB_RANDOM_SEED in 3.2 (a different seed moves them).
LIGAND_1NC3 = os.path.join(DATA, "1nc3_ligand.mol2")
LIGAND_2D = os.path.join(DATA, "ligand_2d.sdf")

needs_obabel = pytest.mark.skipif(tools.settings().obabel is None, reason="obabel not available")
needs_rdkit = pytest.mark.skipif(not conf.rdkit_available(), reason="RDKit not available")


@pytest.fixture(autouse=True)
def _single_thread(monkeypatch):
    monkeypatch.setenv("OMP_NUM_THREADS", "1")


def _prepare_1nc3(out, seed=None):
    args = ["-imol2", LIGAND_1NC3, "-opdbqt", "-O", out, "-h",
            "--partialcharge", "gasteiger"]
    hyd.run_obabel_deterministic(args, out, kind="ligand", seed=seed)
    return open(out, "rb").read()


@needs_obabel
def test_1nc3_ligand_preparation_is_byte_identical(tmp_path):
    outputs = {_prepare_1nc3(str(tmp_path / f"lig{i}.pdbqt")) for i in range(10)}
    assert len(outputs) == 1
    assert not list(tmp_path.glob("*probe*"))


@needs_obabel
def test_1nc3_ligand_does_not_depend_on_the_seed(tmp_path):
    """The amino hydrogens move with OB_RANDOM_SEED (Open Babel >= 3.2), so
    the probes with other seeds expose them and they are re-placed from the
    heavy atoms: every seed gives the same PDBQT."""
    outputs = {_prepare_1nc3(str(tmp_path / f"lig{seed}.pdbqt"), seed=seed)
               for seed in (1, 42, 1000, 2 ** 31 - 1)}
    assert len(outputs) == 1


@needs_obabel
def test_1nc3_torsion_tree_and_typing_unchanged(tmp_path):
    raw = str(tmp_path / "raw.pdbqt")
    utils.run([tools.settings().obabel, "-imol2", LIGAND_1NC3, "-opdbqt", "-O", raw, "-h",
               "--partialcharge", "gasteiger"], capture=True)
    det = _prepare_1nc3(str(tmp_path / "det.pdbqt")).decode().splitlines(True)
    raw_lines = open(raw, encoding="utf-8").read().splitlines(True)
    assert len(det) == len(raw_lines)
    for a, b in zip(det, raw_lines):
        if a.startswith(("ATOM", "HETATM")):
            assert a[:30] == b[:30] and a[54:] == b[54:]
            if a[77:79].strip() not in hyd.PDBQT_HYDROGEN_TYPES:
                assert a[30:54] == b[30:54]
        else:
            assert a == b  # ROOT/BRANCH/ENDBRANCH/TORSDOF/REMARK untouched
    # The amino hydrogens stay bonded (Vina's rule) to their N only, in the
    # plane of the aromatic ring (sp2).
    atoms = [(l[77:79].strip(), tuple(float(l[i:i + 8]) for i in (30, 38, 46)), l)
             for l in det if l.startswith("ATOM")]
    heavy = [a for a in atoms if a[0] not in hyd.PDBQT_HYDROGEN_TYPES]
    for t, xyz, line in atoms:
        if t == "HD":
            bonded = [h for h in heavy
                      if math.dist(xyz, h[1]) < hyd.vina_bond_cutoff(hyd.pdbqt_element(h[0]))]
            assert len(bonded) == 1, line


def _sdf(atoms, bonds=()):
    lines = ["mol", "  test", "", f"{len(atoms):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000"]
    for el, (x, y, z) in atoms:
        lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {el:<3s} 0  0  0  0  0  0  0  0  0  0  0  0")
    for a, b in bonds:
        lines.append(f"{a:3d}{b:3d}  1  0  0  0  0")
    return "\n".join(lines + ["M  END", "$$$$", ""])


def test_stabilize_sdf_replaces_only_differing_hydrogens(tmp_path):
    heavy = [("C", (0.0, 0.0, 0.0)), ("C", (1.5, 0.0, 0.0)), ("O", (2.2, 1.2, 0.0))]
    a = tmp_path / "a.sdf"
    b = tmp_path / "b.sdf"
    a.write_text(_sdf(heavy + [("H", (3.1, 1.0, 0.3))], [(1, 2), (2, 3), (3, 4)]))
    b.write_text(_sdf(heavy + [("H", (2.0, 2.1, 0.2))], [(1, 2), (2, 3), (3, 4)]))
    before = a.read_text().splitlines()
    assert hyd.stabilize_hydrogens(str(a), str(b), kind="ligand") == 1
    after = a.read_text().splitlines()
    assert [l for l in after if " H " not in l] == [l for l in before if " H " not in l]
    h = [l for l in after if " H " in l][0]
    xyz = (float(h[0:10]), float(h[10:20]), float(h[20:30]))
    assert math.dist(xyz, (2.2, 1.2, 0.0)) == pytest.approx(0.97, abs=1e-3)


def test_run_obabel_deterministic_rejects_heavy_atom_differences(tmp_path):
    calls = []

    def fake_obabel(args, outputs=(), seed=None):
        out = args[args.index("-O") + 1]
        calls.append(out)
        x = 0.0 if len(calls) == 1 else 0.5  # a --gen3d-like heavy-atom change
        with open(out, "w", encoding="utf-8") as f:
            f.write(_sdf([("C", (x, 0.0, 0.0)), ("O", (x + 1.4, 0.0, 0.0))], [(1, 2)]))

    out = str(tmp_path / "lig.sdf")
    with pytest.raises(RuntimeError, match="different heavy atoms"):
        hyd.run_obabel_deterministic(["-isdf", "in.sdf", "-osdf", "-O", out], out,
                                     run_obabel=fake_obabel)
    assert not os.path.exists(out)
    assert not list(tmp_path.glob("*probe*"))


def test_run_obabel_deterministic_identical_runs_untouched(tmp_path):
    text = _sdf([("C", (0.0, 0.0, 0.0)), ("O", (1.4, 0.0, 0.0))], [(1, 2)])

    def fake_obabel(args, outputs=(), seed=None):
        with open(args[args.index("-O") + 1], "w", encoding="utf-8") as f:
            f.write(text)

    out = str(tmp_path / "lig.sdf")
    assert hyd.run_obabel_deterministic(["-O", out], out, run_obabel=fake_obabel) == 0
    assert open(out, encoding="utf-8").read() == text


def test_run_obabel_deterministic_needs_output_in_args(tmp_path):
    with pytest.raises(ValueError):
        hyd.run_obabel_deterministic(["-O", "other.pdbqt"], str(tmp_path / "x.pdbqt"),
                                     run_obabel=lambda *a, **k: None)


@needs_rdkit
def test_rdkit_3d_from_2d_is_byte_identical(tmp_path):
    outputs = set()
    for i in range(5):
        out = str(tmp_path / f"lig{i}.sdf")
        assert conf.generate_3d(LIGAND_2D, out) == conf.METHOD_RDKIT
        outputs.add(open(out, "rb").read())
    assert len(outputs) == 1
    from rdkit import Chem
    mol = next(Chem.SDMolSupplier(out, removeHs=False))
    assert mol.GetNumConformers() == 1 and mol.GetConformer().Is3D()
    assert any(abs(p[2]) > 0.1 for p in mol.GetConformer().GetPositions())


def test_generate_3d_falls_back_to_obabel_with_warning(tmp_path, monkeypatch, caplog):
    monkeypatch.setattr(conf, "rdkit_available", lambda: False)
    seen = []

    def fake_obabel(args, outputs=(), seed=None):
        seen.append((args, seed))
        shutil.copy(LIGAND_2D, args[args.index("-O") + 1])

    out = str(tmp_path / "lig.sdf")
    with caplog.at_level(logging.WARNING):
        assert conf.generate_3d(LIGAND_2D, out, seed=5,
                                run_obabel=fake_obabel) == conf.METHOD_OBABEL
    assert "--gen3d" in seen[0][0] and "NOT reproducible" in caplog.text
    assert seen[0][1] == 5  # OB_RANDOM_SEED = the run seed (Open Babel >= 3.2)
