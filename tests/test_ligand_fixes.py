"""Regression tests for the ligand preparation review findings (L1-L8)."""
import json
import logging
import os
import shutil
import subprocess

import pytest

import vspipeline.library as library
from vspipeline import admet, net, tools, utils

DATA = os.path.join(os.path.dirname(__file__), "data", "prep")

needs_obabel = pytest.mark.skipif(tools.settings().obabel is None, reason="obabel not available")
ENV = {**os.environ, "OMP_NUM_THREADS": "1"}


def _atoms(path):
    return [line for line in open(path, encoding="utf-8")
            if line.startswith(("ATOM", "HETATM"))]


# ---------------------------------------------------------------------------
# Ligands
# ---------------------------------------------------------------------------

def _prepare(tmp_path, files, name="run", **kw):
    indir = tmp_path / f"{name}_in"
    indir.mkdir(exist_ok=True)
    for f in files:
        shutil.copy(os.path.join(DATA, f) if not os.path.isabs(f) else f, indir)
    lm = library.LibraryManager(str(tmp_path / f"{name}_out"), str(indir), **kw)
    try:
        out = lm.prepare_local(apply_lipinski=False)
    except RuntimeError:
        out = []
    return lm, out


def _reasons(lm):
    return {os.path.basename(p): r for p, r in lm.failed}


def test_ligand_pdbqt_with_two_roots_is_rejected(tmp_path):
    p = tmp_path / "two.pdbqt"
    tree = "ROOT\nATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  0.00    +0.000 C\nENDROOT\nTORSDOF 0\n"
    p.write_text(tree + tree)
    assert "2 ROOT blocks" in library._ligand_pdbqt_problem(str(p))
    p.write_text(tree)
    assert library._ligand_pdbqt_problem(str(p)) is None


@needs_obabel
def test_salts_keep_the_largest_fragment_real(tmp_path):
    """L1: salts were prepared with two ROOT blocks and failed in Vina."""
    lm, out = _prepare(tmp_path, ["salt3d_sdf.sdf", "salt3d.mol2", "salt3d_pdb.pdb"],
                       minimize=False)
    assert len(out) == 3, lm.failed
    for path in out:
        assert library._ligand_pdbqt_problem(path) is None
        assert [line.strip() for line in open(path, encoding="utf-8")].count("ROOT") == 1
    assert lm.metadata["salt3d_sdf"]["fragments_removed"] == ["[Cl-]"]
    assert lm.metadata["salt3d"]["fragments_removed"] == ["[Cl-]"]
    assert lm.metadata["salt3d_pdb"]["fragments_removed"] == ["[Na+]"]
    # The input title survives the fragment split; the PDB name record is
    # the input file name.
    fragment = os.path.join(lm.outdir, lm.metadata["salt3d_sdf"]["prepared_sdf"])
    assert fragment.endswith("salt3d_sdf_fragment.sdf")
    with open(os.path.join(DATA, "salt3d_sdf.sdf"), encoding="utf-8") as f:
        title = f.readline()
    assert open(fragment, encoding="utf-8").readline() == title
    mol2_fragment = os.path.join(lm.minimized_dir, "salt3d_fragment.mol2")
    assert open(mol2_fragment, encoding="utf-8").read().splitlines()[1] == "*****"
    pdb_out = [p for p in out if p.endswith("salt3d_pdb.pdbqt")][0]
    assert "REMARK  Name = salt3d_pdb.pdb\n" in open(pdb_out, encoding="utf-8").read()


@needs_obabel
def test_equal_fragments_and_multi_model_pdb_are_rejected_real(tmp_path):
    lm, out = _prepare(tmp_path, ["dimer_sdf.sdf", "twomodel.pdb"], minimize=False)
    assert out == []
    reasons = _reasons(lm)
    assert "share the largest size" in reasons["dimer_sdf.sdf"]
    assert "CCO, CCO" in reasons["dimer_sdf.sdf"]
    assert "2 models" in reasons["twomodel.pdb"]


@needs_obabel
@pytest.mark.parametrize("name", ["flat_ibu.mol2", "octane2d.mol2"])
def test_2d_mol2_and_pdb_get_seeded_3d_coordinates_real(tmp_path, monkeypatch, name):
    """L2: 2D MOL2/PDB inputs were docked flat, and a flat hydrocarbon went
    through a non-reproducible --gen3d charge rebuild."""
    calls = []
    real_run = utils.run

    def spy(cmd, *a, **kw):
        calls.append(list(cmd))
        return real_run(cmd, *a, **kw)

    monkeypatch.setattr(utils, "run", spy)
    outputs = []
    for i in range(3):
        lm, out = _prepare(tmp_path, [name], name=f"r{i}")
        assert len(out) == 1, lm.failed
        outputs.append(open(out[0], "rb").read())
        stem = os.path.splitext(name)[0]
        assert lm.metadata[stem]["coordinates"] == "gen3d"
        assert lm.metadata[stem]["coordinates_method"] == "rdkit-etkdgv3-mmff94"
        assert min(library._pdbqt_coordinate_spread(out[0])) > 0.5
    assert outputs[0] == outputs[1] == outputs[2]
    assert not any("--gen3d" in c for c in calls)


@needs_obabel
def test_2d_pdb_without_readable_bond_orders_is_rejected_real(tmp_path):
    """L2: a flat PDB (no bond orders) was docked flat. Open Babel perceives
    nonsense bond orders from this depiction (C=C=C), so it is rejected
    with the reason instead of being docked."""
    lm, out = _prepare(tmp_path, ["flatpdb_ibu.pdb"])
    assert out == []
    reason = _reasons(lm)["flatpdb_ibu.pdb"]
    assert reason.startswith("2D PDB input: 3D coordinate generation failed")
    assert "bond orders perceived from a flat PDB are unreliable" in reason
    assert "/tmp/autodocker_" not in reason


def test_earlier_output_directory_inside_ligand_dir_is_skipped(tmp_path, monkeypatch, caplog):
    """L3: an earlier run's -o inside -l was read back as ligands."""
    monkeypatch.setattr(library.LibraryManager, "_prepare_local_file",
                        lambda self, fmt, path, lip: path)
    indir = tmp_path / "in"
    (indir / "old_run" / "ligands").mkdir(parents=True)
    (indir / "old_run" / "ligand_metadata.json").write_text("{}")
    (indir / "old_run" / "receptor.pdbqt").write_text("x")
    (indir / "old_run" / "ligands" / "b.pdbqt").write_text("x")
    (indir / "interrupted" / "ligands").mkdir(parents=True)
    (indir / "interrupted" / "ligands_minimized").mkdir()
    (indir / "interrupted" / "ligands" / "c.pdbqt").write_text("x")
    (indir / "sub").mkdir()
    (indir / "sub" / "d.sdf").write_text("x")
    (indir / "a.sdf").write_text("x")
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    with caplog.at_level(logging.WARNING):
        files = [os.path.relpath(p, indir) for _, p in lm._collect_local_files()]
    assert files == ["a.sdf", os.path.join("sub", "d.sdf")]
    assert "old_run inside the ligand directory is the output directory" in caplog.text
    assert "interrupted inside" in caplog.text
    with pytest.raises(ValueError, match="is an AutoDocker output directory"):
        library.LibraryManager(str(tmp_path / "out2"),
                               str(indir / "old_run"))._collect_local_files()


@needs_obabel
def test_pdb_ligand_output_does_not_depend_on_its_folder_real(tmp_path):
    """L4: the PDBQT embedded the input's absolute path (--title src)."""
    src = os.path.join(DATA, "salt3d_pdb.pdb")
    outputs = []
    for folder in ("one", "two/deeper"):
        d = tmp_path / folder
        d.mkdir(parents=True)
        shutil.copy(src, d / "lig.pdb")
        lm = library.LibraryManager(str(tmp_path / f"out_{folder.replace('/', '_')}"), str(d))
        (out,) = lm.prepare_local(apply_lipinski=False)
        outputs.append(open(out, "rb").read())
    assert outputs[0] == outputs[1]
    assert b"REMARK  Name = lig.pdb\n" in outputs[0]
    assert str(tmp_path).encode() not in outputs[0]


@needs_obabel
@pytest.mark.parametrize("fmt", ["sdf", "mol2"])
def test_added_polar_hydrogens_get_charges_real(tmp_path, fmt):
    """L5: hydrogens added by -h to a heavy-atom-only SDF (--no-minimize) or
    MOL2 kept 0.000; only the charge column changes."""
    indir = tmp_path / "in"
    indir.mkdir()
    src = indir / f"acetamide.{fmt}"
    subprocess.run(["obabel", "-:CC(=O)NCCO", "--gen3d", "-d", f"-o{fmt}", "-O", str(src)],
                   check=True, capture_output=True, env=ENV)
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir), minimize=False)
    (out,) = lm.prepare_local(apply_lipinski=False)
    hd = [line for line in _atoms(out) if line[77:79].strip() == "HD"]
    assert len(hd) == 2 and all(float(line[70:76]) > 0.1 for line in hd)

    one_step = tmp_path / "one_step.pdbqt"
    extra = ["--partialcharge", "gasteiger"] if fmt == "sdf" else []
    subprocess.run(["obabel", f"-i{fmt}", str(src), "-opdbqt", "-O", str(one_step), "-h",
                    *extra], check=True, capture_output=True, env=ENV)

    def mask(path):
        return [line[:70] + line[76:] if line.startswith(("ATOM", "HETATM")) else line
                for line in open(path, encoding="utf-8")]

    # Hydrogen positions of the reference run are random; compare the rest.
    def heavy(lines):
        return [line for line in lines if line[77:79].strip() not in ("HD", "H")]

    assert heavy(mask(out)) == heavy(mask(str(one_step)))
    assert len(mask(out)) == len(mask(str(one_step)))


@needs_obabel
def test_tiny_and_small_planar_molecules_are_prepared_real(tmp_path):
    """L6: methane (92 bytes) and water/CO2 (embedded in z = 0 by RDKit) were
    rejected."""
    lm, out = _prepare(tmp_path, ["z2d_methane.sdf", "z2d_water.sdf", "z2d_co2.sdf"])
    assert sorted(os.path.basename(p) for p in out) == [
        "z2d_co2.pdbqt", "z2d_methane.pdbqt", "z2d_water.pdbqt"], lm.failed
    assert os.path.getsize(os.path.join(DATA, "z2d_methane.sdf")) < 150
    props = admet.ADMETFilter.parse_sdf_properties(os.path.join(DATA, "z2d_methane.sdf"))
    assert props is not None and props["mw"] == pytest.approx(16.04, abs=0.01)


@needs_obabel
def test_boron_ligand_failure_names_the_element_real(tmp_path):
    """L7: 'Open Babel produced an empty PDBQT' did not say why."""
    lm, out = _prepare(tmp_path, ["h3d_boron.sdf"], minimize=False)
    assert out == []
    reason = _reasons(lm)["h3d_boron.sdf"]
    assert "element(s) B have no Gasteiger charge parameters" in reason
    assert "converted 0 molecules" in reason and reason.count("molecules converted") == 0


def test_pubchem_source_sdf_is_relative(tmp_path, monkeypatch):
    """L8: source_sdf was absolute for PubChem compounds."""
    def fake_run(cmd, capture=False, extra_env=None, **kw):
        if "-O" in cmd:
            out = cmd[cmd.index("-O") + 1]
            with open(out, "w", encoding="utf-8") as f:
                f.write("ROOT\nATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  "
                        "0.00    -0.100 C\nATOM      2  O   UNL     1       1.200   0.300   "
                        "0.400  0.00  0.00    -0.300 OA\nENDROOT\nTORSDOF 0\n"
                        if out.endswith(".pdbqt") else "x")
        return ("", "", 0) if capture else None

    molblock = ("m\n\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
                "    0.0000    0.0000    0.5000 C   0  0\n"
                "    1.2000    0.0000    0.0000 O   0  0\n  1  2  1  0\nM  END\n$$$$\n")
    monkeypatch.setattr(utils, "run", fake_run)
    monkeypatch.setattr(net, "_http_get_bytes", lambda url: molblock.encode())
    lm = library.LibraryManager(str(tmp_path / "out"), None, minimize=False)
    lm._prepare_pubchem_compound("702", "ethanol", {})
    meta = lm.metadata["ethanol_702"]
    assert meta["source_sdf"] == os.path.join("ligands", "ethanol_702.sdf")
    json.dumps(meta)


def test_converted_ligand_pdbqt_is_validated(tmp_path, monkeypatch):
    """L1: every prepared PDBQT goes through _ligand_pdbqt_problem."""
    tree = ("ROOT\nATOM      1  C   UNL     1       0.000   0.000   0.500  0.00  0.00    "
            "-0.100 C\nENDROOT\nTORSDOF 0\n")

    def fake_run(cmd, capture=False, extra_env=None, **kw):
        if "-O" in cmd:
            with open(cmd[cmd.index("-O") + 1], "w", encoding="utf-8") as f:
                f.write(tree + tree)
        return ("", "", 0) if capture else None

    monkeypatch.setattr(utils, "run", fake_run)
    src = tmp_path / "x.mol2"
    src.write_text("@<TRIPOS>MOLECULE\nx\n")
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path))
    with pytest.raises(library.LigandRejected, match="2 ROOT blocks"):
        lm._convert_to_pdbqt("mol2", str(src), str(tmp_path / "x.pdbqt"), ("-h",))
    assert not (tmp_path / "x.pdbqt").exists()
