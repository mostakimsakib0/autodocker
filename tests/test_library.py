import os

import runner
import vspipeline.library as library


def _fake_http(url):
    if "cids/JSON" in url or "cids" in url:
        return b'{"IdentifierList":{"CID":[12345]}}'
    if "property" in url:
        return b'{"PropertyTable":{"Properties":[{"MolecularWeight":"300","XLogP":"2"}]}}'
    if "SDF" in url:
        return b"sdf-bytes"
    return b""


def test_pubchem_cid_for_name(monkeypatch):
    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)
    lm = library.LibraryManager("/tmp/libtest", "/tmp/libtest/in")
    assert lm._pubchem_cid_for_name("aspirin") == ["12345"]


def test_pubchem_cid_for_name_error(monkeypatch):
    monkeypatch.setattr(runner, "_http_get_bytes",
                        lambda u: (_ for _ in ()).throw(RuntimeError("x")))
    lm = library.LibraryManager("/tmp/libtest", "/tmp/libtest/in")
    assert lm._pubchem_cid_for_name("aspirin") == []


def test_pubchem_properties(monkeypatch):
    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)
    lm = library.LibraryManager("/tmp/libtest", "/tmp/libtest/in")
    props = lm._pubchem_properties("123")
    assert props["MolecularWeight"] == "300"


def test_pubchem_download_sdf(monkeypatch, tmp_path):
    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    out = lm._pubchem_download_sdf("123", str(tmp_path / "c.sdf"))
    assert open(out, "rb").read() == b"sdf-bytes"


def test_save_metadata(tmp_path):
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    lm.metadata["x_1"] = {"a": 1}
    lm._save_metadata()
    assert os.path.exists(lm.metadata_file)
    import json
    assert json.load(open(lm.metadata_file))["x_1"]["a"] == 1


_ATOM = "ATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  0.00    -0.100 C\n"
_LIGAND_PDBQT = "ROOT\n" + _ATOM + "ENDROOT\nTORSDOF 0\n"


def test_prepare_local_sdf_pdbqt_direct(tmp_path):
    indir = tmp_path / "in"
    indir.mkdir()
    lig = indir / "mol1.pdbqt"
    lig.write_text(_LIGAND_PDBQT)
    lm = library.LibraryManager(str(tmp_path), str(indir))
    out = lm.prepare_local(apply_admet=False)
    assert out == [str(lig)]


def test_prepare_local_sdf_missing_dir(tmp_path):
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "nope"))
    try:
        lm.prepare_local()
        assert False, "expected FileNotFoundError"
    except FileNotFoundError:
        pass


def test_prepare_local_sdf_no_ligands(tmp_path):
    indir = tmp_path / "in"
    indir.mkdir()
    lm = library.LibraryManager(str(tmp_path), str(indir))
    try:
        lm.prepare_local()
        assert False, "expected FileNotFoundError"
    except FileNotFoundError as e:
        assert "No ligands found" in str(e)


def test_create_fda_library(monkeypatch, tmp_path):
    indir = tmp_path / "in"
    indir.mkdir()

    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)

    def fake_run(cmd, **kw):
        if "-O" in cmd:
            out = cmd[cmd.index("-O") + 1]
            with open(out, "w") as f:
                f.write("ATOM      1  C   ALA A   1      0.0 0.0 0.0  0.0 0.0    -0.1 C\n")

    monkeypatch.setattr(runner, "run", fake_run)
    monkeypatch.setattr(runner, "_ensure_pdbqt_has_charges", lambda p: True)
    monkeypatch.setattr(
        runner.ADMETFilter, "parse_sdf_properties",
        staticmethod(lambda f: {"mw": 300, "logp": 2, "hba": 2,
                                "hbd": 1, "tpsa": 60, "rotors": 3}))
    monkeypatch.setattr(
        runner.ADMETFilter, "check_lipinski",
        staticmethod(lambda p: (True, [])))

    lm = library.LibraryManager(str(tmp_path), str(indir))
    out = lm.create_fda_library(apply_admet=True)
    assert len(out) >= 1
    assert os.path.exists(lm.metadata_file)


def test_prepare_local_rejects_receptor_style_pdbqt(tmp_path):
    """A PDBQT without a torsion tree (ROOT/TORSDOF) is rejected up front
    instead of failing later inside Vina's parser."""
    indir = tmp_path / "in"
    indir.mkdir()
    good = indir / "good.pdbqt"
    good.write_text(_LIGAND_PDBQT)
    (indir / "c1.pdbqt").write_text(_ATOM)
    lm = library.LibraryManager(str(tmp_path), str(indir))
    assert lm.prepare_local(apply_admet=False) == [str(good)]
    assert "ROOT" in library._ligand_pdbqt_problem(str(indir / "c1.pdbqt"))
    assert library._ligand_pdbqt_problem(str(good)) is None


def test_prepare_local_mixed_formats(tmp_path, monkeypatch):
    """A folder with PDBQT and PDB ligands prepares both, not just the first
    format found."""
    indir = tmp_path / "in"
    indir.mkdir()
    lig = indir / "lib1.pdbqt"
    lig.write_text(_LIGAND_PDBQT)
    (indir / "control.pdb").write_text("HETATM    1  C   UNL     1       0.000   0.000   0.000  1.00  0.00           C\n")

    def fake_run(cmd, *a, **k):
        with open(cmd[cmd.index("-O") + 1], "w") as f:
            f.write(_LIGAND_PDBQT)
    monkeypatch.setattr(runner, "run", fake_run)
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    out = lm.prepare_local(apply_admet=False)
    assert sorted(os.path.basename(p) for p in out) == ["control.pdbqt", "lib1.pdbqt"]


def test_prepare_local_duplicate_names_skipped(tmp_path, monkeypatch):
    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "x.pdbqt").write_text(_LIGAND_PDBQT)
    (indir / "x.pdb").write_text("HETATM    1  C   UNL     1       0.000   0.000   0.000  1.00  0.00           C\n")
    monkeypatch.setattr(runner, "run", lambda *a, **k: None)
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    assert lm.prepare_local(apply_admet=False) == [str(indir / "x.pdbqt")]
