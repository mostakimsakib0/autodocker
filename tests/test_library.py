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


def test_prepare_local_sdf_pdbqt_direct(tmp_path):
    indir = tmp_path / "in"
    indir.mkdir()
    lig = indir / "mol1.pdbqt"
    lig.write_text("ATOM      1  C   ALA A   1      0.0 0.0 0.0  0.0 0.0    -0.1 C\n")
    lm = library.LibraryManager(str(tmp_path), str(indir))
    out = lm._prepare_local_sdf(apply_admet=False)
    assert out == [str(lig)]


def test_prepare_local_sdf_missing_dir(tmp_path):
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "nope"))
    try:
        lm._prepare_local_sdf()
        assert False, "expected FileNotFoundError"
    except FileNotFoundError:
        pass


def test_prepare_local_sdf_no_ligands(tmp_path):
    indir = tmp_path / "in"
    indir.mkdir()
    lm = library.LibraryManager(str(tmp_path), str(indir))
    try:
        lm._prepare_local_sdf()
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
