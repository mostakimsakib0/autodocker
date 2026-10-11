import json
import logging
import os
import shutil
import subprocess

import pytest

import runner
import vspipeline.library as library
import vspipeline.net as net

DATA = os.path.join(os.path.dirname(__file__), "data", "runner")
HAVE_OBABEL = shutil.which("obabel") is not None
needs_obabel = pytest.mark.skipif(not HAVE_OBABEL, reason="obabel not installed")


@pytest.fixture(autouse=True)
def _obabel_3d_unless_real(request, monkeypatch):
    """Unit tests fake Open Babel with placeholder files RDKit cannot read, so
    2D inputs take the Open Babel --gen3d path of conformers.generate_3d
    (through the faked runner.run). Tests with real tools keep RDKit."""
    if "real" not in request.node.name:
        from vspipeline import conformers
        monkeypatch.setattr(conformers, "rdkit_available", lambda: False)


def _fake_http(url):
    if "cids/JSON" in url or "cids" in url:
        return b'{"IdentifierList":{"CID":[12345]}}'
    if "SDF" in url:
        return b"sdf-bytes"
    return b""


def _molblock(z, z2=0.0):
    return ("m\n  test\n\n"
            "  3  2  0  0  0  0  0  0  0  0999 V2000\n"
            f"    0.0000    0.0000    {z:.4f} C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            f"    1.5000    0.0000    {z2:.4f} O   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "    0.0000    1.5000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "  1  2  1  0\n  1  3  1  0\nM  END\n$$$$\n")


_ATOM = "ATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  0.00    -0.100 C\n"
_LIGAND_PDBQT = "ROOT\n" + _ATOM + "ENDROOT\nTORSDOF 0\n"
# A 3-atom ligand PDBQT with extent along every axis.
_LIGAND_PDBQT_3D = (
    "ROOT\n"
    "ATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  0.00    -0.100 C\n"
    "ATOM      2  O   UNL     1       1.200   0.300   0.400  0.00  0.00    -0.300 OA\n"
    "ATOM      3  N   UNL     1      -0.500   1.100  -0.700  0.00  0.00     0.200 NA\n"
    "ENDROOT\nTORSDOF 0\n")
_LIGAND_PDBQT_FLAT = _LIGAND_PDBQT_3D.replace("   0.400", "   0.000").replace("  -0.700", "   0.000")


def _fake_obabel(calls=None, pdbqt=_LIGAND_PDBQT_3D, fail_on=None):
    """A runner.run replacement writing plausible Open Babel outputs."""
    def fake_run(cmd, capture=False, extra_env=None, **kw):
        if calls is not None:
            calls.append((list(cmd), extra_env))
        if fail_on and fail_on in cmd:
            raise subprocess.CalledProcessError(1, cmd, output="", stderr="boom")
        if "-O" in cmd:
            out = cmd[cmd.index("-O") + 1]
            with open(out, "w") as f:
                if out.endswith(".sdf"):
                    f.write(_molblock(0.8, 0.3))
                elif out.endswith(".pdb"):
                    f.write("HETATM    1  C   UNL     1       0.000   0.000   0.000  1.00  0.00           C\n")
                else:
                    f.write(pdbqt)
        return ("", "", 0) if capture else None
    return fake_run


def _lipinski_ok(monkeypatch):
    monkeypatch.setattr(
        runner.ADMETFilter, "parse_sdf_properties",
        staticmethod(lambda f: {"mw": 300.0, "logp": 2.0, "hba": 2.0, "hbd": 1.0,
                                "tpsa": 60.0, "rotors": 3.0,
                                "properties_source": "obabel_descriptors",
                                "logp_model": "Open Babel logP"}))


# ---------------------------------------------------------------------------
# PubChem
# ---------------------------------------------------------------------------

def test_pubchem_cid_for_name(monkeypatch, tmp_path):
    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    assert lm._pubchem_cid_for_name("aspirin") == ["12345"]


def test_pubchem_cid_for_name_error_is_logged_as_warning(monkeypatch, tmp_path, caplog):
    monkeypatch.setattr(runner, "_http_get_bytes",
                        lambda u: (_ for _ in ()).throw(ValueError("bad json")))
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    with caplog.at_level(logging.WARNING):
        assert lm._pubchem_cid_for_name("aspirin") == []
    assert "name lookup failed" in caplog.text


def test_pubchem_download_sdf_records_provenance(monkeypatch, tmp_path):
    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    out = str(tmp_path / "c.sdf")
    prov = lm._pubchem_download_sdf("123", out)
    assert open(out, "rb").read() == b"sdf-bytes"
    assert prov["download_url"].endswith("/compound/cid/123/SDF?record_type=3d")
    assert prov["pubchem_record_type"] == "3d"
    assert prov["retrieved_at"].endswith("+00:00")


def test_save_metadata(tmp_path):
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    lm.metadata["x_1"] = {"a": 1}
    lm._save_metadata()
    assert json.load(open(lm.metadata_file))["x_1"]["a"] == 1


def test_create_fda_library_records_pubchem_provenance(monkeypatch, tmp_path):
    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)
    monkeypatch.setattr(runner, "run", _fake_obabel())
    _lipinski_ok(monkeypatch)
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    out = lm.create_fda_library(apply_lipinski=True, max_compounds=2)
    assert len(out) == 2  # --max-compounds is honoured
    meta = json.load(open(lm.metadata_file))
    entry = meta["aspirin_12345"]
    assert entry["pubchem_cid"] == "12345"
    assert entry["pubchem_query"] == "aspirin"
    assert entry["source_url"] == "https://pubchem.ncbi.nlm.nih.gov/compound/12345"
    assert entry["download_url"].endswith("/compound/cid/12345/SDF?record_type=3d")
    assert entry["retrieved_at"]
    assert entry["logp_model"] == "Open Babel logP"


def test_fda_list_has_no_withdrawn_ranitidine():
    assert "ranitidine" not in library.LibraryManager.FDA_DRUG_NAMES
    assert len(library.LibraryManager.FDA_DRUG_NAMES) == 74
    assert library.LibraryManager.MAX_COMPOUNDS == 20


def test_pubchem_does_not_fall_back_to_local_ligands(monkeypatch, tmp_path):
    """No silent substitution of the user's -l ligands for the PubChem library."""
    class Resp:
        status_code = 404

    class NotFound(OSError):
        response = Resp()

    def not_found(url):
        raise NotFound("404")

    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "local.pdbqt").write_text(_LIGAND_PDBQT)
    monkeypatch.setattr(runner, "_http_get_bytes", not_found)
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    with pytest.raises(RuntimeError, match="could be prepared from PubChem"):
        lm.create_fda_library()
    with pytest.raises(RuntimeError, match="not used in --library custom"):
        lm.create_custom_library()


def test_pubchem_unreachable_stops_early(monkeypatch, tmp_path):
    calls = []

    def offline(url):
        calls.append(url)
        raise ConnectionRefusedError("refused")

    monkeypatch.setattr(runner, "_http_get_bytes", offline)
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path / "in"))
    with pytest.raises(RuntimeError, match="unreachable"):
        lm.create_fda_library()
    assert len(calls) == library.LibraryManager.MAX_NETWORK_FAILURES


def test_persistent_429_counts_as_network_failure(monkeypatch, tmp_path):
    def busy(url):
        raise net._RetryableHTTPError("HTTP 429")

    monkeypatch.setattr(runner, "_http_get_bytes", busy)
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path / "in"))
    with pytest.raises(library.PubChemUnreachable):
        lm.create_fda_library()


def test_pubchem_unknown_name_is_not_a_network_failure(monkeypatch, tmp_path):
    class Resp:
        status_code = 404

    class NotFound(OSError):
        response = Resp()

    def not_found(url):
        raise NotFound("404")

    monkeypatch.setattr(runner, "_http_get_bytes", not_found)
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path / "in"))
    for _ in range(5):
        assert lm._pubchem_cid_for_name("nosuchdrug") == []


def test_pubchem_preparation_failures_are_warnings_with_summary(monkeypatch, tmp_path, caplog):
    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)
    monkeypatch.setattr(runner, "run", _fake_obabel(fail_on="-opdbqt"))
    _lipinski_ok(monkeypatch)
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    with caplog.at_level(logging.WARNING), pytest.raises(RuntimeError):
        lm.create_fda_library(max_compounds=1)
    assert "could not be prepared" in caplog.text
    assert len(lm.failed) == len(lm.FDA_DRUG_NAMES)
    # The failed conversion leaves no partial or placeholder PDBQT behind.
    assert not [f for f in os.listdir(lm.lib_dir) if f.endswith(".pdbqt")]


def _custom_props(mw, logp, hbd=1.0, hba=2.0):
    return staticmethod(lambda f: {"mw": mw, "logp": logp, "hbd": hbd, "hba": hba,
                                   "tpsa": 1.0, "rotors": 1.0,
                                   "properties_source": "obabel_descriptors",
                                   "logp_model": "Open Babel logP"})


def test_custom_thresholds_replace_lipinski(monkeypatch, tmp_path):
    """--mw-max 700 admits a MW 600 compound that Lipinski would reject."""
    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)
    monkeypatch.setattr(runner, "run", _fake_obabel())
    monkeypatch.setattr(runner.ADMETFilter, "parse_sdf_properties", _custom_props(600.0, 6.5))
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    out = lm.create_custom_library(mw_min=200, mw_max=700, logp_max=7, max_compounds=1)
    assert len(out) == 1


def test_custom_thresholds_reject_outside_window(monkeypatch, tmp_path):
    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)
    monkeypatch.setattr(runner, "run", _fake_obabel())
    monkeypatch.setattr(runner.ADMETFilter, "parse_sdf_properties", _custom_props(150.0, 1.0))
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    with pytest.raises(RuntimeError, match="rejected by the property filter"):
        lm.create_custom_library(mw_min=200, mw_max=500, logp_max=5)


def test_custom_without_lipinski_skips_hbd_hba_only(monkeypatch, tmp_path):
    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)
    monkeypatch.setattr(runner, "run", _fake_obabel())
    monkeypatch.setattr(runner.ADMETFilter, "parse_sdf_properties",
                        _custom_props(300.0, 1.0, hbd=9, hba=20))
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    with pytest.raises(RuntimeError):
        lm.create_custom_library(apply_lipinski=True)
    assert len(lm.create_custom_library(apply_lipinski=False, max_compounds=1)) == 1


def test_pubchem_parses_properties_once(monkeypatch, tmp_path):
    monkeypatch.setattr(runner, "_http_get_bytes", _fake_http)
    monkeypatch.setattr(runner, "run", _fake_obabel())
    calls = []

    def props(f):
        calls.append(f)
        return {"mw": 300.0, "logp": 2.0, "hbd": 1.0, "hba": 2.0}

    monkeypatch.setattr(runner.ADMETFilter, "parse_sdf_properties", staticmethod(props))
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "in"))
    lm.create_fda_library(max_compounds=1)
    assert len(calls) == 1


# ---------------------------------------------------------------------------
# Local ligands: collection, names, failures
# ---------------------------------------------------------------------------

def test_prepare_local_pdbqt_direct(tmp_path):
    indir = tmp_path / "in"
    indir.mkdir()
    lig = indir / "mol1.pdbqt"
    lig.write_text(_LIGAND_PDBQT)
    lm = library.LibraryManager(str(tmp_path), str(indir))
    assert lm.prepare_local(apply_lipinski=False) == [str(lig)]
    # The deprecated keyword still works.
    assert lm.prepare_local(apply_admet=False) == [str(lig)]


def test_prepare_local_missing_dir(tmp_path):
    lm = library.LibraryManager(str(tmp_path), str(tmp_path / "nope"))
    with pytest.raises(FileNotFoundError):
        lm.prepare_local()


def test_prepare_local_no_ligands(tmp_path):
    indir = tmp_path / "in"
    indir.mkdir()
    lm = library.LibraryManager(str(tmp_path), str(indir))
    with pytest.raises(FileNotFoundError, match="No ligands found"):
        lm.prepare_local()


def test_prepare_local_rejects_receptor_style_pdbqt(tmp_path):
    indir = tmp_path / "in"
    indir.mkdir()
    good = indir / "good.pdbqt"
    good.write_text(_LIGAND_PDBQT)
    (indir / "c1.pdbqt").write_text(_ATOM)
    lm = library.LibraryManager(str(tmp_path), str(indir))
    assert lm.prepare_local(apply_lipinski=False) == [str(good)]
    assert "ROOT" in library._ligand_pdbqt_problem(str(indir / "c1.pdbqt"))
    assert library._ligand_pdbqt_problem(str(good)) is None


def test_prepare_local_mixed_formats(tmp_path, monkeypatch):
    indir = tmp_path / "in"
    indir.mkdir()
    lig = indir / "lib1.pdbqt"
    lig.write_text(_LIGAND_PDBQT)
    (indir / "control.pdb").write_text(
        "HETATM    1  C   UNL     1       0.000   0.000   0.000  1.00  0.00           C\n")
    monkeypatch.setattr(runner, "run", _fake_obabel())
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    out = lm.prepare_local(apply_lipinski=False)
    assert sorted(os.path.basename(p) for p in out) == ["control.pdbqt", "lib1.pdbqt"]
    assert lm.metadata["control"]["input_format"] == "pdb"
    assert lm.metadata["control"]["minimized"] is False


def test_same_name_other_format_is_superseded_not_failed(tmp_path, monkeypatch, caplog):
    """DUD-E folders hold X.mol2 next to X.pdbqt: the PDBQT is docked and the
    MOL2 is reported as superseded, outside the success denominator."""
    indir = tmp_path / "in"
    indir.mkdir()
    for name in ("a", "b"):
        (indir / f"{name}.pdbqt").write_text(_LIGAND_PDBQT)
        (indir / f"{name}.mol2").write_text("@<TRIPOS>MOLECULE\nx\n")
    monkeypatch.setattr(runner, "run", _fake_obabel())
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    with caplog.at_level(logging.INFO):
        out = lm.prepare_local(apply_lipinski=False)
    assert out == [str(indir / "a.pdbqt"), str(indir / "b.pdbqt")]
    assert lm.failed == []
    assert [os.path.basename(p) for p, _ in lm.superseded] == ["a.mol2", "b.mol2"]
    assert "Success: 2/2" in caplog.text
    assert "superseded (same name as a.pdbqt)" in caplog.text


def test_name_is_reserved_on_first_sight(tmp_path, monkeypatch):
    """lig.sdf rejected by Lipinski must not let lig.mol2 be docked unfiltered."""
    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "lig.sdf").write_text(_molblock(0.8))
    (indir / "lig.mol2").write_text("@<TRIPOS>MOLECULE\nx\n")
    (indir / "ok.pdbqt").write_text(_LIGAND_PDBQT)
    monkeypatch.setattr(runner, "run", _fake_obabel())
    monkeypatch.setattr(runner.ADMETFilter, "parse_sdf_properties",
                        _custom_props(900.0, 9.0))
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    assert lm.prepare_local(apply_lipinski=True) == [str(indir / "ok.pdbqt")]
    assert [os.path.basename(p) for p, _ in lm.failed] == ["lig.sdf"]
    assert "Lipinski" in lm.failed[0][1]
    assert [os.path.basename(p) for p, _ in lm.superseded] == ["lig.mol2"]


def test_invalid_file_still_reserves_its_name(tmp_path, monkeypatch):
    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "x.sdf").write_text("not an sdf\n")
    (indir / "x.pdbqt").write_text(_LIGAND_PDBQT)
    monkeypatch.setattr(runner, "run", _fake_obabel())
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    with pytest.raises(RuntimeError, match="No valid ligands"):
        lm.prepare_local(apply_lipinski=True)


def test_same_name_same_format_in_subdirs_is_a_failure(tmp_path):
    indir = tmp_path / "in"
    (indir / "a").mkdir(parents=True)
    (indir / "b").mkdir()
    (indir / "a" / "x.pdbqt").write_text(_LIGAND_PDBQT)
    (indir / "b" / "x.pdbqt").write_text(_LIGAND_PDBQT)
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    assert lm.prepare_local(apply_lipinski=False) == [str(indir / "a" / "x.pdbqt")]
    assert "Duplicate ligand name" in lm.failed[0][1]


def test_output_dir_inside_ligand_dir_is_skipped(tmp_path, caplog):
    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "x.pdbqt").write_text(_LIGAND_PDBQT)
    outdir = indir / "results"
    (outdir / "docked").mkdir(parents=True)
    (outdir / "docked" / "y.pdbqt").write_text(_LIGAND_PDBQT)  # a previous run
    lm = library.LibraryManager(str(outdir), str(indir))
    with caplog.at_level(logging.WARNING):
        assert lm.prepare_local(apply_lipinski=False) == [str(indir / "x.pdbqt")]
    assert "inside the ligand directory" in caplog.text


def test_output_dir_equal_to_ligand_dir_is_rejected(tmp_path):
    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "x.pdbqt").write_text(_LIGAND_PDBQT)
    lm = library.LibraryManager(str(indir), str(indir))
    with pytest.raises(ValueError, match="is the ligand directory"):
        lm.prepare_local(apply_lipinski=False)


def test_multi_model_pdbqt_rejected(tmp_path):
    p = tmp_path / "vina_out.pdbqt"
    p.write_text("MODEL 1\n" + _LIGAND_PDBQT + "ENDMDL\nMODEL 2\n"
                 + _LIGAND_PDBQT + "ENDMDL\n")
    assert "several models" in library._ligand_pdbqt_problem(str(p))


def test_multi_molecule_sdf_rejected(tmp_path, monkeypatch):
    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "two.sdf").write_text(_molblock(0.8) * 2)
    (indir / "ok.pdbqt").write_text(_LIGAND_PDBQT)
    monkeypatch.setattr(runner, "run", _fake_obabel())
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    assert lm.prepare_local(apply_lipinski=False) == [str(indir / "ok.pdbqt")]


def test_local_sdf_properties_parsed_once(tmp_path, monkeypatch):
    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "a.sdf").write_text(_molblock(0.8))
    calls = []

    def props(f):
        calls.append(f)
        return {"mw": 300.0, "logp": 2.0, "hbd": 1.0, "hba": 2.0}

    monkeypatch.setattr(runner.ADMETFilter, "parse_sdf_properties", staticmethod(props))
    monkeypatch.setattr(runner, "run", _fake_obabel())
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    lm.prepare_local(apply_lipinski=True)
    assert len(calls) == 1


def test_missing_descriptors_are_none_not_zero(tmp_path, monkeypatch):
    """Missing descriptors are stored as None (blank in ranking.csv), never 0."""
    monkeypatch.setattr(runner, "run", _fake_obabel())
    monkeypatch.setattr(runner.ADMETFilter, "parse_sdf_properties", staticmethod(lambda f: None))
    sdf = tmp_path / "a.sdf"
    sdf.write_text(_molblock(0.8))
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path))
    lm._prepare_sdf_to_pdbqt(str(sdf), str(tmp_path / "a.pdbqt"), "a")
    entry = lm.metadata["a"]
    assert all(entry[k] is None for k in ("mw", "logp", "tpsa", "rotors", "hbd", "hba"))
    assert entry["properties_source"] == "unavailable"
    lm._save_metadata()

    analyzer = runner.ResultsAnalyzer(str(tmp_path / "out"))
    analyzer.save_ranking([("a", -6.5)], mode="extended",
                          metrics_dict={"a": {"status": "OK"}})
    import csv
    row = next(csv.DictReader(open(analyzer.csv_file)))
    assert row["MW"] == "" and row["LogP"] == "" and row["HBA"] == ""


# ---------------------------------------------------------------------------
# SDF preparation: 3D detection, minimization, Open Babel calls
# ---------------------------------------------------------------------------

def test_sdf_has_3d_coords(tmp_path):
    flat = tmp_path / "flat.sdf"
    flat.write_text(_molblock(0.0))
    solid = tmp_path / "solid.sdf"
    solid.write_text(_molblock(0.8))
    assert library._sdf_has_3d_coords(str(solid))
    assert not library._sdf_has_3d_coords(str(flat))
    assert not library._sdf_has_3d_coords(str(tmp_path / "missing.sdf"))


def test_sdf_has_3d_coords_reads_v3000():
    assert library._sdf_has_3d_coords(os.path.join(DATA, "acetaminophen_v3000_3d.sdf"))
    assert not library._sdf_has_3d_coords(os.path.join(DATA, "acetaminophen_v3000_2d.sdf"))
    assert not library._sdf_has_3d_coords(os.path.join(DATA, "acetaminophen_2d.sdf"))


def test_v3000_3d_input_is_not_regenerated(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(runner, "run", _fake_obabel(calls))
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path), minimize=False)
    lm._prepare_sdf_to_pdbqt(os.path.join(DATA, "acetaminophen_v3000_3d.sdf"),
                             str(tmp_path / "a.pdbqt"), "a")
    assert not any("--gen3d" in cmd for cmd, _ in calls)


def test_minimize_keeps_3d_input_and_runs_single_threaded(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(runner, "run", _fake_obabel(calls))
    monkeypatch.setattr(runner, "_ensure_pdbqt_has_charges", lambda p: True)
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path))
    for name, z, expect_gen3d in (("solid", 0.8, False), ("flat", 0.0, True)):
        sdf = tmp_path / f"{name}.sdf"
        sdf.write_text(_molblock(z))
        calls.clear()
        lm._prepare_sdf_to_pdbqt(str(sdf), str(tmp_path / f"{name}.pdbqt"), name)
        obabel_calls = [(c, env) for c, env in calls if c[0] == runner.OBABEL]
        minimize = [c for c, _ in obabel_calls if "--minimize" in c]
        assert len(minimize) == 1
        # 3D coordinates of a 2D input are generated first (separate call);
        # the minimization then starts from that 3D file.
        assert "--gen3d" not in minimize[0]
        gen3d = [c for c, _ in obabel_calls if "--gen3d" in c]
        assert bool(gen3d) == expect_gen3d
        if expect_gen3d:
            assert minimize[0][minimize[0].index("-isdf") + 1].endswith(f"{name}_gen3d.sdf")
        # Every Open Babel call (descriptors, minimization, conversion) is
        # single-threaded.
        assert obabel_calls and all(env == {"OMP_NUM_THREADS": "1"} for _, env in obabel_calls)


def test_sdf_conversion_requests_gasteiger_charges(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(runner, "run", _fake_obabel(calls))
    sdf = tmp_path / "a.sdf"
    sdf.write_text(_molblock(0.8))
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path), minimize=False)
    lm._prepare_sdf_to_pdbqt(str(sdf), str(tmp_path / "a.pdbqt"), "a")
    # The conversion runs twice (deterministic hydrogen check), identically.
    converts = [c for c, _ in calls if "-opdbqt" in c]
    assert len(converts) == 2
    assert all(c[-3:] == ["-h", "--partialcharge", "gasteiger"] for c in converts)


def test_failed_minimization_falls_back_to_input(tmp_path, monkeypatch, caplog):
    """Open Babel exits 0 but writes nothing when MMFF94 cannot be set up."""
    calls = []

    def fake_run(cmd, capture=False, extra_env=None, **kw):
        calls.append(cmd)
        if "-O" in cmd:
            out = cmd[cmd.index("-O") + 1]
            if "--minimize" in cmd:
                open(out, "w").close()
            else:
                with open(out, "w") as f:
                    f.write(_LIGAND_PDBQT_3D)
        return ("", "", 0) if capture else None

    monkeypatch.setattr(runner, "run", fake_run)
    sdf = tmp_path / "boron.sdf"
    sdf.write_text(_molblock(0.8))
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path))
    with caplog.at_level(logging.WARNING):
        lm._prepare_sdf_to_pdbqt(str(sdf), str(tmp_path / "boron.pdbqt"), "boron")
    convert = [c for c in calls if "-opdbqt" in c][0]
    assert convert[convert.index("-isdf") + 1] == str(sdf)
    assert lm.metadata["boron"]["minimized"] is False
    assert "MMFF94" in lm.metadata["boron"]["minimization_error"]
    assert "MMFF94 minimization failed" in caplog.text
    assert not os.path.exists(os.path.join(lm.minimized_dir, "boron.sdf"))


def test_minimization_timeout_is_a_logged_fallback(tmp_path, monkeypatch, caplog):
    def fake_run(cmd, capture=False, extra_env=None, **kw):
        if "--minimize" in cmd:
            raise subprocess.TimeoutExpired(cmd, 5)
        return _fake_obabel()(cmd, capture=capture, extra_env=extra_env)

    monkeypatch.setattr(runner, "run", fake_run)
    sdf = tmp_path / "slow.sdf"
    sdf.write_text(_molblock(0.8))
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path))
    with caplog.at_level(logging.WARNING):
        lm._prepare_sdf_to_pdbqt(str(sdf), str(tmp_path / "slow.pdbqt"), "slow")
    assert "timed out" in lm.metadata["slow"]["minimization_error"]
    assert "timed out" in caplog.text


def test_unminimized_2d_sdf_gets_3d_intermediate(tmp_path, monkeypatch):
    """A 2D SDF with --no-minimize: --gen3d writes a 3D SDF first, and that
    3D file (never the 2D input) is converted and used for charges."""
    calls = []
    monkeypatch.setattr(runner, "run", _fake_obabel(calls))
    fixes = []
    monkeypatch.setattr(runner, "_ensure_pdbqt_has_charges", lambda p: False)
    monkeypatch.setattr(runner, "_fix_pdbqt_charges", lambda p, src: fixes.append(src))
    sdf = tmp_path / "flat.sdf"
    sdf.write_text(_molblock(0.0))
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path), minimize=False)
    lm._prepare_sdf_to_pdbqt(str(sdf), str(tmp_path / "flat.pdbqt"), "flat")
    obabel_cmds = [c for c, _ in calls if "--append" not in c]
    gen3d, convert, probe = obabel_cmds  # conversion + its determinism probe
    assert "--gen3d" in gen3d and gen3d[gen3d.index("-O") + 1].endswith("flat_gen3d.sdf")
    src = convert[convert.index("-isdf") + 1]
    assert src.endswith("flat_gen3d.sdf") and "-h" in convert
    assert fixes == [src]
    assert lm.metadata["flat"]["coordinates"] == "gen3d"
    assert lm.metadata["flat"]["coordinates_method"] == "obabel-gen3d"
    assert "h-probe" in probe[probe.index("-O") + 1]


def test_failed_minimization_of_2d_input_still_gets_3d(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(runner, "run", _fake_obabel(calls, fail_on="--minimize"))
    sdf = tmp_path / "flat.sdf"
    sdf.write_text(_molblock(0.0))
    lm = library.LibraryManager(str(tmp_path / "out"), str(tmp_path))
    lm._prepare_sdf_to_pdbqt(str(sdf), str(tmp_path / "flat.pdbqt"), "flat")
    convert = [c for c, _ in calls if "-opdbqt" in c][0]
    assert convert[convert.index("-isdf") + 1].endswith("flat_gen3d.sdf")


def test_planar_result_of_2d_input_is_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "run", _fake_obabel(pdbqt=_LIGAND_PDBQT_FLAT))
    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "flat.sdf").write_text(_molblock(0.0))
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir), minimize=False)
    with pytest.raises(RuntimeError, match="No valid ligands"):
        lm.prepare_local(apply_lipinski=False)
    assert "planar" in lm.failed[0][1]
    assert not os.path.exists(os.path.join(lm.lib_dir, "flat.pdbqt"))


def test_failed_conversion_removes_partial_output(tmp_path, monkeypatch):
    def fake_run(cmd, capture=False, extra_env=None, **kw):
        out = cmd[cmd.index("-O") + 1]
        with open(out, "w") as f:
            f.write("partial")
        raise subprocess.CalledProcessError(1, cmd, output="", stderr="kekulize error")

    monkeypatch.setattr(runner, "run", fake_run)
    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "bad.mol2").write_text("@<TRIPOS>MOLECULE\nx\n")
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    with pytest.raises(RuntimeError):
        lm.prepare_local(apply_lipinski=False)
    assert "kekulize error" in lm.failed[0][1]
    assert os.listdir(lm.lib_dir) == []


def test_every_obabel_call_is_single_threaded(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(runner, "run", _fake_obabel(calls))
    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "a.sdf").write_text(_molblock(0.0))
    (indir / "b.mol2").write_text("@<TRIPOS>MOLECULE\nx\n")
    (indir / "c.pdb").write_text("HETATM    1  C   UNL     1       0.000   0.000   0.000  1.00  0.00           C\n")
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    lm.prepare_local(apply_lipinski=False)
    obabel_calls = [(c, env) for c, env in calls if c[0] == runner.OBABEL]
    # SDF: descriptors, gen3d, MMFF94, conversion; MOL2: conversion; PDB:
    # protonation and conversion. Conversions that add hydrogens run twice
    # (run_obabel_deterministic), so 3 + 2 + 2 + 2 + 2 = 11.
    assert len(obabel_calls) == 11
    assert all(env == {"OMP_NUM_THREADS": "1"} for _, env in obabel_calls)


def test_minimization_scope_is_logged_once(tmp_path, monkeypatch, caplog):
    indir = tmp_path / "in"
    indir.mkdir()
    (indir / "a.pdbqt").write_text(_LIGAND_PDBQT)
    (indir / "b.pdbqt").write_text(_LIGAND_PDBQT)
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    with caplog.at_level(logging.INFO):
        lm.prepare_local(apply_lipinski=False)
    assert caplog.text.count("minimization applies to SDF inputs only") == 1


# ---------------------------------------------------------------------------
# Real Open Babel
# ---------------------------------------------------------------------------

def _spread(path):
    return library._pdbqt_coordinate_spread(path)


@pytest.mark.slow
@needs_obabel
def test_minimization_is_repeatable_with_real_obabel(tmp_path):
    sdf = tmp_path / "aspirin.sdf"
    subprocess.run(["obabel", "-:CC(=O)Oc1ccccc1C(=O)O", "--gen3d", "-osdf", "-O", str(sdf)],
                   check=True, capture_output=True, env={**os.environ, "OMP_NUM_THREADS": "1"})
    outputs = []
    for i in range(3):
        lm = library.LibraryManager(str(tmp_path / f"run{i}"), str(tmp_path))
        pdbqt = tmp_path / f"run{i}.pdbqt"
        lm._prepare_sdf_to_pdbqt(str(sdf), str(pdbqt), "aspirin")
        outputs.append(pdbqt.read_text())
    assert outputs[0] == outputs[1] == outputs[2]


@needs_obabel
@pytest.mark.parametrize("minimize", [False, True])
def test_2d_sdf_is_not_docked_flat_with_real_obabel(tmp_path, minimize):
    """Regression: a 2D SDF prepared with --no-minimize was docked flat
    (every z = 0.000) because the charge rebuild re-read the 2D input."""
    indir = tmp_path / "in"
    indir.mkdir()
    shutil.copy(os.path.join(DATA, "acetaminophen_2d.sdf"), indir)
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir), minimize=minimize,
                                minimize_steps=50)
    (out,) = lm.prepare_local(apply_lipinski=False)
    assert min(_spread(out)) > 0.5
    assert runner._ensure_pdbqt_has_charges(out)


@needs_obabel
def test_prepared_ligands_have_3d_extent_with_real_obabel(tmp_path):
    """Every prepared ligand (SDF 2D/3D/V3000, MOL2, PDB) spreads along x, y and z."""
    indir = tmp_path / "in"
    indir.mkdir()
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    smiles = "CC(=O)NCCO"
    for fmt in ("mol2", "pdb"):
        subprocess.run(["obabel", f"-:{smiles}", "--gen3d", "-o" + fmt,
                        "-O", str(indir / f"m_{fmt}.{fmt}")], check=True,
                       capture_output=True, env=env)
    shutil.copy(os.path.join(DATA, "acetaminophen_2d.sdf"), indir)
    shutil.copy(os.path.join(DATA, "acetaminophen_v3000_3d.sdf"), indir)
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir), minimize=False)
    out = lm.prepare_local(apply_lipinski=True)
    assert len(out) == 4, lm.failed
    for path in out:
        assert min(_spread(path)) > 0.5, path


@pytest.mark.slow
@needs_obabel
def test_pdb_ligand_gets_charged_polar_hydrogens_with_real_obabel(tmp_path):
    """A heavy-atom-only PDB ligand is docked with its polar hydrogens, which
    carry Gasteiger charges (hydrogens are added before the charges)."""
    indir = tmp_path / "in"
    indir.mkdir()
    pdb = indir / "acetamide.pdb"
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    subprocess.run(["obabel", "-:CC(=O)NCCO", "--gen3d", "-d", "-opdb",
                    "-O", str(pdb)], check=True, capture_output=True, env=env)
    lm = library.LibraryManager(str(tmp_path / "out"), str(indir))
    (out,) = lm.prepare_local(apply_lipinski=False)
    atoms = [line for line in open(out) if line.startswith(("ATOM", "HETATM"))]
    hd = [line for line in atoms if line[77:79].strip() == "HD"]
    assert len(hd) == 2  # amide N-H and hydroxyl O-H
    assert all(float(line[70:76]) > 0.1 for line in hd)

    # Coordinates, atom types, torsion tree and name record are those of the
    # former one-step conversion; only the charge column differs.
    single = tmp_path / "single.pdbqt"
    subprocess.run(["obabel", "-ipdb", str(pdb), "-opdbqt", "-O", str(single), "-h",
                    "--partialcharge", "gasteiger"], check=True, capture_output=True, env=env)

    def mask(path):
        return [line[:70] + line[76:] if line.startswith(("ATOM", "HETATM")) else line
                for line in open(path)]

    assert mask(out) == mask(str(single))


def test_metadata_paths_inside_output_are_relative(tmp_path):
    from vspipeline.library import LibraryManager
    lm = LibraryManager(str(tmp_path / "out"), str(tmp_path))
    inside = tmp_path / "out" / "ligands" / "a.pdbqt"
    outside = tmp_path / "given.pdbqt"
    assert lm._out_relpath(str(inside)) == os.path.join("ligands", "a.pdbqt")
    assert lm._out_relpath(str(outside)) == str(outside)
