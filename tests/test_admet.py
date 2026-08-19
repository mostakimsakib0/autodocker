import os
import shutil

import pytest

import runner


class TestADMETFilter:
    def test_check_lipinski_ok(self):
        props = {"mw": 350.0, "logp": 3.0, "hbd": 2, "hba": 4}
        ok, violations = runner.ADMETFilter.check_lipinski(props)
        assert ok is True
        assert violations == []

    def test_check_lipinski_violations(self):
        props = {"mw": 600.0, "logp": 6.0, "hbd": 6, "hba": 12}
        ok, violations = runner.ADMETFilter.check_lipinski(props)
        assert ok is False
        assert len(violations) == 4

    def test_check_lipinski_none(self):
        ok, violations = runner.ADMETFilter.check_lipinski(None)
        assert ok is False

    def test_check_lipinski_missing_keys(self):
        ok, violations = runner.ADMETFilter.check_lipinski({"mw": 300.0})
        assert ok is False

    def test_parse_sdf_properties_requires_valid_sdf(self, tmp_path):
        p = tmp_path / "bad.sdf"
        p.write_text("garbage\n")
        assert runner.ADMETFilter.parse_sdf_properties(str(p)) is None


@pytest.fixture
def ligands_dir(tmp_path, small_sdf):
    d = tmp_path / "ligands"
    d.mkdir()
    shutil.copy(small_sdf, str(d / "methane.sdf"))
    return str(d)


class TestNoAdmetNeverDrops:
    """Regression for ROADMAP A7: --no-admet must not drop ligands
    even when descriptor parsing fails."""

    def test_apply_admet_false_keeps_ligand(self, ligands_dir, tmp_path, monkeypatch):
        monkeypatch.setattr(runner.ADMETFilter, "parse_sdf_properties",
                            staticmethod(lambda sdf: None))
        monkeypatch.setattr(runner.LibraryManager, "_prepare_sdf_to_pdbqt",
                            lambda self, sdf, pdbqt, lid: pdbqt)
        lm = runner.LibraryManager(str(tmp_path / "out"), ligands_dir)
        out = lm._prepare_local_sdf(apply_admet=False)
        assert len(out) == 1

    def test_apply_admet_true_drops_when_descriptors_fail(self, ligands_dir, tmp_path, monkeypatch):
        monkeypatch.setattr(runner.ADMETFilter, "parse_sdf_properties",
                            staticmethod(lambda sdf: None))
        monkeypatch.setattr(runner.LibraryManager, "_prepare_sdf_to_pdbqt",
                            lambda self, sdf, pdbqt, lid: pdbqt)
        lm = runner.LibraryManager(str(tmp_path / "out"), ligands_dir)
        with pytest.raises(RuntimeError, match="No valid ligands prepared"):
            lm._prepare_local_sdf(apply_admet=True)

    def test_apply_admet_true_drops_violating(self, ligands_dir, tmp_path, monkeypatch):
        monkeypatch.setattr(
            runner.ADMETFilter, "parse_sdf_properties",
            staticmethod(lambda sdf: {"mw": 700.0, "logp": 7.0,
                                      "hbd": 6, "hba": 12}))
        monkeypatch.setattr(runner.LibraryManager, "_prepare_sdf_to_pdbqt",
                            lambda self, sdf, pdbqt, lid: pdbqt)
        lm = runner.LibraryManager(str(tmp_path / "out"), ligands_dir)
        with pytest.raises(RuntimeError, match="No valid ligands prepared"):
            lm._prepare_local_sdf(apply_admet=True)

    def test_apply_admet_false_keeps_violating(self, ligands_dir, tmp_path, monkeypatch):
        monkeypatch.setattr(
            runner.ADMETFilter, "parse_sdf_properties",
            staticmethod(lambda sdf: {"mw": 700.0, "logp": 7.0,
                                      "hbd": 6, "hba": 12}))
        monkeypatch.setattr(runner.LibraryManager, "_prepare_sdf_to_pdbqt",
                            lambda self, sdf, pdbqt, lid: pdbqt)
        lm = runner.LibraryManager(str(tmp_path / "out"), ligands_dir)
        out = lm._prepare_local_sdf(apply_admet=False)
        assert len(out) == 1