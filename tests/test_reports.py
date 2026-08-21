import csv
import os

import pytest

import runner


def _write(path, content):
    with open(path, "w") as f:
        f.write(content)
    return path


def _pose_pdbqt(ligand, serial=5, x=5.0, y=5.0, z=5.0):
    return f"""MODEL 1
ROOT
ATOM  {serial:5d}  C1  LIG L   1       {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00      C+0
ENDROOT
ENDMDL
"""


def _receptor_pdbqt():
    return (
        "REMARK synthetic receptor\n"
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00      N+0\n"
        "END\n"
    )


def _ranking_csv(outdir, rows):
    path = os.path.join(outdir, "ranking.csv")
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Ligand", "Binding_Affinity", "SimScore", "Binding_Modes",
            "Mode_Coverage", "MW", "LogP", "HBA", "HBD", "TPSA",
            "Rotors", "ZINC_URL", "Docked_File", "Status"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path


class TestRankingCsv:
    def test_columns(self, outdir):
        analyzer = runner.ResultsAnalyzer(outdir)
        results = [("ligA", -7.5), ("ligB", None)]
        analyzer.save_ranking(results, mode="extended",
                              metrics_dict={"ligA": {"status": "OK"},
                                            "ligB": {}})
        with open(analyzer.csv_file) as f:
            reader = csv.DictReader(f)
            assert reader.fieldnames[0] == "Ligand"
            assert "Binding_Affinity" in reader.fieldnames
            rows = list(reader)
        by_name = {r["Ligand"]: r for r in rows}
        assert by_name["ligA"]["Binding_Affinity"] == "-7.5"
        assert by_name["ligB"]["Binding_Affinity"] == "FAILED"


class TestHtmlReport:
    def test_html_generated_on_fresh_run(self, tmp_path):
        """Regression for H5: HTML report must be generated even when the
        output directory is brand-new (i.e. ranking.csv exists by then)."""
        outdir = str(tmp_path / "fresh")
        os.makedirs(outdir)
        csv_path = _ranking_csv(outdir, [
            {"Ligand": "ligA", "Binding_Affinity": "-7.5",
             "SimScore": "0.9", "Binding_Modes": "9", "Status": "OK"},
        ])
        html_path = os.path.join(outdir, "results_report.html")
        runner.generate_html_report(csv_path, os.path.join(outdir, "docked"),
                                    html_path)
        assert os.path.exists(html_path)
        with open(html_path) as f:
            content = f.read()
        assert "Virtual Screening Results Report" in content
        assert "ligA" in content

    def test_html_missing_csv_warns(self, outdir, caplog):
        runner.generate_html_report(
            os.path.join(outdir, "no.csv"), os.path.join(outdir, "docked"),
            os.path.join(outdir, "results_report.html"))
        assert not os.path.exists(os.path.join(outdir, "results_report.html"))

    def test_html_escapes_ligand_names(self, tmp_path):
        """Regression for H6: ligand names must be HTML-escaped."""
        outdir = str(tmp_path / "escape")
        os.makedirs(outdir)
        malicious = '<script>alert("x")</script>'
        csv_path = _ranking_csv(outdir, [
            {"Ligand": malicious, "Binding_Affinity": "-7.5",
             "SimScore": "0.9", "Binding_Modes": "9", "Status": "OK"},
        ])
        html_path = os.path.join(outdir, "results_report.html")
        runner.generate_html_report(csv_path, os.path.join(outdir, "docked"),
                                    html_path)
        with open(html_path) as f:
            content = f.read()
        assert "<script>alert" not in content
        assert "&lt;script&gt;alert" in content


class TestNglViewer:
    """The interactive NGL 3D viewer embedded in the HTML report.

    The viewer is fully self-contained: the NGL library and the receptor/pose
    PDB texts are inlined into the HTML, so the report opens from ``file://``
    with a double-click (no ``viewer/`` directory, no web server, no CDN).
    """

    @staticmethod
    def _prepare(tmp_path, monkeypatch, ngl_bundle=b"/* fake ngl */"):
        monkeypatch.setattr(runner, "_fetch_ngljs", lambda: ngl_bundle)
        outdir = str(tmp_path / "out")
        docked = os.path.join(outdir, "docked")
        os.makedirs(docked, exist_ok=True)
        _write(os.path.join(outdir, "receptor_b.pdbqt"), _receptor_pdbqt())
        _write(os.path.join(docked, "ligA_out.pdbqt"), _pose_pdbqt("ligA"))
        _write(os.path.join(docked, "ligB_out.pdbqt"), _pose_pdbqt("ligB",
                                                                    serial=7,
                                                                    x=7.0, y=7.0, z=7.0))
        csv_path = _ranking_csv(outdir, [
            {"Ligand": "ligA", "Binding_Affinity": "-7.5",
             "SimScore": "0.9", "Binding_Modes": "9", "Status": "OK"},
            {"Ligand": "ligB", "Binding_Affinity": "-6.2",
             "SimScore": "0.8", "Binding_Modes": "7", "Status": "OK"},
        ])
        html = os.path.join(outdir, "results_report.html")
        runner.generate_html_report(
            csv_path, docked, html,
            receptor_pdbqt=os.path.join(outdir, "receptor_b.pdbqt"))
        with open(html) as f:
            content = f.read()
        return outdir, content

    def test_viewer_section_present(self, tmp_path, monkeypatch):
        outdir, content = self._prepare(tmp_path, monkeypatch)
        assert "3D Structure Viewer" in content
        assert "viewer-container" in content
        assert "ligA" in content
        assert "ligB" in content
        # Data is inlined, never referenced as external files:
        assert "viewer/ligA_pose.pdb" not in content
        assert runner.NGL_CDN_URL not in content
        # Inlined PDB text carries the receptor/pose atoms:
        assert "LIG L" in content
        assert "ALA A" in content
        # The NGL library itself is inlined into the report:
        assert "/* fake ngl */" in content

    def test_no_viewer_files_written(self, tmp_path, monkeypatch):
        outdir, content = self._prepare(tmp_path, monkeypatch)
        assert not os.path.exists(os.path.join(outdir, "viewer"))
        assert "viewer/" not in content

    def test_viewer_no_docked_notice(self, tmp_path, monkeypatch):
        """With no docked poses the viewer is replaced by a clear notice."""
        monkeypatch.setattr(runner, "_fetch_ngljs", lambda: b"/* fake ngl */")
        outdir = str(tmp_path / "out")
        os.makedirs(outdir, exist_ok=True)
        csv_path = _ranking_csv(outdir, [
            {"Ligand": "ligA", "Binding_Affinity": "-7.5",
             "SimScore": "0.9", "Binding_Modes": "9", "Status": "OK"},
        ])
        html = os.path.join(outdir, "results_report.html")
        runner.generate_html_report(
            csv_path, os.path.join(outdir, "docked"), html, receptor_pdbqt=None)
        with open(html) as f:
            content = f.read()
        assert "No docked structures were available" in content
        assert "viewer-container" not in content

    def test_ngl_inlined_not_cdn(self, tmp_path, monkeypatch):
        outdir, content = self._prepare(tmp_path, monkeypatch)
        assert 'src="viewer/ngl.min.js"' not in content
        assert runner.NGL_CDN_URL not in content
        assert "/* fake ngl */" in content

    def test_ngl_unavailable_notice(self, tmp_path, monkeypatch):
        """If the NGL library cannot be bundled, show a notice (no CDN)."""
        monkeypatch.setattr(runner, "_fetch_ngljs", lambda: None)
        outdir = str(tmp_path / "out")
        docked = os.path.join(outdir, "docked")
        os.makedirs(docked, exist_ok=True)
        _write(os.path.join(docked, "ligA_out.pdbqt"), _pose_pdbqt("ligA"))
        csv_path = _ranking_csv(outdir, [
            {"Ligand": "ligA", "Binding_Affinity": "-7.5",
             "SimScore": "0.9", "Binding_Modes": "9", "Status": "OK"},
        ])
        html = os.path.join(outdir, "results_report.html")
        runner.generate_html_report(
            csv_path, docked, html, receptor_pdbqt=None)
        with open(html) as f:
            content = f.read()
        assert "3D viewer unavailable" in content
        assert runner.NGL_CDN_URL not in content
        assert "viewer-container" not in content