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

def test_cluster_poses_reads_vina_output_and_skips_other_ligands(tmp_path):
    from vspipeline import report
    docked = tmp_path / "docked"
    docked.mkdir()

    def model(i, dx, aff):
        return (f"MODEL {i}\nREMARK VINA RESULT:    {aff}      0.000      0.000\n"
                + "".join(f"ATOM  {k:5d}  C   UNL     1    {k + dx:8.3f}   0.000   0.000"
                          "  0.00  0.00    +0.000 C \n" for k in range(1, 4))
                + "BEGIN_RES SER A  45\nATOM      9  CB  SER A  45      50.000  50.000  50.000"
                  "  0.00  0.00    +0.000 C \nEND_RES SER A  45\nENDMDL\n")

    (docked / "lig_out.pdbqt").write_text(model(1, 0.0, -8.0) + model(2, 0.5, -7.5)
                                          + model(3, 6.0, -7.0))
    (docked / "old_out.pdbqt").write_text(model(1, 0.0, -6.0))
    report._cluster_cache.clear()
    clusters = report.cluster_poses(str(docked), rmsd_threshold=2.0, ligands=["lig"])
    assert list(clusters) == ["lig"]
    # Models 1 and 2 are 0.5 A apart (flex atoms ignored); model 3 is separate.
    assert [len(c["poses"]) for c in clusters["lig"]] == [2, 1]


def test_cluster_rmsd_uses_heavy_atoms_only(tmp_path):
    from vspipeline import report
    docked = tmp_path / "docked"
    docked.mkdir()

    def model(i, h_x):
        return (f"MODEL {i}\nREMARK VINA RESULT:    -7.0      0.000      0.000\n"
                "ATOM      1  C1  UNL     1       0.000   0.000   0.000  0.00  0.00    +0.000 C \n"
                "ATOM      2  C2  UNL     1       1.500   0.000   0.000  0.00  0.00    +0.000 C \n"
                f"ATOM      3  H1  UNL     1    {h_x:8.3f}   0.000   0.000  0.00  0.00    +0.000 HD\n"
                "ENDMDL\n")

    # Only the polar hydrogen moves (by 9 A): same heavy-atom pose.
    (docked / "lig_out.pdbqt").write_text(model(1, -1.0) + model(2, 8.0))
    report._cluster_cache.clear()
    clusters = report.cluster_poses(str(docked), rmsd_threshold=0.5, ligands=["lig"])
    assert [len(c["poses"]) for c in clusters["lig"]] == [2]


class TestReportContent:
    @staticmethod
    def _render(tmp_path, monkeypatch, rows, ngl=b"/* fake ngl */", poses=None):
        monkeypatch.setattr(runner, "_fetch_ngljs", lambda: ngl)
        outdir = str(tmp_path / "out")
        docked = os.path.join(outdir, "docked")
        os.makedirs(docked, exist_ok=True)
        for name, text in (poses or {}).items():
            _write(os.path.join(docked, name), text)
        csv_path = _ranking_csv(outdir, rows)
        html_path = os.path.join(outdir, "results_report.html")
        assert runner.generate_html_report(csv_path, docked, html_path)
        with open(html_path) as f:
            return f.read()

    ROWS = [
        {"Ligand": "ligA", "Binding_Affinity": "-9.1", "Binding_Modes": "9",
         "Docked_File": "docked/ligA_out.pdbqt", "Status": "OK"},
        {"Ligand": "ligB", "Binding_Affinity": "-5.0", "Binding_Modes": "9", "Status": "OK"},
        {"Ligand": "ligOOW", "Binding_Affinity": "-25.0", "Binding_Modes": "1",
         "Status": "OUT_OF_WINDOW"},
        {"Ligand": "ligF", "Binding_Affinity": "FAILED", "Status": "FAILED"},
    ]

    def test_tables_have_thead_and_tbody(self, tmp_path, monkeypatch):
        html = self._render(tmp_path, monkeypatch, self.ROWS)
        table = html[html.index('<table id="all-results">'):]
        table = table[:table.index("</table>")]
        assert table.index("<thead>") < table.index("<th class") < table.index("</thead>")
        assert table.index("</thead>") < table.index("<tbody>") < table.index("ligA")
        assert html.count("<thead>") == 2

    def test_out_of_window_is_not_a_hit(self, tmp_path, monkeypatch):
        html = self._render(tmp_path, monkeypatch, self.ROWS)
        top = html[html.index('id="top"'):html.index('id="dist"')]
        assert "ligA" in top and "ligOOW" not in top
        assert "-25.00" not in html  # never shown as a best score

    def test_relative_docked_file_feeds_viewer(self, tmp_path, monkeypatch):
        html = self._render(tmp_path, monkeypatch, self.ROWS,
                            poses={"ligA_out.pdbqt": _pose_pdbqt("ligA")})
        assert "viewer-container" in html and '"name": "ligA"' in html

    def test_overlay_closure_captures_index(self, tmp_path, monkeypatch):
        html = self._render(tmp_path, monkeypatch, self.ROWS,
                            poses={"ligA_out.pdbqt": _pose_pdbqt("ligA")})
        assert '"ov" + idx' in html and '"ov" + i ' not in html

    def test_guide_matches_score_bands(self, tmp_path, monkeypatch):
        from vspipeline import report
        html = self._render(tmp_path, monkeypatch, self.ROWS)
        guide = html[html.index('id="guide"'):]
        for upper, _css, label in report.AFFINITY_BANDS:
            assert label in guide
        assert "&lt; -8 kcal/mol:</strong> Excellent" in guide
        assert report._confidence(-8.5)[1] == "Excellent"
        assert report._confidence(-7.9)[1] == "Good"
        assert report._confidence(-4.0)[1] == "Poor"
        assert "-10" not in guide

    def test_product_name_and_version(self, tmp_path, monkeypatch):
        html = self._render(tmp_path, monkeypatch, self.ROWS)
        assert f"AutoDocker {runner.__version__}" in html
        assert "Toolkit" not in html and "professional" not in html.lower()

    def test_no_emoji(self, tmp_path, monkeypatch):
        import re
        html = self._render(tmp_path, monkeypatch, self.ROWS,
                            poses={"ligA_out.pdbqt": _pose_pdbqt("ligA")})
        html = html.replace("/* fake ngl */", "")
        assert not re.search("[\U0001F300-\U0001FAFF☀-⛿✀-➿]", html)


class TestNglIntegrity:
    def test_download_hash_mismatch_fails_closed(self, monkeypatch, caplog):
        from vspipeline import report
        monkeypatch.setattr(report, "_ngl_local_candidates", lambda: [])
        monkeypatch.setattr(runner, "_http_get_bytes", lambda url, timeout=60: b"evil()")
        assert report._fetch_ngljs() is None

    def test_download_with_pinned_hash_accepted(self, monkeypatch):
        import hashlib
        from vspipeline import report
        payload = b"/* pinned */"
        monkeypatch.setattr(report, "_ngl_local_candidates", lambda: [])
        monkeypatch.setattr(report, "NGL_SHA256", hashlib.sha256(payload).hexdigest())
        monkeypatch.setattr(runner, "_http_get_bytes", lambda url, timeout=60: payload)
        assert report._fetch_ngljs() == payload

    def test_pinned_url_is_exact_version(self):
        from vspipeline import report
        assert f"ngl@{report.NGL_VERSION}/dist/ngl.js" in report.NGL_CDN_URL
        assert len(report.NGL_SHA256) == 64

    def test_templates_ship_with_package(self):
        from vspipeline import report
        for name in ("report.css", "report.js", "viewer.html"):
            assert os.path.isfile(os.path.join(report.TEMPLATE_DIR, name))
