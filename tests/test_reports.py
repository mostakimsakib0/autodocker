import csv
import os


import runner
from vspipeline import net
from vspipeline import report as report_mod
from vspipeline import results as results_mod
def _write(path, content):
    with open(path, "w", encoding="utf-8") as f:
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
    with open(path, "w", newline="", encoding="utf-8") as f:
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
        analyzer = results_mod.ResultsAnalyzer(outdir)
        results = [("ligA", -7.5), ("ligB", None)]
        analyzer.save_ranking(results, mode="extended",
                              metrics_dict={"ligA": {"status": "OK"},
                                            "ligB": {}})
        with open(analyzer.csv_file, encoding="utf-8") as f:
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
        report_mod.generate_html_report(csv_path, os.path.join(outdir, "docked"),
                                        html_path)
        assert os.path.exists(html_path)
        with open(html_path, encoding="utf-8") as f:
            content = f.read()
        assert "Virtual Screening Results Report" in content
        assert "ligA" in content

    def test_top_hits_follow_top_n(self, tmp_path):
        """The report highlights --top-n ligands, like Top_hits.txt."""
        outdir = str(tmp_path / "topn")
        os.makedirs(outdir)
        rows = [{"Ligand": f"lig{i}", "Binding_Affinity": f"-{9 - i}.0",
                 "SimScore": "0.9", "Binding_Modes": "9", "Status": "OK"}
                for i in range(5)]
        csv_path = _ranking_csv(outdir, rows)
        html_path = os.path.join(outdir, "results_report.html")
        report_mod.generate_html_report(csv_path, os.path.join(outdir, "docked"),
                                        html_path, top_n=3)
        with open(html_path, encoding="utf-8") as f:
            content = f.read()
        assert "Top 3 Hits" in content
        top_table = content.split("Top 3 Hits", 1)[1].split("</table>", 1)[0]
        assert "lig2" in top_table and "lig3" not in top_table

    def test_html_missing_csv_warns(self, outdir, caplog):
        report_mod.generate_html_report(
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
        report_mod.generate_html_report(csv_path, os.path.join(outdir, "docked"),
                                        html_path)
        with open(html_path, encoding="utf-8") as f:
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
        monkeypatch.setattr(report_mod, "_fetch_ngljs", lambda: ngl_bundle)
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
        report_mod.generate_html_report(
            csv_path, docked, html,
            receptor_pdbqt=os.path.join(outdir, "receptor_b.pdbqt"))
        with open(html, encoding="utf-8") as f:
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
        assert report_mod.NGL_CDN_URL not in content
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
        monkeypatch.setattr(report_mod, "_fetch_ngljs", lambda: b"/* fake ngl */")
        outdir = str(tmp_path / "out")
        os.makedirs(outdir, exist_ok=True)
        csv_path = _ranking_csv(outdir, [
            {"Ligand": "ligA", "Binding_Affinity": "-7.5",
             "SimScore": "0.9", "Binding_Modes": "9", "Status": "OK"},
        ])
        html = os.path.join(outdir, "results_report.html")
        report_mod.generate_html_report(
            csv_path, os.path.join(outdir, "docked"), html, receptor_pdbqt=None)
        with open(html, encoding="utf-8") as f:
            content = f.read()
        assert "No docked structures were available" in content
        assert "viewer-container" not in content

    def test_ngl_inlined_not_cdn(self, tmp_path, monkeypatch):
        outdir, content = self._prepare(tmp_path, monkeypatch)
        assert 'src="viewer/ngl.min.js"' not in content
        assert report_mod.NGL_CDN_URL not in content
        assert "/* fake ngl */" in content

    def test_ngl_unavailable_notice(self, tmp_path, monkeypatch):
        """If the NGL library cannot be bundled, show a notice (no CDN)."""
        monkeypatch.setattr(report_mod, "_fetch_ngljs", lambda: None)
        outdir = str(tmp_path / "out")
        docked = os.path.join(outdir, "docked")
        os.makedirs(docked, exist_ok=True)
        _write(os.path.join(docked, "ligA_out.pdbqt"), _pose_pdbqt("ligA"))
        csv_path = _ranking_csv(outdir, [
            {"Ligand": "ligA", "Binding_Affinity": "-7.5",
             "SimScore": "0.9", "Binding_Modes": "9", "Status": "OK"},
        ])
        html = os.path.join(outdir, "results_report.html")
        report_mod.generate_html_report(
            csv_path, docked, html, receptor_pdbqt=None)
        with open(html, encoding="utf-8") as f:
            content = f.read()
        assert "3D viewer unavailable" in content
        assert report_mod.NGL_CDN_URL not in content
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
        monkeypatch.setattr(report_mod, "_fetch_ngljs", lambda: ngl)
        outdir = str(tmp_path / "out")
        docked = os.path.join(outdir, "docked")
        os.makedirs(docked, exist_ok=True)
        for name, text in (poses or {}).items():
            _write(os.path.join(docked, name), text)
        csv_path = _ranking_csv(outdir, rows)
        html_path = os.path.join(outdir, "results_report.html")
        assert report_mod.generate_html_report(csv_path, docked, html_path)
        with open(html_path, encoding="utf-8") as f:
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
        monkeypatch.setattr(net, "_http_get_bytes", lambda url, timeout=60: b"evil()")
        assert report._fetch_ngljs() is None

    def test_download_with_pinned_hash_accepted(self, monkeypatch):
        import hashlib
        from vspipeline import report
        payload = b"/* pinned */"
        monkeypatch.setattr(report, "_ngl_local_candidates", lambda: [])
        monkeypatch.setattr(report, "NGL_SHA256", hashlib.sha256(payload).hexdigest())
        monkeypatch.setattr(net, "_http_get_bytes", lambda url, timeout=60: payload)
        assert report._fetch_ngljs() == payload

    def test_pinned_url_is_exact_version(self):
        from vspipeline import report
        assert f"ngl@{report.NGL_VERSION}/dist/ngl.js" in report.NGL_CDN_URL
        assert len(report.NGL_SHA256) == 64

    def test_templates_ship_with_package(self):
        from vspipeline import report
        for name in ("report.css", "report.html", "report.js", "viewer.html"):
            assert os.path.isfile(os.path.join(report.TEMPLATE_DIR, name))


# ---------------------------------------------------------------------------
# Reports without any in-window score; stale reports; viewer navigation
# ---------------------------------------------------------------------------

def _failed_rows():
    return [{"Ligand": "ligA", "Binding_Affinity": "FAILED", "Status": "FAILED"},
            {"Ligand": "ligB", "Binding_Affinity": "-0.5", "Status": "OUT_OF_WINDOW"}]


def test_report_without_scores_shows_na_and_no_viewer_link(tmp_path, monkeypatch):
    monkeypatch.setattr(report_mod, "_fetch_ngljs", lambda: b"/* fake ngl */")
    outdir = str(tmp_path)
    csv_path = _ranking_csv(outdir, _failed_rows())
    report = os.path.join(outdir, "report.html")
    assert report_mod.generate_html_report(csv_path, os.path.join(outdir, "docked"), report)
    content = open(report, encoding="utf-8").read()
    assert '<div class="stat-value">n/a</div>' in content
    assert "0.00" not in content  # no invented best/average affinity
    assert 'href="#viewer"' not in content and 'id="viewer"' not in content


def test_report_with_scores_links_the_viewer(tmp_path, monkeypatch):
    monkeypatch.setattr(report_mod, "_fetch_ngljs", lambda: b"/* fake ngl */")
    outdir = str(tmp_path)
    os.makedirs(os.path.join(outdir, "docked"))
    _write(os.path.join(outdir, "docked", "ligA_out.pdbqt"), _pose_pdbqt("ligA"))
    csv_path = _ranking_csv(outdir, [
        {"Ligand": "ligA", "Binding_Affinity": "-7.25", "Status": "OK",
         "Docked_File": "docked/ligA_out.pdbqt"}])
    report = os.path.join(outdir, "report.html")
    assert report_mod.generate_html_report(csv_path, os.path.join(outdir, "docked"), report)
    content = open(report, encoding="utf-8").read()
    assert 'href="#viewer"' in content and 'id="viewer"' in content
    assert '<div class="stat-value">-7.25</div>' in content


def test_report_early_return_removes_stale_report(tmp_path):
    outdir = str(tmp_path)
    report = _write(os.path.join(outdir, "report.html"), "stale report")
    csv_path = _ranking_csv(outdir, [])  # header only: nothing to report
    assert report_mod.generate_html_report(csv_path, os.path.join(outdir, "docked"),
                                           report) is False
    assert not os.path.exists(report)
    _write(report, "stale report")
    assert report_mod.generate_html_report(os.path.join(outdir, "missing.csv"),
                                           os.path.join(outdir, "docked"), report) is False
    assert not os.path.exists(report)


_VIEWER_DOM_STUB = r"""
function el(id) {
  return { id: id, textContent: "", style: {}, checked: true,
           addEventListener: function () {}, scrollIntoView: function () {},
           appendChild: function (o) { this.options.push(o); } };
}
var elements = {};
["viewer-select", "viewer-message", "viewer-container", "viewer-log",
 "viewer-error", "vs-box"].forEach(function (id) { elements[id] = el(id); });
var sel = elements["viewer-select"];
sel.options = [];
sel.selectedIndex = -1;
Object.defineProperty(sel, "value", {
  get: function () { return this.selectedIndex >= 0 ? this.options[this.selectedIndex].value : ""; },
  set: function (v) { for (var i = 0; i < this.options.length; i++)
                        if (this.options[i].value === v) this.selectedIndex = i; } });
var loadHandlers = [];
var document = { getElementById: function (id) { return elements[id]; },
                 createElement: function () { return {}; } };
var window = { addEventListener: function (ev, f) { if (ev === "load") loadHandlers.push(f); } };
"""

_VIEWER_DOM_CHECK = r"""
loadHandlers.forEach(function (f) { f(); });
var msg = elements["viewer-message"];
var shown = [];
// NGL is not defined here: every (re)load of a hit ends with this message.
function clicked(name) { msg.textContent = ""; window.selectHit(name); shown.push([sel.value, msg.textContent]); }
clicked("ligB");      // a top hit: selected and (re)loaded
clicked("outsider");  // not in the viewer: no reload of the current hit
console.log(JSON.stringify(shown));
"""


def test_viewer_click_on_ligand_outside_top_hits_does_not_reload_current(tmp_path):
    """Clicking a scatter point of a ligand that is not in the viewer (outside
    --top-n) must not silently re-show the currently selected hit."""
    import json
    import re
    import shutil
    import subprocess

    import pytest

    node = shutil.which("node") or shutil.which("nodejs")
    if node is None:
        pytest.skip("node not available")
    template = os.path.join(os.path.dirname(runner.__file__), "vspipeline", "templates",
                            "viewer.html")
    text = open(template, encoding="utf-8").read()
    script = re.findall(r"<script>(.*?)</script>", text, re.S)[-1]
    data = {"receptor": "", "gridBox": None,
            "hits": [{"name": "ligA", "affinity": -8.0, "pose": "ATOM"},
                     {"name": "ligB", "affinity": -7.0, "pose": "ATOM"}]}
    script = script.replace("__VS_DATA_JSON__", json.dumps(data))
    js = tmp_path / "viewer_test.js"
    js.write_text(_VIEWER_DOM_STUB + script + _VIEWER_DOM_CHECK, encoding="utf-8")
    proc = subprocess.run([node, str(js)], capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr
    shown = json.loads(proc.stdout.strip().splitlines()[-1])
    assert shown[0] == ["ligB", "NGL viewer failed to load."]
    assert shown[1][0] == "ligB"
    assert "outsider is not among the top hits" in shown[1][1]


def test_reports_of_identical_inputs_are_byte_identical(tmp_path, monkeypatch):
    """No wall-clock timestamp and no output path in the report: identical
    runs (here in two output directories) write identical bytes."""
    monkeypatch.setattr(report_mod, "_fetch_ngljs", lambda: b"/* fake ngl */")
    reports = []
    for run in ("a", "b"):
        outdir = str(tmp_path / run)
        os.makedirs(os.path.join(outdir, "docked"))
        _write(os.path.join(outdir, "docked", "ligA_out.pdbqt"), _pose_pdbqt("ligA"))
        csv_path = _ranking_csv(outdir, [
            {"Ligand": "ligA", "Binding_Affinity": "-7.25", "Status": "OK",
             "Docked_File": "docked/ligA_out.pdbqt"},
            {"Ligand": "ligB", "Binding_Affinity": "FAILED", "Status": "FAILED"}])
        report = os.path.join(outdir, "results_report.html")
        assert report_mod.generate_html_report(
            csv_path, os.path.join(outdir, "docked"), report,
            meta={"Seed": 42, "Exhaustiveness": 8})
        with open(report, "rb") as f:
            reports.append(f.read())
    assert reports[0] == reports[1]
    assert b"Generated:" not in reports[0]
    assert runner.__version__.encode() in reports[0]
