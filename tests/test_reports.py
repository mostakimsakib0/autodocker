import csv
import os

import runner


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