import os


import vspipeline.results as res


def test_score_sort_key_failed_last():
    data = [("a", -7.0), ("b", None), ("c", -8.5)]
    data.sort(key=res._score_sort_key)
    names = [n for n, _ in data]
    # best affinity first, failed at the end
    assert names == ["c", "a", "b"]


def test_score_csv():
    assert res._score_csv(None) == "FAILED"
    assert res._score_csv(-7.25) == -7.25


def test_score_text():
    assert res._score_text(None) == "FAILED"
    assert res._score_text(-7.25) == "-7.25"


def test_save_ranking_simple(tmp_path):
    ra = res.ResultsAnalyzer(str(tmp_path))
    results = [("ligA", -7.0), ("ligB", None), ("ligC", -8.5)]
    ra.save_ranking(results)
    csv_path = tmp_path / "ranking.csv"
    assert csv_path.exists()
    lines = csv_path.read_text().strip().splitlines()
    assert lines[0] == "Ligand,Binding_Affinity,Status"
    # sorted: ligC, ligA, ligB(failed)
    assert "ligC" in lines[1] and "ligB" in lines[3]
    assert (tmp_path / "Results_full.txt").exists()
    assert (tmp_path / "Top_hits.txt").exists()


def test_save_ranking_with_metrics(tmp_path):
    ra = res.ResultsAnalyzer(str(tmp_path))
    results = [("ligA", -7.0), ("ligC", -8.5)]
    metrics = {
        "ligA": {"simscore": 0.8, "binding_modes": 3, "mode_coverage": 1.0,
                 "dock_file": "/x/ligA.pdbqt", "status": "OK"},
        "ligC": {"simscore": 0.6, "binding_modes": 1, "mode_coverage": 0.5,
                 "dock_file": "/x/ligC.pdbqt", "status": "OK"},
    }
    # provide metadata
    (tmp_path / "ligand_metadata.json").write_text(
        '{"ligA": {"mw": 300.0, "logp": 2.1, "source_url": "http://x"}}')
    ra.metadata = ra._load_metadata()
    ra.save_ranking(results, metrics_dict=metrics)
    lines = (tmp_path / "ranking.csv").read_text().strip().splitlines()
    assert "SimScore" in lines[0] and "MW" in lines[0]
    assert "ligC" in lines[1]  # best affinity (-8.5) sorts first
    assert "ligA" in lines[2]
    top = (tmp_path / "Top_hits.txt").read_text()
    assert "MW: 300.00" in top


def test_save_metrics_report(tmp_path):
    ra = res.ResultsAnalyzer(str(tmp_path))
    results = [("ligA", -7.0), ("ligB", None), ("ligC", -8.5)]
    metrics = {"ligA": {"simscore": 0.8, "binding_modes": 3}}
    ra.save_metrics_report(results, metrics)
    txt = (tmp_path / "metrics.txt").read_text()
    assert "METRICS SUMMARY" in txt
    assert "Successful dockings: 2" in txt
    assert "Out-of-window results: 0" in txt
    assert "Failed dockings: 1" in txt
    assert "Best affinity: -8.50" in txt


def _rows(path):
    import csv
    with open(path) as f:
        return list(csv.DictReader(f))


MIXED = [("good", -8.0), ("weak", None), ("bad", None), ("ok2", -6.0)]
MIXED_METRICS = {
    "good": {"status": "OK", "simscore": 0.5, "binding_modes": 9, "engine": "Vina"},
    "ok2": {"status": "OK", "simscore": None, "binding_modes": 1, "engine": "QuickVina"},
    "weak": {"status": "OUT_OF_WINDOW", "raw_score": -0.5, "binding_modes": 9,
             "simscore": 0.25, "engine": "Vina"},
    "bad": {"status": "FAILED", "error": "boom"},
}


def test_out_of_window_is_reported_ranked_last_and_counted(tmp_path):
    ra = res.ResultsAnalyzer(str(tmp_path))
    ra.save_ranking(MIXED, metrics_dict=MIXED_METRICS)
    ra.save_metrics_report(MIXED, MIXED_METRICS)
    rows = _rows(tmp_path / "ranking.csv")
    assert [r["Ligand"] for r in rows] == ["good", "ok2", "weak", "bad"]
    weak = rows[2]
    assert weak["Binding_Affinity"] == "-0.5" and weak["Status"] == "OUT_OF_WINDOW"
    assert rows[3]["Binding_Affinity"] == "FAILED" and rows[3]["Status"] == "FAILED"
    txt = (tmp_path / "metrics.txt").read_text()
    assert "Successful dockings: 2" in txt
    assert "Out-of-window results: 1" in txt
    assert "Failed dockings: 1" in txt
    assert "Worst affinity: -6.00" in txt  # out-of-window not in the summary
    top = (tmp_path / "Top_hits.txt").read_text()
    assert "Status: OUT_OF_WINDOW" in top and "-0.50 kcal/mol" in top


def test_legacy_rejected_status_maps_to_out_of_window(tmp_path):
    ra = res.ResultsAnalyzer(str(tmp_path))
    ra.save_ranking([("w", None)], metrics_dict={
        "w": {"status": "REJECTED (-0.50)", "raw_score": -0.5}})
    assert _rows(tmp_path / "ranking.csv")[0]["Status"] == "OUT_OF_WINDOW"


def test_failed_rows_have_blank_metrics_and_engine_column(tmp_path):
    ra = res.ResultsAnalyzer(str(tmp_path))
    ra.save_ranking(MIXED, metrics_dict=MIXED_METRICS)
    rows = {r["Ligand"]: r for r in _rows(tmp_path / "ranking.csv")}
    assert rows["bad"]["SimScore"] == "" and rows["bad"]["Binding_Modes"] == ""
    assert rows["ok2"]["SimScore"] == ""          # undefined SimScore stays blank
    assert rows["good"]["SimScore"] == "0.5"
    assert rows["good"]["Engine"] == "Vina" and rows["ok2"]["Engine"] == "QuickVina"


def test_docked_file_relative_and_outputs_reproducible(tmp_path):
    outputs = []
    for run in ("a", "b"):
        out = tmp_path / run
        out.mkdir()
        metrics = {k: dict(v) for k, v in MIXED_METRICS.items()}
        metrics["good"]["dock_file"] = str(out / "docked" / "good_out.pdbqt")
        ra = res.ResultsAnalyzer(str(out))
        ra.save_ranking(MIXED, metrics_dict=metrics)
        ra.save_metrics_report(MIXED, metrics)
        assert {r["Ligand"]: r["Docked_File"] for r in _rows(out / "ranking.csv")}[
            "good"] == os.path.join("docked", "good_out.pdbqt")
        outputs.append({n: (out / n).read_bytes() for n in
                        ("ranking.csv", "Results_full.txt", "Top_hits.txt", "metrics.txt")})
    assert outputs[0] == outputs[1]   # no timestamps, no absolute paths
    for blob in outputs[0].values():
        assert b"Date:" not in blob and b"Generated:" not in blob


def test_top_hits_wording_is_neutral(tmp_path):
    ra = res.ResultsAnalyzer(str(tmp_path))
    ra.save_ranking(MIXED, metrics_dict=MIXED_METRICS)
    top = (tmp_path / "Top_hits.txt").read_text()
    assert "PRIORITY" not in top and "Next steps" not in top
    assert "strong predicted affinity" in top
