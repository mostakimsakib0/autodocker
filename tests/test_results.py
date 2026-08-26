import os

import pytest

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
    assert "ADVANCED METRICS ANALYSIS" in txt
    assert "Successful dockings: 2" in txt
    assert "Failed dockings: 1" in txt
    assert "Best affinity: -8.50" in txt
