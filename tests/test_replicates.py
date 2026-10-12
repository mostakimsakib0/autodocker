"""Seed replicates (vspipeline.replicates): seeds, names, checkpoint keys,
aggregation, output rows and summary statistics; plus real Vina docking of
2 replicates of 2 small ligands through DockJob / dock_ligand."""
import csv
import json
import os
import re
import shutil

import pytest

from vspipeline import checkpoint, outputs, replicates as rp

E2E = os.path.join(os.path.dirname(__file__), "data", "e2e")
VINA = shutil.which("vina")


def _ok(score, k, dock_dir="/run/docked", name="lig", status="OK", **extra):
    stem = rp.replicate_stem(name, k)
    return score, {"status": status, "raw_score": score,
                   "dock_file": f"{dock_dir}/{stem}_out.pdbqt",
                   "log_file": f"{dock_dir}/{stem}_vina.log",
                   "simscore": 0.5 + k, "binding_modes": 9, "engine": "vina", **extra}


def _failed(k, name="lig"):
    return None, {"status": "FAILED", "error": f"boom {k}",
                  "dock_file": f"/run/docked/{name}_rep{k}_out.pdbqt",
                  "log_file": f"/run/docked/{name}_rep{k}_vina.log"}


def _out_of_window(raw, k, name="lig"):
    score, m = _ok(raw, k, name=name, status="OUT_OF_WINDOW")
    return None, dict(m, error="outside the window")


# ---------------------------------------------------------------------------
# Seeds, names, keys
# ---------------------------------------------------------------------------

def test_seeds_are_seed_plus_k():
    assert rp.replicate_seeds(42, 3) == [42, 43, 44]
    assert rp.replicate_seed(42, 0) == 42
    with pytest.raises(ValueError):
        rp.replicate_seed(42, -1)


@pytest.mark.parametrize("n,seed,message", [
    (0, 42, "between 1 and"),
    (rp.MAX_REPLICATES + 1, 42, "between 1 and"),
    (3, 2**31 - 2, "exceeds"),
    (True, 42, "integer"),
    (2.0, 42, "integer"),
])
def test_check_replicates_rejects(n, seed, message):
    with pytest.raises(ValueError, match=message):
        rp.check_replicates(n, seed)


def test_check_replicates_accepts():
    rp.check_replicates(1, 2**31 - 1)
    rp.check_replicates(2, 2**31 - 2)
    rp.check_replicates(rp.MAX_REPLICATES, 42)


def test_stems_and_file_names():
    assert rp.replicate_stem("aspirin", 0) == "aspirin_rep0"
    assert rp.checkpoint_key("aspirin", 12) == "aspirin_rep12"
    assert rp.pose_file("d", "aspirin", 1) == os.path.join("d", "aspirin_rep1_out.pdbqt")
    assert rp.log_file("d", "aspirin", 1) == os.path.join("d", "aspirin_rep1_vina.log")
    assert rp.pose_file("d", "x", 2).endswith(outputs.POSE_SUFFIX)


@pytest.mark.parametrize("name,k", [("lig", 0), ("lig_rep1", 3), ("a_rep", 10),
                                    ("rep", 1), ("x_rep01y", 2)])
def test_stem_round_trip(name, k):
    assert rp.parse_replicate_stem(rp.replicate_stem(name, k)) == (name, k)


@pytest.mark.parametrize("stem", ["lig", "lig_rep", "lig_rep01", "lig_repx", "_rep1"])
def test_non_replicate_stems(stem):
    assert rp.parse_replicate_stem(stem) is None


def test_stems_are_unique_across_ligands_and_replicates():
    names = ["a", "a_rep1", "a_rep", "a_rep1_rep0", "b"]
    stems = [rp.replicate_stem(n, k) for n in names for k in range(12)]
    assert len(stems) == len(set(stems))


def test_replicate_params_and_fingerprint():
    params = {"seed": 42, "exhaustiveness": 8}
    assert rp.replicate_params(params, None) is params
    assert rp.replicate_params(params, 2) == {"seed": 44, "exhaustiveness": 8}
    assert params["seed"] == 42                          # not modified
    assert rp.fingerprint_params(params, 1) is params    # N = 1: fingerprint unchanged
    assert rp.fingerprint_params(params, 3) == {**params, "replicates": 3}


def test_expand_ligands(tmp_path):
    lig = tmp_path / "lib" / "aspirin.pdbqt"
    lig.parent.mkdir()
    lig.write_text("ATOM ...\n", encoding="utf-8")
    dock_dir = tmp_path / "docked"
    assert rp.expand_ligands([str(lig)], 1, str(dock_dir)) == [(str(lig), None)]
    assert not dock_dir.exists()                         # N = 1 stages nothing
    units = rp.expand_ligands([str(lig)], 2, str(dock_dir))
    assert [(os.path.basename(p), k) for p, k in units] == [
        ("aspirin_rep0.pdbqt", 0), ("aspirin_rep1.pdbqt", 1)]
    assert all(os.path.dirname(p) == str(dock_dir / rp.REPLICATE_INPUTS_DIR) for p, _ in units)


def test_staged_copy_is_byte_identical_and_follows_the_source(tmp_path):
    lig = tmp_path / "x.pdbqt"
    lig.write_bytes(b"ATOM 1\n")
    staged = rp.stage_ligand(str(lig), 1, str(tmp_path / "docked"))
    assert open(staged, "rb").read() == b"ATOM 1\n"
    mtime = os.stat(staged).st_mtime_ns
    rp.stage_ligand(str(lig), 1, str(tmp_path / "docked"))
    assert os.stat(staged).st_mtime_ns == mtime          # unchanged source: not rewritten
    lig.write_bytes(b"ATOM 2\n")
    rp.stage_ligand(str(lig), 1, str(tmp_path / "docked"))
    assert open(staged, "rb").read() == b"ATOM 2\n"


def test_replicate_job(tmp_path):
    from vspipeline.docking import DockJob, EngineConfig
    lig = tmp_path / "lib" / "x.pdbqt"
    lig.parent.mkdir()
    lig.write_text("ATOM\n", encoding="utf-8")
    job = DockJob("r.pdbqt", str(lig), "grid.conf", str(tmp_path / "docked"),
                  {"seed": 7, "exhaustiveness": 1}, None, EngineConfig("vina"))
    rjob = rp.replicate_job(job, 3)
    assert rjob.vina_params == {"seed": 10, "exhaustiveness": 1}
    assert os.path.basename(rjob.ligand) == "x_rep3.pdbqt"
    assert (rjob.receptor, rjob.grid_file, rjob.dock_dir, rjob.engine) == (
        job.receptor, job.grid_file, job.dock_dir, job.engine)
    assert job.vina_params["seed"] == 7


def test_checkpoint_resumes_per_replicate(tmp_path):
    """The checkpoint keys are the replicate stems; the staged copy has the
    source's hash, so a resumed run finds every finished (ligand, replicate)
    and only those."""
    out = tmp_path / "out"
    out.mkdir()
    lig = tmp_path / "x.pdbqt"
    lig.write_text("ATOM\n", encoding="utf-8")
    staged = [p for p, _ in rp.expand_ligands([str(lig)], 3, str(out / "docked"))]
    from vspipeline.docking import _file_sha256
    assert {_file_sha256(p) for p in staged} == {_file_sha256(str(lig))}
    ckpt = checkpoint.DockingCheckpoint(str(out), fingerprint="f")
    for k in (0, 2):
        ckpt.save_result(rp.checkpoint_key("x", k), -7.0 - k, {"status": "OK"},
                         ligand_hash=_file_sha256(staged[k]))
    ckpt.close()
    ckpt = checkpoint.DockingCheckpoint(str(out), fingerprint="f")
    done = [k for k in range(3) if ckpt.is_completed(rp.checkpoint_key("x", k),
                                                     _file_sha256(staged[k]))]
    ckpt.close()
    assert done == [0, 2]


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def test_best_replicate_is_reported_with_its_pose():
    score, m = rp.aggregate_ligand("lig", {0: _ok(-7.0, 0), 1: _ok(-7.5, 1), 2: _ok(-6.9, 2)},
                                   3, 42)
    assert score == -7.5
    assert m["best_replicate"] == 1 and m["best_seed"] == 43
    assert m["dock_file"].endswith("lig_rep1_out.pdbqt")
    assert m["log_file"].endswith("lig_rep1_vina.log")
    assert m["simscore"] == 1.5                          # the reported replicate's metrics
    assert m["n_replicates"] == 3 and m["replicates_requested"] == 3
    assert m["score_min"] == -7.5 and m["score_max"] == -6.9
    assert m["score_mean"] == pytest.approx(-7.1333333)
    assert m["score_sd"] == pytest.approx(0.3214550)    # sample SD (n - 1)
    assert [r["seed"] for r in m["replicates"]] == [42, 43, 44]
    assert all("_metrics" not in r for r in m["replicates"])


def test_ties_go_to_the_lowest_replicate():
    _, m = rp.aggregate_ligand("lig", {0: _ok(-7.0, 0), 1: _ok(-7.2, 1), 2: _ok(-7.2, 2)},
                               3, 1)
    assert m["best_replicate"] == 1


def test_partial_failure_scores_from_the_others():
    score, m = rp.aggregate_ligand("lig", {0: _failed(0), 1: _ok(-6.0, 1),
                                           2: _out_of_window(0.5, 2)}, 3, 42)
    assert score == -6.0 and m["status"] == "OK"
    assert m["n_replicates"] == 1 and m["score_sd"] is None
    assert m["n_replicates_failed"] == 1 and m["n_replicates_out_of_window"] == 1
    assert m["score_mean"] == m["score_min"] == m["score_max"] == -6.0


def test_missing_replicates_count_as_failed():
    score, m = rp.aggregate_ligand("lig", {1: _ok(-6.0, 1)}, 2, 42)
    assert score == -6.0 and m["n_replicates_failed"] == 1
    assert m["replicates"][0]["status"] == "FAILED"


def test_all_out_of_window():
    score, m = rp.aggregate_ligand("lig", {0: _out_of_window(0.3, 0),
                                           1: _out_of_window(0.1, 1)}, 2, 42)
    assert score is None and m["status"] == "OUT_OF_WINDOW"
    assert m["raw_score"] == 0.3                         # replicate 0's
    assert m["n_replicates"] == 0 and m["score_mean"] is None
    assert m["best_replicate"] is None


def test_all_failed_lists_every_error():
    score, m = rp.aggregate_ligand("lig", {0: _failed(0), 1: _failed(1)}, 2, 42)
    assert score is None and m["status"] == "FAILED"
    assert "replicate 0 (seed 42): boom 0" in m["error"]
    assert "replicate 1 (seed 43): boom 1" in m["error"]


def test_aggregate_in_library_order():
    scores, metrics = {}, {}
    for name, values in (("b", (-5.0, -6.0)), ("a", (-8.0, -7.0))):
        for k, v in enumerate(values):
            s, m = _ok(v, k, name=name)
            scores[rp.replicate_stem(name, k)] = s
            metrics[rp.replicate_stem(name, k)] = m
    results, md = rp.aggregate(["b", "a", "c"], scores, metrics, 2, 42)
    assert results == [("b", -6.0), ("a", -8.0), ("c", None)]
    assert md["c"]["status"] == "FAILED" and md["a"]["best_replicate"] == 0


def test_ranking_values_and_columns():
    _, m = rp.aggregate_ligand("lig", {0: _ok(-7.0, 0), 1: _ok(-7.5, 1)}, 2, 42)
    assert rp.RANKING_REPLICATE_COLUMNS == ("n_replicates", "score_mean", "score_sd",
                                            "score_min", "score_max")
    assert rp.ranking_values(m) == [2, -7.25, 0.354, -7.5, -7.0]
    _, single = rp.aggregate_ligand("lig", {0: _ok(-7.0, 0), 1: _failed(1)}, 2, 42)
    assert rp.ranking_values(single) == [1, -7.0, "", -7.0, -7.0]
    _, none = rp.aggregate_ligand("lig", {}, 2, 42)
    assert rp.ranking_values(none) == [0, "", "", "", ""]
    assert not rp.has_replicates({"status": "OK"}) and rp.has_replicates(m)


def test_top_hits_line():
    _, m = rp.aggregate_ligand("lig", {0: _ok(-7.0, 0), 1: _ok(-7.5, 1)}, 2, 42)
    assert rp.top_hits_line(m) == ("Replicates: 2/2 scored; mean -7.25, SD 0.35, range "
                                   "-7.50 to -7.00 (reported: replicate 1, seed 43)")
    assert rp.top_hits_line({"status": "OK"}) is None


def test_replicate_scores_csv(tmp_path):
    out = tmp_path / "run"
    dock = str(out / "docked")
    _, m = rp.aggregate_ligand("lig", {0: _ok(-7.0, 0, dock_dir=dock),
                                       1: _ok(-7.5, 1, dock_dir=dock)}, 2, 42)
    path = tmp_path / "replicate_scores.csv"
    rp.write_replicate_scores(str(path), str(out), ["lig"], {"lig": m})
    rows = list(csv.reader(open(path, encoding="utf-8")))
    assert rows[0] == list(rp.REPLICATE_SCORES_COLUMNS)
    assert rows[1] == ["lig", "0", "42", "OK", "-7.0", "-7.0", "no",
                       "docked/lig_rep0_out.pdbqt", "docked/lig_rep0_vina.log"]
    assert rows[2][:7] == ["lig", "1", "43", "OK", "-7.5", "-7.5", "yes"]


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

def _screen(values):
    """metrics_dict of a screen: {name: [score of replicate k or None]}."""
    md = {}
    for name, scores in values.items():
        per = {k: (_ok(s, k, name=name) if s is not None else _failed(k, name))
               for k, s in enumerate(scores)}
        md[name] = rp.aggregate_ligand(name, per, len(scores), 42)[1]
    return md


def test_summary_statistics():
    md = _screen({"a": [-8.0, -8.2, -8.1], "b": [-7.0, -7.0, -7.6],
                  "c": [-6.0, -6.5, -6.1], "d": [-5.0, None, None]})
    s = rp.summary(md, 3, 42)
    assert s["seeds"] == [42, 43, 44] and s["ligands"] == 4
    assert s["ligands_with_sd"] == 3
    # per-ligand SDs: a 0.1, b 0.3464, c 0.2646 -> median 0.2646
    assert s["median_sd"] == pytest.approx(0.264575, abs=1e-6)
    assert s["max_sd"] == pytest.approx(0.34641, abs=1e-5)
    assert s["spearman_rep0_rep1"] == {"rho": pytest.approx(1.0), "n": 3}
    assert s["best_replicate_counts"] == [1, 2, 1]       # d: 0; a: 1; b: 2; c: 1
    assert s["ligands_all_replicates_scored"] == 3
    assert s["ligands_with_failed_replicates"] == ["d"]
    json.dumps(s, allow_nan=False)


def test_spearman_uses_the_consensus_implementation():
    from vspipeline.consensus import scorer_agreement_spearman
    md = _screen({"a": [-8.0, -6.0], "b": [-7.0, -7.0], "c": [-6.0, -8.0], "d": [-5.0, -5.5]})
    s = rp.summary(md, 2, 42)
    names = sorted(md)
    expected = scorer_agreement_spearman(rp.replicate_scores(md, 0, names),
                                         rp.replicate_scores(md, 1, names))
    assert (s["spearman_rep0_rep1"]["rho"], s["spearman_rep0_rep1"]["n"]) == expected
    assert s["spearman_rep0_rep1"]["rho"] == pytest.approx(0.2)   # 1 - 6*8/(4*15)


def test_summary_undefined_statistics():
    s = rp.summary(_screen({"a": [-8.0, None]}), 2, 42)
    assert s["median_sd"] is None and s["spearman_rep0_rep1"] == {"rho": None, "n": 0}
    lines = rp.log_lines(s)
    assert "seeds 42..43" in lines[0] and any("failed replicates" in t for t in lines)


def test_summary_file_is_deterministic(tmp_path):
    md = _screen({"b": [-7.0, -7.4], "a": [-8.0, -8.2]})
    p1, p2 = tmp_path / "1.json", tmp_path / "2.json"
    rp.write_summary(str(p1), rp.summary(md, 2, 42))
    rp.write_summary(str(p2), rp.summary(dict(reversed(list(md.items()))), 2, 42))
    assert p1.read_bytes() == p2.read_bytes()


# ---------------------------------------------------------------------------
# Real docking: 2 replicates of 2 ligands
# ---------------------------------------------------------------------------

@pytest.mark.skipif(VINA is None, reason="vina not available on PATH")
def test_two_replicates_dock_with_different_seeds(tmp_path):
    from vspipeline import docking
    with open(os.path.join(E2E, "box.json"), encoding="utf-8") as f:
        box = json.load(f)
    (cx, cy, cz), s = box["center"], box["size"]
    grid = tmp_path / "grid.conf"
    grid.write_text(f"center_x = {cx}\ncenter_y = {cy}\ncenter_z = {cz}\n"
                    f"size_x = {s}\nsize_y = {s}\nsize_z = {s}\n", encoding="utf-8")
    dock_dir = tmp_path / "out" / outputs.DOCKED_DIR
    dock_dir.mkdir(parents=True)
    params = {"exhaustiveness": 1, "binding_modes": 3, "energy_range": 3.0, "seed": 42,
              "min_valid_affinity": -1.0, "threads": 1}
    engine = docking.EngineConfig(vina_primary=VINA, timeout=300)
    names = ["lig_paracetamol", "lig_ibuprofen"]
    scores, metrics = {}, {}
    for name in names:
        job = docking.DockJob(os.path.join(E2E, "receptor.pdbqt"),
                              os.path.join(E2E, name + ".pdbqt"), str(grid), str(dock_dir),
                              params, None, engine)
        for k, (stem, score, m) in enumerate(rp.dock_replicates(job, 2)):
            assert stem == rp.replicate_stem(name, k)
            assert m["status"] in ("OK", "WEAK"), m.get("error")
            assert m["dock_file"] == rp.pose_file(str(dock_dir), name, k)
            scores[stem], metrics[stem] = score, m

    for name in names:
        for k in range(2):
            log = open(rp.log_file(str(dock_dir), name, k), encoding="utf-8").read()
            command = next(line for line in log.splitlines() if line.startswith("# Command:"))
            assert re.search(rf"--seed {42 + k}(\s|$)", command)
            assert f"(random seed: {42 + k})" in log            # Vina's own echo
            assert os.path.isfile(rp.pose_file(str(dock_dir), name, k))
        # no un-suffixed outputs: N > 1 writes replicate files only
        assert not os.path.exists(dock_dir / (name + outputs.POSE_SUFFIX))

    results, md = rp.aggregate(names, scores, metrics, 2, 42)
    for name, score in results:
        rep = [docking.parse_pdbqt_results(rp.pose_file(str(dock_dir), name, k))[0][0]
               for k in range(2)]
        assert score == min(rep)
        best = min(range(2), key=lambda k: (rep[k], k))
        assert md[name]["best_replicate"] == best
        assert md[name]["dock_file"] == rp.pose_file(str(dock_dir), name, best)
        assert md[name]["score_min"] == min(rep) and md[name]["score_max"] == max(rep)
        assert md[name]["n_replicates"] == 2
    summary = rp.summary(md, 2, 42)
    assert summary["ligands_with_sd"] == 2
    assert summary["spearman_rep0_rep1"]["n"] == 2
