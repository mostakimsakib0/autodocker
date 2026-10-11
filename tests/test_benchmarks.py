"""Unit tests for the benchmark scripts (hand-computed expectations).

The scripts live in scripts/ and are loaded by path; they import the shared
helpers as ``_benchlib`` from their own directory.
"""
import csv
import importlib.util
import itertools
import json
import os
import random

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "scripts")
# Optional real-data check: point AUTODOCKER_AKT1_RANKING at an archived DUD-E
# akt1 ranking.csv (150 ligands, some failed) to run it; skipped otherwise.
AKT1_RANKING = os.environ.get("AUTODOCKER_AKT1_RANKING", "")


def _load(name):
    spec = importlib.util.spec_from_file_location(
        f"bench_{name}", os.path.join(SCRIPTS, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


be = _load("benchmark_enrichment")
br = _load("benchmark_redock")
bp = _load("benchmark_vina_parity")
bperf = _load("benchmark_performance")
fc = _load("fetch_casf_rcsb")
bl = _load("_benchlib")


def _ranked(labels):
    """Scores 1, 2, 3, ... in the given label order (lower = better)."""
    return [float(i) for i in range(1, len(labels) + 1)], list(labels)


# ---------------------------------------------------------------------------
# Enrichment factor
# ---------------------------------------------------------------------------

class TestEnrichmentFactor:
    # A D A D D D D D D D  (N = 10, n = 2, n/N = 0.2)
    LABELS = [True, False, True] + [False] * 7

    def test_top_one(self):
        scores, labels = _ranked(self.LABELS)
        ef = be.enrichment_factor(scores, labels, 0.10)
        assert ef["n_top"] == 1 and ef["tp"] == 1
        assert ef["ef"] == pytest.approx(5.0)       # (1/1) / 0.2
        assert ef["ef_max"] == pytest.approx(5.0)

    def test_top_two_and_three(self):
        scores, labels = _ranked(self.LABELS)
        ef2 = be.enrichment_factor(scores, labels, 0.20)
        assert ef2["ef"] == pytest.approx(2.5)      # (1/2) / 0.2
        assert ef2["ef_max"] == pytest.approx(5.0)  # (2/2) / 0.2
        ef3 = be.enrichment_factor(scores, labels, 0.30)
        assert ef3["ef"] == pytest.approx(10 / 3)   # (2/3) / 0.2
        assert ef3["ef_max"] == pytest.approx(10 / 3)

    def test_ceil_is_not_fooled_by_float_error(self):
        # 0.05 * 100 = 5.000000000000001 must give n_top = 5, not 6
        scores, labels = _ranked([True] * 10 + [False] * 90)
        assert be.enrichment_factor(scores, labels, 0.05)["n_top"] == 5

    def test_tie_at_cutoff_uses_expected_value(self):
        # A and D tied at the top; n_top = 1 -> E[TP] = 0.5
        scores = [-10.0, -10.0] + [float(i) for i in range(8)]
        labels = [True, False] + [True] + [False] * 7
        ef = be.enrichment_factor(scores, labels, 0.10)
        assert ef["tp"] == pytest.approx(0.5)
        assert ef["ef"] == pytest.approx(2.5)
        # input order does not matter
        ef_rev = be.enrichment_factor(scores[::-1], labels[::-1], 0.10)
        assert ef_rev["ef"] == pytest.approx(ef["ef"])

    def test_never_exceeds_max(self):
        # 1.1.x bug: (tp/n_pos)/fraction with n_top = ceil(fraction*N)
        # reported EF1% = 4.0 for 50 actives / 150 when the max is 3.0.
        rng = random.Random(1)
        for _ in range(200):
            n = rng.randint(5, 200)
            labels = [rng.random() < 0.3 for _ in range(n)]
            if not any(labels):
                labels[0] = True
            scores = [rng.choice([None, -9.0, -8.0]) if rng.random() < 0.1
                      else rng.uniform(-12, -2) for _ in range(n)]
            for frac in (0.01, 0.05, 0.2):
                ef = be.enrichment_factor(scores, labels, frac)
                assert ef["ef"] <= ef["ef_max"] + 1e-12

    def test_akt1_like_case_is_capped(self):
        scores, labels = _ranked([True] * 2 + [False] * 98 + [True] * 48 + [False] * 2)
        ef = be.enrichment_factor(scores, labels, 0.01)   # N=150, n_top=2
        assert ef["n_top"] == 2
        assert ef["ef"] == pytest.approx(3.0) == ef["ef_max"]


# ---------------------------------------------------------------------------
# AUC, BEDROC, ROC enrichment, Mann-Whitney
# ---------------------------------------------------------------------------

class TestRankingMetrics:
    def test_auc_with_ties(self):
        # actives 1, 2; decoys 2, 3: wins 1v2, 1v3, 2v3, tie 2v2 -> 3.5 / 4
        assert be.roc_auc([1.0, 2.0, 2.0, 3.0], [True, True, False, False]) == 0.875

    def test_auc_failed_ranked_last(self):
        assert be.roc_auc([None, -5.0], [True, False]) == 0.0
        assert be.roc_auc([None, None], [True, False]) == 0.5
        assert be.roc_auc([-5.0, None], [True, False]) == 1.0

    def test_auc_degenerate(self):
        assert be.roc_auc([1.0, 2.0], [False, False]) is None

    def test_bedroc_extremes(self):
        labels = [True] * 3 + [False] * 47
        scores = [float(i) for i in range(50)]
        assert be.bedroc(scores, labels) == pytest.approx(1.0, abs=1e-9)
        assert be.bedroc(scores[::-1], labels) == pytest.approx(0.0, abs=1e-6)

    def test_bedroc_matches_rdkit(self):
        scoring = pytest.importorskip("rdkit.ML.Scoring.Scoring")
        rng = random.Random(7)
        for _ in range(20):
            n = rng.randint(20, 120)
            labels = [rng.random() < 0.1 for _ in range(n)]
            labels[rng.randrange(n)] = True
            labels[rng.randrange(n)] = False
            scores = rng.sample(range(n * 10), n)  # distinct: no ties
            order = sorted(range(n), key=lambda i: scores[i])
            ref = scoring.CalcBEDROC([[labels[i]] for i in order], 0, 80.5)
            assert be.bedroc([float(s) for s in scores], labels) == pytest.approx(ref, abs=1e-12)

    def test_bedroc_ties_are_mean_over_tie_breaks(self):
        # BEDROC is affine in the exponential sum, so the tie rule must equal
        # the mean over every ordering of the tied block.
        n, alpha = 6, 20.0
        labels = [True, True, False, False, False, False]
        tied = be.bedroc([0.0] * n, labels, alpha)
        values = []
        for pos in itertools.combinations(range(n), 2):
            lab = [i in pos for i in range(n)]
            values.append(be.bedroc([float(i) for i in range(n)], lab, alpha))
        assert tied == pytest.approx(sum(values) / len(values), abs=1e-12)

    def test_roc_enrichment(self):
        scores, labels = _ranked([True, False, True, False, False, False])
        assert be.roc_enrichment(scores, labels, 0.25) == pytest.approx(4.0)   # TPR 1.0
        assert be.roc_enrichment(scores, labels, 0.125) == pytest.approx(4.0)  # TPR 0.5
        # A and D tied at the top: diagonal ROC step -> TPR(0.125) = 0.25
        tied = [1.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        assert be.roc_enrichment(tied, [True, False, True, False, False, False],
                                 0.125) == pytest.approx(2.0)

    def test_mannwhitney_small_p_not_rounded_to_zero(self):
        _u, p1, p2 = be.mannwhitney([-10.0 - i for i in range(50)],
                                    [-5.0 + i * 0.01 for i in range(100)])
        assert 0.0 < p1 < 1e-10
        assert p2 == pytest.approx(2 * p1, rel=1e-6)

    def test_mannwhitney_ties_and_failed(self):
        u1, p1, _ = be.mannwhitney([-1.0, -1.0], [-1.0, -1.0])
        assert u1 is not None and p1 is not None
        # failed actives rank last: actives look worse, p_one large
        _u, p1, _ = be.mannwhitney([None, None, None], [-5.0, -6.0, -7.0])
        assert p1 > 0.9

    def test_degenerate_returns_none(self):
        # no actives -> every metric None
        m = be.metric_set([-5.0 + i for i in range(10)], [False] * 10)
        assert m["auc"] is None and m["ef1"] is None and m["bedroc"] is None


# ---------------------------------------------------------------------------
# Bootstrap
# ---------------------------------------------------------------------------

def _screen(seed=3, n_act=30, n_dec=90):
    rng = random.Random(seed)
    names = [f"A_{i}" for i in range(n_act)] + [f"D_{i}" for i in range(n_dec)]
    scores = ([rng.gauss(-8.0, 1.0) for _ in range(n_act)]
              + [rng.gauss(-6.5, 1.0) for _ in range(n_dec)])
    scores[0] = scores[n_act] = None  # one failed active, one failed decoy
    return names, scores


class TestBootstrap:
    def test_point_estimates_inside_ci(self):
        names, scores = _screen()
        row = be.compute_stats(names, scores, seed=11, n_boot=400, ci=0.95)
        for m in ("auc", "bedroc", "ef5", "roce5"):
            assert row[f"{m}_lo"] <= row[m] <= row[f"{m}_hi"], m

    def test_deterministic_for_seed(self):
        names, scores = _screen()
        labels = [be.is_active(n) for n in names]
        a = be.bootstrap_ci(scores, labels, 100, seed=5)
        b = be.bootstrap_ci(scores, labels, 100, seed=5)
        c = be.bootstrap_ci(scores, labels, 100, seed=6)
        assert a == b
        assert a != c

    def test_ci_level_is_used(self):
        names, scores = _screen()
        labels = [be.is_active(n) for n in names]
        wide = be.bootstrap_ci(scores, labels, 300, seed=5, ci=0.99)["auc"]
        narrow = be.bootstrap_ci(scores, labels, 300, seed=5, ci=0.50)["auc"]
        assert wide[0] < narrow[0] < narrow[1] < wide[1]
        row = be.compute_stats(names, scores, 5, 50, ci=0.8)
        assert row["ci_level"] == 0.8

    def test_resamples_keep_size_and_duplicates(self, monkeypatch):
        seen = []
        real = be.metric_set

        def spy(scores, labels):
            seen.append((len(scores), sum(labels)))
            return real(scores, labels)

        monkeypatch.setattr(be, "metric_set", spy)
        scores = [-9.0, -8.0, -7.0, -6.0, -5.0]
        labels = [True, True, False, False, False]
        be.bootstrap_ci(scores, labels, 50, seed=1)
        assert seen and all(s == (5, 2) for s in seen)  # stratified, by index

    def test_failed_counts_reported(self):
        names, scores = _screen()
        row = be.compute_stats(names, scores, 1, 0)
        assert row["total"] == 120
        assert row["n_failed_actives"] == 1 and row["n_failed_decoys"] == 1
        assert row["auc_lo"] is None  # no bootstrap requested

    def test_small_n_top_warns(self):
        names, scores = _screen()
        warnings = be.compute_stats(names, scores, 1, 0)["warnings"]
        assert any("EF1%" in w for w in warnings)


# ---------------------------------------------------------------------------
# ranking.csv input (old and new formats, failed ligands)
# ---------------------------------------------------------------------------

def _write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    return str(path)


OLD_HEADER = ["Ligand", "Binding_Affinity", "SimScore", "Binding_Modes", "Mode_Coverage",
              "MW", "LogP", "HBA", "HBD", "TPSA", "Rotors", "ZINC_URL", "Docked_File",
              "Status"]


class TestRankingInput:
    def test_old_format_failed_rows_stay_in_n(self, tmp_path):
        path = _write_csv(tmp_path / "r.csv", OLD_HEADER, [
            ["A_0_x", "-9.1", 0, 1, 1, "", "", "", "", "", "", "", "/abs/A_0_x_out.pdbqt", "OK"],
            ["D_0_y", "-7.0", 0, 1, 1, "", "", "", "", "", "", "", "/abs/D_0_y_out.pdbqt", "OK"],
            ["A_1_z", "FAILED", 0, 1, 1, "", "", "", "", "", "", "", "", "FAILED"],
            ["D_1_w", "-0.5", 0, 1, 1, "", "", "", "", "", "", "", "", "REJECTED (-0.5)"],
        ])
        names, scores, info = be.collect_scores(
            path, expected=["A_0_x", "D_0_y", "A_1_z", "D_1_w", "D_2_missing"])
        assert len(names) == 5
        by = dict(zip(names, scores))
        assert by["A_0_x"] == -9.1 and by["A_1_z"] is None
        assert by["D_1_w"] is None and by["D_2_missing"] is None
        assert info["status"]["D_2_missing"] == "MISSING"
        row = be.compute_stats(names, scores, 1, 0)
        assert row["n_failed_actives"] == 1 and row["n_failed_decoys"] == 2

    def test_new_format_out_of_window_and_engine(self, tmp_path):
        header = OLD_HEADER[:-1] + ["Engine", "Status"]
        path = _write_csv(tmp_path / "r.csv", header, [
            ["A_0", "-9.123", 0, 1, 1, "", "", "", "", "", "", "", "docked/A_0_out.pdbqt", "vina", "OK"],
            ["D_0", "-0.4", 0, 1, 1, "", "", "", "", "", "", "", "docked/D_0_out.pdbqt", "vina",
             "OUT_OF_WINDOW"],
            ["D_1", "nan", 0, 1, 1, "", "", "", "", "", "", "", "", "vina", "OK"],
        ])
        names, scores, info = be.collect_scores(path, out_of_window="score")
        assert dict(zip(names, scores)) == {"A_0": -9.123, "D_0": -0.4, "D_1": None}
        _, scores_fail, _ = be.collect_scores(path, out_of_window="fail")
        assert dict(zip(names, scores_fail))["D_0"] is None
        rows = bl.read_ranking(path)
        assert rows[0]["engine"] == "vina"
        assert bl.resolve_docked_file(path, rows[0]["docked_file"]) == str(
            tmp_path / "docked" / "A_0_out.pdbqt")

    def test_unexpected_rows_ignored(self, tmp_path):
        path = _write_csv(tmp_path / "r.csv", ["Ligand", "Affinity"],
                          [["A_0.pdbqt", "-8"], ["stale", "-12"]])
        names, scores, info = be.collect_scores(path, expected=["A_0"])
        assert names == ["A_0"] and scores == [-8.0]
        assert info["unexpected"] == ["stale"]

    def test_missing_columns_raise(self, tmp_path):
        path = _write_csv(tmp_path / "r.csv", ["Name", "Score"], [["x", "1"]])
        with pytest.raises(ValueError):
            bl.read_ranking(path)

    @pytest.mark.skipif(not os.path.exists(AKT1_RANKING),
                        reason="set AUTODOCKER_AKT1_RANKING to run")
    def test_real_akt1_ranking_keeps_failures(self):
        names, scores, _ = be.collect_scores(AKT1_RANKING)
        row = be.compute_stats(names, scores, 1, 0)
        assert row["total"] == 150
        assert row["n_failed_actives"] + row["n_failed_decoys"] == sum(
            s is None for s in scores) > 0
        assert row["ef1"] <= row["ef1_max"]


class TestLigandPrep:
    def test_intermediates_kept_out_of_ligand_dir(self, tmp_path, monkeypatch):
        import gzip
        blocks = "".join(f"@<TRIPOS>MOLECULE\nmol{i}\n 1 0 0 0 0\nSMALL\n" for i in range(4))
        gz = tmp_path / "act.mol2.gz"
        with gzip.open(gz, "wt") as f:
            f.write(blocks)

        def fake_convert(mol2, pdbqt):
            if "mol2" in open(mol2, encoding="utf-8").read():   # 'mol2' name fails
                raise RuntimeError("boom")
            open(pdbqt, "w", encoding="utf-8").write("ATOM\n")

        monkeypatch.setattr(be, "mol2_to_pdbqt", fake_convert)
        lig_dir, scratch = tmp_path / "ligands", tmp_path / "scratch"
        lig_dir.mkdir()
        (lig_dir / "A_9_old.mol2").write_text("stale")
        paths, failures = be.sample_ligands(str(gz), str(lig_dir), "A_", 3, 1, str(scratch))
        be.clean_ligand_dir(str(lig_dir), paths)
        assert len(paths) == 3 and len(failures) == 1
        assert sorted(os.listdir(lig_dir)) == sorted(os.path.basename(p) for p in paths)
        assert all(p.endswith(".pdbqt") for p in os.listdir(lig_dir))
        assert any(p.endswith(".mol2") for p in os.listdir(scratch))


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

POSE = """MODEL 1
REMARK VINA RESULT:    -7.512      0.000      0.000
ROOT
ATOM      1  C1  UNL     1       0.000   0.000   0.000  0.00  0.00    +0.000 C
ATOM      2  H1  UNL     1       5.000   0.000   0.000  0.00  0.00    +0.000 HD
ENDROOT
BEGIN_RES SER A  45
ATOM      3  CB  SER A  45      50.000  50.000  50.000  0.00  0.00    +0.000 C
END_RES SER A  45
ENDMDL
MODEL 2
REMARK VINA RESULT:    -7.004      1.100      2.200
ATOM      1  C1  UNL     1       1.000   0.000   0.000  0.00  0.00    +0.000 C
ENDMDL
"""


class TestBenchlib:
    def test_parse_vina_models(self, tmp_path):
        path = tmp_path / "out.pdbqt"
        path.write_text(POSE)
        models = bl.parse_vina_models(str(path))
        assert [m["model"] for m in models] == [1, 2]
        assert [m["affinity_str"] for m in models] == ["-7.512", "-7.004"]
        assert models[0]["heavy_coords"] == [(0.0, 0.0, 0.0)]  # no H, no flex
        assert not any("SER" in ln for ln in models[0]["lines"])
        assert bl.first_vina_result(str(path)) == ("-7.512", -7.512)

    @pytest.mark.parametrize("text", [
        POSE,
        # SMINA docking output: minimizedAffinity only.
        "MODEL 1\nREMARK minimizedAffinity -8.12345\nATOM      1  C1  UNL     1"
        "       0.000   0.000   0.000  0.00  0.00    +0.000 C\nENDMDL\n",
        # smina --score_only copies the input's Vina remark: minimizedAffinity wins.
        "MODEL 1\nREMARK VINA RESULT:    -7.000      0.000      0.000\n"
        "REMARK minimizedAffinity -6.54321\nENDMDL\n",
        # A model without a score, and a file without MODEL records.
        "MODEL 1\nATOM      1  C1  UNL     1       0.000   0.000   0.000  0.00  0.00"
        "    +0.000 C\nENDMDL\nMODEL 2\nREMARK VINA RESULT:    -5.250      1.0      2.0\n"
        "ENDMDL\n",
        "REMARK VINA RESULT:    -9.876      0.000      0.000\n",
    ])
    def test_scores_match_pipeline_parser(self, tmp_path, text):
        """The standalone parser gives the pipeline's scores (same rules)."""
        from vspipeline.docking import parse_pdbqt_results
        path = tmp_path / "out.pdbqt"
        path.write_text(text)
        ours = [m["affinity"] for m in bl.parse_vina_models(str(path))
                if m["affinity_str"] is not None]
        pipeline = [aff for aff, _lb, _ub in parse_pdbqt_results(str(path))]
        assert ours == pipeline and pipeline
        first = bl.first_vina_result(str(path))[1]
        assert first == pipeline[0]

    def test_ligand_box_rule(self):
        center, size = bl.ligand_box([(0, 0, 0), (10, 2, 2)])
        assert center == (5.0, 1.0, 1.0) and size == 22.5     # max(22.5, 10 + 12)
        _, size = bl.ligand_box([(0, 0, 0), (15, 0, 0)])
        assert size == 27.0

    def test_frame_offset(self, tmp_path):
        (tmp_path / "frame_offset.txt").write_text(
            "# comment = ignored\noffset_x = 1.5\noffset_y = -2\noffset_z = 3.25\n")
        assert bl.read_frame_offset(str(tmp_path)) == (1.5, -2.0, 3.25)
        rf = _load("restore_frame")
        assert rf.read_offset(str(tmp_path)) == (1.5, -2.0, 3.25)

    def test_select_sample(self):
        ids = [f"{i}abc" for i in range(100)]
        a = bl.select_sample(ids, 10, 42)
        assert a == bl.select_sample(list(reversed(ids)), 10, 42)
        assert a == sorted(a) and len(set(a)) == 10
        assert a != bl.select_sample(ids, 10, 43)
        assert bl.select_sample(ids, 0, 1) == sorted(ids)

    def test_mol2_shift(self, tmp_path, mol2_file):
        out = tmp_path / "shifted.mol2"
        bl.shift_mol2(mol2_file, (1.0, 2.0, 3.0), str(out))
        atoms = bl.parse_mol2_atoms(str(out))
        assert (atoms[0]["x"], atoms[0]["y"], atoms[0]["z"]) == (-1.0, -2.0, -3.0)
        assert bl.mol2_heavy_coords(str(out)) == [(-1.0, -2.0, -3.0)]


# ---------------------------------------------------------------------------
# Re-docking
# ---------------------------------------------------------------------------

THREE_MODES = "".join(
    f"MODEL {i}\nREMARK VINA RESULT:    {-8 + i}.000      0.000      0.000\n"
    f"ATOM      1  C1  UNL     1    {i:8.3f}   0.000   0.000  0.00  0.00    +0.000 C\nENDMDL\n"
    for i in (1, 2, 3))


def _atom_x(path):
    line = next(ln for ln in open(path, encoding="utf-8") if ln.startswith("ATOM"))
    return float(line.split()[5])


class TestRedock:
    def test_mode1_failure_is_not_replaced_by_mode2(self, tmp_path):
        out = tmp_path / "out.pdbqt"
        out.write_text(THREE_MODES)

        def fake_rmsd(ref, query):
            x = _atom_x(query)
            if x == 1.0:
                raise RuntimeError("obrms crashed")
            return x  # mode 2 -> 2.0, mode 3 -> 3.0

        rmsds, errors = br.mode_rmsds("ref.mol2", str(out), rmsd_fn=fake_rmsd)
        assert rmsds == [None, 2.0, 3.0] and len(errors) == 1
        rec = br.apply_rmsds({"error": None}, rmsds, errors)
        assert rec["rmsd_top1"] is None
        assert rec["rmsd_best"] == 2.0
        assert rec["status"] == "rmsd_failed"

    def test_modes_kept_in_file_order(self, tmp_path):
        out = tmp_path / "out.pdbqt"
        out.write_text(THREE_MODES)
        rmsds, _ = br.mode_rmsds("ref", str(out),
                                 rmsd_fn=lambda r, q: _atom_x(q))
        rec = br.apply_rmsds({"error": None}, rmsds, [])
        assert rec["rmsd_top1"] == 1.0 and rec["status"] == "ok"

    def test_summary_denominators_and_median(self):
        results = [
            {"pdbid": "a", "status": "ok", "rmsd_top1": 1.0, "rmsd_best": 0.5, "error": None},
            {"pdbid": "b", "status": "ok", "rmsd_top1": 3.0, "rmsd_best": 1.5, "error": None},
            {"pdbid": "c", "status": "ok", "rmsd_top1": 1.5, "rmsd_best": 1.5, "error": None},
            {"pdbid": "d", "status": "ok", "rmsd_top1": 5.0, "rmsd_best": 4.0, "error": None},
            {"pdbid": "e", "status": "error", "rmsd_top1": None, "rmsd_best": None,
             "error": "docking failed"},
            {"pdbid": "f", "status": "missing", "rmsd_top1": None, "rmsd_best": None,
             "error": "fetch skipped"},
            {"pdbid": "g", "status": "skipped", "rmsd_top1": None, "rmsd_best": None,
             "error": "too large"},
        ]
        s = br.summarize(results, 2.0)
        assert (s["n_selected"], s["n_attempted"], s["n_evaluated"]) == (7, 5, 4)
        assert s["success_top1"] == 2
        assert s["success_top1_frac_evaluated"] == 0.5
        assert s["success_top1_frac_attempted"] == 0.4
        assert s["success_best"] == 3 and s["success_best_frac_attempted"] == 0.6
        assert s["median_rmsd_top1"] == 2.25   # mean of 1.5 and 3.0
        assert {e["pdbid"] for e in s["excluded"]} == {"e", "f", "g"}

    def test_resume_requires_same_params(self, tmp_path):
        params = br.docking_params({"seed": 42, "exhaustiveness": 8})
        assert not br.params_match(str(tmp_path), params)
        (tmp_path / "params.json").write_text(json.dumps({"docking": params}))
        assert br.params_match(str(tmp_path), params)
        changed = dict(params, exhaustiveness=16)
        assert not br.params_match(str(tmp_path), changed)

    def test_rdkit_conformer_is_deterministic(self, tmp_path, mol2_file):
        pytest.importorskip("rdkit")
        if not __import__("shutil").which("obabel"):
            pytest.skip("obabel not available")
        a, b = tmp_path / "a.mol2", tmp_path / "b.mol2"
        meta = bl.rdkit_start_conformer(mol2_file, str(a), seed=42)
        bl.rdkit_start_conformer(mol2_file, str(b), seed=42)
        assert meta["seed"] == 42 and meta["method"].startswith("ETKDGv3")
        strip = lambda p: [ln for ln in open(p, encoding="utf-8") if not ln.startswith("#")]
        assert strip(a) == strip(b)


# ---------------------------------------------------------------------------
# Parity
# ---------------------------------------------------------------------------

class TestParity:
    def _workdir(self, tmp_path, ranking_value="-7.512"):
        wd = tmp_path / "1abc"
        (wd / "docked").mkdir(parents=True)
        out = wd / "docked" / "1abc_ligand_out.pdbqt"
        out.write_text(POSE)
        (wd / "docked" / "1abc_ligand_vina.log").write_text(
            "# Attempt: Vina (primary)\n"
            f"# Command: /nonexistent/vina --receptor r --ligand l --out {out} --seed 42\n")
        _write_csv(wd / "ranking.csv", ["Ligand", "Binding_Affinity", "Status"],
                   [["1abc_ligand", ranking_value, "OK"]])
        return wd

    @staticmethod
    def _fake_vina(text):
        def run(cmd, timeout):
            open(cmd[cmd.index("--out") + 1], "w", encoding="utf-8").write(text)
        return run

    def test_pass(self, tmp_path):
        wd = self._workdir(tmp_path)
        rec = bp.check_complex(str(wd), "1abc", str(tmp_path), 10, runner=self._fake_vina(POSE))
        assert rec["status"] == "ok" and rec["modes"] == rec["matched"] == 2
        assert rec["ranking_string_match"] and rec["output_identical"]
        assert bp.verdict([rec])

    def test_ranking_rounding_is_a_mismatch(self, tmp_path):
        wd = self._workdir(tmp_path, ranking_value="-7.51")
        rec = bp.check_complex(str(wd), "1abc", str(tmp_path), 10, runner=self._fake_vina(POSE))
        assert rec["status"] == "mismatch"
        assert rec["ranking_numeric_match"] is False

    def test_errors_fail_the_run(self, tmp_path):
        wd = self._workdir(tmp_path)

        def boom(cmd, timeout):
            raise RuntimeError("vina not found")

        rec = bp.check_complex(str(wd), "1abc", str(tmp_path), 10, runner=boom)
        assert rec["status"] == "error"
        assert not bp.verdict([rec])
        assert not bp.verdict([])

    def test_string_compare_at_printed_precision(self):
        assert bp.ranking_matches("-7.1", "-7.100") == (True, True)
        assert bp.ranking_matches("-7.1484", "-7.148") == (False, True)
        assert bp.ranking_matches("FAILED", "-7.148") == (False, False)

    def test_last_logged_command_is_used(self, tmp_path):
        log = tmp_path / "x_vina.log"
        log.write_text("# Command: qvina --out a\n# Command: vina --out b\n")
        assert bp.extract_command(str(log)) == ["vina", "--out", "b"]


# ---------------------------------------------------------------------------
# fetch_casf_rcsb
# ---------------------------------------------------------------------------

def _pdb(record, serial, name, resname, chain, resnum, x, alt=" ", occ=1.0, element="C"):
    return (f"{record:<6s}{serial:5d} {name:<4s}{alt}{resname:>3s} {chain}{resnum:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{occ:6.2f}{20.0:6.2f}          {element:>2s}")


class TestFetch:
    def test_altloc_a_kept(self):
        lines = [
            _pdb("ATOM", 1, "CA", "SER", "A", 10, 1.0),
            _pdb("ATOM", 2, "OG", "SER", "A", 10, 2.0, alt="A", occ=0.4),
            _pdb("ATOM", 3, "OG", "SER", "A", 10, 3.0, alt="B", occ=0.6),
        ]
        out = fc.filter_altlocs(lines)
        assert len(out) == 2
        assert out[1][16] == " " and float(out[1][30:38]) == 2.0

    def test_altloc_highest_occupancy_without_a(self):
        lines = [
            _pdb("HETATM", 1, "C1", "LIG", "A", 900, 1.0, alt="B", occ=0.3),
            _pdb("HETATM", 2, "C2", "LIG", "A", 900, 2.0, alt="B", occ=0.3),
            _pdb("HETATM", 3, "C1", "LIG", "A", 900, 5.0, alt="C", occ=0.7),
            _pdb("HETATM", 4, "C2", "LIG", "A", 900, 6.0, alt="C", occ=0.7),
        ]
        out = fc.filter_altlocs(lines)
        assert [float(ln[30:38]) for ln in out] == [5.0, 6.0]
        # one conformer -> the formula counts the ligand once
        assert fc.het_groups(out)[0][1] == {"C": 2}

    def test_receptor_keeps_metals_not_ligand_or_water(self):
        recs = [
            _pdb("ATOM", 1, "CA", "ALA", "A", 1, 0.0),
            "TER",
            _pdb("HETATM", 2, "ZN", "ZN", "A", 301, 1.0, element="ZN"),
            _pdb("HETATM", 3, "O", "HOH", "A", 401, 2.0, element="O"),
            _pdb("HETATM", 4, "C1", "LIG", "A", 501, 3.0),
            _pdb("HETATM", 5, "C1", "NAG", "A", 601, 4.0),
        ]
        metals = fc.receptor_records(recs, "metals", ("LIG", "A", "501"))
        assert [ln[17:20].strip() for ln in metals if ln.startswith("HETATM")] == ["ZN"]
        every = fc.receptor_records(recs, "all", ("LIG", "A", "501"))
        assert [ln[17:20].strip() for ln in every if ln.startswith("HETATM")] == ["ZN", "NAG"]
        assert not any(ln.startswith("HETATM") for ln in fc.receptor_records(recs, "none"))

    def test_skeleton_ignores_charges_and_bond_orders(self):
        def sdf(elements, bonds):
            lines = ["x", "", "", f"{len(elements):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000"]
            lines += [f"    0.0000    0.0000    0.0000 {e:<3s} 0  0" for e in elements]
            lines += [f"{a:3d}{b:3d}{o:3d}  0" for a, b, o in bonds]
            return "\n".join(lines + ["M  END", "$$$$"])

        # acetate as C-C(=O)O vs a different bond-order/charge pattern
        a = fc.skeleton_hash(*fc.sdf_heavy_graph(sdf("CCOO", [(1, 2, 1), (2, 3, 2), (2, 4, 1)])))
        b = fc.skeleton_hash(*fc.sdf_heavy_graph(sdf("CCOOH", [(1, 2, 1), (2, 3, 1), (2, 4, 2),
                                                               (4, 5, 1)])))
        c = fc.skeleton_hash(*fc.sdf_heavy_graph(sdf("COCO", [(1, 2, 1), (2, 3, 1), (3, 4, 1)])))
        assert a == b
        assert a != c


# ---------------------------------------------------------------------------
# Performance
# ---------------------------------------------------------------------------

class TestPerformance:
    def _runs(self):
        runs = []
        for size, p, walls in ((10, 1, [100, 110, 90]), (10, 2, [60, 50, 55]),
                               (20, 1, [200, 220, 210]), (20, 4, [70, 80, 75])):
            for rep, w in enumerate(walls):
                runs.append({"n_ligands": size, "processes": p, "repeat": rep, "wall_s": w,
                             "exit": 0, "prep_s": 5, "dock_s": w - 6, "post_s": 1,
                             "peak_tree_rss_kb": 1000 * p, "peak_tree_pss_kb": None,
                             "max_single_process_rss_kb": 500, "vina_threads": 8 // p,
                             "docked_ok": size})
        return runs

    def test_speedup_vs_p1_same_size(self):
        rows = {(r["n_ligands"], r["processes"]): r for r in bperf.aggregate(self._runs())}
        assert rows[(10, 1)]["wall_s_median"] == 100
        assert (rows[(10, 2)]["wall_s_min"], rows[(10, 2)]["wall_s_max"]) == (50, 60)
        assert rows[(10, 2)]["speedup_vs_p1"] == pytest.approx(100 / 55)
        assert rows[(10, 2)]["efficiency"] == pytest.approx(100 / 55 / 2)
        assert rows[(20, 4)]["speedup_vs_p1"] == pytest.approx(210 / 75)
        assert rows[(10, 1)]["speedup_vs_p1"] == 1.0
        assert rows[(10, 2)]["throughput_ligands_per_min"] == pytest.approx(10 / 55 * 60)

    def test_no_p1_means_no_speedup(self):
        rows = bperf.aggregate([r for r in self._runs() if r["processes"] != 1])
        assert all(r["speedup_vs_p1"] is None for r in rows)

    def test_failed_repeat_excluded_from_median(self):
        runs = self._runs()
        runs[0]["exit"] = 1
        row = next(r for r in bperf.aggregate(runs)
                   if (r["n_ligands"], r["processes"]) == (10, 1))
        assert row["repeats_ok"] == 2 and row["wall_s_median"] == 100.0

    @pytest.mark.skipif(not os.path.isdir("/proc"), reason="needs /proc")
    def test_tree_memory_counts_children(self):
        import subprocess
        import sys
        import time
        proc = subprocess.Popen([sys.executable, "-c",
                                 "import subprocess,sys,time;"
                                 "subprocess.Popen([sys.executable,'-c','import time;time.sleep(5)']);"
                                 "time.sleep(5)"])
        try:
            time.sleep(1.5)
            pids = bperf._proc_tree_pids(proc.pid)
            rss, _pss = bperf.tree_memory_kb(proc.pid)
            assert len(pids) >= 2
            assert rss > 0
        finally:
            for pid in bperf._proc_tree_pids(proc.pid)[::-1]:
                try:
                    os.kill(pid, 9)
                except OSError:
                    pass
            proc.wait()

    def test_build_library_requires_pdbqt(self, tmp_path):
        with pytest.raises(RuntimeError, match="no \\*.pdbqt"):
            bperf.build_library(str(tmp_path), 5, 42, str(tmp_path))


def test_script_output_names_match_the_pipeline():
    """scripts/_benchlib.py mirrors vspipeline/outputs.py; they must agree."""
    from vspipeline import outputs
    assert bl.RANKING_CSV == outputs.RANKING_CSV
    assert bl.DOCKED_DIR == outputs.DOCKED_DIR
    assert bl.FRAME_OFFSET_TXT == outputs.FRAME_OFFSET_TXT
    for name in ("COL_LIGAND", "COL_AFFINITY", "COL_STATUS", "COL_DOCKED_FILE", "COL_ENGINE"):
        assert getattr(bl, name) == getattr(outputs, name), name
