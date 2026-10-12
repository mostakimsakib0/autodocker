"""scripts/merge_shards.py: sharded benchmark outputs merge into exactly the
files a single run writes.

The shard and single-run outputs are written by the benchmark scripts' own
main() functions; only the docking step (process_complex /
process_target / check_complex) is replaced by fake records, so the
tables, summaries and protocols are the scripts' real output.
"""
import csv
import importlib.util
import json
import os

import pytest

SCRIPTS = os.path.join(os.path.dirname(__file__), "..", "scripts")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(SCRIPTS, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ms = _load("merge_shards")
br, be, bp, bl = ms.br, ms.be, ms.bp, ms._benchlib


class _SerialPool:
    """multiprocessing.Pool stand-in (the fakes need not be picklable)."""

    def __init__(self, processes=None):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def imap_unordered(self, func, jobs):
        # Reverse order: the scripts must not depend on completion order.
        return [func(j) for j in reversed(list(jobs))]


class _FakeMultiprocessing:
    Pool = _SerialPool


def _read(path):
    with open(path, "rb") as f:
        return f.read()


def _tree(root):
    """{relative path: bytes} of every file below root."""
    out = {}
    for dirpath, _, files in os.walk(root):
        for name in files:
            path = os.path.join(dirpath, name)
            out[os.path.relpath(path, root)] = _read(path)
    return out


@pytest.fixture(autouse=True)
def _no_external_tools(monkeypatch):
    monkeypatch.setattr(bl, "runner_version", lambda *a, **k: "AutoDocker 2.0.0")
    monkeypatch.setattr(bl, "tool_version", lambda *a, **k: "AutoDock Vina v1.2.7")
    monkeypatch.setattr(bl, "rdkit_version", lambda: "2023.09.3")


# ---------------------------------------------------------------------------
# J1 re-docking
# ---------------------------------------------------------------------------

def _redock_record(pdbid, status, rmsds=None, error=None):
    rec = {"pdbid": pdbid, "status": status, "affinity": None, "affinity_mode1": None,
           "rmsd_top1": None, "rmsd_best": None, "rmsd_modes": [], "reused": False,
           "error": error}
    if rmsds is not None:
        rec["affinity"] = -7.123
        rec["affinity_mode1"] = -7.123
        br.apply_rmsds(rec, rmsds, [])
        rec["start_conformer"] = {"method": "ETKDGv3+MMFF94", "seed": 42}
    return rec


# Full-precision RMSDs: 2.0004 and 1.9996 both print as 2.000 in
# results.csv but fall on different sides of the 2 A cutoff.
FAKE_REDOCK = {
    "1aaa": _redock_record("1aaa", "ok", [1.23456789, 0.9, 3.5]),
    "1bbb": _redock_record("1bbb", "ok", [2.0004, 1.9996]),
    "1ccc": _redock_record("1ccc", "ok", [1.9996, 4.0]),
    "1ddd": _redock_record("1ddd", "error", error='docking failed (rc=1): "x", y'),
    "1eee": _redock_record("1eee", "rmsd_failed", [None, 1.5], error="mode 1: obrms"),
    "1fff": _redock_record("1fff", "missing", error="fetch skipped: peptide"),
    "1ggg": _redock_record("1ggg", "ok", [5.4321, 0.7654321]),
    "1hhh": _redock_record("1hhh", "ok", [0.3333333, 0.3333333]),
}
ALL_IDS = sorted(FAKE_REDOCK)


@pytest.fixture
def fake_redock(monkeypatch):
    monkeypatch.setattr(br, "multiprocessing", _FakeMultiprocessing)
    monkeypatch.setattr(br, "process_complex",
                        lambda job: json.loads(json.dumps(FAKE_REDOCK[job[0]])))


def _run_redock(tmp_path, name, ids, *extra):
    outdir = str(tmp_path / name)
    status = ms.main(["capture", "redock", "--", "--data-dir", str(tmp_path / "data"),
                      "--pdbids", ",".join(ids), "--outdir", outdir,
                      "--processes", "4", "--threads", "1", "--seed", "42", *extra])
    assert status == 0
    return outdir


def _redock_shards(tmp_path, n=3, ids=ALL_IDS):
    return [_run_redock(tmp_path, f"shard{k}", ids[k::n]) for k in range(n)]


def test_redock_merge_equals_single_run(tmp_path, fake_redock):
    single = _run_redock(tmp_path, "single", ALL_IDS)
    merged = str(tmp_path / "merged")
    assert ms.main(["redock", "--shards", *_redock_shards(tmp_path), "--outdir", merged,
                    "--expected-ids", ",".join(ALL_IDS)]) == 0
    for name in ("results.csv", "summary.json"):
        assert _read(os.path.join(merged, name)) == _read(os.path.join(single, name)), name
    summary = json.loads(_read(os.path.join(merged, "summary.json")))
    # Exact values, not the 3-decimal results.csv ones: 1bbb fails, 1ccc passes.
    assert summary["success_top1"] == 3
    assert summary["median_rmsd_top1"] == 1.9996
    assert summary["protocol"]["pdbids"] == ALL_IDS


def test_redock_capture_writes_full_precision_records(tmp_path, fake_redock):
    outdir = _run_redock(tmp_path, "single", ALL_IDS)
    records = json.loads(_read(os.path.join(outdir, ms.RECORDS_JSON)))
    assert [r["pdbid"] for r in records] == ALL_IDS
    assert records[1]["rmsd_top1"] == 2.0004
    with open(os.path.join(outdir, "results.csv"), newline="", encoding="utf-8") as f:
        rows = {r["pdbid"]: r for r in csv.DictReader(f)}
    assert rows["1bbb"]["rmsd_top1"] == "2.000"


def test_redock_merge_is_order_independent(tmp_path, fake_redock):
    shards = _redock_shards(tmp_path, n=4)
    outs = []
    for k, order in enumerate((shards, shards[::-1], shards[1:] + shards[:1])):
        out = str(tmp_path / f"merged{k}")
        ms.merge_redock(order, out)
        outs.append(_tree(out))
    assert outs[0] == outs[1] == outs[2]


def test_redock_protocol_mismatch_is_refused(tmp_path, fake_redock, capsys):
    shards = [_run_redock(tmp_path, "a", ALL_IDS[:4]),
              _run_redock(tmp_path, "b", ALL_IDS[4:], "--exhaustiveness", "4")]
    with pytest.raises(ms.MergeError, match="exhaustiveness"):
        ms.merge_redock(shards, str(tmp_path / "m1"))
    assert ms.main(["redock", "--shards", *shards, "--outdir", str(tmp_path / "m2")]) == 1
    assert "different protocols" in capsys.readouterr().err


def test_redock_missing_shard_is_reported(tmp_path, fake_redock):
    shards = _redock_shards(tmp_path, n=3)
    report = str(tmp_path / "report.json")
    status = ms.main(["redock", "--shards", *shards[:2], "--outdir", str(tmp_path / "m"),
                      "--expected-ids", ",".join(ALL_IDS), "--report", report])
    assert status == ms.EXIT_MISSING
    rep = json.loads(_read(report))
    assert rep["missing"] == ALL_IDS[2::3]
    assert rep["complete"] is False
    # The shards that ran are still merged.
    summary = json.loads(_read(str(tmp_path / "m" / "summary.json")))
    assert summary["n_selected"] == len(ALL_IDS) - len(ALL_IDS[2::3])
    assert ms.main(["redock", "--shards", *shards[:2], "--outdir", str(tmp_path / "m2"),
                    "--expected-ids", ",".join(ALL_IDS), "--allow-missing"]) == 0


def test_redock_duplicate_and_unexpected_ids_are_refused(tmp_path, fake_redock):
    a = _run_redock(tmp_path, "a", ALL_IDS[:3])
    b = _run_redock(tmp_path, "b", ALL_IDS[2:5])
    with pytest.raises(ms.MergeError, match="two shards"):
        ms.merge_redock([a, b], str(tmp_path / "m1"))
    with pytest.raises(ms.MergeError, match="not in --expected-ids"):
        ms.merge_redock([a], str(tmp_path / "m2"), expected=ALL_IDS[1:])


def test_redock_records_must_reproduce_the_summary(tmp_path, fake_redock):
    shard = _run_redock(tmp_path, "a", ALL_IDS)
    path = os.path.join(shard, ms.RECORDS_JSON)
    records = json.loads(_read(path))
    records[0]["rmsd_top1"] = 9.0
    with open(path, "w", encoding="utf-8") as f:
        json.dump(records, f)
    with pytest.raises(ms.MergeError, match="does not reproduce"):
        ms.merge_redock([shard], str(tmp_path / "m"))
    os.remove(path)
    with pytest.raises(ms.MergeError, match="capture"):
        ms.merge_redock([shard], str(tmp_path / "m2"))


def test_redock_selection_override_and_nonempty_outdir(tmp_path, fake_redock):
    shards = _redock_shards(tmp_path, n=2)
    out = str(tmp_path / "m")
    ms.merge_redock(shards, out, selection="sample 40 seed 42 of casf2016_core_set.csv")
    summary = json.loads(_read(os.path.join(out, "summary.json")))
    assert summary["protocol"]["selection"] == "sample 40 seed 42 of casf2016_core_set.csv"
    with pytest.raises(ms.MergeError, match="not empty"):
        ms.merge_redock(shards, out)


def test_unavailable_engine_summary_is_kept(tmp_path):
    shard = tmp_path / "easydock"
    shard.mkdir()
    ms.write_json(str(shard / "summary.json"),
                  {"engine": "easydock", "available": False, "reason": "not installed",
                   "tools": {}})
    report = ms.merge_redock([str(shard)], str(tmp_path / "m"), suite="baseline-redock")
    assert report["unavailable"] == [{"shard": str(shard), "reason": "not installed"}]
    assert _read(str(tmp_path / "m" / "summary.json")) == _read(str(shard / "summary.json"))


# ---------------------------------------------------------------------------
# J2 enrichment
# ---------------------------------------------------------------------------

TARGETS = ["akt1", "braf", "cdk2", "hivpr"]


def _enrichment_record(target):
    seed = sum(map(ord, target))
    return {"target": target, "family": be.TARGET_FAMILY.get(target, "?"), "status": "ok",
            "n_actives": 50, "n_decoys": 500, "box_source": "crystal", "error": None,
            "box_center": [1.5, -2.25, 3.125], "box_size": 22.5,
            "auc": 0.5 + seed / 10007.0, "auc_lo": 0.4, "auc_hi": 0.9,
            "ef5": 3.0 + seed / 997.0, "ef5_max": 11.0, "p_mw_one": 1.234e-7,
            "warnings": [f"{target}: EF1% uses only n_top=6 ligands"]}


@pytest.fixture
def fake_enrichment(monkeypatch):
    def process_target(target, workdir, opts):
        d = os.path.join(opts["outdir"], target)
        os.makedirs(d, exist_ok=True)
        for name in ("ligands.csv", "metrics.csv", "ranking.csv"):
            with open(os.path.join(d, name), "w", encoding="utf-8") as f:
                f.write(f"{target},{name}\n")
        return _enrichment_record(target)
    monkeypatch.setattr(be, "process_target", process_target)


def _run_enrichment(tmp_path, name, targets, *extra):
    outdir = str(tmp_path / name)
    assert be.main(["--targets", ",".join(targets), "--workdir", str(tmp_path / "dude"),
                    "--outdir", outdir, "--processes", "4", "--exhaustiveness", "4",
                    "--seed", "42", *extra]) == 0
    return outdir


def test_enrichment_merge_equals_single_run_in_any_order(tmp_path, fake_enrichment):
    single = _run_enrichment(tmp_path, "single", TARGETS)
    shards = [_run_enrichment(tmp_path, f"s_{t}", [t]) for t in TARGETS]
    trees = []
    for k, order in enumerate((shards, shards[::-1])):
        out = str(tmp_path / f"merged{k}")
        assert ms.main(["enrichment", "--shards", *order, "--outdir", out,
                        "--expected-ids", ",".join(TARGETS)]) == 0
        trees.append(_tree(out))
    assert trees[0] == trees[1] == _tree(single)


def test_enrichment_protocol_mismatch_and_missing_target(tmp_path, fake_enrichment):
    a = _run_enrichment(tmp_path, "a", ["akt1"])
    b = _run_enrichment(tmp_path, "b", ["braf"], "--max-decoys", "100")
    with pytest.raises(ms.MergeError, match="max_decoys"):
        ms.merge_enrichment([a, b], str(tmp_path / "m1"))
    report = ms.merge_enrichment([a], str(tmp_path / "m2"), expected=["akt1", "braf"])
    assert report["missing"] == ["braf"] and report["exit_status"] == ms.EXIT_MISSING


def test_enrichment_processes_may_differ(tmp_path, fake_enrichment):
    a = _run_enrichment(tmp_path, "a", ["akt1"])
    b = _run_enrichment(tmp_path, "b", ["braf"], "--processes", "2")
    ms.merge_enrichment([a, b], str(tmp_path / "m"))
    summary = json.loads(_read(str(tmp_path / "m" / "summary.json")))
    assert summary["protocol"]["processes"] == [2, 4]
    assert summary["protocol"]["targets"] == ["akt1", "braf"]


# ---------------------------------------------------------------------------
# J3 parity
# ---------------------------------------------------------------------------

def _parity_record(pdbid):
    k = int(pdbid[1:], 16)
    status = "error" if k % 5 == 0 else ("mismatch" if k % 7 == 0 else "ok")
    rec = {"pdbid": pdbid, "status": status, "modes": 9, "matched": 9 if status == "ok" else 8,
           "max_diff": 0.0, "ranking_affinity": "-7.1", "raw_mode1": "-7.1",
           "ranking_numeric_match": True, "ranking_string_match": True,
           "output_identical": k % 2 == 0, "binary": "/usr/local/bin/vina",
           "error": None if status == "ok" else "per-mode affinities, differ"}
    if status == "error":
        rec.update({"modes": 0, "matched": 0, "binary": None, "error": "no *_vina.log found"})
    return rec


PARITY_IDS = [f"{k:04x}" for k in range(0x1a0, 0x1a0 + 12)]


@pytest.fixture
def fake_parity(monkeypatch, tmp_path):
    work = tmp_path / "work"
    for p in PARITY_IDS:
        (work / p).mkdir(parents=True)

    def check_complex(workdir, pdbid, tmpdir, timeout, vina_bin=None):
        with open(os.path.join(tmpdir, f"{pdbid}_raw.pdbqt"), "w", encoding="utf-8") as f:
            f.write(f"MODEL 1\nREMARK VINA RESULT: -7.1 0 0\n{pdbid}\nENDMDL\n")
        return _parity_record(pdbid)
    monkeypatch.setattr(bp, "check_complex", check_complex)
    return str(work)


def _run_parity(tmp_path, work, name, ids):
    outdir = str(tmp_path / name)
    bp.main(["--workdir", work, "--pdbids", ",".join(ids), "--outdir", outdir])
    return outdir


def test_parity_merge_equals_single_run_in_any_order(tmp_path, fake_parity):
    single = _run_parity(tmp_path, fake_parity, "single", PARITY_IDS)
    shards = [_run_parity(tmp_path, fake_parity, f"s{k}", PARITY_IDS[k::5]) for k in range(5)]
    trees = []
    for k, order in enumerate((shards, shards[::-1])):
        out = str(tmp_path / f"m{k}")
        ms.merge_parity(order, out, expected=PARITY_IDS)
        trees.append(_tree(out))
    assert trees[0] == trees[1] == _tree(single)
    summary = json.loads(trees[0]["summary.json"])
    assert summary["result"] == "FAIL" and summary["checked"] == len(PARITY_IDS)


# ---------------------------------------------------------------------------
# Split
# ---------------------------------------------------------------------------

def test_split_is_sorted_round_robin_without_excluded_ids():
    ids = [f"{k}x{k:02d}" for k in range(1, 10)] + ["1ABC", " 1abc "]
    got, sample, shards = ms.split_ids(ids, "full", ["2x02"], ["3x03"], 4, sample=3, seed=42)
    pool = sorted({i.strip().lower() for i in ids} - {"2x02", "3x03"})
    assert got == pool
    assert shards == [pool[k::4] for k in range(4)]
    assert sorted(sum(shards, [])) == pool
    # The sample is drawn from every ID, so it equals the protocol's sample.
    assert set(sample) <= set(pool)
    assert ms.split_ids(list(reversed(ids)), "full", ["2x02"], ["3x03"], 4,
                        sample=3, seed=42) == (got, sample, shards)
    # More shards than IDs: one ID per shard.
    assert len(ms.split_ids(ids, "full", [], [], 100)[2]) == len(set(i.strip().lower()
                                                                    for i in ids))


def test_split_checks_the_fetch_manifest(tmp_path):
    ids_csv = tmp_path / "core.csv"
    ids_csv.write_text("pdbid,smiles\n" + "".join(f"{k}aa{k},C\n" for k in range(1, 9)),
                       encoding="utf-8")
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("pdbid,status,reason\n1aa1,ok,\n2aa2,skipped,peptide\n3aa3,ok,\n",
                        encoding="utf-8")
    base = ["split", "--ids-csv", str(ids_csv), "--set", "full", "--exclude", "3aa3",
            "--manifest", str(manifest), "--shards", "3"]
    out = tmp_path / "plan"
    assert ms.main(base + ["--skipped", "2aa2", "--outdir", str(out)]) == 0
    plan = json.loads((out / "plan.json").read_text(encoding="utf-8"))
    assert plan["n_ids"] == 6 and plan["shards"] == ["000", "001", "002"]
    assert (out / "shard_000.txt").read_text(encoding="utf-8").split() == ["1aa1", "6aa6"]
    assert ms.main(base + ["--skipped", "4aa4", "--outdir", str(tmp_path / "p2")]) == 1


def test_report_collects_shards_and_hashes(tmp_path, fake_redock):
    shards_root = tmp_path / "artifacts"
    shard = shards_root / "casf-shard-000"
    shard.mkdir(parents=True)
    ms.write_json(str(shard / "shard-info.json"), {"cpu_model": "AMD EPYC 7763", "nproc": 4})
    shards = _redock_shards(tmp_path, n=2)
    rep = str(tmp_path / "redock.json")
    assert ms.main(["redock", "--shards", *shards, "--outdir", str(tmp_path / "res"),
                    "--report", rep]) == 0
    out = str(tmp_path / "benchmarks_summary.json")
    assert ms.main(["report", "--out", out, "--shards-root", str(shards_root),
                    "--merge-report", rep, "--meta", "git_sha=abc"]) == 0
    summary = json.loads(_read(out))
    assert summary["run"] == {"git_sha": "abc"}
    assert summary["runners"] == {"cpu_models": ["AMD EPYC 7763"], "nproc": [4]}
    hashes = summary["suites"]["redock"]["outputs_sha256"]
    assert set(hashes) == {"results.csv", "summary.json"}
    assert summary["complete"] is True
