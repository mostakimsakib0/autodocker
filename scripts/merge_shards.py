#!/usr/bin/env python3
"""Split the benchmark sets into shards and merge sharded benchmark outputs.

The GitHub Actions workflow ``.github/workflows/benchmarks.yml`` runs the
J1/J2/J3/J5 benchmarks in shards (CASF complexes round-robin over N jobs,
one job per DUD-E target). This script is both ends of that:

``split``
    The CASF-2016 complexes of a run (``--set full``: every core-set ID,
    ``--set sample``: the protocol's seeded sample of
    ``BENCH_SAMPLE_SIZE`` IDs with ``BENCH_SEED``), minus ``--exclude``
    (documented faulty entries) and ``--skipped`` (IDs fetch_casf_rcsb.py
    cannot build). The IDs are sorted and dealt round-robin: shard k of n
    gets ``ids[k::n]``. Writes ``ids.txt``, ``sample.txt`` (the seeded
    sample, the J1 crystal-start / J3 / J5 set), ``shard_KKK.txt`` and
    ``plan.json``. With ``--manifest`` (fetch_casf_rcsb.py's manifest.csv)
    the IDs the fetch skipped must be exactly ``--skipped``.

``capture {redock,baseline} -- ARGS``
    Runs ``benchmark_redock.py ARGS`` (or ``benchmark_baseline.py ARGS``)
    unchanged and also writes ``OUTDIR/records.json``: the per-complex
    records the script's ``summarize`` received, at full precision.
    ``results.csv`` rounds RMSDs to 3 decimals, and a merged
    ``summary.json`` must be computed from the exact values (medians,
    ``rmsd <= cutoff``) to equal a single run's.

``redock | baseline-redock | enrichment | baseline-enrichment | parity``
    ``--shards DIR... --outdir OUT``: merge shard output directories (each
    one benchmark ``--outdir``) into the files a single run over the union
    of the shards writes, using the scripts' own code
    (``benchmark_redock.summarize``, ``benchmark_enrichment.save_tables``,
    ``benchmark_vina_parity.verdict``). The shards must share the protocol
    (every protocol field equal except the ID/target list, which is
    merged, and ``processes``, which does not change results); an ID or
    target in two shards, a shard whose records do not reproduce its own
    summary, or a protocol difference is an error (exit 1). Rows and
    records are ordered by ID (enrichment: by ``--targets`` order, default
    benchmark_enrichment.DEFAULT_TARGETS), never by shard order, so any
    order of ``--shards`` gives byte-identical files. ``--expected-ids``
    names the full set: IDs no shard produced are reported and the exit
    status is 3 (outputs are still written) unless ``--allow-missing``.

``report``
    Writes the top-level ``benchmarks_summary.json`` of a workflow run:
    run metadata, the image record of the build job, the input checksum
    list, per-shard runner records and the merge reports with the SHA-256
    of every merged file.

Exit status: 0 ok, 1 inconsistent shards, 2 usage error, 3 IDs missing.
"""

import argparse
import csv
import datetime
import glob
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import _benchlib  # noqa: E402
import benchmark_enrichment as be  # noqa: E402
import benchmark_redock as br  # noqa: E402
import benchmark_vina_parity as bp  # noqa: E402

EXIT_OK = 0
EXIT_INCONSISTENT = 1
EXIT_USAGE = 2
EXIT_MISSING = 3

RECORDS_JSON = "records.json"
SUMMARY_JSON = "summary.json"
RESULTS_CSV = "results.csv"
PARITY_CSV = "parity.csv"
ENRICHMENT_CSV = "enrichment.csv"
WORK_DIR = "work"
PLAN_JSON = "plan.json"
# Protocol fields that differ between shards by construction: merged.
REDOCK_LIST_KEY = "pdbids"
ENRICHMENT_LIST_KEY = "targets"
# Protocol fields that do not change results (parallel workers); kept when
# all shards agree, else the sorted list of the values.
FREE_KEYS = ("processes",)
REDOCK_SUITES = ("redock", "baseline-redock")
ENRICHMENT_SUITES = ("enrichment", "baseline-enrichment")
SUITES = REDOCK_SUITES + ENRICHMENT_SUITES + ("parity",)


class MergeError(Exception):
    """The shards cannot be merged into one consistent result."""


# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------

def read_json(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def write_json(path, obj):
    """As the benchmark scripts write summary.json (indent 2, no newline)."""
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)


def read_csv_rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    if not rows:
        raise MergeError(f"{path}: empty CSV (no header)")
    return rows[0], rows[1:]


def write_csv_rows(path, header, rows):
    """csv.writer defaults, as the scripts write their tables."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def read_ids(spec):
    """IDs from a file (whitespace/comma separated) or a comma list."""
    if spec is None:
        return None
    text = spec
    if os.path.exists(spec):
        with open(spec, encoding="utf-8") as f:
            text = f.read()
    return sorted({_benchlib.normalize_pdbid(t) for t in text.replace(",", " ").split()} - {""})


def write_ids(path, ids):
    with open(path, "w", encoding="utf-8") as f:
        f.write("".join(f"{i}\n" for i in ids))


def file_sha256(path):
    return _benchlib.file_sha256(path)


def prepare_outdir(outdir):
    if os.path.isdir(outdir) and os.listdir(outdir):
        raise MergeError(f"--outdir {outdir} is not empty")
    os.makedirs(outdir, exist_ok=True)


def copy_tree(src, dst):
    if os.path.exists(dst):
        raise MergeError(f"{dst} exists (the same directory in two shards?)")
    shutil.copytree(src, dst, symlinks=True)


# ---------------------------------------------------------------------------
# Protocols
# ---------------------------------------------------------------------------

def merge_protocols(protocols, list_key, order=None):
    """One protocol from the shards' protocols ``[(shard, protocol)]``.

    Every field must be equal across the shards except ``list_key`` (the
    union, ordered by ``order`` or sorted) and FREE_KEYS. Raises MergeError
    naming the differing fields and shards."""
    if not protocols:
        return None
    keys = []
    for _, proto in protocols:
        for k in proto:
            if k not in keys:
                keys.append(k)
    merged, problems = {}, []
    for k in keys:
        values = [(shard, proto.get(k, "<absent>")) for shard, proto in protocols]
        if k == list_key:
            union = set()
            for _, v in values:
                union.update(v or [])
            merged[k] = order(union) if order else sorted(union)
            continue
        distinct = []
        for _, v in values:
            if v not in distinct:
                distinct.append(v)
        if len(distinct) == 1:
            merged[k] = distinct[0]
        elif k in FREE_KEYS:
            merged[k] = sorted(distinct, key=lambda v: json.dumps(v, sort_keys=True))
        else:
            shown = sorted({json.dumps(v, sort_keys=True)[:120] for _, v in values})
            problems.append(f"{k}: {' | '.join(shown)}")
    if problems:
        raise MergeError("shards ran different protocols: " + "; ".join(problems))
    return merged


def _expected_report(suite, got, expected, allow_missing):
    report = {"suite": suite, "n_merged": len(got)}
    if expected is not None:
        unexpected = sorted(set(got) - set(expected))
        if unexpected:
            raise MergeError(f"IDs not in --expected-ids: {', '.join(unexpected)}")
        missing = sorted(set(expected) - set(got))
        report.update({"n_expected": len(expected), "missing": missing,
                       "complete": not missing})
    else:
        report.update({"missing": [], "complete": None})
    report["exit_status"] = (EXIT_MISSING if report["missing"] and not allow_missing
                             else EXIT_OK)
    return report


def _unavailable(summary):
    """benchmark_baseline.py's summary when its engine could not run."""
    return isinstance(summary, dict) and summary.get("available") is False


def keep_unavailable_summary(unavailable, outdir):
    """When no shard could run its engine, the merged summary.json is one
    shard's (the same one whatever the shard order)."""
    shard = min(unavailable, key=lambda su: json.dumps(su[1], sort_keys=True))[0]
    shutil.copyfile(os.path.join(shard, SUMMARY_JSON), os.path.join(outdir, SUMMARY_JSON))


def _shard_order(shards):
    """Shards in a fixed order (only error messages depend on it)."""
    return sorted(shards, key=lambda s: os.path.abspath(s))


# ---------------------------------------------------------------------------
# J1 / J5 re-docking
# ---------------------------------------------------------------------------

def merge_redock(shards, outdir, expected=None, selection=None, with_work=False,
                 allow_missing=False, suite="redock"):
    """Merge benchmark_redock.py (or benchmark_baseline.py redock) shard
    outdirs; see the module docstring."""
    prepare_outdir(outdir)
    header, rows, records, protocols, cutoffs = None, {}, {}, [], set()
    owner, unavailable = {}, []
    for shard in _shard_order(shards):
        summary = read_json(os.path.join(shard, SUMMARY_JSON))
        if _unavailable(summary):
            unavailable.append((shard, summary))
            continue
        shard_header, shard_rows = read_csv_rows(os.path.join(shard, RESULTS_CSV))
        if header is None:
            header = shard_header
        elif shard_header != header:
            raise MergeError(f"{shard}: results.csv columns differ from the other shards")
        rec_path = os.path.join(shard, RECORDS_JSON)
        if not os.path.exists(rec_path):
            raise MergeError(f"{shard}: no {RECORDS_JSON} (run the shard through "
                             "'merge_shards.py capture')")
        shard_records = read_json(rec_path)
        cutoff = summary.get("rmsd_cutoff")
        expected_summary = br.summarize(shard_records, cutoff)
        if {k: v for k, v in summary.items() if k != "protocol"} != expected_summary:
            raise MergeError(f"{shard}: {RECORDS_JSON} does not reproduce {SUMMARY_JSON}")
        row_ids = [r[0] for r in shard_rows]
        rec_ids = [r["pdbid"] for r in shard_records]
        proto_ids = list(summary.get("protocol", {}).get(REDOCK_LIST_KEY) or [])
        if sorted(row_ids) != sorted(rec_ids) or sorted(rec_ids) != sorted(proto_ids):
            raise MergeError(f"{shard}: results.csv, {RECORDS_JSON} and the protocol's "
                             "pdbids list different complexes")
        for pdbid, row in zip(row_ids, shard_rows):
            if pdbid in owner:
                raise MergeError(f"{pdbid} is in two shards: {owner[pdbid]} and {shard}")
            owner[pdbid] = shard
            rows[pdbid] = row
        for rec in shard_records:
            records[rec["pdbid"]] = rec
        protocols.append((shard, summary.get("protocol", {})))
        cutoffs.add(json.dumps(cutoff))

    report = {"unavailable": [{"shard": s, "reason": u.get("reason")}
                              for s, u in unavailable]}
    if not protocols:
        if not unavailable:
            raise MergeError("no shard to merge")
        # Every shard says why its engine could not run: keep that summary.
        keep_unavailable_summary(unavailable, outdir)
        report.update(_expected_report(suite, [], expected, True))
        report["exit_status"] = EXIT_OK
        return report
    if len(cutoffs) != 1:
        raise MergeError(f"shards used different RMSD cutoffs: {sorted(cutoffs)}")
    protocol = merge_protocols(protocols, REDOCK_LIST_KEY)
    if selection is not None:
        protocol["selection"] = selection
    ids = sorted(records)
    write_csv_rows(os.path.join(outdir, RESULTS_CSV), header, [rows[i] for i in ids])
    ordered = [records[i] for i in ids]
    summary = br.summarize(ordered, json.loads(cutoffs.pop()))
    summary["protocol"] = protocol
    write_json(os.path.join(outdir, SUMMARY_JSON), summary)
    if with_work:
        for pdbid in ids:
            src = os.path.join(owner[pdbid], WORK_DIR, pdbid)
            if os.path.isdir(src):
                copy_tree(src, os.path.join(outdir, WORK_DIR, pdbid))
    report.update(_expected_report(suite, ids, expected, allow_missing))
    return report


# ---------------------------------------------------------------------------
# J2 / J5 enrichment
# ---------------------------------------------------------------------------

def target_order(order=None):
    """Sort key function: targets in ``order`` first (its order), others
    after them by name."""
    order = list(order or be.DEFAULT_TARGETS)
    pos = {t: i for i, t in enumerate(order)}

    def key(t):
        return (pos.get(t, len(pos)), t)
    return lambda targets: sorted(targets, key=key)


def merge_enrichment(shards, outdir, expected=None, order=None, with_work=False,
                     allow_missing=False, suite="enrichment"):
    """Merge benchmark_enrichment.py (or benchmark_baseline.py enrichment)
    shard outdirs; see the module docstring."""
    prepare_outdir(outdir)
    sort_targets = target_order(order)
    records, owner, protocols, unavailable = {}, {}, [], []
    constants = set()
    for shard in _shard_order(shards):
        summary = read_json(os.path.join(shard, SUMMARY_JSON))
        if _unavailable(summary):
            unavailable.append((shard, summary))
            continue
        proto = summary.get("protocol", {})
        shard_targets = [r["target"] for r in summary.get("targets", [])]
        if sorted(shard_targets) != sorted(proto.get(ENRICHMENT_LIST_KEY) or []):
            raise MergeError(f"{shard}: summary.json targets differ from its protocol")
        for rec in summary["targets"]:
            t = rec["target"]
            if t in owner:
                raise MergeError(f"target {t} is in two shards: {owner[t]} and {shard}")
            owner[t] = shard
            records[t] = rec
        protocols.append((shard, proto))
        constants.add(json.dumps([summary.get("primary_metrics"),
                                  summary.get("secondary_metrics_caveats")], sort_keys=True))
    report = {"unavailable": [{"shard": s, "reason": u.get("reason")}
                              for s, u in unavailable]}
    if not protocols:
        if not unavailable:
            raise MergeError("no shard to merge")
        keep_unavailable_summary(unavailable, outdir)
        report.update(_expected_report(suite, [], expected, True))
        report["exit_status"] = EXIT_OK
        return report
    expected_constants = json.dumps([list(be.PRIMARY_METRICS), be.SECONDARY_METRIC_CAVEATS],
                                    sort_keys=True)
    if constants != {expected_constants}:
        raise MergeError("shards were scored with metric definitions other than "
                         "this checkout's benchmark_enrichment.py")
    protocol = merge_protocols(protocols, ENRICHMENT_LIST_KEY, order=sort_targets)
    targets = sort_targets(records)
    for t in targets:
        src = os.path.join(owner[t], t)
        if os.path.isdir(src):
            copy_tree(src, os.path.join(outdir, t))
        work = os.path.join(owner[t], WORK_DIR, t)
        if with_work and os.path.isdir(work):
            copy_tree(work, os.path.join(outdir, WORK_DIR, t))
    be.save_tables([records[t] for t in targets], outdir, protocol)
    report.update(_expected_report(suite, targets, expected, allow_missing))
    return report


# ---------------------------------------------------------------------------
# J3 parity
# ---------------------------------------------------------------------------

def parity_record(header, row):
    """The fields of a parity.csv row that its summary counts."""
    rec = dict(zip(header, row))

    def as_int(v):
        return int(v) if v not in ("", None) else 0
    return {"pdbid": rec["pdbid"], "status": rec["status"],
            "modes": as_int(rec.get("modes")), "matched": as_int(rec.get("matched")),
            "output_identical": rec.get("output_identical") == "True",
            "binary": rec.get("binary") or None}


def parity_summary(records, binaries):
    """benchmark_vina_parity.py's summary.json for ``records``; ``binaries``
    maps each binary path to the version its shard recorded."""
    used = sorted({r["binary"] for r in records if r.get("binary")})
    return {
        "checked": len(records),
        "pass": sum(1 for r in records if r["status"] == "ok"),
        "mismatch": sum(1 for r in records if r["status"] == "mismatch"),
        "error": sum(1 for r in records if r["status"] == "error"),
        "modes_compared": sum(r["modes"] for r in records),
        "modes_matched": sum(r["matched"] for r in records),
        "outputs_byte_identical": sum(1 for r in records if r.get("output_identical")),
        "binaries": {b: binaries.get(b) for b in used},
        "result": "PASS" if bp.verdict(records) else "FAIL",
    }


def merge_parity(shards, outdir, expected=None, allow_missing=False, suite="parity"):
    """Merge benchmark_vina_parity.py shard outdirs."""
    prepare_outdir(outdir)
    header, rows, records, owner, binaries = None, {}, {}, {}, {}
    for shard in _shard_order(shards):
        summary = read_json(os.path.join(shard, SUMMARY_JSON))
        shard_header, shard_rows = read_csv_rows(os.path.join(shard, PARITY_CSV))
        if header is None:
            header = shard_header
        elif shard_header != header:
            raise MergeError(f"{shard}: parity.csv columns differ from the other shards")
        shard_records = [parity_record(header, r) for r in shard_rows]
        for b, version in (summary.get("binaries") or {}).items():
            if b in binaries and binaries[b] != version:
                raise MergeError(f"{b}: version {binaries[b]!r} in one shard, "
                                 f"{version!r} in {shard}")
            binaries[b] = version
        if parity_summary(shard_records, binaries) != summary:
            raise MergeError(f"{shard}: parity.csv does not reproduce {SUMMARY_JSON}")
        for row, rec in zip(shard_rows, shard_records):
            if rec["pdbid"] in owner:
                raise MergeError(f"{rec['pdbid']} is in two shards: "
                                 f"{owner[rec['pdbid']]} and {shard}")
            owner[rec["pdbid"]] = shard
            rows[rec["pdbid"]] = row
            records[rec["pdbid"]] = rec
    if header is None:
        raise MergeError("no shard to merge")
    ids = sorted(records)
    write_csv_rows(os.path.join(outdir, PARITY_CSV), header, [rows[i] for i in ids])
    write_json(os.path.join(outdir, SUMMARY_JSON),
               parity_summary([records[i] for i in ids], binaries))
    for pdbid in ids:
        raw = os.path.join(owner[pdbid], f"{pdbid}_raw.pdbqt")
        if os.path.exists(raw):
            shutil.copy2(raw, os.path.join(outdir, os.path.basename(raw)))
    return {**_expected_report(suite, ids, expected, allow_missing), "unavailable": []}


# ---------------------------------------------------------------------------
# Split
# ---------------------------------------------------------------------------

def manifest_skipped(manifest_csv):
    """IDs fetch_casf_rcsb.py could not build (status != ok)."""
    with open(manifest_csv, newline="", encoding="utf-8") as f:
        return sorted({_benchlib.normalize_pdbid(r["pdbid"]) for r in csv.DictReader(f)
                       if r.get("status") != "ok"})


def split_ids(all_ids, which, exclude, skipped, n_shards,
              sample=_benchlib.BENCH_SAMPLE_SIZE, seed=_benchlib.BENCH_SEED):
    """``(ids, sample_ids, shards)`` of a run (see the module docstring)."""
    pool = sorted({_benchlib.normalize_pdbid(i) for i in all_ids} - {""})
    drop = set(exclude) | set(skipped)
    sample_ids = [i for i in _benchlib.select_complexes(pool, sample, seed) if i not in drop]
    ids = sample_ids if which == "sample" else [i for i in pool if i not in drop]
    if n_shards < 1:
        raise MergeError("--shards must be >= 1")
    n = max(1, min(n_shards, len(ids)))
    return ids, sample_ids, [ids[k::n] for k in range(n)]


def run_split(args):
    rows = _benchlib.read_pdbid_rows(args.ids_csv)
    exclude, skipped = read_ids(args.exclude) or [], read_ids(args.skipped) or []
    ids, sample_ids, shards = split_ids(rows, args.set, exclude, skipped, args.shards,
                                        args.sample, args.sample_seed)
    if args.manifest:
        selected = set(ids) | set(skipped) | set(exclude)
        if args.set == "sample":
            selected = set(_benchlib.select_complexes(rows, args.sample, args.sample_seed))
        got = sorted(set(manifest_skipped(args.manifest)) & selected)
        want = sorted(set(skipped) & selected)
        if got != want:
            raise MergeError(f"{args.manifest}: fetch skipped {got or 'nothing'}, "
                             f"expected {want or 'nothing'} (--skipped)")
    os.makedirs(args.outdir, exist_ok=True)
    write_ids(os.path.join(args.outdir, "ids.txt"), ids)
    write_ids(os.path.join(args.outdir, "sample.txt"), sample_ids)
    names = []
    for k, shard in enumerate(shards):
        names.append(f"{k:03d}")
        write_ids(os.path.join(args.outdir, f"shard_{k:03d}.txt"), shard)
    plan = {"set": args.set, "ids_csv": os.path.basename(args.ids_csv),
            "sample": args.sample, "sample_seed": args.sample_seed,
            "excluded": exclude, "skipped": skipped, "n_ids": len(ids),
            "n_sample": len(sample_ids), "n_shards": len(shards), "shards": names,
            "selection": (f"sample {args.sample} seed {args.sample_seed} of "
                          f"{os.path.basename(args.ids_csv)}")}
    write_json(os.path.join(args.outdir, PLAN_JSON), plan)
    print(json.dumps(plan, separators=(",", ":")))
    return EXIT_OK


# ---------------------------------------------------------------------------
# Capture
# ---------------------------------------------------------------------------

def _argv_outdir(argv, default):
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--outdir", default=default)
    known, _ = ap.parse_known_args(argv)
    return known.outdir


def capture_redock(argv, module):
    """Run ``module.main(argv)`` (benchmark_redock.py, or the
    benchmark_redock module that benchmark_baseline.py uses) with its
    summarize() wrapped, and write the records it received."""
    captured = []
    original = module.summarize

    def summarize(results, cutoff):
        captured.append([dict(r) for r in results])
        return original(results, cutoff)

    module.summarize = summarize
    return captured, original


def run_capture(script, argv):
    if argv and argv[0] == "--":
        argv = argv[1:]
    if script == "redock":
        target, main, default = br, br.main, "benchmarks/results/redock"
    else:
        import benchmark_baseline as bb
        target, main, default = bb.br, bb.main, None
    captured, original = capture_redock(argv, target)
    try:
        status = main(argv)
    finally:
        target.summarize = original
    if captured:
        outdir = _argv_outdir(argv, default)
        with open(os.path.join(outdir, RECORDS_JSON), "w", encoding="utf-8") as f:
            json.dump(captured[-1], f, indent=1)
    return status


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def output_hashes(root):
    """``{relative path: sha256}`` of every file below ``root``, sorted."""
    out = {}
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        for name in sorted(filenames):
            path = os.path.join(dirpath, name)
            out[os.path.relpath(path, root).replace(os.sep, "/")] = file_sha256(path)
    return dict(sorted(out.items()))


def run_report(args):
    def maybe_json(path):
        return read_json(path) if path and os.path.exists(path) else None

    meta = {}
    for item in args.meta:
        key, _, value = item.partition("=")
        meta[key] = value
    shards = []
    for path in sorted(glob.glob(os.path.join(args.shards_root, "*", "shard-info.json"))):
        info = read_json(path)
        info.setdefault("artifact", os.path.basename(os.path.dirname(path)))
        shards.append(info)
    shards.sort(key=lambda s: s.get("artifact", ""))
    suites = {}
    for path in sorted(args.merge_report):
        rep = read_json(path)
        outdir = rep.get("outdir")
        if outdir and os.path.isdir(outdir):
            rep["outputs_sha256"] = output_hashes(outdir)
        suites[rep.get("name") or rep.get("suite")] = rep
    inputs = None
    if args.inputs_sums and os.path.exists(args.inputs_sums):
        with open(args.inputs_sums, encoding="utf-8") as f:
            n_files = sum(1 for line in f if line.strip())
        inputs = {"file": os.path.basename(args.inputs_sums), "n_files": n_files,
                  "sha256": file_sha256(args.inputs_sums)}
    summary = {
        "generated_at": datetime.datetime.now(datetime.timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"),
        "run": meta,
        "image": maybe_json(args.build_info),
        "inputs": inputs,
        "plan": maybe_json(args.plan),
        "complete": all(s.get("complete") is not False and not s.get("error")
                        for s in suites.values()),
        "suites": suites,
        "runners": {
            "cpu_models": sorted({s.get("cpu_model") for s in shards if s.get("cpu_model")}),
            "nproc": sorted({s.get("nproc") for s in shards if s.get("nproc") is not None}),
        },
        "shards": shards,
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")
    return EXIT_OK


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for suite in SUITES:
        p = sub.add_parser(suite, help=f"merge {suite} shard outdirs")
        p.add_argument("--shards", nargs="+", required=True, metavar="DIR",
                       help="Shard output directories (one benchmark --outdir each)")
        p.add_argument("--outdir", required=True, help="Merged output (new or empty)")
        p.add_argument("--expected-ids", default=None, metavar="FILE|LIST",
                       help="Every ID (targets for enrichment) the shards should cover")
        p.add_argument("--allow-missing", action="store_true",
                       help="Exit 0 even when --expected-ids are missing (still reported)")
        p.add_argument("--report", default=None, metavar="FILE",
                       help="Write the merge report (JSON) here")
        p.add_argument("--name", default=None,
                       help="Name of the merged result in the report (default: the suite)")
        if suite in REDOCK_SUITES:
            p.add_argument("--selection", default=None,
                           help="protocol.selection of the merged summary (the selection "
                                "a single run would record, e.g. 'sample 40 seed 42 of "
                                "casf2016_core_set.csv'); default: the shards' value")
        if suite in ENRICHMENT_SUITES:
            p.add_argument("--targets", default=",".join(be.DEFAULT_TARGETS),
                           help="Target order of the merged tables (default: the "
                                "benchmark's default order)")
        if suite != "parity":
            p.add_argument("--with-work", action="store_true",
                           help="Also copy the shards' work/<id> directories")

    sp = sub.add_parser("split", help="CASF IDs of a run, dealt round-robin into shards")
    sp.add_argument("--ids-csv", default=br.DEFAULT_IDS)
    sp.add_argument("--set", choices=["full", "sample"], default="full")
    sp.add_argument("--sample", type=int, default=_benchlib.BENCH_SAMPLE_SIZE)
    sp.add_argument("--sample-seed", type=int, default=_benchlib.BENCH_SEED)
    sp.add_argument("--exclude", default="", help="IDs left out on purpose (comma list)")
    sp.add_argument("--skipped", default="",
                    help="IDs fetch_casf_rcsb.py cannot build (comma list)")
    sp.add_argument("--manifest", default=None,
                    help="fetch_casf_rcsb.py manifest.csv: its skipped IDs must equal --skipped")
    sp.add_argument("--shards", type=int, default=20)
    sp.add_argument("--outdir", required=True)

    cp = sub.add_parser("capture", help="run a re-docking benchmark and keep its records")
    cp.add_argument("script", choices=["redock", "baseline"])
    cp.add_argument("args", nargs=argparse.REMAINDER,
                    help="-- then the script's own arguments")

    rp = sub.add_parser("report", help="write benchmarks_summary.json")
    rp.add_argument("--out", required=True)
    rp.add_argument("--shards-root", required=True,
                    help="Directory with one sub-directory per downloaded shard artifact")
    rp.add_argument("--merge-report", nargs="*", default=[])
    rp.add_argument("--build-info", default=None)
    rp.add_argument("--inputs-sums", default=None)
    rp.add_argument("--plan", default=None)
    rp.add_argument("--meta", action="append", default=[], metavar="KEY=VALUE")
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.cmd == "split":
            return run_split(args)
        if args.cmd == "capture":
            return run_capture(args.script, args.args)
        if args.cmd == "report":
            return run_report(args)
        expected = read_ids(args.expected_ids)
        common = {"expected": expected, "allow_missing": args.allow_missing,
                  "suite": args.cmd}
        if args.cmd in REDOCK_SUITES:
            report = merge_redock(args.shards, args.outdir, selection=args.selection,
                                  with_work=args.with_work, **common)
        elif args.cmd in ENRICHMENT_SUITES:
            report = merge_enrichment(args.shards, args.outdir,
                                      order=[t for t in args.targets.split(",") if t],
                                      with_work=args.with_work, **common)
        else:
            report = merge_parity(args.shards, args.outdir, **common)
    except MergeError as e:
        print(f"[x] {e}", file=sys.stderr)
        if getattr(args, "report", None):
            write_json(args.report, {"suite": args.cmd, "name": args.name or args.cmd,
                                     "error": str(e), "complete": False,
                                     "exit_status": EXIT_INCONSISTENT})
        return EXIT_INCONSISTENT
    report["name"] = args.name or args.cmd
    report["outdir"] = args.outdir
    report["n_shards"] = len(args.shards)
    if report["missing"]:
        print(f"[!] {len(report['missing'])} expected ID(s) missing: "
              f"{', '.join(report['missing'])}", file=sys.stderr)
    if args.report:
        write_json(args.report, report)
    print(json.dumps({k: v for k, v in report.items() if k != "missing"}, sort_keys=True))
    return report["exit_status"]


if __name__ == "__main__":
    sys.exit(main())
