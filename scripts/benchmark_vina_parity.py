#!/usr/bin/env python3
"""J3: Vina parity check.

Checks that the autodocker pipeline reports exactly what Vina computed:

  1. For each complex docked by ``benchmark_redock.py``, read the Vina
     command the pipeline logged (last ``# Command:`` line of *_vina.log,
     i.e. the attempt that produced the output).
  2. Re-run that command with the raw docking binary, writing to a
     separate file.
  3. Compare mode by mode: every ``REMARK VINA RESULT`` affinity of the
     pipeline output must equal the raw run's (numeric equality).
  4. Compare what the user sees: the affinity in ``ranking.csv`` must equal
     the raw run's mode-1 ``REMARK VINA RESULT`` value, both as numbers and
     as strings at the precision Vina printed.

A complex PASSes only when all checks hold; a complex that could not be
checked (missing log, binary failure, ...) counts as FAIL, and a run with no
complexes checked FAILs. Whether the raw output file is byte-identical to the
pipeline's is recorded as well (informational).

The binary is the one in the logged command unless ``--vina-bin`` is given;
its ``--version`` is recorded.
"""

import argparse
import csv
import filecmp
import glob
import json
import os
import shlex
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _benchlib  # noqa: E402


def parse_vina_results(pdbqt_path):
    """Per-mode ``(token, value)`` pairs of a Vina output file."""
    return [(m["affinity_str"], m["affinity"])
            for m in _benchlib.parse_vina_models(pdbqt_path)]


def extract_command(log_path):
    """The last ``# Command:`` recorded in a *_vina.log (the attempt that
    wrote the output when a fallback binary was used)."""
    cmd = None
    with open(log_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("# Command:"):
                # The pipeline logs shlex.join(cmd); shlex.split inverts it
                # (paths with spaces or quotes survive).
                cmd = shlex.split(line.split(":", 1)[1].strip())
    return cmd


def run_vina(cmd, timeout=_benchlib.DOCKING_RUN_TIMEOUT):
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError(
            f"vina failed (rc={proc.returncode}): "
            f"{(proc.stdout or '')[-500:]} {(proc.stderr or '')[-500:]}")


def ranking_matches(ranking_text, remark_token):
    """``(numeric_equal, string_equal)`` for a ranking.csv affinity and a
    REMARK VINA RESULT token. The string check formats the ranking value at
    the number of decimals Vina printed (``-7.1`` vs ``-7.100`` is equal,
    ``-7.1484`` vs ``-7.148`` is not, numerically)."""
    try:
        rank_val, raw_val = float(ranking_text), float(remark_token)
    except (TypeError, ValueError):
        return False, False
    decimals = len(remark_token.split(".")[1]) if "." in remark_token else 0
    return rank_val == raw_val, f"{rank_val:.{decimals}f}" == remark_token


def check_complex(workdir, pdbid, tmpdir, timeout, vina_bin=None, runner=run_vina):
    rec = {"pdbid": pdbid, "status": "error", "modes": 0, "matched": 0,
           "max_diff": None, "ranking_affinity": None, "raw_mode1": None,
           "ranking_numeric_match": None, "ranking_string_match": None,
           "output_identical": None, "binary": None, "error": None}
    logs = sorted(glob.glob(os.path.join(workdir, _benchlib.DOCKED_DIR, "*" + _benchlib.VINA_LOG_SUFFIX)))
    if not logs:
        rec["error"] = "no *_vina.log found"
        return rec
    log = logs[0]
    ligand = os.path.basename(log)[:-len(_benchlib.VINA_LOG_SUFFIX)]
    cmd = extract_command(log)
    if not cmd:
        rec["error"] = "no command header in log"
        return rec
    if vina_bin:
        cmd = [vina_bin] + cmd[1:]
    rec["binary"] = cmd[0]

    pipeline_out = None
    for i, c in enumerate(cmd):
        if c == "--out" and i + 1 < len(cmd):
            pipeline_out = cmd[i + 1]
    local_out = os.path.join(workdir, _benchlib.DOCKED_DIR, ligand + _benchlib.POSE_SUFFIX)
    if not pipeline_out or not os.path.exists(pipeline_out):
        pipeline_out = local_out
    if not os.path.exists(pipeline_out):
        rec["error"] = "pipeline output missing"
        return rec

    raw_out = os.path.join(tmpdir, f"{pdbid}_raw.pdbqt")
    if os.path.exists(raw_out):
        os.remove(raw_out)
    raw_cmd = list(cmd)
    for i, c in enumerate(raw_cmd):
        if c == "--out" and i + 1 < len(raw_cmd):
            raw_cmd[i + 1] = raw_out
    try:
        runner(raw_cmd, timeout)
    except Exception as e:
        rec["error"] = str(e)
        return rec
    if not os.path.exists(raw_out):
        rec["error"] = "raw run wrote no output"
        return rec

    pipeline_aff = parse_vina_results(pipeline_out)
    raw_aff = parse_vina_results(raw_out)
    if not raw_aff or raw_aff[0][0] is None:
        rec["error"] = "no REMARK VINA RESULT in raw output"
        return rec
    if len(pipeline_aff) != len(raw_aff):
        rec["error"] = (f"mode count mismatch pipeline={len(pipeline_aff)} "
                        f"raw={len(raw_aff)}")
        rec["status"] = "mismatch"
        return rec

    rec["modes"] = len(raw_aff)
    rec["matched"] = sum(1 for a, b in zip(pipeline_aff, raw_aff)
                         if a[1] is not None and a[1] == b[1])
    diffs = [abs(a[1] - b[1]) for a, b in zip(pipeline_aff, raw_aff)
             if a[1] is not None and b[1] is not None]
    rec["max_diff"] = max(diffs) if diffs else None
    rec["raw_mode1"] = raw_aff[0][0]
    rec["output_identical"] = filecmp.cmp(pipeline_out, raw_out, shallow=False)

    ranking_csv = os.path.join(workdir, _benchlib.RANKING_CSV)
    row = None
    if os.path.exists(ranking_csv):
        row = next((r for r in _benchlib.read_ranking(ranking_csv)
                    if r["name"] == ligand), None)
    if row is None:
        rec["error"] = f"{ligand} not found in ranking.csv"
        return rec
    rec["ranking_affinity"] = row["affinity_str"]
    num_ok, str_ok = ranking_matches(row["affinity_str"], raw_aff[0][0])
    rec["ranking_numeric_match"], rec["ranking_string_match"] = num_ok, str_ok

    all_modes = rec["matched"] == rec["modes"] and rec["modes"] > 0
    rec["status"] = "ok" if (all_modes and num_ok and str_ok) else "mismatch"
    if rec["status"] == "mismatch" and not rec["error"]:
        rec["error"] = ("ranking.csv affinity differs from raw mode 1"
                        if all_modes else "per-mode affinities differ")
    return rec


def verdict(results):
    """PASS only when at least one complex was checked and all are ok."""
    return bool(results) and all(r["status"] == "ok" for r in results)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", required=True,
                    help="benchmark redock work dir (one subdir per complex)")
    ap.add_argument("--pdbids", help="Comma-separated IDs to check (else all)")
    ap.add_argument("--outdir", default="benchmarks/results/parity")
    ap.add_argument("--vina-bin", default=None,
                    help="Docking binary for the raw runs (default: the logged one)")
    ap.add_argument("--timeout", type=int, default=_benchlib.DOCKING_RUN_TIMEOUT)
    args = ap.parse_args(argv)

    if args.pdbids:
        pdbids = [p.strip() for p in args.pdbids.split(",") if p.strip()]
    else:
        pdbids = [d for d in sorted(os.listdir(args.workdir))
                  if os.path.isdir(os.path.join(args.workdir, d))]

    os.makedirs(args.outdir, exist_ok=True)
    results = []
    for p in pdbids:
        wd = os.path.join(args.workdir, p)
        if not os.path.isdir(wd):
            rec = {"pdbid": p, "status": "error", "modes": 0, "matched": 0,
                   "max_diff": None, "error": "no work dir"}
        else:
            rec = check_complex(wd, p, args.outdir, args.timeout, args.vina_bin)
        results.append(rec)
        print(f"[{p}] {rec['status']} modes={rec['modes']} matched={rec['matched']} "
              f"max_diff={rec['max_diff']} ranking={rec.get('ranking_affinity')} "
              f"raw1={rec.get('raw_mode1')} err={rec['error']}", flush=True)

    cols = ["pdbid", "status", "modes", "matched", "max_diff", "ranking_affinity",
            "raw_mode1", "ranking_numeric_match", "ranking_string_match",
            "output_identical", "binary", "error"]
    with open(os.path.join(args.outdir, "parity.csv"), "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    binaries = sorted({r["binary"] for r in results if r.get("binary")})
    ok = verdict(results)
    summary = {
        "checked": len(results),
        "pass": sum(1 for r in results if r["status"] == "ok"),
        "mismatch": sum(1 for r in results if r["status"] == "mismatch"),
        "error": sum(1 for r in results if r["status"] == "error"),
        "modes_compared": sum(r["modes"] for r in results),
        "modes_matched": sum(r["matched"] for r in results),
        "outputs_byte_identical": sum(1 for r in results if r.get("output_identical")),
        "binaries": {b: _benchlib.tool_version(b) for b in binaries},
        "result": "PASS" if ok else "FAIL",
    }
    with open(os.path.join(args.outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print("\n=== VINA PARITY ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
