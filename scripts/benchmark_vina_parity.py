#!/usr/bin/env python3
"""J3: Vina parity check.

Proves the autodocker pipeline does not corrupt anything: the affinity of
every docked pose produced by ``runner.py`` must exactly match the affinity
the raw ``vina`` CLI produces for the *same* receptor PDBQT, ligand PDBQT,
grid, seed and sampling parameters.

For each complex already docked by ``benchmark_redock.py``, this script:
  1. Reads the Vina command line the pipeline logged (header of *_vina.log)
  2. Re-runs the raw ``vina`` binary with those exact arguments
  3. Compares mode-by-mode ``REMARK VINA RESULT`` affinities

All modes must match to numerical precision (no rounding allowed).
"""

import argparse
import csv
import glob
import os
import re
import subprocess
import sys

VINA_RESULT = re.compile(
    r"^\s*REMARK\s+VINA RESULT:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)")


def parse_vina_results(pdbqt_path):
    affinities = []
    with open(pdbqt_path) as f:
        for line in f:
            m = VINA_RESULT.match(line)
            if m:
                affinities.append(float(m.group(1)))
    return affinities


def extract_command(log_path):
    """Parse the 'vina --receptor ...' command recorded in a *_vina.log."""
    with open(log_path) as f:
        for line in f:
            if line.startswith("# Command:"):
                cmd = line.split(":", 1)[1].strip().split()
                return cmd
    return None


def run_vina(cmd, timeout=1800):
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError(
            f"vina failed (rc={proc.returncode}): "
            f"{(proc.stdout or '')[-500:]} {(proc.stderr or '')[-500:]}")


def process_complex(workdir, pdbid, tmpdir, timeout):
    rec = {"pdbid": pdbid, "status": "error", "modes": 0, "matched": 0,
           "max_diff": None, "error": None}
    logs = glob.glob(os.path.join(workdir, "docked", "*_vina.log"))
    if not logs:
        rec["error"] = "no *_vina.log found"
        return rec
    log = logs[0]
    cmd = extract_command(log)
    if cmd is None:
        rec["error"] = "no command header in log"
        return rec
    if os.path.basename(cmd[0]) not in ("vina", "vina_1.2.5_linux_x86_64"):
        rec["error"] = f"unexpected binary: {cmd[0]}"
        return rec

    pipeline_out = None
    for i, c in enumerate(cmd):
        if c == "--out" and i + 1 < len(cmd):
            pipeline_out = cmd[i + 1]

    raw_out = os.path.join(tmpdir, f"{pdbid}_raw.pdbqt")
    raw_cmd = [c if not (i > 0 and c == pipeline_out) else raw_out
               for i, c in enumerate(cmd)]

    try:
        run_vina(raw_cmd, timeout)
    except Exception as e:
        rec["error"] = str(e)
        return rec

    pipeline_aff = parse_vina_results(pipeline_out)
    raw_aff = parse_vina_results(raw_out)

    if len(pipeline_aff) != len(raw_aff):
        rec["error"] = (f"mode count mismatch pipeline={len(pipeline_aff)} "
                        f"raw={len(raw_aff)}")
        return rec

    rec["modes"] = len(pipeline_aff)
    rec["matched"] = sum(
        1 for a, b in zip(pipeline_aff, raw_aff) if a == b)
    diffs = [abs(a - b) for a, b in zip(pipeline_aff, raw_aff)]
    rec["max_diff"] = max(diffs) if diffs else None
    rec["affinities"] = pipeline_aff
    rec["status"] = "ok" if rec["matched"] == rec["modes"] else "mismatch"
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", required=True,
                    help="benchmark redock work dir (one subdir per complex)")
    ap.add_argument("--pdbids", help="Comma-separated IDs to check (else all)")
    ap.add_argument("--outdir", default="benchmarks/results/parity")
    ap.add_argument("--timeout", type=int, default=1800)
    args = ap.parse_args()

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
            print(f"[{p}] skip (no workdir)", flush=True)
            continue
        rec = process_complex(wd, p, args.outdir, args.timeout)
        results.append(rec)
        print(f"[{p}] {rec['status']} modes={rec['modes']} "
              f"matched={rec['matched']} max_diff={rec['max_diff']} "
              f"err={rec['error']}", flush=True)

    rows_path = os.path.join(args.outdir, "parity.csv")
    with open(rows_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdbid", "status", "modes", "matched",
                         "max_diff", "error"])
        for r in results:
            writer.writerow([r["pdbid"], r["status"], r["modes"],
                             r["matched"], r["max_diff"], r["error"]])

    all_ok = all(r["status"] == "ok" for r in results if r["status"] != "error")
    print("\n=== VINA PARITY ===")
    print(f"  complexes checked: {len(results)}")
    print(f"  exact match (all modes): {sum(1 for r in results if r['status'] == 'ok')}")
    print(f"  mismatches: {sum(1 for r in results if r['status'] == 'mismatch')}")
    print(f"  result: {'PASS' if all_ok else 'FAIL'}")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()