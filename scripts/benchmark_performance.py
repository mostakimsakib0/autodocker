#!/usr/bin/env python3
"""Performance & scalability benchmark for the autodocker pipeline.

Runs the ``runner.py`` pipeline on a fixed receptor with increasing
ligand-library sizes and process counts, and records:

  * wall-clock time per run
  * peak resident set size (process group)
  * throughput (ligands / min)
  * speedup relative to the smallest size for the same process count

Ligands are a deterministic, seeded subsampling (with replacement) of a
source PDBQT directory, so results are reproducible for a given seed and
hardware. This benchmark measures throughput/scaling behaviour, not docking
accuracy (see BENCHMARKS.md for accuracy).

Usage
-----
python3 scripts/benchmark_performance.py \\
    --receptor ensemble/1xgi.pdbqt \\
    --ligands benchmarks/data/perf_ligs \\
    --sizes 10,25,50,100 \\
    --processes 1,2,4,8 \\
    --exhaustiveness 4 \\
    --outdir benchmarks/results/performance
"""

import argparse
import csv
import glob
import json
import os
import random
import shutil
import subprocess
import sys
import time


def build_library(src_ligands, size, seed, outdir):
    """Deterministically sample ``size`` ligands (with replacement)."""
    src = sorted(glob.glob(os.path.join(src_ligands, "*.pdbqt")))
    if not src:
        raise RuntimeError(f"no *.pdbqt files in {src_ligands}")
    rng = random.Random(seed)
    chosen = [rng.choice(src) for _ in range(size)]
    lib = os.path.join(outdir, "ligs")
    os.makedirs(lib, exist_ok=True)
    for i, path in enumerate(chosen):
        shutil.copy(path, os.path.join(lib, f"P{size:04d}_{i:05d}.pdbqt"))
    return lib


def time_run(receptor, lig_dir, outdir, processes, exhaustive, seed,
             runner):
    dock_out = os.path.join(outdir, "dock")
    cmd = [
        sys.executable, runner,
        "-r", receptor, "-l", lig_dir, "-o", dock_out,
        "--seed", str(seed),
        "--exhaustiveness", str(exhaustive),
        "--binding-modes", "1",
        "--energy-range", "1.0",
        "--processes", str(processes),
    ]
    t0 = time.monotonic()
    # Route output to /dev/null: the runner/obabel can emit huge stderr
    # (e.g. per-atom PDB warnings), which would fill a pipe buffer and
    # deadlock this process which never reads it. wait4() still gives us
    # the exit status and peak RSS.
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
    _, status, usage = os.wait4(proc.pid, 0)
    wall = time.monotonic() - t0
    rc = os.waitstatus_to_exitcode(status)
    peak_rss_kb = getattr(usage, "ru_maxrss", None)
    return {"wall_s": wall, "exit": rc, "peak_rss_kb": peak_rss_kb}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--receptor", required=True, help="Receptor PDB/PDBQT")
    ap.add_argument("--ligands", required=True,
                    help="Source dir of *.pdbqt ligands")
    ap.add_argument("--sizes", default="10,25,50,100",
                    help="Comma-separated library sizes")
    ap.add_argument("--processes", default="1,2,4,8",
                    help="Comma-separated process counts")
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default="benchmarks/results/performance")
    ap.add_argument("--runner", default="autodocker/runner.py")
    args = ap.parse_args()

    sizes = [int(s) for s in args.sizes.split(",")]
    procs = [int(p) for p in args.processes.split(",")]
    os.makedirs(args.outdir, exist_ok=True)
    work = os.path.join(args.outdir, "work")
    os.makedirs(work, exist_ok=True)

    results = []
    for size in sizes:
        for nproc in procs:
            run_dir = os.path.join(work, f"s{size}_p{nproc}")
            if os.path.exists(run_dir):
                shutil.rmtree(run_dir)
            os.makedirs(run_dir)
            lib = build_library(args.ligands, size, args.seed, run_dir)
            r = time_run(args.receptor, lib, run_dir, nproc,
                         args.exhaustiveness, args.seed,
                         os.path.abspath(args.runner))
            top = len(glob.glob(os.path.join(run_dir, "dock", "docked",
                                              "*_out.pdbqt")))
            row = {
                "n_ligands": size, "processes": nproc,
                "wall_s": round(r["wall_s"], 2),
                "peak_rss_kb": r["peak_rss_kb"],
                "throughput_ligands_per_min":
                    round(size / r["wall_s"] * 60, 2) if r["wall_s"] else 0.0,
                "docked_ok": top, "exit": r["exit"],
            }
            results.append(row)
            print(f"size={size:5d} procs={nproc} wall={r['wall_s']:7.2f}s "
                  f"rss={r['peak_rss_kb']}KB tput={row['throughput_ligands_per_min']:7.2f}/min "
                  f"docked={top} exit={r['exit']}", flush=True)

    # Speedup columns relative to the smallest size at the same process count.
    base = {p: next((r for r in results
                     if r["processes"] == p and r["n_ligands"] == min(sizes)),
                    None) for p in procs}
    for r in results:
        b = base.get(r["processes"])
        r["speedup_vs_min"] = round(b["wall_s"] / r["wall_s"],
                                    2) if b and r["wall_s"] else None

    rows_path = os.path.join(args.outdir, "performance.csv")
    with open(rows_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)
    with open(os.path.join(args.outdir, "summary.json"), "w") as f:
        json.dump({"seed": args.seed, "exhaustiveness": args.exhaustiveness,
                   "results": results}, f, indent=2)
    print(f"\nResults written to {rows_path}")


if __name__ == "__main__":
    main()