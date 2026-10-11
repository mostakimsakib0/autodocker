#!/usr/bin/env python3
"""Performance & scalability benchmark for the autodocker pipeline.

Runs ``runner.py`` on a fixed receptor with increasing ligand-library sizes
and process counts, each cell ``--repeats`` times (default 3, interleaved so
slow drift of the host does not favour one cell), and records per run:

  * wall-clock time, split into phases from the runner's step markers
    (time-stamped as the log lines arrive): ``prep`` = start until
    "[3/5] Docking" (library + receptor preparation, grid), ``dock`` =
    "[3/5]" until "[4/5]" (docking incl. optional rescoring), ``post`` =
    the rest (ranking, reports)
  * peak memory of the whole process tree (runner + obabel/vina children),
    sampled every ``--mem-interval`` seconds: summed RSS and, where
    readable, summed PSS (shared pages split between processes). Both are
    sampled, so a peak shorter than the interval can be missed. Without
    psutil the tree is read from /proc (Linux).
    ``max_single_process_rss_kb`` (``ru_maxrss`` from wait4) is the largest
    single process, not the group.
  * Vina threads per process (``--cpu`` value in the logged Vina command;
    None when the runner passed no ``--cpu``, i.e. Vina used all cores,
    which is the runner default for ``-p 1``)

and per cell the median and range over repeats, throughput from the median,
and the speedup vs. ``-p 1`` at the same library size (only when p=1 is in
``--processes``), with parallel efficiency = speedup / p.

Ligands are a deterministic, seeded subsample (with replacement) of a
source PDBQT directory. This measures throughput, not docking accuracy.

Usage
-----
python3 scripts/benchmark_performance.py \\
    --receptor receptor.pdb --ligands ligs/ \\
    --sizes 10,25,50,100 --processes 1,2,4,8 --repeats 3 \\
    --exhaustiveness 4 --outdir benchmarks/results/performance
"""

import argparse
import csv
import glob
import json
import os
import platform
import random
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _benchlib  # noqa: E402

try:
    import psutil
except ImportError:  # the AutoDocker image does not ship psutil
    psutil = None

PHASE_MARKERS = {"[3/5]": "dock_start", "[4/5]": "post_start"}


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


# ---------------------------------------------------------------------------
# Process-tree memory
# ---------------------------------------------------------------------------

def _proc_tree_pids(root):
    """PIDs of ``root`` and all its descendants (via /proc)."""
    children = {}
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open(f"/proc/{entry}/stat", encoding="utf-8") as f:
                stat = f.read()
            ppid = int(stat.rsplit(")", 1)[1].split()[1])
        except (OSError, ValueError, IndexError):
            continue
        children.setdefault(ppid, []).append(int(entry))
    out, todo = [], [root]
    while todo:
        pid = todo.pop()
        out.append(pid)
        todo.extend(children.get(pid, []))
    return out


def _proc_kb(path, field):
    try:
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.startswith(field + ":"):
                    return int(line.split()[1])
    except (OSError, ValueError, IndexError):
        return None
    return None


def tree_memory_kb(root):
    """``(rss_sum_kb, pss_sum_kb or None)`` over the process tree of ``root``.

    PIDs come from psutil when installed, else from /proc; RSS/PSS are read
    from /proc/<pid>/status and smaps_rollup (Linux). PSS is None when it is
    unreadable for any process of the tree."""
    if psutil is not None:
        try:
            top = psutil.Process(root)
            pids = [root] + [c.pid for c in top.children(recursive=True)]
        except psutil.Error:
            return 0, None
    else:
        pids = _proc_tree_pids(root)
    rss, pss, pss_ok = 0, 0, True
    for pid in pids:
        value = _proc_kb(f"/proc/{pid}/status", "VmRSS")
        if value is None and psutil is not None:
            try:
                value = psutil.Process(pid).memory_info().rss // 1024
            except psutil.Error:
                value = None
        rss += value or 0
        value = _proc_kb(f"/proc/{pid}/smaps_rollup", "Pss")
        if value is None:
            pss_ok = False
        else:
            pss += value
    return rss, (pss if pss_ok else None)


# ---------------------------------------------------------------------------
# One run
# ---------------------------------------------------------------------------

def time_run(receptor, lig_dir, outdir, processes, exhaustive, seed, runner,
             threads=None, mem_interval=0.5):
    dock_out = os.path.join(outdir, "dock")
    cmd = [
        sys.executable, runner,
        "-r", receptor, "-l", lig_dir, "-o", dock_out,
        "--seed", str(seed),
        "--exhaustiveness", str(exhaustive),
        "--binding-modes", "1",
        "--energy-range", "1.0",
        "--processes", str(processes),
        "--no-resume",
    ]
    if threads:
        cmd += ["--threads", str(threads)]
    marks = {}
    peak = {"rss": 0, "pss": None}
    done = threading.Event()
    t0 = time.monotonic()
    # Output is read continuously (never left in a full pipe buffer) and
    # stored in runner.log; the arrival time of each step marker gives the
    # phase split.
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    def reader():
        with open(os.path.join(outdir, "runner.log"), "wb") as log:
            for raw in proc.stdout:
                now = time.monotonic()
                log.write(raw)
                text = raw.decode("utf-8", "replace")
                for marker, key in PHASE_MARKERS.items():
                    if marker in text and key not in marks:
                        marks[key] = now - t0

    def sampler():
        while not done.wait(mem_interval):
            rss, pss = tree_memory_kb(proc.pid)
            peak["rss"] = max(peak["rss"], rss)
            if pss is not None:
                peak["pss"] = max(peak["pss"] or 0, pss)

    threads_ = [threading.Thread(target=reader, daemon=True),
                threading.Thread(target=sampler, daemon=True)]
    for t in threads_:
        t.start()
    _, status, usage = os.wait4(proc.pid, 0)
    wall = time.monotonic() - t0
    proc.returncode = os.waitstatus_to_exitcode(status)
    done.set()
    for t in threads_:
        t.join(timeout=10)
    phases = {"prep_s": marks.get("dock_start"),
              "dock_s": (marks["post_start"] - marks["dock_start"]
                         if "dock_start" in marks and "post_start" in marks else None),
              "post_s": (wall - marks["post_start"] if "post_start" in marks else None)}
    return {"wall_s": wall, "exit": proc.returncode, **phases,
            "peak_tree_rss_kb": peak["rss"], "peak_tree_pss_kb": peak["pss"],
            "max_single_process_rss_kb": getattr(usage, "ru_maxrss", None)}


def vina_threads(dock_dir):
    """``--cpu`` value of the first logged Vina command (None if absent)."""
    for log in sorted(glob.glob(os.path.join(dock_dir, _benchlib.DOCKED_DIR, "*_vina.log"))):
        with open(log, encoding="utf-8") as f:
            for line in f:
                m = re.search(r"--cpu\s+(\d+)", line) if line.startswith("# Command:") else None
                if m:
                    return int(m.group(1))
    return None


def count_docked(dock_dir):
    ranking = os.path.join(dock_dir, _benchlib.RANKING_CSV)
    if not os.path.exists(ranking):
        return 0
    return sum(1 for r in _benchlib.read_ranking(ranking) if r["affinity"] is not None)


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate(runs):
    """Per-(size, processes) median and range over repeats."""
    cells = {}
    for r in runs:
        cells.setdefault((r["n_ligands"], r["processes"]), []).append(r)
    rows = []
    for (size, procs), rs in sorted(cells.items()):
        ok = [r for r in rs if r["exit"] == 0]
        walls = [r["wall_s"] for r in ok]

        def med(key):
            vals = [r[key] for r in ok if r.get(key) is not None]
            return statistics.median(vals) if vals else None

        wall = statistics.median(walls) if walls else None
        rows.append({
            "n_ligands": size, "processes": procs,
            "repeats": len(rs), "repeats_ok": len(ok),
            "wall_s_median": wall,
            "wall_s_min": min(walls) if walls else None,
            "wall_s_max": max(walls) if walls else None,
            "prep_s_median": med("prep_s"), "dock_s_median": med("dock_s"),
            "post_s_median": med("post_s"),
            "throughput_ligands_per_min": (size / wall * 60.0) if wall else None,
            "peak_tree_rss_kb_max": max((r["peak_tree_rss_kb"] for r in ok), default=None),
            "peak_tree_pss_kb_max": max((r["peak_tree_pss_kb"] for r in ok
                                         if r["peak_tree_pss_kb"] is not None), default=None),
            "max_single_process_rss_kb": max((r["max_single_process_rss_kb"] for r in ok
                                              if r["max_single_process_rss_kb"]), default=None),
            "vina_threads_per_process": next((r["vina_threads"] for r in ok
                                              if r.get("vina_threads")), None),
            "docked_ok_min": min((r["docked_ok"] for r in rs), default=0),
        })
    return add_speedup(rows)


def add_speedup(rows):
    """Speedup vs. p=1 at the same library size, and parallel efficiency."""
    base = {r["n_ligands"]: r["wall_s_median"] for r in rows
            if r["processes"] == 1 and r["wall_s_median"]}
    for r in rows:
        b = base.get(r["n_ligands"])
        if b and r["wall_s_median"]:
            r["speedup_vs_p1"] = b / r["wall_s_median"]
            r["efficiency"] = r["speedup_vs_p1"] / r["processes"]
        else:
            r["speedup_vs_p1"] = r["efficiency"] = None
    return rows


def host_info(vina_bin=None):
    info = {"logical_cpus": os.cpu_count(),
            "usable_cpus": (len(os.sched_getaffinity(0))
                            if hasattr(os, "sched_getaffinity") else None),
            "physical_cores": None, "cpu_model": None,
            "platform": platform.platform(), "python": platform.python_version(),
            "psutil": getattr(psutil, "__version__", None)}
    try:
        out = subprocess.run(["lscpu", "-p=CORE,SOCKET"], capture_output=True,
                             text=True, timeout=30).stdout
        cores = {ln for ln in out.splitlines() if ln and not ln.startswith("#")}
        info["physical_cores"] = len(cores) or None
        model = subprocess.run(["lscpu"], capture_output=True, text=True,
                               timeout=30).stdout
        m = re.search(r"^Model name:\s*(.+)$", model, re.M)
        info["cpu_model"] = m.group(1).strip() if m else None
    except (OSError, subprocess.TimeoutExpired):
        pass
    if info["cpu_model"] is None:
        try:
            with open("/proc/cpuinfo", encoding="utf-8") as f:
                m = re.search(r"^model name\s*:\s*(.+)$", f.read(), re.M)
            info["cpu_model"] = m.group(1).strip() if m else None
        except OSError:
            pass
    info["vina_version"] = _benchlib.tool_version(vina_bin or "vina")
    return info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--receptor", required=True, help="Receptor PDB")
    ap.add_argument("--ligands", required=True, help="Source dir of *.pdbqt ligands")
    ap.add_argument("--sizes", default="10,25,50,100", help="Comma-separated library sizes")
    ap.add_argument("--processes", default="1,2,4,8", help="Comma-separated process counts")
    ap.add_argument("--repeats", type=int, default=3, help="Runs per cell (default 3)")
    ap.add_argument("--threads", type=int, default=None,
                    help="Vina threads per process (default: runner default, cores / -p)")
    ap.add_argument("--mem-interval", type=float, default=0.5,
                    help="Memory sampling interval in seconds (default 0.5)")
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default="benchmarks/results/performance")
    ap.add_argument("--runner", default="autodocker/runner.py")
    args = ap.parse_args(argv)
    if args.repeats < 1:
        ap.error("--repeats must be >= 1")

    sizes = [int(s) for s in args.sizes.split(",")]
    procs = [int(p) for p in args.processes.split(",")]
    runner = os.path.abspath(args.runner)
    os.makedirs(args.outdir, exist_ok=True)
    work = os.path.join(args.outdir, "work")
    os.makedirs(work, exist_ok=True)
    libs = {}
    for size in sizes:
        lib_root = os.path.join(work, f"lib_s{size}")
        shutil.rmtree(lib_root, ignore_errors=True)
        libs[size] = build_library(args.ligands, size, args.seed, lib_root)

    runs = []
    for rep in range(args.repeats):
        for size in sizes:
            for nproc in procs:
                run_dir = os.path.join(work, f"s{size}_p{nproc}_r{rep}")
                shutil.rmtree(run_dir, ignore_errors=True)
                os.makedirs(run_dir)
                r = time_run(args.receptor, libs[size], run_dir, nproc,
                             args.exhaustiveness, args.seed, runner,
                             threads=args.threads, mem_interval=args.mem_interval)
                dock_dir = os.path.join(run_dir, "dock")
                r.update({"n_ligands": size, "processes": nproc, "repeat": rep,
                          "docked_ok": count_docked(dock_dir),
                          "vina_threads": vina_threads(dock_dir)})
                runs.append(r)
                print(f"rep={rep} size={size:5d} procs={nproc} wall={r['wall_s']:8.2f}s "
                      f"prep={r['prep_s']} dock={r['dock_s']} "
                      f"rss_tree={r['peak_tree_rss_kb']}KB docked={r['docked_ok']} "
                      f"threads={r['vina_threads']} exit={r['exit']}", flush=True)

    rows = aggregate(runs)
    with open(os.path.join(args.outdir, "runs.csv"), "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(runs[0].keys()))
        writer.writeheader()
        writer.writerows(runs)
    rows_path = os.path.join(args.outdir, "performance.csv")
    with open(rows_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    with open(os.path.join(args.outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump({"seed": args.seed, "exhaustiveness": args.exhaustiveness,
                   "repeats": args.repeats, "threads_option": args.threads,
                   "mem_interval_s": args.mem_interval,
                   "runner_version": _benchlib.runner_version(runner),
                   "host": host_info(), "cells": rows, "runs": runs}, f, indent=2)
    print(f"\nResults written to {rows_path}")
    return 0 if all(r["exit"] == 0 for r in runs) else 1


if __name__ == "__main__":
    sys.exit(main())
