#!/usr/bin/env python3
"""J1: CASF-2016 core-set re-docking (docking power) benchmark.

Protocol
--------
  * Complexes: a seeded random sample of the CASF-2016 core set
               (``--sample 40 --sample-seed 42``, IDs from
               ``benchmarks/data/casf2016_core_set.csv``), or ``--pdbids``.
               Sampled complexes missing from ``--data-dir`` (e.g. peptide
               ligands that ``fetch_casf_rcsb.py`` cannot build) are listed
               with the reason recorded in ``manifest.csv``.
  * Receptor : ``{pdbid}_protein.pdb``, all chains; metal ions are kept
               (``--receptor-hetatm metals`` -> runner ``--detect-metals``).
               No receptor-size filter unless ``--max-receptor-atoms`` is set.
  * Ligand   : start conformer from ``--start-conformer``:
               ``rdkit`` (default) = one ETKDGv3 conformer embedded from the
               CCD-derived molecule (fixed randomSeed, single thread, then
               MMFF94 with a fixed iteration cap; deterministic), so the
               docked ring conformations, bond lengths and angles do not
               come from the crystal pose; ``crystal`` = the crystal-ligand
               geometry itself (1.1.x protocol, optimistic).
  * Box      : cubic, centred on the crystal-ligand heavy-atom centroid,
               edge = max(22.5, ligand extent + 2*padding) A
               (``_benchlib.ligand_box``), passed as --grid-center/--grid-size
               in input-PDB coordinates.
  * Docking  : autodocker ``runner.py`` (Vina), fixed seed, exhaustiveness,
               num_modes and energy range.
  * Metric   : heavy-atom, symmetry-corrected RMSD of each output mode vs.
               the crystal pose (``obrms``, no superposition). Top-1 is
               always mode 1 (the first MODEL); if its RMSD cannot be
               computed the complex is reported as ``rmsd_failed``.
  * Success  : RMSD <= --rmsd-cutoff (default 2.0 A). Success rates are
               reported over complexes evaluated and over complexes
               attempted (failures count as misses).

Resume: a complex is reused only when ``params.json`` in its work directory
matches the current docking parameters; otherwise it is re-docked.

Results: ``outdir/results.csv`` and ``outdir/summary.json``.
"""

import argparse
import csv
import json
import multiprocessing
import os
import shutil
import statistics
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _benchlib  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_IDS = os.path.join(HERE, "..", "benchmarks", "data", "casf2016_core_set.csv")
HETATM_FLAGS = {"none": [], "metals": ["--detect-metals"], "all": ["--keep-hetero"]}


def ligand_grid(mol2_path, padding=_benchlib.DEFAULT_PADDING,
                min_size=_benchlib.MIN_BOX_EDGE):
    """Box around the crystal ligand: ``(cx, cy, cz, edge)``."""
    center, size = _benchlib.ligand_box(_benchlib.mol2_heavy_coords(mol2_path),
                                        padding=padding, min_size=min_size)
    return center[0], center[1], center[2], size


def ranking_affinity(ranking_csv, ligand):
    """Affinity (float or None) and the text written for ``ligand``."""
    if not os.path.exists(ranking_csv):
        return None, None
    for row in _benchlib.read_ranking(ranking_csv):
        if row["name"] == ligand:
            return row["affinity"], row["affinity_str"]
    return None, None


def mode_rmsds(crystal_ref, out_pdbqt, rmsd_fn=None):
    """RMSD of every output mode vs. the crystal reference, in file order.

    Returns ``(rmsds, errors)``: ``rmsds[i]`` is the RMSD of mode i+1 or
    None when it could not be computed; ``errors`` lists the reasons.
    """
    rmsd_fn = rmsd_fn or _benchlib.obrms
    models = _benchlib.parse_vina_models(out_pdbqt)
    rmsds, errors = [], []
    with tempfile.TemporaryDirectory() as tmp:
        for m in models:
            path = _benchlib.write_model(m, os.path.join(tmp, f"mode_{m['model']}.pdbqt"))
            try:
                rmsds.append(rmsd_fn(crystal_ref, path))
            except Exception as e:
                rmsds.append(None)
                errors.append(f"mode {m['model']}: {e}")
    return rmsds, errors


def apply_rmsds(rec, rmsds, errors):
    """Fill rmsd_top1 / rmsd_best / status from per-mode RMSDs."""
    rec["rmsd_modes"] = rmsds
    rec["rmsd_top1"] = rmsds[0] if rmsds else None
    valid = [r for r in rmsds if r is not None]
    rec["rmsd_best"] = min(valid) if valid else None
    rec["n_modes"] = len(rmsds)
    rec["n_modes_rmsd_failed"] = len(rmsds) - len(valid)
    if errors:
        rec["error"] = "; ".join(errors)[:500]
    if not rmsds:
        rec["status"] = "rmsd_failed"
        rec["error"] = rec.get("error") or "no MODEL records in the docked output"
    elif rec["rmsd_top1"] is None:
        rec["status"] = "rmsd_failed"
    else:
        rec["status"] = "ok"
    return rec


def docking_params(opts):
    """Everything that changes the docked poses (stored per complex)."""
    keys = ("seed", "exhaustiveness", "binding_modes", "energy_range",
            "padding", "min_box", "threads", "start_conformer",
            "conformer_seed", "receptor_hetatm", "runner_version")
    return {k: opts.get(k) for k in keys}


def params_match(work, params):
    try:
        with open(os.path.join(work, "params.json"), encoding="utf-8") as f:
            return json.load(f).get("docking") == params
    except (OSError, ValueError):
        return False


def prepare_start_ligand(pdbid, complex_dir, work, opts):
    """Ligand file handed to the runner, plus metadata about its origin."""
    crystal = os.path.join(complex_dir, f"{pdbid}_ligand.mol2")
    lig_dir = os.path.join(work, "input")
    os.makedirs(lig_dir, exist_ok=True)
    out = os.path.join(lig_dir, f"{pdbid}_ligand.mol2")
    if opts["start_conformer"] == "crystal":
        shutil.copy(crystal, out)
        return out, {"method": "crystal"}
    sdf = os.path.join(complex_dir, f"{pdbid}_ligand.sdf")
    pre = os.path.join(complex_dir, f"{pdbid}_ligand_rdkit.mol2")
    pre_meta = pre[:-len(".mol2")] + ".json"
    try:
        meta = _benchlib.rdkit_start_conformer(
            sdf if os.path.exists(sdf) else crystal, out, opts["conformer_seed"])
        return out, meta
    except ImportError:
        pass
    # No RDKit here (e.g. inside the AutoDocker image): use the conformer
    # fetch_casf_rcsb.py generated, if it was made with the same seed.
    if os.path.exists(pre) and os.path.exists(pre_meta):
        with open(pre_meta, encoding="utf-8") as f:
            meta = json.load(f)
        if meta.get("seed") == opts["conformer_seed"]:
            shutil.copy(pre, out)
            meta["source"] = "precomputed"
            return out, meta
    raise RuntimeError(
        "--start-conformer rdkit needs RDKit, or a precomputed "
        f"{os.path.basename(pre)} (+ .json, same seed) from fetch_casf_rcsb.py")


def run_docking(pdbid, complex_dir, work, opts):
    protein = os.path.join(complex_dir, f"{pdbid}_protein.pdb")
    crystal = os.path.join(complex_dir, f"{pdbid}_ligand.mol2")
    cx, cy, cz, size = ligand_grid(crystal, opts["padding"], opts["min_box"])
    ligand, conf_meta = prepare_start_ligand(pdbid, complex_dir, work, opts)
    cmd = [
        sys.executable, opts["runner"],
        "-r", protein, "-l", ligand, "-o", work,
        "--chain=all", "--no-resume", "--no-html-report",
        *_benchlib.box_cli_args((cx, cy, cz), size),
        *HETATM_FLAGS[opts["receptor_hetatm"]],
        "--seed", str(opts["seed"]),
        "--exhaustiveness", str(opts["exhaustiveness"]),
        "--binding-modes", str(opts["binding_modes"]),
        "--energy-range", str(opts["energy_range"]),
        "--timeout", str(opts["timeout"]),
    ]
    if opts.get("threads"):
        cmd.extend(["--threads", str(opts["threads"])])
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          timeout=opts["timeout"])
    with open(os.path.join(work, "runner.log"), "w", encoding="utf-8") as f:
        f.write(" ".join(cmd) + "\n")
        f.write(proc.stdout or "")
        f.write(proc.stderr or "")
    out_pdbqt = os.path.join(work, _benchlib.DOCKED_DIR, f"{pdbid}_ligand_out.pdbqt")
    if proc.returncode != 0 or not os.path.exists(out_pdbqt):
        raise RuntimeError(
            f"docking failed (rc={proc.returncode}): "
            f"{(proc.stdout or '')[-300:]} {(proc.stderr or '')[-300:]}")
    with open(os.path.join(work, "params.json"), "w", encoding="utf-8") as f:
        json.dump({"docking": docking_params(opts), "start_conformer": conf_meta,
                   "box": {"center_input_frame": [cx, cy, cz], "edge": size}},
                  f, indent=2)
    return {"box_edge": size, "start_conformer": conf_meta}


def process_complex(job):
    pdbid, data_dir, workroot, opts = job
    rec = {"pdbid": pdbid, "status": "error", "affinity": None,
           "affinity_mode1": None, "rmsd_top1": None, "rmsd_best": None,
           "rmsd_modes": [], "reused": False, "error": None}
    complex_dir = os.path.join(data_dir, pdbid)
    protein = os.path.join(complex_dir, f"{pdbid}_protein.pdb")
    crystal = os.path.join(complex_dir, f"{pdbid}_ligand.mol2")
    work = os.path.join(workroot, pdbid)
    out_pdbqt = os.path.join(work, _benchlib.DOCKED_DIR, f"{pdbid}_ligand_out.pdbqt")
    try:
        if not (os.path.exists(protein) and os.path.exists(crystal)):
            rec["status"] = "missing"
            rec["error"] = opts["manifest"].get(pdbid) or "not in --data-dir"
            return rec
        if opts["max_receptor_atoms"]:
            with open(protein, encoding="utf-8") as f:
                n_atom = sum(1 for ln in f if ln.startswith("ATOM"))
            if n_atom > opts["max_receptor_atoms"]:
                rec["status"] = "skipped"
                rec["error"] = (f"receptor too large ({n_atom} ATOM > "
                                f"{opts['max_receptor_atoms']})")
                return rec
        params = docking_params(opts)
        if os.path.exists(out_pdbqt) and params_match(work, params):
            rec["reused"] = True
        else:
            shutil.rmtree(work, ignore_errors=True)
            os.makedirs(work)
            rec.update(run_docking(pdbid, complex_dir, work, opts))
        # Poses are in the docking frame; move the crystal ligand into it
        # with the offset the runner recorded.
        shift = _benchlib.read_frame_offset(work)
        crystal_ref = os.path.join(work, f"{pdbid}_crystal_centered.mol2")
        _benchlib.shift_mol2(crystal, shift, crystal_ref)
        rec["affinity"], _ = ranking_affinity(os.path.join(work, _benchlib.RANKING_CSV),
                                              f"{pdbid}_ligand")
        _, rec["affinity_mode1"] = _benchlib.first_vina_result(out_pdbqt)
        apply_rmsds(rec, *mode_rmsds(crystal_ref, out_pdbqt))
    except Exception as e:
        rec["status"] = "error"
        rec["error"] = str(e)[:500]
    return rec


def _frac(k, n):
    return k / n if n else None


def summarize(results, cutoff):
    """Success rates with explicit denominators, medians, exclusions."""
    evaluated = [r for r in results if r["status"] == "ok"]
    attempted = [r for r in results if r["status"] not in ("missing", "skipped")]
    best = [r for r in attempted if r.get("rmsd_best") is not None]
    top1_ok = sum(1 for r in evaluated if r["rmsd_top1"] <= cutoff)
    best_ok = sum(1 for r in best if r["rmsd_best"] <= cutoff)
    return {
        "n_selected": len(results),
        "n_missing": sum(1 for r in results if r["status"] == "missing"),
        "n_skipped": sum(1 for r in results if r["status"] == "skipped"),
        "n_attempted": len(attempted),
        "n_docking_failed": sum(1 for r in results if r["status"] == "error"),
        "n_rmsd_failed": sum(1 for r in results if r["status"] == "rmsd_failed"),
        "n_evaluated": len(evaluated),
        "success_top1": top1_ok,
        "success_top1_frac_evaluated": _frac(top1_ok, len(evaluated)),
        "success_top1_frac_attempted": _frac(top1_ok, len(attempted)),
        "success_best": best_ok,
        "success_best_frac_evaluated": _frac(best_ok, len(best)),
        "success_best_frac_attempted": _frac(best_ok, len(attempted)),
        "median_rmsd_top1": (statistics.median(r["rmsd_top1"] for r in evaluated)
                             if evaluated else None),
        "median_rmsd_best": (statistics.median(r["rmsd_best"] for r in best)
                             if best else None),
        "rmsd_cutoff": cutoff,
        "excluded": [{"pdbid": r["pdbid"], "status": r["status"],
                      "reason": r["error"]}
                     for r in results if r["status"] != "ok"],
    }


def read_manifest(data_dir):
    """pdbid -> reason for complexes fetch_casf_rcsb.py could not build."""
    path = os.path.join(data_dir, "manifest.csv")
    reasons = {}
    if os.path.exists(path):
        with open(path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                if row.get("status") != "ok":
                    reasons[row["pdbid"]] = f"fetch skipped: {row.get('reason', '')}"
    return reasons


def select_ids(args):
    if args.pdbids:
        return sorted({p.strip() for p in args.pdbids.split(",") if p.strip()})
    if args.first:
        return sorted(d for d in os.listdir(args.data_dir)
                      if os.path.isdir(os.path.join(args.data_dir, d)))[:args.first]
    with open(args.ids_csv, newline="", encoding="utf-8") as f:
        ids = [row["pdbid"].strip().lower() for row in csv.DictReader(f)]
    return _benchlib.select_sample(ids, args.sample, args.sample_seed)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True,
                    help="Dir with one sub-dir per PDB ID (fetch_casf_rcsb.py output "
                         "or extracted PDBbind core set)")
    ap.add_argument("--pdbids", help="Comma-separated PDB IDs (overrides --sample)")
    ap.add_argument("--sample", type=int, default=40,
                    help="Seeded random sample size from --ids-csv (default 40; "
                         "0 = the whole core set)")
    ap.add_argument("--sample-seed", type=int, default=42,
                    help="Seed of the complex sample (default 42)")
    ap.add_argument("--ids-csv", default=DEFAULT_IDS,
                    help="CSV with a pdbid column (default: CASF-2016 core set)")
    ap.add_argument("--first", "--max-complexes", dest="first", type=int, default=None,
                    help="Legacy 1.1.x selection: first N sorted IDs in --data-dir")
    ap.add_argument("--outdir", default="benchmarks/results/redock")
    ap.add_argument("--runner", default="autodocker/runner.py")
    ap.add_argument("--processes", type=int, default=_benchlib.available_cpus(),
                    help="Complexes docked in parallel (default: available CPUs)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--binding-modes", type=int, default=9)
    ap.add_argument("--energy-range", type=float, default=3.0)
    ap.add_argument("--padding", type=float, default=_benchlib.DEFAULT_PADDING)
    ap.add_argument("--min-box", type=float, default=_benchlib.MIN_BOX_EDGE)
    ap.add_argument("--threads", type=int, default=0,
                    help="Vina threads per complex (0 = runner default)")
    ap.add_argument("--start-conformer", choices=["rdkit", "crystal"], default="rdkit",
                    help="Ligand start geometry (default rdkit: seeded ETKDGv3 + MMFF)")
    ap.add_argument("--conformer-seed", type=int, default=42,
                    help="ETKDG randomSeed for --start-conformer rdkit (default 42)")
    ap.add_argument("--receptor-hetatm", choices=sorted(HETATM_FLAGS), default="metals",
                    help="Receptor HETATM records kept by the runner: metals "
                         "(--detect-metals, default), all (--keep-hetero) or none")
    ap.add_argument("--max-receptor-atoms", type=int, default=0,
                    help="Skip receptors with more ATOM records (0 = no filter, default)")
    ap.add_argument("--rmsd-cutoff", type=float, default=2.0)
    ap.add_argument("--timeout", type=int, default=1800)
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    pdbids = select_ids(args)
    runner = os.path.abspath(args.runner)
    opts = {"padding": args.padding, "min_box": args.min_box, "seed": args.seed,
            "exhaustiveness": args.exhaustiveness,
            "binding_modes": args.binding_modes,
            "energy_range": args.energy_range, "threads": args.threads,
            "start_conformer": args.start_conformer,
            "conformer_seed": (args.conformer_seed
                               if args.start_conformer == "rdkit" else None),
            "receptor_hetatm": args.receptor_hetatm,
            "runner": runner, "runner_version": _benchlib.runner_version(runner),
            "max_receptor_atoms": args.max_receptor_atoms,
            "timeout": args.timeout, "manifest": read_manifest(args.data_dir)}
    workroot = os.path.join(args.outdir, "work")
    jobs = [(p, args.data_dir, workroot, opts) for p in pdbids]

    results = []
    with multiprocessing.Pool(processes=max(1, args.processes)) as pool:
        for rec in pool.imap_unordered(process_complex, jobs):
            results.append(rec)
            print(f"[{rec['pdbid']}] {rec['status']}{' (reused)' if rec['reused'] else ''} "
                  f"aff={rec['affinity']} rmsd_top1={rec['rmsd_top1']} "
                  f"rmsd_best={rec['rmsd_best']} err={rec['error']}", flush=True)

    results.sort(key=lambda r: r["pdbid"])
    with open(os.path.join(args.outdir, "results.csv"), "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["pdbid", "status", "affinity", "affinity_mode1",
                         "rmsd_top1", "rmsd_best", "rmsd_modes",
                         "n_modes_rmsd_failed", "reused", "error"])
        for r in results:
            writer.writerow([
                r["pdbid"], r["status"],
                "" if r["affinity"] is None else r["affinity"],
                "" if r["affinity_mode1"] is None else r["affinity_mode1"],
                "" if r["rmsd_top1"] is None else f"{r['rmsd_top1']:.3f}",
                "" if r["rmsd_best"] is None else f"{r['rmsd_best']:.3f}",
                " ".join("NA" if x is None else f"{x:.3f}" for x in r["rmsd_modes"]),
                r.get("n_modes_rmsd_failed", ""), r["reused"], r["error"] or ""])

    summary = summarize(results, args.rmsd_cutoff)
    summary["protocol"] = {
        "selection": ("pdbids" if args.pdbids else
                      f"first {args.first}" if args.first else
                      f"sample {args.sample} seed {args.sample_seed} of "
                      f"{os.path.basename(args.ids_csv)}"),
        "pdbids": pdbids,
        **docking_params(opts),
        "max_receptor_atoms": args.max_receptor_atoms,
        "rmsd": "obrms heavy-atom, symmetry-corrected, no superposition",
        "top1": "mode 1 (first MODEL)",
        "vina_version": _benchlib.tool_version("vina"),
    }
    with open(os.path.join(args.outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n=== SUMMARY ===")
    for k, v in summary.items():
        if k not in ("excluded", "protocol"):
            print(f"  {k}: {v}")
    for e in summary["excluded"]:
        print(f"  excluded {e['pdbid']}: {e['status']} ({e['reason']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
