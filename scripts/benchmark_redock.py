#!/usr/bin/env python3
"""J1: PDBbind core-set re-docking (docking power) benchmark.

Protocol (standard CASF "docking power" evaluation):

  * Receptor : PDBbind prepared ``{pdbid}_protein.pdb``
  * Ligand   : PDBbind ``{pdbid}_ligand.mol2`` (crystal-bound geometry)
  * Grid     : box centered on crystal-ligand heavy-atom centroid,
               size = max(22.5, ligand_extent + 2*padding)
  * Docking  : autodocker ``runner.py`` (Vina) with fixed seed, exhaustiveness,
               num_modes and energy_range (crystal-ligand-centered grid via
               --grid-center/--grid-size overrides)
  * Metric   : heavy-atom RMSD of docked pose vs crystal pose (``obrms``)
  * Success  : RMSD <= --rmsd-cutoff (default 2.0 A)

Results are written to ``outdir/results.csv`` and ``outdir/summary.json``.
"""

import argparse
import csv
import json
import multiprocessing
import os
import shutil
import subprocess
import sys
import tempfile


def parse_mol2_heavy_coords(mol2_path):
    coords, in_atoms = [], False
    with open(mol2_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>"):
                in_atoms = False
                continue
            if in_atoms and line:
                parts = line.split()
                if len(parts) >= 6 and not parts[5].startswith("H"):
                    coords.append((float(parts[2]), float(parts[3]),
                                   float(parts[4])))
    if not coords:
        raise ValueError("no heavy atoms in ligand")
    return coords


def ligand_grid(mol2_path, padding=6.0, min_size=22.5):
    coords = parse_mol2_heavy_coords(mol2_path)
    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)
    extent = max(
        max(c[0] for c in coords) - min(c[0] for c in coords),
        max(c[1] for c in coords) - min(c[1] for c in coords),
        max(c[2] for c in coords) - min(c[2] for c in coords),
    )
    size = max(min_size, extent + 2 * padding)
    return cx, cy, cz, size


def receptor_centroid(pdb_path):
    """Mean of all ATOM-record coordinates (matches the obabel -c centering
    the runner applies to the receptor PDBQT)."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                coords.append((float(line[30:38]), float(line[38:46]),
                               float(line[46:54])))
    if not coords:
        raise ValueError("no ATOM records in receptor PDB")
    return (sum(c[0] for c in coords) / len(coords),
            sum(c[1] for c in coords) / len(coords),
            sum(c[2] for c in coords) / len(coords))


def shift_mol2(mol2_path, shift, out_path):
    """Write a copy of a mol2 with all ATOM-section coordinates translated by
    ``shift`` (used to bring the crystal pose into the centered frame)."""
    with open(mol2_path) as src, open(out_path, "w") as dst:
        in_atoms = False
        for line in src:
            stripped = line.strip()
            if stripped.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                dst.write(line)
                continue
            if stripped.startswith("@<TRIPOS>"):
                in_atoms = False
                dst.write(line)
                continue
            if in_atoms and line.strip():
                parts = line.split()
                x = float(parts[2]) - shift[0]
                y = float(parts[3]) - shift[1]
                z = float(parts[4]) - shift[2]
                parts[2], parts[3], parts[4] = f"{x:.3f}", f"{y:.3f}", f"{z:.3f}"
                dst.write(" ".join(parts) + "\n")
            else:
                dst.write(line)


def read_ranking(ranking_csv):
    affinity = None
    if os.path.exists(ranking_csv):
        with open(ranking_csv) as f:
            for row in csv.DictReader(f):
                key = "Binding_Affinity" if "Binding_Affinity" in row else "Affinity"
                affinity = row.get(key)
                break
    return affinity


def split_vina_models(pdbqt_path, tmpdir):
    """Split all Vina MODEL blocks into per-mode PDBQT files. Returns list of
    (mode_number, affinity, file_path) sorted by affinity (mode 1 first)."""
    models = []
    current, model_no = [], 0
    affinity = None
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                current, model_no = [], model_no + 1
            elif line.startswith("ENDMDL"):
                if current:
                    path = os.path.join(tmpdir, f"mode_{model_no}.pdbqt")
                    with open(path, "w") as out:
                        out.write("\n".join(current) + "\n")
                    models.append((model_no, affinity, path))
            elif line.startswith("REMARK VINA RESULT:"):
                affinity = float(line.split()[-3])
            else:
                current.append(line.rstrip())
    models.sort(key=lambda m: m[1] if m[1] is not None else 0.0)
    return models


def obrms(ref, query):
    proc = subprocess.run(
        ["obrms", ref, query], capture_output=True, text=True, timeout=120)
    if proc.returncode != 0:
        raise RuntimeError(f"obrms failed: {proc.stderr.strip()}")
    for tok in proc.stdout.split():
        try:
            return float(tok)
        except ValueError:
            continue
    raise RuntimeError(f"no RMSD in obrms output: {proc.stdout!r}")


def process_complex(job):
    pdbid, data_dir, workdir, runner_script, opts = job
    rec = {
        "pdbid": pdbid,
        "status": "error",
        "affinity": None,
        "rmsd_top1": None,
        "rmsd_best": None,
        "error": None,
    }
    complex_dir = os.path.join(data_dir, pdbid)
    protein = os.path.join(complex_dir, f"{pdbid}_protein.pdb")
    ligand = os.path.join(complex_dir, f"{pdbid}_ligand.mol2")
    outdir = os.path.join(workdir, pdbid)
    try:
        if opts["max_receptor_atoms"]:
            n_atom = sum(1 for l in open(protein) if l.startswith("ATOM"))
            if n_atom > opts["max_receptor_atoms"]:
                rec["status"] = "skipped"
                rec["error"] = (f"receptor too large ({n_atom} ATOM > "
                                f"{opts['max_receptor_atoms']})")
                return rec
        cx, cy, cz, size = ligand_grid(ligand, padding=opts["padding"])
        rc = receptor_centroid(protein)
        # The runner centers the receptor at the origin (obabel -c); express the
        # crystal-ligand-centered box in that centered frame.
        cx, cy, cz = cx - rc[0], cy - rc[1], cz - rc[2]
        crystal_shifted = os.path.join(outdir, f"{pdbid}_crystal_centered.mol2")
        os.makedirs(outdir, exist_ok=True)
        shift_mol2(ligand, rc, crystal_shifted)
        cmd = [
            sys.executable, runner_script,
            "-r", protein, "-l", ligand, "-o", outdir,
            "--chain=all",
            "--no-fpocket",
            f"--grid-center={cx:.3f},{cy:.3f},{cz:.3f}",
            f"--grid-size={size:.1f},{size:.1f},{size:.1f}",
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
        best_pdbqt = os.path.join(outdir, "docked",
                                  f"{pdbid}_ligand_out.pdbqt")
        if proc.returncode != 0 or not os.path.exists(best_pdbqt):
            raise RuntimeError(
                f"docking failed (rc={proc.returncode}): "
                f"{(proc.stdout or '')[-400:]} {(proc.stderr or '')[-400:]}")
        rec["affinity"] = read_ranking(
            os.path.join(outdir, "ranking.csv"))
        with tempfile.TemporaryDirectory() as tmp:
            modes = split_vina_models(best_pdbqt, tmp)
            rmsds = []
            for _, _, mode_path in modes:
                try:
                    rmsds.append(obrms(crystal_shifted, mode_path))
                except Exception as e:
                    rec["error"] = f"obrms: {e}"
            if rmsds:
                rec["rmsd_top1"] = rmsds[0]
                rec["rmsd_best"] = min(rmsds)
        rec["status"] = "ok"
    except Exception as e:
        rec["error"] = str(e)
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True,
                    help="Extracted PDBbind core-set dir (dirs per pdbid)")
    ap.add_argument("--pdbids", help="Comma-separated PDB IDs (else default subset)")
    ap.add_argument("--max-complexes", type=int, default=40,
                    help="Use first N sorted PDB IDs when --pdbids is empty")
    ap.add_argument("--outdir", default="benchmarks/results/redock",
                    help="Results directory")
    ap.add_argument("--runner", default="autodocker/runner.py",
                    help="Path to runner.py")
    ap.add_argument("--processes", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--binding-modes", type=int, default=9)
    ap.add_argument("--energy-range", type=float, default=3.0)
    ap.add_argument("--padding", type=float, default=6.0)
    ap.add_argument("--threads", type=int, default=0,
                    help="Vina threads per complex (0 = let Vina decide; "
                         "only honored by Vina builds that support --threads)")
    ap.add_argument("--max-receptor-atoms", type=int, default=6000,
                    help="Skip receptors larger than this many ATOM records "
                         "(large oligomeric assemblies dominate runtime)")
    ap.add_argument("--rmsd-cutoff", type=float, default=2.0)
    ap.add_argument("--timeout", type=int, default=1800)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if args.pdbids:
        pdbids = [p.strip() for p in args.pdbids.split(",") if p.strip()]
        missing = [p for p in pdbids if not os.path.isdir(os.path.join(args.data_dir, p))]
        if missing:
            raise SystemExit(f"missing complexes in data dir: {missing}")
    else:
        pdbids = sorted(os.listdir(args.data_dir))[:args.max_complexes]

    # Resume: skip complexes already docked successfully; drop partial dirs.
    workroot = os.path.join(args.outdir, "work")
    todo, done = [], []
    for p in pdbids:
        out_pdbqt = os.path.join(workroot, p, "docked", f"{p}_ligand_out.pdbqt")
        if os.path.exists(out_pdbqt):
            done.append(p)
        else:
            shutil.rmtree(os.path.join(workroot, p), ignore_errors=True)
            todo.append(p)
    if done:
        print(f"[*] Skipping {len(done)} already-docked complexes: {', '.join(done)}")

    jobs = [(p, args.data_dir, workroot,
             os.path.abspath(args.runner),
             {"padding": args.padding, "seed": args.seed,
              "exhaustiveness": args.exhaustiveness,
              "binding_modes": args.binding_modes,
              "energy_range": args.energy_range,
              "threads": args.threads,
              "max_receptor_atoms": args.max_receptor_atoms,
              "timeout": args.timeout})
            for p in todo]

    results = []
    with multiprocessing.Pool(processes=args.processes) as pool:
        for rec in pool.imap_unordered(process_complex, jobs):
            results.append(rec)
            print(f"[{rec['pdbid']}] {rec['status']} "
                  f"aff={rec['affinity']} rmsd_top1={rec['rmsd_top1']} "
                  f"rmsd_best={rec['rmsd_best']} err={rec['error']}",
                  flush=True)

    # Merge previously-completed results back in for the final report.
    for p in done:
        results.append({"pdbid": p, "status": "ok", "affinity": None,
                        "rmsd_top1": None, "rmsd_best": None, "error": None})
        best_pdbqt = os.path.join(workroot, p, "docked", f"{p}_ligand_out.pdbqt")
        try:
            rc = receptor_centroid(
                os.path.join(args.data_dir, p, f"{p}_protein.pdb"))
            crystal_ref = os.path.join(workroot, p, f"{p}_crystal_centered.mol2")
            here = os.path.join(args.data_dir, p, f"{p}_ligand.mol2")
            if not os.path.exists(crystal_ref):
                shift_mol2(here, rc, crystal_ref)
            with tempfile.TemporaryDirectory() as tmp:
                modes = split_vina_models(best_pdbqt, tmp)
                rmsds = [obrms(crystal_ref, m[2]) for m in modes]
            results[-1]["rmsd_top1"] = rmsds[0]
            results[-1]["rmsd_best"] = min(rmsds)
            results[-1]["affinity"] = read_ranking(
                os.path.join(workroot, p, "ranking.csv"))
        except Exception as e:
            results[-1]["error"] = f"resume: {e}"

    results.sort(key=lambda r: r["pdbid"])
    rows_path = os.path.join(args.outdir, "results.csv")
    with open(rows_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pdbid", "status", "affinity",
                         "rmsd_top1", "rmsd_best", "error"])
        for r in results:
            writer.writerow([r["pdbid"], r["status"], r["affinity"],
                             (f"{r['rmsd_top1']:.3f}" if r["rmsd_top1"] is not None else ""),
                             (f"{r['rmsd_best']:.3f}" if r["rmsd_best"] is not None else ""),
                             r["error"]])

    ok = [r for r in results if r["rmsd_top1"] is not None]
    best = [r for r in results if r["rmsd_best"] is not None]
    summary = {
        "total": len(results),
        "skipped": len([r for r in results if r["status"] == "skipped"]),
        "failed": len([r for r in results if r["status"] == "error"]),
        "docked": len(ok),
        "success_top1": len([r for r in ok if r["rmsd_top1"] <= args.rmsd_cutoff]),
        "success_top1_frac": (len([r for r in ok if r["rmsd_top1"] <= args.rmsd_cutoff]) / len(ok)) if ok else 0.0,
        "success_best": len([r for r in best if r["rmsd_best"] <= args.rmsd_cutoff]),
        "success_best_frac": (len([r for r in best if r["rmsd_best"] <= args.rmsd_cutoff]) / len(best)) if best else 0.0,
        "median_rmsd_top1": sorted(r["rmsd_top1"] for r in ok)[len(ok) // 2] if ok else None,
        "median_rmsd_best": sorted(r["rmsd_best"] for r in best)[len(best) // 2] if best else None,
        "rmsd_cutoff": args.rmsd_cutoff,
        "seed": args.seed,
        "exhaustiveness": args.exhaustiveness,
        "binding_modes": args.binding_modes,
        "energy_range": args.energy_range,
    }
    with open(os.path.join(args.outdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== SUMMARY ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()