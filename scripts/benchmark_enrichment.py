#!/usr/bin/env python3
"""J2: DUD-E enrichment test.

Protocol (screening-power subset evaluation):

  * Targets : small curated set of DUD-E targets (default akt1, braf, aldr)
  * Ligands : sampled actives + decoys from DUD-E (``actives_final.mol2.gz``,
              ``decoys_final.mol2.gz``), converted to PDBQT
  * Receptor: DUD-E ``receptor.pdb`` (pocket grid via fpocket, runner default)
  * Docking : autodocker ``runner.py`` (Vina), fixed seed, exhaustiveness 4,
              single binding mode
  * Metrics : AUC (rank-based, ties averaged) and enrichment factor EF1%
              (and EF5%) computed from predicted affinities vs true labels.

Actives and decoys are labelled by filename prefix ``A_`` / ``D_``.
"""

import argparse
import csv
import gzip
import json
import math
import os
import shutil
import subprocess
import sys
import urllib.request

DUD_E_BASE = "https://dude.docking.org/targets/{t}/{f}"


def download(url, dest):
    if os.path.exists(dest):
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "autodocker-bench"})
    with urllib.request.urlopen(req) as r, open(dest, "wb") as out:
        shutil.copyfileobj(r, out)


def mol2_to_pdbqt(mol2_path, pdbqt_path):
    proc = subprocess.run(
        ["obabel", "-imol2", mol2_path, "-opdbqt", "-O", pdbqt_path],
        capture_output=True, text=True, timeout=300)
    if proc.returncode != 0 or not os.path.exists(pdbqt_path):
        raise RuntimeError(f"obabel mol2->pdbqt failed for {mol2_path}: {proc.stderr[-300:]}")


def sample_ligands(mol2_gz, out_dir, prefix, max_n):
    """Decompress a *_final.mol2.gz and split into per-ligand mol2 files.
    Returns list of output PDBQT paths."""
    os.makedirs(out_dir, exist_ok=True)
    with gzip.open(mol2_gz, "rt", errors="replace") as f:
        text = f.read()
    molecules = text.split("@<TRIPOS>MOLECULE")
    paths = []
    count = 0
    for block in molecules[1:]:
        lines = block.splitlines()
        name_line = next((l.strip() for l in lines if l.strip()), str(count))
        name = name_line.split()[0]
        mol_text = "@<TRIPOS>MOLECULE" + block
        mol2_path = os.path.join(out_dir, f"{prefix}{count}_{name}.mol2")
        with open(mol2_path, "w") as out:
            out.write(mol_text)
        pdbqt_path = mol2_path.replace(".mol2", ".pdbqt")
        try:
            mol2_to_pdbqt(mol2_path, pdbqt_path)
        except Exception:
            continue
        paths.append(pdbqt_path)
        count += 1
        if max_n and count >= max_n:
            break
    return paths


def read_ranking(ranking_csv):
    scores = {}
    with open(ranking_csv) as f:
        for row in csv.DictReader(f):
            name = row.get("Ligand", "").replace(".pdbqt", "")
            key = "Binding_Affinity" if "Binding_Affinity" in row else "Affinity"
            try:
                scores[name] = float(row[key])
            except (TypeError, ValueError):
                continue
    return scores


def auc_and_ef(scores, labels, percentile=0.01):
    """scores: dict name->affinity (lower=better); labels: dict name->bool(active).

    AUC = P(random active ranks better than random decoy).
    1.0 = perfect enrichment, 0.5 = random, 0.0 = inverted.
    Ties contribute 0.5 to the pairwise comparison.
    """
    names = sorted(scores)
    n_pos = sum(1 for n in names if labels.get(n))
    n_neg = len(names) - n_pos
    if not n_pos or not n_neg:
        return None, None, None

    by_score = {}
    for n in names:
        by_score.setdefault(scores[n], []).append(n)

    pairs = 0.0
    dec_before = 0
    for grp in sorted(by_score):
        group = by_score[grp]
        act_in_grp = sum(1 for n in group if labels.get(n))
        dec_in_grp = len(group) - act_in_grp
        dec_after = n_neg - dec_before - dec_in_grp
        pairs += act_in_grp * dec_after + 0.5 * act_in_grp * dec_in_grp
        dec_before += dec_in_grp
    auc = pairs / (n_pos * n_neg)

    ordered = sorted(names, key=lambda n: scores[n])
    n_top = max(1, math.ceil(percentile * len(ordered)))
    top_names = ordered[:n_top]
    tp = sum(1 for n in top_names if labels.get(n))
    ef = (tp / n_pos) / percentile if n_pos else None
    return auc, ef, len(ordered)


def process_target(target, workdir, opts):
    rec = {"target": target, "status": "error", "n_actives": 0,
           "n_decoys": 0, "auc": None, "ef1": None, "ef5": None,
           "error": None}
    tdir = os.path.join(workdir, target)
    try:
        receptor = os.path.join(tdir, "receptor.pdb")
        actives_gz = os.path.join(tdir, "actives_final.mol2.gz")
        decoys_gz = os.path.join(tdir, "decoys_final.mol2.gz")
        download(DUD_E_BASE.format(t=target, f="receptor.pdb"), receptor)
        download(DUD_E_BASE.format(t=target, f="actives_final.mol2.gz"), actives_gz)
        download(DUD_E_BASE.format(t=target, f="decoys_final.mol2.gz"), decoys_gz)

        lig_dir = os.path.join(tdir, "ligands")
        actives = sample_ligands(actives_gz, lig_dir, "A_", opts["max_actives"])
        decoys = sample_ligands(decoys_gz, lig_dir, "D_", opts["max_decoys"])
        rec["n_actives"], rec["n_decoys"] = len(actives), len(decoys)
        if not actives or not decoys:
            raise RuntimeError("no usable actives/decoys after conversion")

        outdir = os.path.join(tdir, "dock")
        cmd = [
            sys.executable, opts["runner"],
            "-r", receptor, "-l", lig_dir, "-o", outdir,
            "--seed", str(opts["seed"]),
            "--exhaustiveness", str(opts["exhaustiveness"]),
            "--binding-modes", "1",
            "--energy-range", "1.0",
            "--processes", str(opts["processes"]),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=opts["timeout"])
        ranking = os.path.join(outdir, "ranking.csv")
        if proc.returncode != 0 or not os.path.exists(ranking):
            raise RuntimeError(
                f"docking failed (rc={proc.returncode}): "
                f"{(proc.stdout or '')[-300:]}{(proc.stderr or '')[-300:]}")
        scores = read_ranking(ranking)
        labels = {n: n.startswith("A_") for n in scores}
        auc, ef1, _ = auc_and_ef(scores, labels, percentile=0.01)
        _, ef5, _ = auc_and_ef(scores, labels, percentile=0.05)
        rec.update({"auc": auc, "ef1": ef1, "ef5": ef5})
        rec["status"] = "ok"
    except Exception as e:
        rec["error"] = str(e)
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--targets", default="akt1,braf,aldr",
                    help="Comma-separated DUD-E target codes")
    ap.add_argument("--max-actives", type=int, default=50)
    ap.add_argument("--max-decoys", type=int, default=100)
    ap.add_argument("--workdir", default="/tmp/autodocker_dude")
    ap.add_argument("--outdir", default="benchmarks/results/enrichment")
    ap.add_argument("--runner", default="autodocker/runner.py")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--exhaustiveness", type=int, default=4)
    ap.add_argument("--processes", type=int, default=8)
    ap.add_argument("--timeout", type=int, default=7200)
    args = ap.parse_args()

    targets = [t.strip() for t in args.targets.split(",") if t.strip()]
    os.makedirs(args.outdir, exist_ok=True)
    results = []
    for t in targets:
        rec = process_target(t, args.workdir, {
            "max_actives": args.max_actives, "max_decoys": args.max_decoys,
            "seed": args.seed, "exhaustiveness": args.exhaustiveness,
            "processes": args.processes,
            "timeout": args.timeout, "runner": os.path.abspath(args.runner)})
        results.append(rec)
        print(f"[{t}] {rec['status']} actives={rec['n_actives']} "
              f"decoys={rec['n_decoys']} auc={rec['auc']} "
              f"ef1={rec['ef1']} ef5={rec['ef5']} err={rec['error']}",
              flush=True)

    rows_path = os.path.join(args.outdir, "enrichment.csv")
    with open(rows_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["target", "status", "n_actives", "n_decoys",
                         "auc", "ef1", "ef5", "error"])
        for r in results:
            writer.writerow([r["target"], r["status"], r["n_actives"],
                             r["n_decoys"], r["auc"], r["ef1"], r["ef5"],
                             r["error"]])
    with open(os.path.join(args.outdir, "summary.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {rows_path}")


if __name__ == "__main__":
    main()