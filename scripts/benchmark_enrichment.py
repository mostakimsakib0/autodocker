#!/usr/bin/env python3
"""J2: DUD-E enrichment test (screening power).

Protocol (screening-power evaluation):

  * Targets : curated set of DUD-E targets spanning several protein
              families — protein kinases (akt1, braf, cdk2, egfr), serine
              protease (fa10), HIV-1 aspartic protease (hivpr), nuclear
              receptor (andr), reductases (aldr, hmdh) and the esterase
              superfamily (aces)
  * Ligands : actives + decoys from DUD-E (``actives_final.mol2.gz``,
              ``decoys_final.mol2.gz``), seeded-sampled then converted to
              PDBQT with Open Babel
  * Receptor: DUD-E ``receptor.pdb`` (pocket grid via fpocket, runner
              default)
  * Docking : autodocker ``runner.py`` (Vina), fixed seed, exhaustiveness 4,
              single binding mode
  * Metrics : AUC (rank-based, tie-aware), enrichment factors EF1%/EF5%,
              95% bootstrap confidence intervals (stratified resampling of
              actives/decoys) and a one-sided Mann-Whitney U test asking
              whether actives score better than decoys.

Actives and decoys are labelled by filename prefix ``A_`` / ``D_``.

Modes
-----
* default: full pipeline (download + prep + dock with ``runner.py``)
* ``--prepare-only``: stop after ligand prep; useful for handing the exact
  same prepared ligands to a competitor tool
* ``--ingest DIR``: skip docking; convert an external ranking.csv (columns
  Ligand, Binding_Affinity|Affinity) found under ``DIR/<target>/`` into the
  standard per-target metrics.csv so results from other pipelines can be
  compared on the identical ligand set
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import random
import shutil
import subprocess
import sys
import urllib.request

DUD_E_BASE = "https://dude.docking.org/targets/{t}/{f}"

# Functional diversity across protein families so screening power is not
# claimed on a single receptor class.
DEFAULT_TARGETS = [
    # kinases
    "akt1", "braf", "cdk2", "egfr",
    # serine protease
    "fa10",
    # aspartic protease
    "hivpr",
    # nuclear hormone receptor
    "andr",
    # reductases
    "aldr", "hmdh",
    # esterase superfamily
    "aces",
]

TARGET_FAMILY = {
    "akt1": "kinase", "braf": "kinase", "cdk2": "kinase", "egfr": "kinase",
    "fa10": "serine protease", "hivpr": "aspartic protease",
    "andr": "nuclear receptor", "aldr": "reductase", "hmdh": "reductase",
    "aces": "esterase",
}


def download(url, dest):
    if os.path.exists(dest):
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "autodocker-bench"})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as out:
        shutil.copyfileobj(r, out)


def mol2_to_pdbqt(mol2_path, pdbqt_path):
    proc = subprocess.run(
        ["obabel", "-imol2", mol2_path, "-opdbqt", "-O", pdbqt_path],
        capture_output=True, text=True, timeout=300)
    if proc.returncode != 0 or not os.path.exists(pdbqt_path):
        raise RuntimeError(f"obabel mol2->pdbqt failed for {mol2_path}: {proc.stderr[-300:]}")


def read_molecules(mol2_gz):
    """Return the list of @<TRIPOS>MOLECULE blocks from a compound file."""
    with gzip.open(mol2_gz, "rt", errors="replace") as f:
        text = f.read()
    blocks = text.split("@<TRIPOS>MOLECULE")
    return [("@<TRIPOS>MOLECULE" + b) for b in blocks[1:] if b.strip()]


def sample_ligands(mol2_gz, out_dir, prefix, max_n, seed, skip_existing=True):
    """Seeded-sample up to ``max_n`` ligands and convert them to PDBQT.

    Returns the list of output PDBQT paths. Conversion failures are logged
    and skipped (never fatal).
    """
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    blocks = read_molecules(mol2_gz)
    rng = random.Random(seed)
    shuffled = list(range(len(blocks)))
    rng.shuffle(shuffled)
    count = 0
    for idx in shuffled:
        if len(paths) >= max_n:
            break
        block = blocks[idx]
        name_line = next(
            (l.strip() for l in block.splitlines()
             if l.strip() and not l.strip().startswith("@<TRIPOS>")),
            str(idx))
        name = (name_line.split()[0] or "_").replace("/", "_")
        mol2_path = os.path.join(out_dir, f"{prefix}{count}_{name}.mol2")
        pdbqt_path = mol2_path.replace(".mol2", ".pdbqt")
        if skip_existing and os.path.exists(pdbqt_path):
            paths.append(pdbqt_path)
            count += 1
            continue
        try:
            with open(mol2_path, "w") as out:
                out.write(block)
            mol2_to_pdbqt(mol2_path, pdbqt_path)
        except Exception as e:
            sys.stderr.write(f"  [!] convert failed [{prefix}{count}]: {e}\n")
            continue
        paths.append(pdbqt_path)
        count += 1
    return paths


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

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


def _ranked_with_ties(values):
    """Average ranks for a list of values (ties share the mean rank)."""
    if not values:
        return []
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while (j + 1 < len(order)
               and values[order[j + 1]] == values[order[i]]):
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def mannwhitney(active_scores, decoy_scores):
    """One- and two-sided Mann-Whitney U with normal approximation.

    Scores are binding affinities (lower = better). H1 (one-sided) =
    actives rank significantly better than decoys. Ties are corrected.
    Returns (U1, p_one_sided, p_two_sided).
    """
    n1, n2 = len(active_scores), len(decoy_scores)
    if n1 == 0 or n2 == 0:
        return None, None, None
    allv = list(active_scores) + list(decoy_scores)
    ranks = _ranked_with_ties(allv)
    r1 = sum(ranks[:n1])
    u1 = r1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    n = n1 + n2
    ties = {}
    for v in allv:
        ties.setdefault(v, 0)
        ties[v] += 1
    tie_adj = sum(t ** 3 - t for t in ties.values() if t > 1)
    var_u = (n1 * n2 / 12.0
             * ((n + 1) - tie_adj / (n * (n - 1.0))))
    if var_u <= 0:
        return u1, 1.0, 1.0
    z = (u1 - mu) / math.sqrt(var_u)

    def _norm_cdf(z):
        return 0.5 * (1.0 + math.erf(z / math.sqrt(2)))

    p_one = _norm_cdf(z)
    p_two = 2.0 * min(p_one, 1.0 - p_one)
    return u1, p_one, p_two


def _percentile(values, p):
    values = sorted(values)
    k = (len(values) - 1) * p
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return values[int(k)]
    d0 = values[int(f)] * (c - k)
    d1 = values[int(c)] * (k - f)
    return d0 + d1


def bootstrap_ci(scores, labels, n_boot, seed, ci=0.95):
    """Stratified bootstrap CIs for AUC, EF1%, EF5%.

    Each replicate resamples the observed actives with replacement and the
    observed decoys with replacement, then recomputes the metrics. The CI is
    the (ci/2, 1-ci/2) percentile band across replicates.
    """
    names = sorted(scores)
    pos = [n for n in names if labels.get(n)]
    neg = [n for n in names if not labels.get(n)]
    if not pos or not neg:
        return None
    rng = random.Random(seed)
    alpha = (1.0 - ci) / 2.0
    if n_boot <= 0:
        return None
    aucs, ef1s, ef5s = [], [], []
    for _ in range(n_boot):
        sample = [rng.choice(pos) for _ in range(len(pos))] \
            + [rng.choice(neg) for _ in range(len(neg))]
        sub = {n: scores[n] for n in sample}
        sub_labels = {n: n.startswith("A_") for n in sub}
        auc, ef1, _ = auc_and_ef(sub, sub_labels, percentile=0.01)
        _, ef5, _ = auc_and_ef(sub, sub_labels, percentile=0.05)
        if auc is not None:
            aucs.append(auc)
            ef1s.append(ef1 if ef1 is not None else 0.0)
            ef5s.append(ef5 if ef5 is not None else 0.0)
    if not aucs:
        return None
    return {
        "auc": (min(aucs), _percentile(aucs, alpha), _percentile(aucs, 1 - alpha)),
        "ef1": (min(ef1s), _percentile(ef1s, alpha), _percentile(ef1s, 1 - alpha)),
        "ef5": (min(ef5s), _percentile(ef5s, alpha), _percentile(ef5s, 1 - alpha)),
    }


# ---------------------------------------------------------------------------
# Per-target pipeline
# ---------------------------------------------------------------------------

def compute_stats(scores, seed, n_boot):
    labels = {n: n.startswith("A_") for n in scores}
    auc, ef1, total = auc_and_ef(scores, labels, percentile=0.01)
    _, ef5, _ = auc_and_ef(scores, labels, percentile=0.05)
    pos_s = [scores[n] for n in sorted(scores) if labels[n]]
    neg_s = [scores[n] for n in sorted(scores) if not labels[n]]
    _u, p1, p2 = mannwhitney(pos_s, neg_s)
    ci = bootstrap_ci(scores, labels, n_boot, seed)

    def _sig(p):
        if p is None:
            return ""
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return ""

    row = {
        "auc": round(auc, 4) if auc is not None else None,
        "ef1": round(ef1, 3) if ef1 is not None else None,
        "ef5": round(ef5, 3) if ef5 is not None else None,
        "auc_lo": round(ci["auc"][1], 4) if ci else None,
        "auc_hi": round(ci["auc"][2], 4) if ci else None,
        "ef1_lo": round(ci["ef1"][1], 3) if ci else None,
        "ef1_hi": round(ci["ef1"][2], 3) if ci else None,
        "ef5_lo": round(ci["ef5"][1], 3) if ci else None,
        "ef5_hi": round(ci["ef5"][2], 3) if ci else None,
        "p_mw_one": round(p1, 6) if p1 is not None else None,
        "p_mw_two": round(p2, 6) if p2 is not None else None,
        "sig": _sig(p1),
        "total": total,
    }
    return row


def write_metrics(target, tdir, outdir, scores):
    os.makedirs(os.path.join(outdir, target), exist_ok=True)
    metrics = os.path.join(outdir, target, "metrics.csv")
    with open(metrics, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["name", "label", "affinity"])
        for name, aff in sorted(scores.items()):
            writer.writerow([name, "active" if name.startswith("A_") else "decoy",
                             f"{aff:.4f}"])
    shutil.copy(os.path.join(tdir, "dock", "ranking.csv"),
                os.path.join(outdir, target, "ranking.csv"))
    return metrics


def process_target(target, workdir, opts):
    rec = {"target": target, "family": TARGET_FAMILY.get(target, "?"),
           "status": "error", "n_actives": 0, "n_decoys": 0,
           "error": None}
    tdir = os.path.join(workdir, target)
    target_seed = int(hashlib.md5(target.encode()).hexdigest()[:8], 16)
    try:
        receptor = os.path.join(tdir, "receptor.pdb")
        actives_gz = os.path.join(tdir, "actives_final.mol2.gz")
        decoys_gz = os.path.join(tdir, "decoys_final.mol2.gz")
        download(DUD_E_BASE.format(t=target, f="receptor.pdb"), receptor)
        download(DUD_E_BASE.format(t=target, f="actives_final.mol2.gz"),
                 actives_gz)
        download(DUD_E_BASE.format(t=target, f="decoys_final.mol2.gz"),
                 decoys_gz)

        lig_dir = os.path.join(tdir, "ligands")
        actives = sample_ligands(actives_gz, lig_dir, "A_",
                                 opts["max_actives"], seed=target_seed)
        decoys = sample_ligands(decoys_gz, lig_dir, "D_",
                                opts["max_decoys"], seed=target_seed ^ 0xDEAD)
        rec["n_actives"], rec["n_decoys"] = len(actives), len(decoys)
        if not actives or not decoys:
            raise RuntimeError("no usable actives/decoys after conversion")

        if opts["prepare_only"]:
            rec["status"] = "prepared"
            return rec

        outdir = os.path.join(tdir, "dock")
        if not opts["ingest"]:
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
                    f"{(proc.stdout or '')[-400:]}{(proc.stderr or '')[-400:]}")
        else:
            ranking = os.path.join(opts["ingest"], target, "ranking.csv")
            if not os.path.exists(ranking):
                raise RuntimeError(f"ingest ranking.csv missing: {ranking}")
            os.makedirs(outdir, exist_ok=True)
            shutil.copy(ranking, os.path.join(outdir, "ranking.csv"))

        scores = read_ranking(ranking)
        metrics = write_metrics(target, tdir, opts["outdir"], scores)
        rec.update(compute_stats(scores, target_seed, opts["bootstrap"]))
        rec["metrics"] = metrics
        rec["status"] = "ok"
    except Exception as e:
        rec["error"] = str(e)
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--targets",
                    default=",".join(DEFAULT_TARGETS),
                    help="Comma-separated DUD-E target codes "
                         f"(default: {','.join(DEFAULT_TARGETS)})")
    ap.add_argument("--max-actives", type=int, default=50,
                    help="Ligands sampled per target (default 50)")
    ap.add_argument("--max-decoys", type=int, default=100,
                    help="Decoys sampled per target (default 100)")
    ap.add_argument("--bootstrap", type=int, default=1000,
                    help="Bootstrap replicates for CI (default 1000)")
    ap.add_argument("--ci", type=float, default=0.95)
    ap.add_argument("--prepare-only", action="store_true",
                    help="Download+prep ligands only (no docking)")
    ap.add_argument("--ingest", default=None, metavar="DIR",
                    help="Do not dock; read DIR/<target>/ranking.csv instead")
    ap.add_argument("--workdir", default="/tmp/opencode/dude")
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
        print(f"\n=== [{t}] ({TARGET_FAMILY.get(t, '?')}) ===", flush=True)
        rec = process_target(t, args.workdir, {
            "max_actives": args.max_actives, "max_decoys": args.max_decoys,
            "seed": args.seed, "exhaustiveness": args.exhaustiveness,
            "processes": args.processes, "timeout": args.timeout,
            "runner": os.path.abspath(args.runner),
            "bootstrap": args.bootstrap,
            "prepare_only": args.prepare_only,
            "ingest": args.ingest, "outdir": args.outdir})
        results.append(rec)
        print(f"[{t}] {rec['status']} actives={rec['n_actives']} "
              f"decoys={rec['n_decoys']} auc={rec.get('auc')}"
              f" ({rec.get('auc_lo')}-{rec.get('auc_hi')}) "
              f"ef1={rec.get('ef1')} p_mw={rec.get('p_mw_one')} "
              f"{rec.get('sig')} err={rec['error']}", flush=True)
        # Incremental save so a partial run still leaves usable results.
        save_tables(results, args.outdir)

    save_tables(results, args.outdir)
    n_ok = [r for r in results if r["status"] == "ok"]
    if n_ok:
        aucs = [r["auc"] for r in n_ok if r["auc"] is not None]
        if aucs:
            mean = sum(aucs) / len(aucs)
            print(f"\nMean AUC over {len(aucs)} targets: {mean:.3f}")


def save_tables(results, outdir):
    rows_path = os.path.join(outdir, "enrichment.csv")
    cols = ["target", "family", "status", "n_actives", "n_decoys", "total",
            "auc", "auc_lo", "auc_hi", "ef1", "ef1_lo", "ef1_hi",
            "ef5", "ef5_lo", "ef5_hi", "p_mw_one", "p_mw_two", "sig", "error"]
    with open(rows_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            writer.writerow(r)
    with open(os.path.join(outdir, "summary.json"), "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {rows_path}")


if __name__ == "__main__":
    main()