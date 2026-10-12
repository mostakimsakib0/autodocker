"""Seed replicates (--replicates N): dock every ligand N times.

AutoDock Vina's search is stochastic; one seed gives one sample of it.
With --replicates N, replicate k (k = 0..N-1) of every ligand is docked with
seed --seed + k, and the screen reports per ligand:

* the best (lowest) in-window score of the N replicates as Binding_Affinity,
  with that replicate's pose file, SimScore and modes (the reported pose is
  the one the reported score belongs to);
* n_replicates (replicates with an in-window score), score_mean, score_sd
  (sample SD, n - 1; blank below two scores), score_min, score_max.

Replicate k of ligand <name> writes docked/<name>_rep<k>_out.pdbqt and
docked/<name>_rep<k>_vina.log and is checkpointed under the key
<name>_rep<k>; a resumed run re-docks only the missing (ligand, replicate)
pairs. dock_ligand derives every output name from the ligand file name, so
replicate k docks a byte-identical copy of the ligand named <name>_rep<k>.pdbqt
in docked/replicate_inputs/ (stage_ligand): the docking worker, its log and
the checkpoint need no replicate logic.

N = 1 is the plain run: nothing here is called, and every output, file
name, log line and the checkpoint fingerprint stay byte-identical to a run
without the option (see the wiring in docs/features_2_0.md).

Imports only modules that do not import docking, so docking.dock_all can use
it; docking and consensus are imported where they are needed.
"""
import csv
import json
import math
import os
import re
import shutil
import statistics
from dataclasses import replace
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from vspipeline import defaults, outputs, utils

REPLICATES = defaults.REPLICATES
# Upper bound of --replicates: a typo guard (each replicate is a full docking).
MAX_REPLICATES = 100

REPLICATE_TAG = outputs.REPLICATE_TAG
REPLICATE_INPUTS_DIR = outputs.REPLICATE_INPUTS_DIR
REPLICATE_SCORES_CSV = outputs.REPLICATE_SCORES_CSV
REPLICATE_SUMMARY_JSON = outputs.REPLICATE_SUMMARY_JSON
RANKING_REPLICATE_COLUMNS = outputs.RANKING_REPLICATE_COLUMNS
REPLICATE_SCORES_COLUMNS = ("Ligand", "Replicate", "Seed", "Status", "Score", "Raw_Score",
                            "Reported", "Docked_File", "Log_File")

# Digits of the mean/SD written to ranking.csv (Vina prints 3 decimals).
STAT_DIGITS = 3

_STEM_RE = re.compile(r"^(?P<name>.+)" + re.escape(REPLICATE_TAG) + r"(?P<k>0|[1-9][0-9]*)$")

# Status values of docking.Status (kept as strings: no docking import here).
_OK, _WEAK, _OUT, _FAILED = "OK", "WEAK", "OUT_OF_WINDOW", "FAILED"


# ---------------------------------------------------------------------------
# Seeds, names, checkpoint keys
# ---------------------------------------------------------------------------

def check_replicates(n: int, seed: int, max_seed: int = 2**31 - 1) -> None:
    """Validate --replicates N against --seed (raises ValueError): 1 <= N <=
    MAX_REPLICATES, and the last seed (seed + N - 1) must stay a valid
    engine seed."""
    if isinstance(n, bool) or not isinstance(n, int):
        raise ValueError(f"--replicates must be an integer, got {n!r}")
    if not 1 <= n <= MAX_REPLICATES:
        raise ValueError(f"--replicates must be between 1 and {MAX_REPLICATES}, got {n}")
    if seed + n - 1 > max_seed:
        raise ValueError(
            f"--replicates {n} with --seed {seed}: the last replicate's seed "
            f"{seed + n - 1} exceeds {max_seed}; use a smaller --seed")


def replicate_seed(seed: int, k: int) -> int:
    """Engine seed of replicate k: seed + k (replicate 0 = the run's seed)."""
    if k < 0:
        raise ValueError(f"replicate index must be >= 0, got {k}")
    return seed + k


def replicate_seeds(seed: int, n: int) -> List[int]:
    return [replicate_seed(seed, k) for k in range(n)]


def replicate_stem(name: str, k: int) -> str:
    """Output stem and checkpoint key of replicate k of ligand ``name``:
    <name>_rep<k>. Injective: the last '_rep<digits>' is the replicate."""
    if k < 0:
        raise ValueError(f"replicate index must be >= 0, got {k}")
    return f"{name}{REPLICATE_TAG}{k}"


checkpoint_key = replicate_stem


def parse_replicate_stem(stem: str) -> Optional[Tuple[str, int]]:
    """(name, k) of a replicate stem, or None if ``stem`` is not one."""
    m = _STEM_RE.match(stem)
    return (m.group("name"), int(m.group("k"))) if m else None


def pose_file(dock_dir: str, name: str, k: int) -> str:
    return os.path.join(dock_dir, replicate_stem(name, k) + outputs.POSE_SUFFIX)


def log_file(dock_dir: str, name: str, k: int) -> str:
    return os.path.join(dock_dir, replicate_stem(name, k) + outputs.VINA_LOG_SUFFIX)


def replicate_params(vina_params: Dict, k: Optional[int]) -> Dict:
    """``vina_params`` for replicate k (a copy with seed + k); the same dict
    when k is None (no replicates)."""
    if k is None:
        return vina_params
    return {**vina_params, "seed": replicate_seed(int(vina_params["seed"]), k)}


def fingerprint_params(params: Dict, n: int) -> Dict:
    """Docking parameters hashed into the run fingerprint: ``params`` itself
    for N = 1 (the fingerprint of a plain run is unchanged), else a copy with
    "replicates": N. Changing N re-docks everything: the replicate files of
    N = 2 are not the files of N = 1, and the reported score depends on N."""
    if n <= 1:
        return params
    return {**params, "replicates": n}


# ---------------------------------------------------------------------------
# Jobs
# ---------------------------------------------------------------------------

def _same_bytes(a: str, b: str) -> bool:
    if not os.path.exists(b) or os.path.getsize(a) != os.path.getsize(b):
        return False
    with open(a, "rb") as fa, open(b, "rb") as fb:
        return fa.read() == fb.read()


def stage_ligand(ligand: str, k: int, dock_dir: str) -> str:
    """Byte-identical copy of ``ligand`` named <name>_rep<k>.pdbqt in
    <dock_dir>/replicate_inputs/ (rewritten only when the source changed, so
    its checkpoint hash matches the source's). Returns its path."""
    stage_dir = os.path.join(dock_dir, REPLICATE_INPUTS_DIR)
    os.makedirs(stage_dir, exist_ok=True)
    staged = os.path.join(stage_dir, replicate_stem(utils.ligand_name(ligand), k) + ".pdbqt")
    if not _same_bytes(ligand, staged):
        tmp = f"{staged}.tmp.{os.getpid()}"
        shutil.copyfile(ligand, tmp)
        os.replace(tmp, staged)
    return staged


def expand_ligands(ligands: Sequence[str], n: int, dock_dir: str
                   ) -> List[Tuple[str, Optional[int]]]:
    """Docking units of a screen: [(ligand file, replicate index)] in library
    order, replicates of a ligand together. N = 1: [(ligand, None)] for every
    ligand, nothing staged. N > 1: the staged copy of every (ligand, k)."""
    if n <= 1:
        return [(lig, None) for lig in ligands]
    return [(stage_ligand(lig, k, dock_dir), k) for lig in ligands for k in range(n)]


def replicate_job(job, k: int):
    """The docking.DockJob of replicate k of ``job``: the staged ligand copy
    (so dock_ligand writes <name>_rep<k>_out.pdbqt/_vina.log) and seed + k."""
    staged = stage_ligand(job.ligand, k, job.dock_dir)
    return replace(job, ligand=staged, vina_params=replicate_params(job.vina_params, k))


def dock_replicates(job, n: int) -> List[Tuple[str, Optional[float], Dict]]:
    """Dock replicates 0..n-1 of one ligand in this process (docking.dock_ligand);
    returns their (stem, score, metrics) results in replicate order."""
    from vspipeline import docking
    return [docking.dock_ligand(replicate_job(job, k)) for k in range(n)]


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def _finite(value) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def _status(score: Optional[float], metrics: Dict) -> str:
    status = (metrics or {}).get("status")
    if status in (_OK, _WEAK, _OUT, _FAILED):
        if status in (_OK, _WEAK) and score is None:
            return _FAILED
        return status
    return _FAILED if score is None else _OK


def score_stats(scores: Sequence[float]) -> Dict[str, Optional[float]]:
    """n, mean, sample SD (None below 2 scores), min and max of ``scores``."""
    scores = [float(s) for s in scores]
    if not scores:
        return {"n": 0, "mean": None, "sd": None, "min": None, "max": None}
    return {"n": len(scores),
            "mean": statistics.fmean(scores),
            "sd": statistics.stdev(scores) if len(scores) > 1 else None,
            "min": min(scores), "max": max(scores)}


def aggregate_ligand(name: str, replicate_results: Dict[int, Tuple[Optional[float], Dict]],
                     n: int, seed: int) -> Tuple[Optional[float], Dict]:
    """Combine the replicates of one ligand into its reported (score, metrics).

    ``replicate_results`` maps k -> (score, metrics) as dock_ligand returned
    them (score None for OUT_OF_WINDOW and FAILED); missing k count as FAILED.

    * Some replicate in-window: the reported score is the lowest in-window
      score (ties: the lowest k); metrics are that replicate's (pose file,
      log, SimScore, modes, engine, status) plus the replicate statistics.
    * None in-window, some OUT_OF_WINDOW: OUT_OF_WINDOW, with the metrics of
      the lowest such k (no statistics: they cover in-window scores only).
    * All FAILED: FAILED, with the error of every replicate.
    """
    rows = []
    for k in range(n):
        score, metrics = replicate_results.get(
            k, (None, {"status": _FAILED, "error": "no result"}))
        metrics = metrics or {}
        rows.append({
            "replicate": k, "seed": replicate_seed(seed, k),
            "status": _status(score, metrics),
            "score": score if _finite(score) else None,
            "raw_score": metrics.get("raw_score"),
            "dock_file": metrics.get("dock_file"),
            "log_file": metrics.get("log_file"),
            "error": metrics.get("error"),
            "_metrics": metrics,
        })
    valid = [r for r in rows if r["status"] in (_OK, _WEAK) and r["score"] is not None]
    out_of_window = [r for r in rows if r["status"] == _OUT]
    failed = [r for r in rows if r["status"] == _FAILED]

    if valid:
        best = min(valid, key=lambda r: (r["score"], r["replicate"]))
        reported_score: Optional[float] = best["score"]
    elif out_of_window:
        best = out_of_window[0]
        reported_score = None
    else:
        best = rows[0]
        reported_score = None

    stats = score_stats([r["score"] for r in valid])
    metrics = dict(best["_metrics"])
    if not valid and not out_of_window:
        metrics["status"] = _FAILED
        metrics["error"] = "; ".join(
            f"replicate {r['replicate']} (seed {r['seed']}): {r['error'] or 'failed'}"
            for r in rows)
    metrics.update({
        "replicates_requested": n,
        "best_replicate": best["replicate"] if valid else None,
        "best_seed": best["seed"] if valid else None,
        "n_replicates": stats["n"],
        "n_replicates_failed": len(failed),
        "n_replicates_out_of_window": len(out_of_window),
        "score_mean": stats["mean"],
        "score_sd": stats["sd"],
        "score_min": stats["min"],
        "score_max": stats["max"],
        "replicates": [{k: v for k, v in r.items() if k != "_metrics"} for r in rows],
    })
    return reported_score, metrics


def aggregate(ligand_names: Sequence[str], scores: Dict[str, Optional[float]],
              metrics: Dict[str, Dict], n: int, seed: int
              ) -> Tuple[List[Tuple[str, Optional[float]]], Dict[str, Dict]]:
    """Per-ligand (results, metrics_dict), in library order, from the
    per-replicate results keyed by replicate stem (what dock_all and the
    checkpoint hold). Replicates absent from both maps count as FAILED."""
    results: List[Tuple[str, Optional[float]]] = []
    metrics_dict: Dict[str, Dict] = {}
    for name in ligand_names:
        per_rep = {}
        for k in range(n):
            stem = replicate_stem(name, k)
            if stem in scores or stem in metrics:
                per_rep[k] = (scores.get(stem), metrics.get(stem, {}))
        score, m = aggregate_ligand(name, per_rep, n, seed)
        results.append((name, score))
        metrics_dict[name] = m
    return results, metrics_dict


# ---------------------------------------------------------------------------
# Output rows
# ---------------------------------------------------------------------------

def _cell(value, digits: Optional[int] = None):
    if value is None:
        return ""
    return round(float(value), digits) if digits is not None else value


def ranking_values(metrics: Dict) -> List:
    """ranking.csv cells for RANKING_REPLICATE_COLUMNS (blank when unknown).
    min/max are engine scores as written; mean and SD are rounded to
    STAT_DIGITS decimals."""
    metrics = metrics or {}
    return [metrics.get("n_replicates", 0),
            _cell(metrics.get("score_mean"), STAT_DIGITS),
            _cell(metrics.get("score_sd"), STAT_DIGITS),
            _cell(metrics.get("score_min")),
            _cell(metrics.get("score_max"))]


def has_replicates(metrics: Optional[Dict]) -> bool:
    return bool(metrics) and "replicates_requested" in metrics


def top_hits_line(metrics: Dict) -> Optional[str]:
    """'Replicates: 3/3 scored; mean -7.10, SD 0.12, range -7.25 to -6.98
    (reported: replicate 1, seed 43)' for Top_hits.txt; None without
    replicates."""
    if not has_replicates(metrics):
        return None
    n_ok, n = metrics["n_replicates"], metrics["replicates_requested"]
    text = f"Replicates: {n_ok}/{n} scored"
    if n_ok:
        text += f"; mean {metrics['score_mean']:.2f}"
        if metrics.get("score_sd") is not None:
            text += f", SD {metrics['score_sd']:.2f}"
        text += f", range {metrics['score_min']:.2f} to {metrics['score_max']:.2f}"
        text += (f" (reported: replicate {metrics['best_replicate']}, "
                 f"seed {metrics['best_seed']})")
    return text


def _relative(path: Optional[str], outdir: str) -> str:
    if not path:
        return ""
    return os.path.relpath(os.path.abspath(path), os.path.abspath(outdir))


def write_replicate_scores(path: str, outdir: str, ligand_names: Iterable[str],
                           metrics_dict: Dict[str, Dict]) -> None:
    """replicate_scores.csv: one row per (ligand, replicate), ligands in
    ``ligand_names`` order (pass them sorted for a stable file)."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(REPLICATE_SCORES_COLUMNS)
        for name in ligand_names:
            m = metrics_dict.get(name) or {}
            for r in m.get("replicates", []):
                writer.writerow([
                    name, r["replicate"], r["seed"], r["status"],
                    _cell(r["score"]), _cell(r["raw_score"]),
                    "yes" if r["replicate"] == m.get("best_replicate") else "no",
                    _relative(r["dock_file"], outdir), _relative(r["log_file"], outdir),
                ])


# ---------------------------------------------------------------------------
# Run-level summary
# ---------------------------------------------------------------------------

def replicate_scores(metrics_dict: Dict[str, Dict], k: int,
                     names: Sequence[str]) -> List[Optional[float]]:
    """In-window score of replicate k of every ligand of ``names`` (None when
    it has none)."""
    out = []
    for name in names:
        reps = (metrics_dict.get(name) or {}).get("replicates", [])
        row = next((r for r in reps if r["replicate"] == k), None)
        ok = row is not None and row["status"] in (_OK, _WEAK) and row["score"] is not None
        out.append(row["score"] if ok else None)
    return out


def summary(metrics_dict: Dict[str, Dict], n: int, seed: int) -> Dict:
    """run_manifest.json / replicate_summary.json entry of a screen with N > 1.

    * median_sd: median over ligands with >= 2 in-window replicates of the
      per-ligand sample SD of their scores (kcal/mol);
    * spearman_rep0_rep1: Spearman rank correlation of replicate 0's and
      replicate 1's scores over the ligands scored in both (consensus.py's
      tie-aware implementation; None when undefined) and its n;
    * best_replicate_counts: how often each replicate gave the reported score.
    """
    from vspipeline.consensus import scorer_agreement_spearman
    names = sorted(metrics_dict)
    sds = [m["score_sd"] for m in (metrics_dict[nm] for nm in names)
           if m.get("score_sd") is not None]
    ranges = [m["score_max"] - m["score_min"] for m in (metrics_dict[nm] for nm in names)
              if m.get("n_replicates", 0) >= 2]
    rho, n_pairs = scorer_agreement_spearman(replicate_scores(metrics_dict, 0, names),
                                             replicate_scores(metrics_dict, 1, names))
    counts = [0] * n
    for nm in names:
        best = metrics_dict[nm].get("best_replicate")
        if best is not None and 0 <= best < n:
            counts[best] += 1
    return {
        "replicates": n,
        "seeds": replicate_seeds(seed, n),
        "reported_score": "best (lowest) in-window score of the replicates; its pose "
                          "file is the reported pose",
        "sd_definition": "sample standard deviation (n - 1) of a ligand's in-window "
                         "replicate scores",
        "ligands": len(names),
        "ligands_with_sd": len(sds),
        "ligands_all_replicates_scored": sum(
            1 for nm in names if metrics_dict[nm].get("n_replicates") == n),
        "ligands_with_failed_replicates": sorted(
            nm for nm in names if metrics_dict[nm].get("n_replicates_failed")),
        "median_sd": statistics.median(sds) if sds else None,
        "max_sd": max(sds) if sds else None,
        "median_range": statistics.median(ranges) if ranges else None,
        "spearman_rep0_rep1": {"rho": rho, "n": n_pairs},
        "best_replicate_counts": counts,
    }


def write_summary(path: str, data: Dict) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")


def log_lines(data: Dict) -> List[str]:
    """INFO lines for the end of the docking stage."""
    rho = data["spearman_rep0_rep1"]
    lines = [f"[*] Replicates: {data['replicates']} per ligand (seeds "
             f"{data['seeds'][0]}..{data['seeds'][-1]}); reported score = best replicate"]
    if data["median_sd"] is not None:
        lines.append(f"[*] Replicate score SD: median {data['median_sd']:.3f}, max "
                     f"{data['max_sd']:.3f} kcal/mol over {data['ligands_with_sd']} ligands")
    if rho["rho"] is not None:
        lines.append(f"[*] Spearman rho, replicate 0 vs 1: {rho['rho']:.3f} (n={rho['n']})")
    if data["ligands_with_failed_replicates"]:
        lines.append(f"[!] {len(data['ligands_with_failed_replicates'])} ligand(s) have "
                     "failed replicates (scored from the others; re-docked on resume)")
    return lines
