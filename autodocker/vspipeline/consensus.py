#!/usr/bin/env python3
"""PHASE 4: CONSENSUS SCORING + SMINA DOCKING.

Split out of runner.py (behavior-preserving refactor). All cross-module
references go through the ``runner`` namespace at call time so the public
entry point and its monkeypatch surface stay unchanged.
"""
import os
import re
import csv
import json
import math
import time
import shlex
import uuid
import shutil
import fnmatch
import glob
import html
import statistics
import subprocess
import logging
from multiprocessing import Pool
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import List, Tuple, Optional, Dict
from urllib.parse import quote
from multiprocessing import Pool, cpu_count as mp_cpu_count

# The runner module owns the tool-path globals (OBABEL, VINA*, SMINA,
# COMMAND_TIMEOUT) and every public name; import it lazily-safe: during the
# initial import of runner this binds the partially-initialized module
# object, which is fine because attributes are only accessed at call time.
import runner  # noqa: E402

logger = logging.getLogger(__name__)

def consensus_rank(vina_affinity: float, smina_affinity: Optional[float]) -> Dict:
    """Calculate a per-ligand consensus score for quick UI coloring/sorting.

    WARNING: the ``agreement`` field here is a *heuristic UI indicator only*,
    ``max(0, 100 - |vina - smina| * 10)`` (1 point per 0.1 kcal/mol of
    disagreement, floored at 10 kcal/mol). It is an arbitrary rule of thumb
    with no statistical basis and MUST NOT be cited as a quantitative result.
    For a defensible, citable measure of scorer concordance use
    :func:`scorer_agreement_spearman` (global Spearman rank correlation).
    """
    if smina_affinity is None:
        return {'vina': vina_affinity, 'smina': None, 'consensus': vina_affinity, 'agreement': 0}

    # Calculate agreement (0-100, higher = better agreement) -- HEURISTIC ONLY
    diff = abs(vina_affinity - smina_affinity)
    agreement = max(0, 100 - (diff * 10))  # Rule of thumb, not a real metric

    # Average the scores
    consensus = (vina_affinity + smina_affinity) / 2

    return {
        'vina': vina_affinity,
        'smina': smina_affinity,
        'consensus': consensus,
        'agreement': agreement
    }


def _rankdata(values: List[float]) -> List[float]:
    """Assign 1-based average ranks, handling ties (e.g. [10, 20, 20, 30]
    -> [1.0, 2.5, 2.5, 4.0]). Used for the Spearman correlation."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    n = len(values)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def _pearson(a: List[float], b: List[float]) -> float:
    n = len(a)
    if n < 2:
        return 0.0
    ma = sum(a) / n
    mb = sum(b) / n
    num = sum((a[i] - ma) * (b[i] - mb) for i in range(n))
    da = sum((a[i] - ma) ** 2 for i in range(n))
    db = sum((b[i] - mb) ** 2 for i in range(n))
    if da == 0 or db == 0:
        return 0.0
    return num / (da ** 0.5 * db ** 0.5)


def scorer_agreement_spearman(vina_scores: List[Optional[float]],
                              smina_scores: List[Optional[float]]):
    """Global, citable scorer-concordance metric (Spearman rank correlation).

    Computes the Spearman rank correlation (rho) between two scoring functions
    across all paired ligands in a screen, together with the number of paired
    points ``n``. This is the statistically grounded measure of agreement and
    is safe to report as a quantitative result.

    Returns ``(rho, n)`` or ``(None, n)`` when undefined -- fewer than 2
    paired values, or either scorer has zero variance (constant ranks).

    Pure-Python (no scipy/numpy dependency required).
    """
    pairs = [(float(v), float(s)) for v, s in zip(vina_scores, smina_scores)
             if v is not None and s is not None]
    n = len(pairs)
    if n < 2:
        return None, n
    rv = _rankdata([p[0] for p in pairs])
    rs = _rankdata([p[1] for p in pairs])
    # Constant-rank guard: zero variance -> correlation undefined.
    if all(r == rv[0] for r in rv) or all(r == rs[0] for r in rs):
        return None, n
    return _pearson(rv, rs), n


def run_smina_docking(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    output_pdbqt: str,
    cx: float, cy: float, cz: float,
    sx: float, sy: float, sz: float,
    exhaustiveness: int = 8,
    seed: Optional[int] = None,
    cpu: Optional[int] = None,
    flex_file: Optional[str] = None,
) -> Optional[float]:
    """Run SMINA docking and return best affinity."""

    if runner.SMINA is None:
        logger.debug("SMINA not available")
        return None

    cmd = [
        runner.SMINA,
        "-r", receptor_pdbqt,
        "-l", ligand_pdbqt,
        "-o", output_pdbqt,
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", str(sx),
        "--size_y", str(sy),
        "--size_z", str(sz),
        "--num_modes", "1",
        "--exhaustiveness", str(exhaustiveness),
    ]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    if cpu:
        cmd += ["--cpu", str(cpu)]
    if flex_file:
        cmd += ["--flex", flex_file]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=runner.COMMAND_TIMEOUT
    )

    # 🔴 Fail loudly
    if result.returncode != 0:
        raise RuntimeError(
            f"SMINA docking failed:\nSTDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
        )

    # 🔴 Validate output file
    if not os.path.exists(output_pdbqt):
        raise RuntimeError("SMINA did not produce output file")

    with open(output_pdbqt) as f:
        if not any(line.startswith("MODEL") for line in f):
            raise RuntimeError("SMINA produced no docking poses")

    # 🔴 Parse properly (mode table)
    aff = _parse_smina_affinity(result.stdout)
    if aff is None:
        raise RuntimeError("Could not extract affinity from SMINA output")
    return aff


def _parse_smina_affinity(stdout: str) -> Optional[float]:
    """Extract the best-mode affinity (kcal/mol) from SMINA stdout.

    SMINA prints the Vina-style table: a ``mode | affinity`` header, a units
    row, a ``-----+-----`` separator, then whitespace-separated rows
    (``1  -7.5  0.000  0.000``); the first row is the best mode. Output
    without that table falls back to the first digit-leading line.
    """
    lines = runner._mode_table_lines(stdout)
    for line in lines:
        toks = line.replace("|", " ").split()
        if len(toks) >= 2 and toks[0].isdigit():
            try:
                return float(toks[1])
            except ValueError:
                continue
    return None


def _smina_score_one(task):
    """Worker: dock a single ligand with SMINA. Picklable for Pool."""
    (receptor, lig, out, cx, cy, cz, sx, sy, sz, exhaustiveness, seed,
     cpu, flex_file, min_valid) = task
    name = runner.ligand_name(lig)
    try:
        score = run_smina_docking(
            receptor, lig, out, cx, cy, cz, sx, sy, sz,
            exhaustiveness=exhaustiveness, seed=seed, cpu=cpu, flex_file=flex_file)
    except Exception as e:
        return name, None, str(e)
    # Same validity window as the Vina scores they are compared with.
    if score is not None and (score >= min_valid or score < -20):
        return name, None, (f"score {score:.2f} kcal/mol outside the accepted "
                            f"window (-20 to <{min_valid})")
    return name, score, None


def _run_smina_scoring(receptor: str, ligands: List[str], outdir: str,
                       grid: Tuple[float, float, float, float, float, float],
                       vina_params: Dict,
                       num_processes: int = 1,
                       flex_file: Optional[str] = None) -> Dict[str, Optional[float]]:
    """Run SMINA docking for every ligand.

    Returns {name: smina_affinity} with None for ligands SMINA could not
    score, so callers can report them instead of dropping them. Per-ligand
    failures are logged, never fatal. Docking is parallelized with a process
    pool (like Vina) when ``num_processes > 1``, with the cores split between
    workers.
    """
    cx, cy, cz, sx, sy, sz = grid
    vina_params = vina_params or {}
    exhaustiveness = int(vina_params.get('exhaustiveness', 8))
    seed = vina_params.get('seed')
    min_valid = float(vina_params.get('min_valid_affinity', -1.0))
    cpu = vina_params.get('threads')
    if not cpu and num_processes > 1:
        cpu = max(1, runner.available_cpus() // num_processes)
    smina_scores = {}
    dock_dir = os.path.join(outdir, "docked")
    os.makedirs(dock_dir, exist_ok=True)

    tasks = [
        (receptor, lig, os.path.join(dock_dir, runner.ligand_name(lig) + "_smina.pdbqt"),
         cx, cy, cz, sx, sy, sz, exhaustiveness, seed, cpu, flex_file, min_valid)
        for lig in ligands
    ]

    def _collect(results):
        for name, score, err in results:
            if err:
                logger.warning(f"  [!] SMINA failed for {name}: {err}")
            smina_scores[name] = score if not err else None

    if num_processes > 1:
        try:
            with Pool(num_processes) as pool:
                _collect(pool.imap_unordered(_smina_score_one, tasks))
        except Exception as e:
            logger.warning(
                f"  [!] SMINA parallel scoring failed ({e}); retrying serially")
            _collect(_smina_score_one(t) for t in tasks)
    else:
        _collect(_smina_score_one(t) for t in tasks)

    logger.info(
        f"[*] SMINA scored {sum(v is not None for v in smina_scores.values())}"
        f"/{len(ligands)} ligands")
    return smina_scores
