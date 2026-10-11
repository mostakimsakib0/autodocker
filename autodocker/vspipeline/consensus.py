#!/usr/bin/env python3
"""SMINA docking (--smina-only) and consensus rescoring (--consensus).

Consensus design: SMINA does *not* re-dock. It rescores the top pose Vina
produced for each ligand (``smina --score_only``) with a different scoring
function, Vinardo (``--scoring vinardo``; ``dkoes_scoring`` when the SMINA
build lacks Vinardo). Re-docking with SMINA's default scoring would apply the
Vina function a second time, so the Vina/SMINA rank correlation would only
measure search noise. Rescoring one fixed pose with a second function
measures agreement between scoring functions.

The global agreement statistic (Spearman rho, n) is written to
``consensus_summary.json`` next to ``consensus_ranking.csv``.

SMINA paths and the timeout are passed to workers explicitly, so they
survive the ``spawn`` start method. Every SMINA command is recorded in
``docked/<name>_smina.log``.
"""
from __future__ import annotations

import csv
import json
import logging
import os
import re
import shlex
import subprocess
import traceback
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from vspipeline.docking import (
    runner,  # lazily imported runner module (see docking._LazyRunner)
    DEFAULT_MIN_VALID_AFFINITY,
    _mode_table_lines,
    _scores_agree,
    _stdout_best_affinity_token,
    _parse_pdbqt_result_tokens,
    affinity_in_window,
    window_text,
)

logger = logging.getLogger(__name__)

# Scoring functions tried, in order, for consensus rescoring.
SMINA_RESCORING_FUNCTIONS = ("vinardo", "dkoes_scoring")

_SMINA_SCORING_CACHE: Dict[str, Tuple[str, ...]] = {}
_SCORE_ONLY_RE = re.compile(r"^Affinity:\s+(\S+)", re.MULTILINE)


def consensus_rank(vina_affinity: float, smina_affinity: Optional[float]) -> Dict:
    """Per-ligand consensus row: the mean of the Vina and SMINA scores.

    The ``agreement`` field is deprecated (kept only because the current
    runner.py still writes it); it is an arbitrary rule of thumb,
    ``max(0, 100 - 10 * |vina - smina|)``, and must not be reported. Use
    :func:`scorer_agreement_spearman` for scorer concordance.
    """
    if smina_affinity is None:
        return {'vina': vina_affinity, 'smina': None, 'consensus': vina_affinity, 'agreement': 0}
    diff = abs(vina_affinity - smina_affinity)
    return {
        'vina': vina_affinity,
        'smina': smina_affinity,
        'consensus': (vina_affinity + smina_affinity) / 2,
        'agreement': max(0, 100 - (diff * 10)),
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
    """Spearman rank correlation between two scoring functions over all
    paired ligands of a screen.

    Returns ``(rho, n)``, or ``(None, n)`` when undefined: fewer than 2
    paired values, or either scorer has zero variance (constant ranks).
    Pure Python (no scipy/numpy).
    """
    pairs = [(float(v), float(s)) for v, s in zip(vina_scores, smina_scores)
             if v is not None and s is not None]
    n = len(pairs)
    if n < 2:
        return None, n
    rv = _rankdata([p[0] for p in pairs])
    rs = _rankdata([p[1] for p in pairs])
    if all(r == rv[0] for r in rv) or all(r == rs[0] for r in rs):
        return None, n
    return _pearson(rv, rs), n


# ---------------------------------------------------------------------------
# Running SMINA
# ---------------------------------------------------------------------------

def _write_log(log_file: Optional[str], text: str, mode: str = "a") -> None:
    if not log_file:
        return
    with open(log_file, mode) as f:
        f.write(text)


def _run_smina_command(cmd: List[str], timeout: int, log_file: Optional[str]
                       ) -> subprocess.CompletedProcess:
    """Run SMINA, recording the command first and its output afterwards.
    Raises RuntimeError on a nonzero exit status."""
    _write_log(log_file, f"# Command: {shlex.join(cmd)}\n", mode="w")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    _write_log(log_file, (result.stdout or "") + "\n")
    if result.stderr and result.stderr.strip():
        _write_log(log_file, "# STDERR:\n" + result.stderr + "\n")
    if result.returncode != 0:
        _write_log(log_file, f"# Error: exit status {result.returncode}\n")
        raise RuntimeError(
            f"SMINA failed (exit {result.returncode}):\nSTDERR:\n{result.stderr}\n"
            f"STDOUT:\n{result.stdout}")
    return result


def _first_pdbqt_affinity(path: str) -> str:
    """Affinity token of the first model of a SMINA output PDBQT
    (``REMARK minimizedAffinity``). Raises RuntimeError when absent."""
    if not os.path.exists(path):
        raise RuntimeError("SMINA did not produce output file")
    with open(path) as f:
        if not any(line.startswith("MODEL") for line in f):
            raise RuntimeError("SMINA produced no docking poses")
    tokens = _parse_pdbqt_result_tokens(path)
    if not tokens:
        raise RuntimeError(f"Could not extract affinity from SMINA output {path}")
    return tokens[0][0]


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
    smina: Optional[str] = None,
    timeout: Optional[int] = None,
    log_file: Optional[str] = None,
) -> Optional[float]:
    """Dock one ligand with SMINA (default scoring) and return the best
    affinity, read from ``REMARK minimizedAffinity`` of the output PDBQT
    (the stdout table, rounded to one decimal, is only a cross-check).

    ``smina``/``timeout`` default to the runner globals. Returns None when no
    SMINA binary is available; raises RuntimeError on any failure.
    """
    smina = smina or runner.SMINA
    if smina is None:
        logger.debug("SMINA not available")
        return None
    timeout = int(timeout or runner.COMMAND_TIMEOUT)

    cmd = [
        smina,
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

    if os.path.exists(output_pdbqt):
        os.remove(output_pdbqt)
    result = _run_smina_command(cmd, timeout, log_file)
    token = _first_pdbqt_affinity(output_pdbqt)
    stdout_token = _stdout_best_affinity_token(result.stdout or "")
    if stdout_token is not None and not _scores_agree(token, stdout_token):
        logger.warning(
            f"[!] SMINA output PDBQT score {token} disagrees with stdout "
            f"({stdout_token}) for {ligand_pdbqt}; using the PDBQT value")
    return float(token)


def _parse_smina_affinity(stdout: str) -> Optional[float]:
    """Best-mode affinity from SMINA stdout (rounded to print precision).

    SMINA prints the Vina-style table: a ``mode | affinity`` header, a units
    row, a ``-----+-----`` separator, then whitespace-separated rows; the
    first row is the best mode. Output without that table falls back to the
    first digit-leading line.
    """
    for line in _mode_table_lines(stdout):
        toks = line.replace("|", " ").split()
        if len(toks) >= 2 and toks[0].isdigit():
            try:
                return float(toks[1])
            except ValueError:
                continue
    return None


# ---------------------------------------------------------------------------
# Worker plumbing
# ---------------------------------------------------------------------------

def _collect_with_retry(worker: Callable, jobs: Sequence, num_processes: int,
                        label: str) -> Dict[str, Tuple[Optional[float], Optional[str]]]:
    """Run ``worker`` over ``jobs`` (in a pool when num_processes > 1).
    If the pool itself breaks, only the jobs without a result are retried
    serially. Returns {name: (score, error)}."""
    collected: Dict[str, Tuple[Optional[float], Optional[str]]] = {}
    if num_processes > 1 and len(jobs) > 1:
        try:
            with Pool(num_processes) as pool:
                for name, score, err in pool.imap_unordered(worker, jobs):
                    collected[name] = (score, err)
        except Exception as e:
            remaining = [j for j in jobs if j.name not in collected]
            logger.warning(
                f"  [!] {label}: parallel run failed ({e}); retrying the "
                f"{len(remaining)} unfinished ligand(s) serially")
            for job in remaining:
                name, score, err = worker(job)
                collected[name] = (score, err)
    else:
        for job in jobs:
            name, score, err = worker(job)
            collected[name] = (score, err)
    return collected


@dataclass(frozen=True)
class SminaDockJob:
    name: str
    receptor: str
    ligand: str
    out: str
    log_file: str
    grid: Tuple[float, float, float, float, float, float]
    exhaustiveness: int
    seed: Optional[int]
    cpu: Optional[int]
    flex_file: Optional[str]
    min_valid: float
    smina: str
    timeout: int


def _smina_score_one(job: SminaDockJob):
    """Worker: dock a single ligand with SMINA. Picklable for Pool."""
    cx, cy, cz, sx, sy, sz = job.grid
    try:
        score = run_smina_docking(
            job.receptor, job.ligand, job.out, cx, cy, cz, sx, sy, sz,
            exhaustiveness=job.exhaustiveness, seed=job.seed, cpu=job.cpu,
            flex_file=job.flex_file, smina=job.smina, timeout=job.timeout,
            log_file=job.log_file)
    except Exception as e:
        try:
            _write_log(job.log_file, f"# Error: {e}\n# Traceback:\n{traceback.format_exc()}\n")
        except OSError:
            pass
        return job.name, None, str(e)
    # Same validity window as Vina docking scores.
    if score is not None and not affinity_in_window(score, job.min_valid):
        return job.name, None, (f"score {score:.3f} kcal/mol outside the accepted "
                                f"window ({window_text(job.min_valid)})")
    return job.name, score, None


def _run_smina_scoring(receptor: str, ligands: List[str], outdir: str,
                       grid: Tuple[float, float, float, float, float, float],
                       vina_params: Dict,
                       num_processes: int = 1,
                       flex_file: Optional[str] = None) -> Dict[str, Optional[float]]:
    """Dock every ligand with SMINA (``--smina-only`` mode).

    Returns {name: smina_affinity} with None for ligands SMINA could not
    score or that scored outside the accepted window. Per-ligand failures are
    logged, never fatal. Uses a process pool when ``num_processes > 1``, with
    the cores split between workers. For consensus scoring use
    :func:`rescore_with_smina` instead.
    """
    vina_params = vina_params or {}
    exhaustiveness = int(vina_params.get('exhaustiveness', 8))
    seed = vina_params.get('seed')
    min_valid = float(vina_params.get('min_valid_affinity', DEFAULT_MIN_VALID_AFFINITY))
    cpu = vina_params.get('threads')
    if not cpu and num_processes > 1:
        cpu = max(1, runner.available_cpus() // num_processes)
    dock_dir = os.path.join(outdir, "docked")
    os.makedirs(dock_dir, exist_ok=True)
    smina = runner.SMINA
    timeout = int(runner.COMMAND_TIMEOUT)

    jobs = []
    for lig in ligands:
        name = runner.ligand_name(lig)
        jobs.append(SminaDockJob(
            name=name, receptor=receptor, ligand=lig,
            out=os.path.join(dock_dir, name + "_smina.pdbqt"),
            log_file=os.path.join(dock_dir, name + "_smina.log"),
            grid=tuple(grid), exhaustiveness=exhaustiveness, seed=seed, cpu=cpu,
            flex_file=flex_file, min_valid=min_valid, smina=smina, timeout=timeout))

    collected = _collect_with_retry(_smina_score_one, jobs, num_processes, "SMINA docking")
    smina_scores: Dict[str, Optional[float]] = {}
    for job in jobs:
        score, err = collected.get(job.name, (None, "no result"))
        if err:
            logger.warning(f"  [!] SMINA failed for {job.name}: {err}")
        smina_scores[job.name] = None if err else score

    logger.info(
        f"[*] SMINA scored {sum(v is not None for v in smina_scores.values())}"
        f"/{len(ligands)} ligands")
    return smina_scores


# ---------------------------------------------------------------------------
# Consensus: rescoring Vina's top pose
# ---------------------------------------------------------------------------

def _smina_scoring_functions(smina: str, receptor: str, ligand: str,
                             timeout: int) -> Tuple[str, ...]:
    """Builtin scoring functions of this SMINA build (cached per binary).

    SMINA lists them when given an invalid ``--scoring`` name. Returns ()
    when the list cannot be read.
    """
    if smina in _SMINA_SCORING_CACHE:
        return _SMINA_SCORING_CACHE[smina]
    names: Tuple[str, ...] = ()
    try:
        result = subprocess.run(
            [smina, "-r", receptor, "-l", ligand, "--score_only",
             "--scoring", "__autodocker_list__", "-q"],
            capture_output=True, text=True, timeout=timeout)
        text = (result.stdout or "") + "\n" + (result.stderr or "")
        if "Options are:" in text:
            tail = text.split("Options are:", 1)[1]
            names = tuple(tok for tok in tail.split() if re.fullmatch(r"[A-Za-z0-9_]+", tok))
    except Exception as e:
        logger.debug(f"Could not list SMINA scoring functions: {e}")
    _SMINA_SCORING_CACHE[smina] = names
    return names


def select_rescoring_function(smina: str, receptor: str, ligand: str,
                              timeout: int) -> str:
    """First of SMINA_RESCORING_FUNCTIONS that this SMINA build offers
    (vinardo when the list cannot be read)."""
    available = _smina_scoring_functions(smina, receptor, ligand, timeout)
    for name in SMINA_RESCORING_FUNCTIONS:
        if not available or name in available:
            return name
    raise RuntimeError(
        f"SMINA {smina} offers none of the rescoring functions "
        f"{', '.join(SMINA_RESCORING_FUNCTIONS)} (available: {', '.join(available)})")


def _extract_top_pose(dock_file: str, ligand_out: str, flex_out: str) -> bool:
    """Write the first MODEL of a Vina output PDBQT to ``ligand_out``; its
    flexible side chains (BEGIN_RES..END_RES blocks of flexible docking), if
    any, go to ``flex_out``. Returns True when flex residues were written."""
    ligand_lines: List[str] = []
    flex_lines: List[str] = []
    in_model = in_res = False
    with open(dock_file) as f:
        for line in f:
            if line.startswith("MODEL"):
                if in_model:
                    break
                in_model = True
                continue
            if line.startswith("ENDMDL"):
                break
            if not in_model:
                continue
            if line.startswith("BEGIN_RES"):
                in_res = True
            (flex_lines if in_res else ligand_lines).append(line)
            if line.startswith("END_RES"):
                in_res = False
    if not ligand_lines:
        raise RuntimeError(f"No pose in {dock_file}")
    with open(ligand_out, "w") as f:
        f.writelines(ligand_lines)
    if flex_lines:
        with open(flex_out, "w") as f:
            f.writelines(flex_lines)
        return True
    return False


@dataclass(frozen=True)
class SminaRescoreJob:
    name: str
    receptor: str
    dock_file: str
    dock_dir: str
    flex_file: Optional[str]
    scoring: str
    smina: str
    timeout: int


def rescoring_command(smina: str, receptor: str, pose: str, out: str, scoring: str,
                      flex: Optional[str] = None) -> List[str]:
    """``smina --score_only`` command that rescores one fixed pose."""
    cmd = [smina, "-r", receptor, "-l", pose, "--score_only",
           "--scoring", scoring, "-o", out]
    if flex:
        cmd += ["--flex", flex]
    return cmd


def _smina_rescore_one(job: SminaRescoreJob):
    """Worker: rescore Vina's top pose of one ligand. Picklable for Pool."""
    log_file = os.path.join(job.dock_dir, job.name + "_smina.log")
    pose = os.path.join(job.dock_dir, job.name + "_vina_top.pdbqt")
    pose_flex = os.path.join(job.dock_dir, job.name + "_vina_top_flex.pdbqt")
    out = os.path.join(job.dock_dir, job.name + "_smina_rescore.pdbqt")
    try:
        if os.path.exists(out):
            os.remove(out)
        has_flex = _extract_top_pose(job.dock_file, pose, pose_flex)
        if job.flex_file and not has_flex:
            raise RuntimeError(
                f"flexible docking output {job.dock_file} has no flexible residues")
        cmd = rescoring_command(job.smina, job.receptor, pose, out, job.scoring,
                                pose_flex if has_flex else None)
        result = _run_smina_command(cmd, job.timeout, log_file)
        token = _first_pdbqt_affinity(out)
        m = _SCORE_ONLY_RE.search(result.stdout or "")
        if m and not _scores_agree(token, m.group(1)):
            logger.warning(
                f"[!] {job.name}: SMINA rescoring PDBQT score {token} disagrees with "
                f"stdout ({m.group(1)}); using the PDBQT value")
        return job.name, float(token), None
    except Exception as e:
        try:
            _write_log(log_file, f"# Error: {e}\n# Traceback:\n{traceback.format_exc()}\n")
        except OSError:
            pass
        return job.name, None, str(e)


def rescore_with_smina(receptor: str, results: List[Tuple[str, Optional[float]]],
                       outdir: str, metrics_dict: Optional[Dict] = None,
                       num_processes: int = 1, flex_file: Optional[str] = None,
                       scoring: Optional[str] = None, smina: Optional[str] = None,
                       timeout: Optional[int] = None
                       ) -> Tuple[Dict[str, Optional[float]], str]:
    """Rescore Vina's top pose of every ligand with an in-window Vina score.

    ``smina --score_only --scoring <scoring>`` (default: the first available
    of SMINA_RESCORING_FUNCTIONS) is run on the first model of
    ``docked/<name>_out.pdbqt`` (or ``metrics_dict[name]["dock_file"]``);
    in flexible docking the docked side chains are passed with ``--flex``.
    No affinity window is applied to the rescored values (they come from a
    different function). Returns ``({name: score or None}, scoring)``.
    """
    smina = smina or runner.SMINA
    if smina is None:
        raise RuntimeError("SMINA binary not available for consensus rescoring")
    timeout = int(timeout or runner.COMMAND_TIMEOUT)
    metrics_dict = metrics_dict or {}
    dock_dir = os.path.join(outdir, "docked")
    os.makedirs(dock_dir, exist_ok=True)

    names = sorted(name for name, score in results if score is not None)
    dock_files = {n: metrics_dict.get(n, {}).get("dock_file")
                  or os.path.join(dock_dir, n + "_out.pdbqt") for n in names}
    scores: Dict[str, Optional[float]] = {name: None for name, _ in results}
    if not names:
        return scores, scoring or SMINA_RESCORING_FUNCTIONS[0]

    if scoring is None:
        probe_pose = os.path.join(dock_dir, names[0] + "_vina_top.pdbqt")
        _extract_top_pose(dock_files[names[0]], probe_pose,
                          os.path.join(dock_dir, names[0] + "_vina_top_flex.pdbqt"))
        scoring = select_rescoring_function(smina, receptor, probe_pose, timeout)
    logger.info(f"[*] Consensus: rescoring Vina top poses with smina --score_only "
                f"--scoring {scoring}")

    jobs = [SminaRescoreJob(n, receptor, dock_files[n], dock_dir, flex_file,
                            scoring, smina, timeout) for n in names]
    collected = _collect_with_retry(_smina_rescore_one, jobs, num_processes,
                                    "SMINA rescoring")
    for name in names:
        score, err = collected.get(name, (None, "no result"))
        if err:
            logger.warning(f"  [!] SMINA rescoring failed for {name}: {err}")
        scores[name] = None if err else score
    logger.info(f"[*] SMINA rescored {sum(v is not None for v in scores.values())}"
                f"/{len(names)} ligands")
    return scores, scoring


def write_consensus_outputs(outdir: str, results: List[Tuple[str, Optional[float]]],
                            smina_scores: Dict[str, Optional[float]],
                            scoring: str) -> Dict:
    """Write consensus_ranking.csv and consensus_summary.json.

    consensus_ranking.csv: Ligand, Vina_Affinity, SMINA_Affinity, Consensus
    (mean of the two), sorted by Consensus then name; ligands one scorer could
    not score follow with NA. consensus_summary.json holds the global
    Spearman rho and n. Returns ``{"rho": rho, "n": n}``.
    """
    paired, unpaired = [], []
    for name, vina_score in results:
        smina_score = smina_scores.get(name)
        if vina_score is None or smina_score is None:
            unpaired.append((name, vina_score, smina_score))
        else:
            paired.append((name, vina_score, smina_score, (vina_score + smina_score) / 2))

    def fmt(v):
        return "FAILED" if v is None else v

    consensus_file = os.path.join(outdir, "consensus_ranking.csv")
    with open(consensus_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Ligand", "Vina_Affinity", "SMINA_Affinity", "Consensus"])
        for name, v, s, c in sorted(paired, key=lambda r: (r[3], r[0])):
            writer.writerow([name, v, s, c])
        for name, v, s in sorted(unpaired, key=lambda r: r[0]):
            writer.writerow([name, fmt(v), fmt(s), "NA"])

    rho, n = scorer_agreement_spearman([r[1] for r in paired], [r[2] for r in paired])
    summary = {
        "metric": "spearman_rho",
        "rho": rho,
        "n": n,
        "n_unpaired": len(unpaired),
        "scorer_a": "AutoDock Vina docking score (top pose)",
        "scorer_b": f"smina --score_only --scoring {scoring} (Vina top pose)",
    }
    with open(os.path.join(outdir, "consensus_summary.json"), "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")
    logger.info(f"[*] Consensus ranking saved: {consensus_file} ({len(paired)} paired ligands; "
                f"Spearman rho = {'undefined' if rho is None else f'{rho:.3f}'}, n={n})")
    return {"rho": rho, "n": n}
