#!/usr/bin/env python3
"""fpocket pockets: running fpocket, reading its pockets and scores, ranking
and selecting them, and pocket_summary.csv.

fpocket is found with ``tools.find_tool`` and run through ``utils.run``,
both called through their module so tests can replace them.
"""
import csv
import logging
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from vspipeline import tools, utils

logger = logging.getLogger(__name__)

# fpocket's documented druggability threshold: "A druggability score at 0.5
# (the threshold) indicates that binding of prodrugs or druglike molecules can
# be possible" (fpocket doc/GETTINGSTARTED.md, "Druggability Score").
FPOCKET_DRUGGABILITY_THRESHOLD = 0.5


def _run_tool(cmd: List[str], tool: str, log_stderr: bool = True) -> str:
    """Run an external tool through ``utils.run``; its stderr goes to the
    log (unless ``log_stderr`` is False). Returns stdout. Raises RuntimeError
    with the tool's stderr on failure."""
    try:
        result = utils.run(cmd, capture=True)
    except subprocess.CalledProcessError as e:
        stderr = (e.stderr or "").strip()
        raise RuntimeError(
            f"{tool} failed (exit code {e.returncode}): {' '.join(cmd)}"
            + (f"\n{stderr[-2000:]}" if stderr else "")) from e
    stdout, stderr = (result or ("", ""))[:2]
    if log_stderr:
        utils.log_tool_stderr(tool, stderr)
    return stdout or ""


def _fpocket_version(fpocket: str) -> Optional[str]:
    """fpocket's version from its usage text (it has no --version flag)."""
    try:
        proc = subprocess.run([fpocket, "-h"], capture_output=True, text=True,
                              timeout=tools.settings().command_timeout)
    except (OSError, subprocess.SubprocessError):
        return None
    match = re.search(r"fpocket\s+(\d[\w.]*)", re.sub(r"\x1b\[[0-9;]*m", "",
                                                     proc.stdout + proc.stderr))
    return match.group(1) if match else None


def is_druggable(pocket: Dict) -> bool:
    """fpocket druggability score at or above the documented 0.5 threshold."""
    return pocket.get("druggability_score", 0) >= FPOCKET_DRUGGABILITY_THRESHOLD


def rank_pockets(pockets: List[Dict]) -> List[Dict]:
    """Pockets best first: fpocket score, then druggability score, then the
    lower pocket number. fpocket's volume is a Monte Carlo estimate seeded
    from the clock, so it is not used."""
    return sorted(pockets, key=lambda p: (-p.get("score", -1),
                                          -p.get("druggability_score", -1),
                                          p["number"]))


def run_fpocket(target: str) -> str:
    """Run fpocket on ``target``; returns its output directory. Raises
    when fpocket is missing or fails."""
    fpocket = tools.find_tool("fpocket")
    if fpocket is None:
        raise RuntimeError(
            "fpocket not found on PATH; install fpocket, or pass --no-fpocket "
            "to dock in a box around the whole receptor, or give the box with "
            "--grid-center/--grid-size")
    version = _fpocket_version(fpocket)
    logger.info(f"[*] fpocket version: {version or '(unknown)'} ({fpocket})")
    pocket_root = Path(target).with_suffix("").as_posix() + "_out"
    shutil.rmtree(pocket_root, ignore_errors=True)
    _run_tool([fpocket, "-f", target], "fpocket")
    return pocket_root


def parse_fpocket_info(info_file: str) -> Dict[int, Dict]:
    """Parse fpocket's *_info.txt file into a pocket metadata map.

    Raises if the file cannot be read.
    """
    info = {}
    current = None
    try:
        with open(info_file, encoding="utf-8") as f:
            for line in f:
                pocket_match = re.match(
                    r"Pocket\s+(\d+)\s*:", line.strip())
                if pocket_match:
                    current = int(pocket_match.group(1))
                    info[current] = {}
                    continue

                if current is None or ":" not in line:
                    continue

                key, value = line.split(":", 1)
                key = key.strip().lower().replace(" ", "_").replace(".", "")
                key = key.replace("-", "_")
                try:
                    info[current][key] = float(value.strip().split()[0])
                except (ValueError, IndexError):
                    continue
    except IOError as e:
        raise RuntimeError(f"Could not read fpocket info file {info_file}: {e}") from e

    return info


def collect_pockets(pocket_root: str, target: str) -> List[Dict]:
    """The pockets fpocket wrote for ``target`` in ``pocket_root``: its
    scores (``*_info.txt``) plus number, file, path and the volume of the
    pocket atoms' bounding box. Raises when there is none."""
    pocket_dir = os.path.join(pocket_root, "pockets")
    if not os.path.isdir(pocket_dir):
        raise RuntimeError(f"fpocket wrote no pocket directory ({pocket_dir})")
    info_file = os.path.join(pocket_root, f"{Path(target).stem}_info.txt")
    if not os.path.exists(info_file):
        raise RuntimeError(f"fpocket wrote no pocket scores ({info_file})")
    pocket_info = parse_fpocket_info(info_file)

    pockets = []
    for f in sorted(os.listdir(pocket_dir)):
        match = re.match(r"pocket(\d+)_atm\.pdb$", f)
        if not match:
            continue
        number = int(match.group(1))
        if number not in pocket_info:
            raise RuntimeError(f"Pocket {number} has no entry in {info_file}")
        p = os.path.join(pocket_dir, f)
        xs, ys, zs = utils.parse_pdb_coords(p)
        data = dict(pocket_info[number])
        data.update({
            "number": number,
            "file": f,
            "path": p,
            "box_volume": (max(xs) - min(xs)) * (max(ys) - min(ys)) * (max(zs) - min(zs)),
        })
        pockets.append(data)
    if not pockets:
        raise RuntimeError(
            f"fpocket found no pockets in {target}; pass --grid-center and "
            "--grid-size, or --no-fpocket for a box around the whole receptor")
    return pockets


def write_pocket_summary(path: str, pockets: List[Dict]) -> None:
    """pocket_summary.csv: every pocket by number (a write error is only
    logged)."""
    try:
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow([
                "Pocket", "Score", "Druggability_Score", "Volume",
                "Box_Volume", "Could_Be_Druggable", "File"
            ])
            for pocket in sorted(pockets, key=lambda p: p["number"]):
                writer.writerow([
                    pocket["number"],
                    pocket.get("score", ""),
                    pocket.get("druggability_score", ""),
                    pocket.get("volume", ""),
                    f"{pocket.get('box_volume', 0):.3f}",
                    "YES" if is_druggable(pocket) else "NO",
                    pocket["file"],
                ])
        logger.info(
            f"[*] Pocket summary saved: {path}")
    except IOError as e:
        logger.warning(f"Could not save pocket summary: {e}")


def log_top_pockets(pockets: List[Dict]) -> None:
    """Log the five best-ranked pockets."""
    logger.info(f"[*] Found {len(pockets)} pockets. Top candidates:")
    for pocket in pockets[:5]:
        logger.info(
            "  pocket%-2s score=%6.3f druggability=%5.3f volume=%7.1f %s"
            % (
                pocket["number"],
                pocket.get("score", 0),
                pocket.get("druggability_score", 0),
                pocket.get("volume", 0),
                "druggable" if is_druggable(pocket) else "low-druggability",
            )
        )


def parse_pocket_selection(pocket_spec: Optional[str], pockets: List[Dict]) -> List[int]:
    """The pocket numbers of --pockets (``pockets`` ranked best first):
    the best pocket for None/"auto", else the listed numbers in order."""
    available = {p["number"] for p in pockets}
    if not pocket_spec or pocket_spec.lower() == "auto":
        return [pockets[0]["number"]]

    selected = []
    for token in re.split(r"[,\s]+", pocket_spec.strip()):
        if not token:
            continue
        try:
            pocket_number = int(token)
        except ValueError as e:
            raise ValueError(
                f"Invalid pocket selection '{token}'. Use numbers like 1 or 1,3.") from e
        if pocket_number not in available:
            raise ValueError(
                f"Pocket {pocket_number} was not found. Available: {sorted(available)}")
        selected.append(pocket_number)

    if not selected:
        raise ValueError("No pockets selected")
    return selected


def select_pockets(pocket_spec: Optional[str], pockets: List[Dict]) -> List[Dict]:
    """The pockets chosen by --pockets, in the order given."""
    selected_numbers = parse_pocket_selection(pocket_spec, pockets)
    selected = [p for p in pockets if p["number"] in selected_numbers]
    selected.sort(key=lambda p: selected_numbers.index(p["number"]))
    return selected
