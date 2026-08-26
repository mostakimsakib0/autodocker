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
    """Calculate consensus score."""
    if smina_affinity is None:
        return {'vina': vina_affinity, 'smina': None, 'consensus': vina_affinity, 'agreement': 0}

    # Calculate agreement (0-100, higher = better agreement)
    diff = abs(vina_affinity - smina_affinity)
    agreement = max(0, 100 - (diff * 10))  # Rule of thumb

    # Average the scores
    consensus = (vina_affinity + smina_affinity) / 2

    return {
        'vina': vina_affinity,
        'smina': smina_affinity,
        'consensus': consensus,
        'agreement': agreement
    }


def run_smina_docking(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    output_pdbqt: str,
    cx: float, cy: float, cz: float,
    sx: float, sy: float, sz: float
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
        "--exhaustiveness", str(os.environ.get('VS_EXHAUSTIVENESS', 8))
    ]

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
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                return float(parts[1])
            except ValueError:
                continue

    raise RuntimeError("Could not extract affinity from SMINA output")


def _run_smina_scoring(receptor: str, ligands: List[str], outdir: str,
                       grid: Tuple[float, float, float, float, float, float],
                       vina_params: Dict) -> Dict[str, float]:
    """Run SMINA scoring for every ligand. Returns {name: smina_affinity}.

    Per-ligand failures are logged and skipped (never fatal).
    """
    cx, cy, cz, sx, sy, sz = grid
    smina_scores = {}
    dock_dir = os.path.join(outdir, "docked")
    os.makedirs(dock_dir, exist_ok=True)

    for lig in ligands:
        name = os.path.basename(lig).replace(".pdbqt", "")
        out = os.path.join(dock_dir, name + "_smina.pdbqt")
        try:
            score = run_smina_docking(
                receptor, lig, out, cx, cy, cz, sx, sy, sz)
        except Exception as e:
            logger.warning(f"  [!] SMINA failed for {name}: {e}")
            continue
        if score is not None:
            smina_scores[name] = score

    logger.info(
        f"[*] SMINA scored {len(smina_scores)}/{len(ligands)} ligands")
    return smina_scores
