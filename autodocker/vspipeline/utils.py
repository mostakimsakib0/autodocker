#!/usr/bin/env python3
"""UTILITIES.

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

def run(cmd_list: List[str], capture: bool = False) -> Optional[Tuple[str, str, int]]:
    """Run a command safely without shell injection risk."""
    try:
        if capture:
            result = subprocess.run(
                cmd_list, capture_output=True, text=True, timeout=runner.COMMAND_TIMEOUT)
            if result.returncode != 0:
                raise subprocess.CalledProcessError(
                    result.returncode, cmd_list, output=result.stdout, stderr=result.stderr
                )
            return result.stdout, result.stderr, result.returncode
        else:
            subprocess.run(cmd_list, check=True, timeout=runner.COMMAND_TIMEOUT)
    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out: {' '.join(cmd_list)}")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(
            f"Command failed with return code {e.returncode}: {' '.join(cmd_list)}")
        raise
    except FileNotFoundError as e:
        logger.error(f"Command not found: {' '.join(cmd_list)}")
        raise


def parse_pdb_coords(pdb_path: str) -> Tuple[List[float], List[float], List[float]]:
    """Extract XYZ coordinates from PDB file."""
    xs, ys, zs = [], [], []
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
                    except (ValueError, IndexError):
                        logger.warning(
                            f"Skipping malformed PDB line: {line.strip()}")
                        continue
    except IOError as e:
        logger.error(f"Cannot read PDB file {pdb_path}: {e}")
        raise

    if not xs:
        raise ValueError(f"No ATOM records found in {pdb_path}")

    return xs, ys, zs


def _safe_ligand_id(name: str) -> str:
    """Sanitize a compound name into a safe file basename token."""
    token = re.sub(r"[^A-Za-z0-9_.-]", "_", name.strip().lower())
    return token or "compound"


def get_chains(pdb_file: str) -> List[str]:
    """Extract unique chains from PDB file."""
    chains = set()
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain = line[21]
                    if chain.strip():
                        chains.add(chain)
    except Exception as e:
        logger.warning(f"Could not parse chains: {e}")
        return ["A"]

    return sorted(list(chains)) if chains else ["A"]


def _parse_chain_selection(chain: str, available: List[str]) -> List[str]:
    choice = (chain or "").strip()
    if not choice or choice.lower() in {"all", "*"}:
        return available

    selected = [c.strip().upper()
                for c in re.split(r"[,\s]+", choice) if c.strip()]
    invalid = [c for c in selected if c not in available]
    if invalid:
        raise ValueError(
            f"Invalid chain(s): {', '.join(invalid)}. Available chains: {', '.join(available)}"
        )
    return selected
