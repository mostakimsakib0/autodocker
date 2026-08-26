#!/usr/bin/env python3
"""PDBQT CHARGE VALIDATION & FIXING.

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

def _extract_pdbqt_charge(line: str) -> Optional[float]:
    """Return the atomic charge from a PDBQT ATOM/HETATM line, or None.

    PDBQT charge placement varies by OpenBabel/Vina build, so the first
    numeric token scanning backwards from the atom-type is used.
    """
    parts = line.split()
    for token in reversed(parts[-3:]):
        try:
            return float(token)
        except ValueError:
            continue
    return None


def _ensure_pdbqt_has_charges(pdbqt_file: str) -> bool:
    """Check if PDBQT file contains at least one atom with a NONZERO charge.

    A file whose atoms all carry 0.000 charge is chemically unusable for
    docking and is treated as 'no charges' (zero-charge receptors would
    otherwise pass validation and silently produce garbage results).
    """
    try:
        saw_atom = False
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    saw_atom = True
                    charge = _extract_pdbqt_charge(line)
                    if charge is not None and charge != 0.0:
                        return True
        if saw_atom:
            return False
        return False
    except Exception as e:
        logger.warning(f"Charge check failed: {e}")
        return False


def _pdbqt_has_atoms(pdbqt_file: str) -> bool:
    """Return True when a PDBQT contains at least one ATOM/HETATM record."""
    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    return True
        return False
    except Exception as e:
        logger.warning(f"Could not inspect atom records in {pdbqt_file}: {e}")
        return False


def _sanitize_receptor_pdbqt(pdbqt_file: str) -> None:
    """Remove ligand-only tags from receptor PDBQT for strict Vina compatibility."""
    allowed_prefixes = ("ATOM", "HETATM", "TER", "END", "REMARK")
    removed = 0
    kept = []

    try:
        with open(pdbqt_file, "r") as f:
            for line in f:
                stripped = line.lstrip()
                if stripped.startswith(("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):
                    removed += 1
                    continue
                if line.startswith(allowed_prefixes):
                    kept.append(line)
                    continue
                # Keep unrecognized lines as REMARK to preserve traceability.
                if line.strip():
                    kept.append("REMARK " + line)

        with open(pdbqt_file, "w") as f:
            f.writelines(kept)

        if removed:
            logger.info(
                f"[*] Sanitized receptor PDBQT: removed {removed} ligand-only tags")
    except Exception as e:
        raise RuntimeError(
            f"Failed to sanitize receptor PDBQT {pdbqt_file}: {e}")


def _fix_pdbqt_charges(pdbqt_file: str, source_file: str) -> None:
    """
    STRICT charge assignment.
    No silent guessing. Only real tools allowed.
    """

    logger.warning(f"[!] Fixing charges for: {pdbqt_file}")

    # Try OpenBabel only (most reliable baseline)
    try:
        runner.run([
            runner.OBABEL,
            "-isdf" if source_file.endswith(".sdf") else "-ipdb",
            source_file,
            "-opdbqt",
            "-O",
            pdbqt_file,
            "--partialcharge",
            "gasteiger"
        ], capture=False)

        if _ensure_pdbqt_has_charges(pdbqt_file):
            logger.info("[✔] Charges assigned via OpenBabel Gasteiger")
            return

    except Exception as e:
        logger.error(f"Charge assignment failed: {e}")

    raise RuntimeError(
        f"CRITICAL: Cannot assign valid charges for {pdbqt_file}. "
        "Install OpenBabel with Gasteiger support."
    )


def _assign_simple_charges(pdbqt_file: str) -> None:
    """Deprecated legacy helper (kept for backward compatibility; not used)."""
    # Simple charge assignment based on common atom types
    charge_map = {
        'C': 0.0,    # Carbon
        'N': -0.3,   # Nitrogen
        'O': -0.5,   # Oxygen
        'S': -0.3,   # Sulfur
        'P': 1.1,    # Phosphorus
        'H': 0.0,    # Hydrogen (connected to C)
        'HD': 0.12,  # Polar H
        'OA': -0.7,  # Acceptor O
        'NA': -0.7,  # Acceptor N
        'SA': -0.7,  # Acceptor S
    }

    lines = []
    modified = False

    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Check if charge is zero (needs assignment)
                    if len(line) > 76:
                        try:
                            charge_str = line[70:76].strip()
                            charge = float(charge_str) if charge_str else 0.0

                            if charge == 0.0:
                                # Get atom type from columns 77-78
                                atom_type = line[77:79].strip() if len(
                                    line) > 78 else ''

                                # Look up charge
                                new_charge = charge_map.get(atom_type, 0.0)

                                if new_charge != 0.0:
                                    # Replace charge in the line
                                    charge_field = f"{new_charge:+7.3f}"
                                    line_list = list(line)
                                    for i, c in enumerate(charge_field):
                                        if 70 + i <= 76 and i < len(charge_field):
                                            if 70 + i < len(line_list):
                                                line_list[70 + i] = c
                                    line = ''.join(
                                        line_list[:81]) + (line[81:] if len(line) > 81 else '')
                                    modified = True
                        except (ValueError, IndexError):
                            pass

                lines.append(line)

        # Write back if modified
        if modified:
            with open(pdbqt_file, 'w') as f:
                f.writelines(lines)
            logger.debug(
                f"Heuristic charges assigned to {os.path.basename(pdbqt_file)}")
    except Exception as e:
        logger.debug(f"Could not assign charges heuristically: {e}")
