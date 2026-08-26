#!/usr/bin/env python3
"""PHASE 3: FLEXIBLE DOCKING.

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

def detect_flexible_residues(pdbqt_file: str, binding_site_coords: Tuple[float, float, float],
                             radius: float = 8.0, max_residues: int = 10) -> List[str]:
    """Auto-detect residues near binding site that should be flexible."""
    flexible = []

    try:
        cx, cy, cz = binding_site_coords
        residues = {}

        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        dist = math.sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2)
                        if dist <= radius:
                            chain = line[21].strip()
                            res_num = line[22:26].strip()
                            res_key = f"{chain}{res_num}"

                            if res_key not in residues:
                                residues[res_key] = dist
                    except (ValueError, IndexError):
                        pass

        # Sort by distance and take closest
        sorted_residues = sorted(residues.items(), key=lambda x: x[1])
        flexible = [res[0] for res in sorted_residues[:max_residues]]

    except Exception as e:
        logger.debug(f"Could not detect flexible residues: {e}")

    return flexible


FLEX_BACKBONE_NAMES = {"N", "CA", "C", "O", "H"}
_BOND_TOLERANCE = 0.45

_ATOMIC_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05, "P": 1.07,
    "F": 0.57, "CL": 1.02, "BR": 1.20, "I": 1.39,
    "ZN": 1.39, "FE": 1.32, "CA": 1.76, "MG": 1.41, "MN": 1.39,
    "CU": 1.32, "NI": 1.24, "CO": 1.26, "HG": 1.32, "CD": 1.44,
    "NA": 1.66, "K": 2.03,
}


def _pdbqt_atom_element(atom_type: str) -> str:
    """Map a PDBQT atom type back to its element for bond detection."""
    if atom_type in ("C", "A"):
        return "C"
    if atom_type in ("N", "NA"):
        return "N"
    if atom_type == "OA":
        return "O"
    if atom_type == "SA":
        return "S"
    if atom_type in ("HD", "H"):
        return "H"
    if atom_type in _ATOMIC_RADII:
        return atom_type
    return "C"


def _parse_residue_key(key: str) -> Optional[Tuple[str, str]]:
    """Parse a residue key like 'A45' into (chain, resnum)."""
    match = re.match(r"^([A-Za-z]+?)(\d+)$", key.strip())
    if not match:
        return None
    return match.group(1).upper(), match.group(2)


def _residue_atoms_from_pdbqt(pdbqt_file: str,
                              key: str) -> List[Dict]:
    """Collect PDBQT atoms belonging to a residue identified by 'A45'."""
    parsed = _parse_residue_key(key)
    if parsed is None:
        logger.warning(f"[!] Invalid residue key '{key}' (expected e.g. A45)")
        return []
    chain, resnum = parsed
    atoms = []
    with open(pdbqt_file, "r") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_chain = line[21].strip().upper()
            atom_resnum = line[22:26].strip()
            if atom_chain != chain:
                continue
            if not _resnum_matches(atom_resnum, resnum):
                continue
            try:
                serial = int(line[6:11])
            except ValueError:
                continue
            name = line[12:16].strip()
            parts = line.split()
            atom_type = parts[-1] if parts else "C"
            atoms.append({
                "serial": serial,
                "name": name,
                "atom_type": atom_type,
                "element": _pdbqt_atom_element(atom_type),
                "coords": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                "line": line if line.endswith("\n") else line + "\n",
            })
    return atoms


def _resnum_matches(atom_resnum: str, key_resnum: str) -> bool:
    """Compare a PDB resnum (possibly '45') with a key resnum ('45')."""
    if atom_resnum.isdigit() and key_resnum.isdigit():
        return int(atom_resnum) == int(key_resnum)
    return atom_resnum == key_resnum


def _build_flex_residue_block(atoms: List[Dict]) -> List[str]:
    """Build one BEGIN_RES ... END_RES block for a flexible residue.

    The tree is rooted at CA; backbone atoms (N/CA/C/O/H) are the immobile
    ROOT, and side-chain branches are emitted as nested BRANCH/ENDBRANCH
    pairs (the Vina 1.2 flexible-receptor format).
    """
    by_name = {a["name"]: a for a in atoms}
    ca = by_name.get("CA")
    if ca is None:
        logger.warning(
            f"[!] Cannot build flex residue: no CA atom found")
        return []

    root_atoms = [a for a in atoms if a["name"] in FLEX_BACKBONE_NAMES]
    side = [a for a in atoms if a["name"] not in FLEX_BACKBONE_NAMES]
    if not side:
        logger.warning(
            f"[!] Residue has no flexible side-chain atoms; skipping")
        return []

    by_serial = {a["serial"]: a for a in atoms}
    nodes = side + [ca]
    node_serials = {a["serial"] for a in nodes}

    def bonded(a: Dict, b: Dict) -> bool:
        radius = (_ATOMIC_RADII.get(a["element"], 0.76)
                  + _ATOMIC_RADII.get(b["element"], 0.76) + _BOND_TOLERANCE)
        dx = a["coords"][0] - b["coords"][0]
        dy = a["coords"][1] - b["coords"][1]
        dz = a["coords"][2] - b["coords"][2]
        return math.sqrt(dx * dx + dy * dy + dz * dz) < radius

    # Grow a tree from CA, visiting each side-chain atom at most once.
    tree = {}  # serial -> list of child serials
    parent_of = {ca["serial"]: None}

    def grow(node_serial: int):
        node = by_serial[node_serial]
        neighbors = [nb for nb in nodes
                     if nb["serial"] in node_serials
                     and nb["serial"] != node_serial
                     and bonded(node, nb)]
        for nb in sorted(neighbors, key=lambda x: x["serial"]):
            if nb["serial"] in parent_of:
                continue
            parent_of[nb["serial"]] = node_serial
            tree.setdefault(node_serial, []).append(nb["serial"])
            grow(nb["serial"])

    grow(ca["serial"])
    if not tree.get(ca["serial"]):
        logger.warning(
            f"[!] No side-chain atoms bonded to CA; skipping residue")
        return []

    lines = ["BEGIN_RES\n"]
    lines.append("ROOT\n")
    for a in sorted(root_atoms, key=lambda x: x["serial"]):
        lines.append(a["line"])
    lines.append("ENDROOT\n")

    def emit(node_serial: int):
        for child_serial in tree.get(node_serial, []):
            child = by_serial[child_serial]
            lines.append(f"BRANCH {node_serial} {child_serial}\n")
            lines.append(child["line"])
            emit(child_serial)
            lines.append(f"ENDBRANCH {node_serial} {child_serial}\n")

    emit(ca["serial"])
    lines.append("END_RES\n")
    return lines


def build_flexible_residue_pdbqt(receptor_pdbqt: str, residues: List[str],
                                 outdir: str) -> Optional[str]:
    """Build a Vina flexible-receptor PDBQT for the given residue keys.

    Residues are given as e.g. ``['A45', 'B102']``. Returns the flex file
    path, or None if no residue could be converted.
    """
    if not residues:
        return None

    flex_file = os.path.join(outdir, "flex.pdbqt")
    written = 0
    with open(flex_file, "w") as dst:
        dst.write("REMARK flexible residues: " + ", ".join(residues) + "\n")
        for key in residues:
            atoms = _residue_atoms_from_pdbqt(receptor_pdbqt, key)
            if not atoms:
                logger.warning(
                    f"[!] Residue '{key}' not found in receptor; skipping")
                continue
            block = _build_flex_residue_block(atoms)
            if not block:
                logger.warning(f"[!] Could not build flex block for '{key}'")
                continue
            dst.writelines(block)
            written += 1

    if written == 0:
        logger.warning(
            f"[!] No flexible residues could be built; flex docking disabled")
        return None

    logger.info(
        f"[✔] Flexible receptor PDBQT written: {flex_file} "
        f"({written} residue(s))")
    return flex_file


def _remove_residues_from_pdbqt(receptor_pdbqt: str, residues: List[str],
                                outpath: str) -> str:
    """Build a rigid receptor PDBQT excluding the flexible residue atoms.

    Vina requires flexible residues to be absent from the rigid receptor;
    including them double-counts atoms and corrupts the scoring.
    """
    keys = []
    for key in residues:
        parsed = _parse_residue_key(key)
        if parsed is not None:
            keys.append(parsed)

    removed = 0
    kept = []
    with open(receptor_pdbqt, "r") as src:
        for line in src:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21].strip().upper()
                resnum = line[22:26].strip()
                if any(c == chain and _resnum_matches(resnum, rn)
                       for c, rn in keys):
                    removed += 1
                    continue
            kept.append(line)

    with open(outpath, "w") as dst:
        dst.writelines(kept)

    if removed:
        logger.info(
            f"[*] Rigid receptor: removed {removed} flexible-residue atoms "
            f"({outpath})")
    return outpath
