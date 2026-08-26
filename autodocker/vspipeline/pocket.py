#!/usr/bin/env python3
"""PHASE 2: SMART PREPROCESSING.

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

WATER_RESNAMES = {"HOH", "WAT", "TIP3", "TIP",
                  "TIP3P", "TIP4", "H2O", "DOD", "SOL"}

METAL_ELEMENTS = {
    "ZN": "Zinc", "CU": "Copper", "FE": "Iron", "MG": "Magnesium",
    "CA": "Calcium", "MN": "Manganese", "NI": "Nickel",
    "CO": "Cobalt", "LI": "Lithium", "K": "Potassium", "NA": "Sodium",
    "AL": "Aluminum", "CD": "Cadmium", "HG": "Mercury",
}

COFACTOR_RESNAMES = {
    "NAD", "NADH", "NADP", "NADPH", "ATP", "ADP", "GTP", "GDP",
    "FAD", "FADH", "FMN", "HEM", "HEME", "SAM", "SAH", "PLP",
    "UDP", "UMP", "UTP", "CMP", "ACP", "COA", "THF", "H4B",
}


def _parse_hetatm(line: str) -> Optional[Dict]:
    """Parse a PDB HETATM line into a dict, or None if malformed.

    Element is taken from the element column with fallback to the atom name
    (many PDB files leave the element column blank).
    """
    if not line.startswith("HETATM"):
        return None
    try:
        coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    except (ValueError, IndexError):
        return None

    atom_name = line[12:16].strip()
    element = line[76:78].strip().upper()
    if not element:
        element = re.sub(r"[^A-Za-z]", "", atom_name).upper()

    return {
        "index": line[6:11].strip(),
        "atom": atom_name,
        "resname": line[17:20].strip().upper(),
        "chain": line[21].strip() or "A",
        "resnum": line[22:26].strip(),
        "element": element,
        "coords": coords,
        "line": line,
    }


def _iter_hetatm_residues(pdb_file: str, predicate) -> List[Dict]:
    """Parse all HETATM residues from a PDB file matching ``predicate``.

    Detection runs on the *original* PDB so HETATM records (waters, metal
    ions, cofactors) are still present.
    """
    residues = []
    seen = set()
    try:
        with open(pdb_file, "r") as f:
            for line in f:
                parsed = _parse_hetatm(line)
                if parsed is None or not predicate(parsed):
                    continue
                key = (parsed["chain"], parsed["resnum"], parsed["resname"])
                if key in seen:
                    continue
                seen.add(key)
                residues.append(parsed)
    except Exception as e:
        logger.debug(f"Could not detect hetero residues in {pdb_file}: {e}")
    return residues


def detect_water_molecules(pdb_file: str, binding_site_coords: Tuple[float, float, float],
                           distance_threshold: float = 4.0) -> List[Dict]:
    """Detect water molecules near the binding site.

    Returns a list of residue dicts (chain/resnum/resname/coords) for waters
    whose oxygen lies within ``distance_threshold`` A of the site center.
    """
    cx, cy, cz = binding_site_coords
    waters = []
    for res in _iter_hetatm_residues(
            pdb_file, lambda r: r["resname"] in WATER_RESNAMES):
        x, y, z = res["coords"]
        if math.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2) <= distance_threshold:
            waters.append(res)
    return waters


def detect_metal_ions(pdb_file: str) -> List[Dict]:
    """Detect metal ions in the protein. Returns list with ion info."""
    metals = []
    for res in _iter_hetatm_residues(
            pdb_file, lambda r: r["element"] in METAL_ELEMENTS):
        metals.append({
            "element": res["element"],
            "name": METAL_ELEMENTS[res["element"]],
            "chain": res["chain"],
            "resnum": res["resnum"],
            "resname": res["resname"],
            "coords": res["coords"],
        })
    return metals


def detect_cofactors(pdb_file: str) -> List[Dict]:
    """Detect cofactors/ligands already present in the protein."""
    cofactors = []
    for res in _iter_hetatm_residues(
            pdb_file, lambda r: r["resname"] in COFACTOR_RESNAMES):
        cofactors.append({
            "name": res["resname"],
            "chain": res["chain"],
            "resnum": res["resnum"],
            "coords": res["coords"],
        })
    return cofactors


VINA_METAL_TYPE_MAP = {
    "ZN": "Zn", "FE": "Fe", "CA": "Ca", "MG": "Mg", "MN": "Mn",
    "CU": "Cu", "NI": "Ni", "CO": "Co", "HG": "Hg", "CD": "Cd",
    "NA": "Na", "K": "K",
}


def _max_pdbqt_serial(pdbqt_file: str) -> int:
    """Return the largest atom serial in a PDBQT file (0 if none)."""
    max_serial = 0
    try:
        with open(pdbqt_file, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        max_serial = max(max_serial, int(line[6:11]))
                    except ValueError:
                        pass
    except IOError:
        pass
    return max_serial


def _hetatm_to_pdbqt_line(parsed: Dict, serial: int) -> Optional[str]:
    """Build a receptor PDBQT ATOM line for a kept HETATM record.

    Water oxygens are typed OA (matching OpenBabel); metals use the element
    symbol for the types Vina accepts; other elements get conservative
    AutoDock types. Returns None for unsupported elements.
    """
    element = parsed["element"]
    if element == "O":
        atom_type = "OA"
    elif element == "N":
        atom_type = "NA"
    elif element == "S":
        atom_type = "SA"
    elif element == "H":
        atom_type = "HD"
    elif element == "C":
        atom_type = "C"
    elif element in VINA_METAL_TYPE_MAP:
        atom_type = VINA_METAL_TYPE_MAP[element]
    else:
        return None

    x, y, z = parsed["coords"]
    name = parsed["atom"] or element
    resname = parsed["resname"]
    chain = parsed["chain"]
    resnum = parsed["resnum"]
    # Column layout matches the AutoDock PDBQT convention Vina parses:
    # coords at 31-38/39-46/47-54, charge at 69-76, atom type at 78-79.
    return (
        f"HETATM{serial:5d} {name:<4s} {resname:3s} {chain} {resnum:>3s}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}    {0.000:+.3f} {atom_type}\n"
    )


def _append_hetatm_to_receptor(pdbqt_file: str, pdb_file: str,
                               keep_residues, selected_chains) -> int:
    """Append PDBQT records for the selected HETATM residues to a receptor.

    Works on the *original* PDB so the HETATM records still exist. Returns
    the number of atoms appended. Existing protein charges are preserved.
    """
    keep = set(keep_residues)
    serial = _max_pdbqt_serial(pdbqt_file) + 1
    appended = 0

    with open(pdb_file, "r") as src, open(pdbqt_file, "a") as dst:
        for line in src:
            if not line.startswith("HETATM"):
                continue
            parsed = _parse_hetatm(line)
            if parsed is None:
                continue
            if parsed["chain"] not in selected_chains:
                continue
            key = (parsed["chain"], parsed["resnum"], parsed["resname"])
            if key not in keep:
                continue
            record = _hetatm_to_pdbqt_line(parsed, serial)
            if record is None:
                logger.warning(
                    f"[!] Skipping unsupported hetero atom {parsed['atom']} "
                    f"({parsed['element']}) in {parsed['resname']} {parsed['resnum']}")
                continue
            dst.write(record)
            serial += 1
            appended += 1

    return appended
