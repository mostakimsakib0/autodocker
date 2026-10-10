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
    """Auto-detect residues near the binding site that can be made flexible.

    Residues are ranked by their closest atom to the site; residues without
    rotatable side-chain bonds (GLY, ALA, PRO) are skipped.
    """
    flexible = []

    try:
        cx, cy, cz = binding_site_coords
        residues = {}

        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    if line[17:20].strip() not in SIDE_CHAIN_TORSIONS:
                        continue
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        dist = math.sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2)
                        if dist <= radius:
                            res_key = (f"{line[21].strip()}{line[22:26].strip()}"
                                       f"{line[26].strip()}")
                            residues[res_key] = min(dist, residues.get(res_key, dist))
                    except (ValueError, IndexError):
                        pass

        # Sort by distance and take closest
        sorted_residues = sorted(residues.items(), key=lambda x: (x[1], x[0]))
        flexible = [res[0] for res in sorted_residues[:max_residues]]

    except Exception as e:
        logger.debug(f"Could not detect flexible residues: {e}")

    return flexible


# Rotatable side-chain bonds (chi angles, plus the hydroxyl/thiol/amine
# rotors) per amino acid, as in AutoDockTools. Ring, amide, carboxylate and
# guanidinium bonds are not rotatable and keep their atoms in one rigid
# fragment. GLY, ALA and PRO have no side-chain torsions.
SIDE_CHAIN_TORSIONS = {
    "SER": [("CA", "CB"), ("CB", "OG")],
    "THR": [("CA", "CB"), ("CB", "OG1")],
    "CYS": [("CA", "CB"), ("CB", "SG")],
    "VAL": [("CA", "CB")],
    "LEU": [("CA", "CB"), ("CB", "CG")],
    "ILE": [("CA", "CB"), ("CB", "CG1")],
    "MET": [("CA", "CB"), ("CB", "CG"), ("CG", "SD")],
    "PHE": [("CA", "CB"), ("CB", "CG")],
    "TYR": [("CA", "CB"), ("CB", "CG"), ("CZ", "OH")],
    "TRP": [("CA", "CB"), ("CB", "CG")],
    "HIS": [("CA", "CB"), ("CB", "CG")],
    "ASP": [("CA", "CB"), ("CB", "CG")],
    "ASN": [("CA", "CB"), ("CB", "CG")],
    "GLU": [("CA", "CB"), ("CB", "CG"), ("CG", "CD")],
    "GLN": [("CA", "CB"), ("CB", "CG"), ("CG", "CD")],
    "LYS": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")],
    "ARG": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE")],
}
for _alias, _base in (("HID", "HIS"), ("HIE", "HIS"), ("HIP", "HIS"), ("CYX", "CYS")):
    SIDE_CHAIN_TORSIONS[_alias] = SIDE_CHAIN_TORSIONS[_base]


FLEX_BACKBONE_NAMES = {"N", "C", "O", "H", "HN", "OXT"}
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
    if atom_type in _ATOMIC_RADII or atom_type.upper() in _ATOMIC_RADII:
        return atom_type if atom_type in _ATOMIC_RADII else atom_type.upper()
    return "C"


def _parse_residue_key(key: str) -> Optional[Tuple[str, str]]:
    """Parse a residue key like 'A45' (or 'A100B' with an insertion code, or
    '45' for a blank chain) into (chain, resnum-with-icode)."""
    match = re.match(r"^([A-Za-z]?)(-?\d+)([A-Za-z]?)$", key.strip())
    if not match:
        return None
    return match.group(1).upper(), match.group(2) + match.group(3).upper()


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
            atom_resnum = line[22:27].strip()
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
                "resname": line[17:20].strip(),
                "atom_type": atom_type,
                "element": _pdbqt_atom_element(atom_type),
                "coords": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                "line": line if line.endswith("\n") else line + "\n",
            })
    return atoms


def _resnum_matches(atom_resnum: str, key_resnum: str) -> bool:
    """Compare a PDB resnum plus insertion code ('45', '100B') with a key's."""
    def split(r):
        m = re.match(r"^(-?\d+)([A-Za-z]?)$", r.strip())
        return (int(m.group(1)), m.group(2).upper()) if m else (None, r.strip())
    return split(atom_resnum) == split(key_resnum)


def _bonded(a: Dict, b: Dict) -> bool:
    """Covalent bond by distance (sum of covalent radii plus a tolerance)."""
    radius = (_ATOMIC_RADII.get(a["element"], 0.76)
              + _ATOMIC_RADII.get(b["element"], 0.76) + _BOND_TOLERANCE)
    dx = a["coords"][0] - b["coords"][0]
    dy = a["coords"][1] - b["coords"][1]
    dz = a["coords"][2] - b["coords"][2]
    return math.sqrt(dx * dx + dy * dy + dz * dz) < radius


def _build_flex_residue_block(atoms: List[Dict]) -> List[str]:
    """Build one BEGIN_RES ... END_RES block for a flexible residue.

    Follows the AutoDock flexible-receptor convention: the ROOT is CA, the
    side chain hangs from it, and the backbone (N, C, O, H) stays in the
    rigid receptor. Only the residue's rotatable side-chain bonds become
    BRANCHes; ring and other rigid groups stay together in one fragment.
    """
    resname = atoms[0].get("resname", "") if atoms else ""
    torsions = SIDE_CHAIN_TORSIONS.get(resname)
    if torsions is None:
        logger.warning(f"[!] {resname or 'Residue'} has no rotatable side chain "
                       "(or is not a standard amino acid); skipping")
        return []
    by_name = {}
    for a in atoms:
        by_name.setdefault(a["name"], a)
    ca = by_name.get("CA")
    if ca is None:
        logger.warning(
            f"[!] Cannot build flex residue: no CA atom found")
        return []

    # Open Babel names every added hydrogen "H", so a hydrogen belongs to the
    # backbone only if it is bonded to the backbone N; side-chain hydrogens
    # (Ser OG-H, Lys NZ-H, ...) move with their heavy atom.
    backbone_n = by_name.get("N")
    nodes = []
    for a in atoms:
        if a["element"] == "H":
            if backbone_n is not None and _bonded(a, backbone_n):
                continue
        elif a["name"] in FLEX_BACKBONE_NAMES:
            continue
        nodes.append(a)
    if len(nodes) < 2:
        logger.warning(
            f"[!] Residue has no flexible side-chain atoms; skipping")
        return []

    bonded = _bonded
    serials = [a["serial"] for a in nodes]
    by_serial = {a["serial"]: a for a in nodes}
    neighbors = {s: [] for s in serials}
    heavy = [a for a in nodes if a["element"] != "H"]
    for i, a in enumerate(heavy):
        for b in heavy[i + 1:]:
            if bonded(a, b):
                neighbors[a["serial"]].append(b["serial"])
                neighbors[b["serial"]].append(a["serial"])
    # A hydrogen bonds to exactly one heavy atom: its nearest bonded one.
    for h in (a for a in nodes if a["element"] == "H"):
        partners = [b for b in heavy if bonded(h, b)]
        if partners:
            near = min(partners, key=lambda b: sum(
                (h["coords"][k] - b["coords"][k]) ** 2 for k in range(3)))
            neighbors[h["serial"]].append(near["serial"])
            neighbors[near["serial"]].append(h["serial"])

    rotatable = set()
    for n1, n2 in torsions:
        a, b = by_name.get(n1), by_name.get(n2)
        if a and b and b["serial"] in neighbors.get(a["serial"], []):
            rotatable.add(frozenset((a["serial"], b["serial"])))

    # Rigid fragments: atoms connected without crossing a rotatable bond.
    fragment_of = {}
    fragments = []
    for start in serials:
        if start in fragment_of:
            continue
        idx = len(fragments)
        stack, members = [start], []
        fragment_of[start] = idx
        while stack:
            cur = stack.pop()
            members.append(cur)
            for nb in neighbors[cur]:
                if nb not in fragment_of and frozenset((cur, nb)) not in rotatable:
                    fragment_of[nb] = idx
                    stack.append(nb)
        fragments.append(sorted(members))

    lines = ["BEGIN_RES\n", "ROOT\n"]
    lines += [by_serial[s]["line"] for s in fragments[fragment_of[ca["serial"]]]]
    lines.append("ENDROOT\n")
    placed = {fragment_of[ca["serial"]]}

    def emit(frag_idx: int):
        # Branch out of every atom of this fragment across rotatable bonds.
        for atom_serial in fragments[frag_idx]:
            for nb in sorted(neighbors[atom_serial]):
                child = fragment_of[nb]
                if child in placed or frozenset((atom_serial, nb)) not in rotatable:
                    continue
                placed.add(child)
                lines.append(f"BRANCH {atom_serial:>3d} {nb:>3d}\n")
                lines.extend(by_serial[s]["line"] for s in fragments[child])
                emit(child)
                lines.append(f"ENDBRANCH {atom_serial:>3d} {nb:>3d}\n")

    emit(fragment_of[ca["serial"]])
    if len(placed) < 2:
        logger.warning(
            f"[!] No side-chain atoms bonded to CA; skipping residue")
        return []
    if len(placed) != len(fragments):
        logger.warning(f"[!] {resname}: side-chain atoms not connected to CA; skipping")
        return []
    lines.append("END_RES\n")
    return lines


def build_flexible_residue_pdbqt(receptor_pdbqt: str, residues: List[str],
                                 outdir: str,
                                 built: Optional[List[str]] = None) -> Optional[str]:
    """Build a Vina flexible-receptor PDBQT for the given residue keys.

    Residues are given as e.g. ``['A45', 'B102']``. Returns the flex file
    path, or None if no residue could be converted. When ``built`` is given,
    the keys of the residues actually written are appended to it, so only
    those are removed from the rigid receptor.
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
            if built is not None:
                built.append(key)

    if written == 0:
        logger.warning(
            f"[!] No flexible residues could be built; flex docking disabled")
        return None

    logger.info(
        f"[✔] Flexible receptor PDBQT written: {flex_file} "
        f"({written} residue(s))")
    return flex_file


def _remove_residues_from_pdbqt(receptor_pdbqt: str, residues: List[str],
                                outpath: str, flex_file: Optional[str] = None) -> str:
    """Build a rigid receptor PDBQT excluding the flexible residue atoms.

    Vina requires the flexible atoms (CA and side chain) to be absent from
    the rigid receptor; including them double-counts atoms and corrupts the
    scoring. The backbone N, C, O and H stay in the rigid receptor.

    With ``flex_file``, exactly the atoms written to it are removed (matched
    by serial, chain and residue), so the two files cannot disagree; this
    matters because Open Babel names every added hydrogen "H".
    """
    keys = []
    for key in residues:
        parsed = _parse_residue_key(key)
        if parsed is not None:
            keys.append(parsed)

    flex_atoms = None
    if flex_file and os.path.exists(flex_file):
        flex_atoms = set()
        with open(flex_file) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    flex_atoms.add((line[6:11].strip(), line[21].strip().upper(),
                                    line[22:27].strip()))

    removed = 0
    kept = []
    with open(receptor_pdbqt, "r") as src:
        for line in src:
            if line.startswith(("ATOM", "HETATM")):
                chain = line[21].strip().upper()
                resnum = line[22:27].strip()
                if flex_atoms is not None:
                    drop = (line[6:11].strip(), chain, resnum) in flex_atoms
                else:
                    drop = (line[12:16].strip() not in FLEX_BACKBONE_NAMES
                            and any(c == chain and _resnum_matches(resnum, rn)
                                    for c, rn in keys))
                if drop:
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
