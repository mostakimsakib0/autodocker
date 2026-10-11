#!/usr/bin/env python3
"""HETATM handling: detection of waters, metal ions and cofactors, and
appending the kept ones to the receptor PDBQT.

Detection and appending read the original PDB through the same record
filter as receptor cleaning (first model, one alternate location per
residue), so they see the same atoms.
"""
import logging
import math
from collections import OrderedDict
from typing import Dict, List, Optional, Tuple

from vspipeline.protein import (
    BLANK_CHAIN_ID,
    OB_BOND_TOLERANCE,
    OB_COVALENT_RADII,
    OB_UNKNOWN_COVALENT_RADIUS,
    _element_from_atom_name,
    read_pdb_records,
)

logger = logging.getLogger(__name__)

__all__ = [
    "COFACTOR_RESNAMES", "METAL_ELEMENTS", "VINA_METAL_TYPE_MAP", "WATER_RESNAMES",
    "VINA_ELEMENT_TYPES", "_append_hetatm_to_receptor", "_element_from_atom_name",
    "_hetatm_to_pdbqt_line", "_iter_hetatm_residues", "_max_pdbqt_serial",
    "_parse_hetatm", "detect_cofactors", "detect_metal_ions", "detect_water_molecules",
]

# Residue names are matched against the three-character PDB resname field
# (columns 18-20), so only names of up to three characters can match.
WATER_RESNAMES = {"HOH", "WAT", "H2O", "DOD", "SOL", "TIP", "T3P", "T4P", "SPC"}

# PDB element symbol -> AutoDock Vina receptor atom type for metal ions. The
# types are Vina's AD types (Mg, Mn, Zn, Ca, Fe) or the metals Vina accepts
# without AD parameters (Cu, Fe, Na, K, Hg, Co, Cd, Ni; atom_constants.h
# non_ad_metal_names).
VINA_METAL_TYPE_MAP = {
    "ZN": "Zn", "FE": "Fe", "CA": "Ca", "MG": "Mg", "MN": "Mn",
    "CU": "Cu", "NI": "Ni", "CO": "Co", "HG": "Hg", "CD": "Cd",
    "NA": "Na", "K": "K",
}

# Metals detected by --detect-metals: exactly those that can be written
# to the receptor PDBQT (see VINA_METAL_TYPE_MAP).
METAL_ELEMENTS = {
    "ZN": "Zinc", "CU": "Copper", "FE": "Iron", "MG": "Magnesium",
    "CA": "Calcium", "MN": "Manganese", "NI": "Nickel",
    "CO": "Cobalt", "K": "Potassium", "NA": "Sodium",
    "CD": "Cadmium", "HG": "Mercury",
}

COFACTOR_RESNAMES = {
    "NAD", "NAP", "NDP", "NAI", "ATP", "ADP", "GTP", "GDP",
    "FAD", "FMN", "HEM", "SAM", "SAH", "PLP",
    "UDP", "UMP", "UTP", "CMP", "ACP", "COA", "THF", "H4B",
}

# Non-metal elements of kept HETATM residues -> Vina receptor atom type.
# N and H depend on the bonded atoms (see _residue_atom_types).
VINA_ELEMENT_TYPES = {
    "C": "C", "O": "OA", "S": "SA", "P": "P", "F": "F", "CL": "Cl",
    "BR": "Br", "I": "I", "N": "NA", "H": "HD",
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
        element = _element_from_atom_name(line[12:16])

    return {
        "index": line[6:11].strip(),
        "atom": atom_name,
        "resname": line[17:20].strip().upper(),
        "chain": line[21].strip() or BLANK_CHAIN_ID,
        "resnum": line[22:26].strip(),
        "icode": line[26:27].strip() if len(line) > 26 else "",
        "element": element,
        "coords": coords,
        "line": line,
    }


def _residue_key(parsed: Dict) -> Tuple[str, str, str]:
    return parsed["chain"], parsed["resnum"], parsed["resname"]


def _hetatm_residues(pdb_file: str) -> "OrderedDict[Tuple[str, str, str], List[Dict]]":
    """All HETATM atoms of ``pdb_file`` (first model, chosen altlocs),
    grouped by (chain, resnum, resname) in file order."""
    lines, _ = read_pdb_records(pdb_file)
    residues = OrderedDict()
    for line in lines:
        parsed = _parse_hetatm(line)
        if parsed is not None:
            residues.setdefault(_residue_key(parsed), []).append(parsed)
    return residues


def _iter_hetatm_residues(pdb_file: str, predicate) -> List[Dict]:
    """The first atom matching ``predicate`` of each HETATM residue.

    Detection runs on the *original* PDB so HETATM records (waters, metal
    ions, cofactors) are still present. Raises if the file cannot be read.
    """
    residues = []
    for atoms in _hetatm_residues(pdb_file).values():
        for parsed in atoms:
            if predicate(parsed):
                residues.append(parsed)
                break
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
            pdb_file, lambda r: r["resname"] in WATER_RESNAMES and r["element"] == "O"):
        x, y, z = res["coords"]
        if math.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2) <= distance_threshold:
            waters.append(res)
    return waters


def detect_metal_ions(pdb_file: str) -> List[Dict]:
    """Detect metal ions: HETATM residues made of a single metal atom.

    A metal inside a larger residue (the iron of HEM, the irons of SF4) is
    not an ion; keeping its residue key would keep the whole residue, so it
    is skipped with a warning (``--detect-cofactors`` keeps known cofactor
    residues whole). Returns one dict per ion, including its atom name.
    """
    metals = []
    for key, atoms in _hetatm_residues(pdb_file).items():
        metal_atoms = [a for a in atoms if a["element"] in METAL_ELEMENTS]
        if not metal_atoms:
            continue
        if len(atoms) > 1:
            chain, resnum, resname = key
            logger.warning(
                f"[!] {', '.join(a['element'] for a in metal_atoms)} in {resname} "
                f"{chain}{resnum} is part of a {len(atoms)}-atom residue, not an "
                "ion; not kept as a metal (use --detect-cofactors or --keep-hetero "
                "to keep the residue)")
            continue
        ion = metal_atoms[0]
        metals.append({
            "element": ion["element"],
            "name": METAL_ELEMENTS[ion["element"]],
            "atom": ion["atom"],
            "chain": ion["chain"],
            "resnum": ion["resnum"],
            "resname": ion["resname"],
            "coords": ion["coords"],
        })
    return metals


def detect_cofactors(pdb_file: str) -> List[Dict]:
    """Detect cofactors/ligands already present in the protein."""
    cofactors = []
    for res in _iter_hetatm_residues(
            pdb_file, lambda r: r["resname"] in COFACTOR_RESNAMES):
        cofactors.append({
            "name": res["resname"],
            "resname": res["resname"],
            "chain": res["chain"],
            "resnum": res["resnum"],
            "coords": res["coords"],
        })
    return cofactors


def _max_pdbqt_serial(pdbqt_file: str) -> int:
    """Return the largest atom serial in a PDBQT file (0 if none)."""
    max_serial = 0
    with open(pdbqt_file, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    max_serial = max(max_serial, int(line[6:11]))
                except ValueError:
                    pass
    return max_serial


def _hetatm_to_pdbqt_line(parsed: Dict, serial: int,
                          atom_type: Optional[str] = None) -> Optional[str]:
    """Build a receptor PDBQT ATOM line for a kept HETATM record.

    ``atom_type`` overrides the element-based type (see
    ``_residue_atom_types``). Without it: water and other oxygens are OA
    (matching Open Babel), N is NA, S is SA, H is HD, C is C, P/F/Cl/Br/I
    keep their symbol, and metals use the Vina metal types. Returns None for
    elements Vina cannot type. The charge is 0.000 (Vina ignores charges).
    """
    element = parsed["element"]
    if atom_type is None:
        atom_type = VINA_ELEMENT_TYPES.get(element) or VINA_METAL_TYPE_MAP.get(element)
    if atom_type is None:
        return None

    x, y, z = parsed["coords"]
    name = parsed["atom"] or element
    resname = parsed["resname"]
    chain = parsed["chain"]
    resnum = parsed["resnum"]
    # Column layout matches the AutoDock PDBQT convention Vina parses:
    # coords at 31-38/39-46/47-54, charge at 69-76, atom type at 78-79.
    icode = parsed.get("icode", "")
    # Fixed PDB columns: serial 7-11 (wrapped at 5 digits), name 13-16,
    # resname 18-20, chain 22, resnum 23-26, insertion code 27.
    return (
        f"HETATM{serial % 100000:5d} {name:<4s} {resname:>3s} {chain:1s}{resnum:>4s}"
        f"{icode:1s}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}    {0.000:+.3f} {atom_type}\n"
    )


def _bonded(a: Dict, b: Dict) -> bool:
    """Covalent bond by Open Babel's distance rule."""
    cutoff = (OB_COVALENT_RADII.get(a["element"], OB_UNKNOWN_COVALENT_RADIUS)
              + OB_COVALENT_RADII.get(b["element"], OB_UNKNOWN_COVALENT_RADIUS)
              + OB_BOND_TOLERANCE)
    return math.dist(a["coords"], b["coords"]) < cutoff


def _residue_atom_types(atoms: List[Dict]) -> List[Optional[str]]:
    """Vina atom type for each atom of one kept HETATM residue (None: drop).

    Without bond orders the typing is approximate: an N bonded to three or
    more heavy atoms of the residue is N (no free lone pair, e.g. the amide
    or pyrrole-type N), any other N is NA (acceptor). Deposited hydrogens on
    N/O/S are kept as HD; hydrogens on carbon are dropped (united-atom
    receptor, as Open Babel's -xr writes it). Hydrogens are not added, so an
    N/O of a kept residue is a Vina donor only if the input PDB contains its
    hydrogen. Raises ValueError if an atom has an element Vina cannot type,
    so a residue is never kept partially.
    """
    heavy = [a for a in atoms if a["element"] not in ("H", "D")]
    unsupported = sorted({a["element"] for a in atoms
                          if a["element"] not in ("H", "D")
                          and a["element"] not in VINA_ELEMENT_TYPES
                          and a["element"] not in VINA_METAL_TYPE_MAP})
    if unsupported:
        a = atoms[0]
        raise ValueError(
            f"Cannot keep {a['resname']} {a['chain']}{a['resnum']}: element(s) "
            f"{', '.join(unsupported)} have no AutoDock Vina atom type")
    types = []
    for atom in atoms:
        el = atom["element"]
        if el in ("H", "D"):
            parent = min(heavy, key=lambda b: math.dist(b["coords"], atom["coords"]),
                         default=None)
            types.append("HD" if parent is not None and parent["element"] in ("N", "O", "S")
                         and _bonded(atom, parent) else None)
        elif el == "N":
            n_heavy = sum(1 for b in heavy if b is not atom
                          and b["element"] not in VINA_METAL_TYPE_MAP and _bonded(atom, b))
            types.append("N" if n_heavy >= 3 else "NA")
        elif el in VINA_ELEMENT_TYPES:
            types.append(VINA_ELEMENT_TYPES[el])
        else:
            types.append(VINA_METAL_TYPE_MAP[el])
    return types


def _pdbqt_residue_atoms(pdbqt_file: str) -> set:
    """(resname, coordinate columns) of every atom of a PDBQT file."""
    with open(pdbqt_file) as f:
        return {(line[17:20].strip(), line[30:54]) for line in f
                if line.startswith(("ATOM", "HETATM"))}


def _append_hetatm_to_receptor(pdbqt_file: str, pdb_file: str,
                               keep_residues, selected_chains,
                               offset: Tuple[float, float, float] = (0.0, 0.0, 0.0)) -> int:
    """Append PDBQT records for the selected HETATM residues to a receptor.

    Works on the *original* PDB so the HETATM records still exist, read
    with the same first-model and altloc filter as receptor cleaning.
    ``offset`` is subtracted from each coordinate to move the atoms into the
    receptor PDBQT's frame (``prepare_receptor`` translates it by the mean of
    the receptor heavy atoms). Residues whose heavy atoms are all already in
    the receptor (same residue name and coordinates, e.g. kept by
    --keep-hetero) are skipped; a residue found
    only partly is an error. Returns the number of atoms appended. Existing
    protein charges are preserved.
    """
    keep = set(keep_residues)
    present = _pdbqt_residue_atoms(pdbqt_file)
    serial = _max_pdbqt_serial(pdbqt_file) + 1
    records = []
    skipped = []

    for key, atoms in _hetatm_residues(pdb_file).items():
        if key[0] not in selected_chains or key not in keep:
            continue
        moved = []
        for parsed in atoms:
            x, y, z = parsed["coords"]
            moved.append(dict(parsed, coords=(x - offset[0], y - offset[1], z - offset[2])))
        heavy = [a for a in moved if a["element"] not in ("H", "D")]
        found = [a for a in heavy
                 if (a["resname"], "{:8.3f}{:8.3f}{:8.3f}".format(*a["coords"])) in present]
        if found:
            if len(found) != len(heavy):
                raise ValueError(
                    f"{key[2]} {key[0]}{key[1]}: {len(found)} of {len(heavy)} heavy atoms "
                    "are already in the receptor; cannot keep the residue consistently")
            skipped.append(key)
            continue
        for parsed, atom_type in zip(moved, _residue_atom_types(moved)):
            if atom_type is None:
                continue
            records.append(_hetatm_to_pdbqt_line(parsed, serial, atom_type))
            serial += 1

    if skipped:
        logger.info(
            f"[*] {len(skipped)} requested HETATM residue(s) are already in the "
            "receptor (e.g. kept by --keep-hetero); not added twice: "
            + ", ".join(f"{n} {c}{r}" for c, r, n in skipped))
    with open(pdbqt_file, "a") as dst:
        dst.writelines(records)
    return len(records)
