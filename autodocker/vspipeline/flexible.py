#!/usr/bin/env python3
"""Flexible receptor residues: selection and the Vina flex PDBQT.

Builds the BEGIN_RES/END_RES blocks Vina reads with --flex and the rigid
receptor without those residues' side chains.
"""
import logging
import math
import os
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from vspipeline import defaults  # noqa: E402
from vspipeline import outputs  # noqa: E402

logger = logging.getLogger(__name__)

# Chain ID a blank PDB chain column is read as (protein.BLANK_CHAIN_ID,
# utils.DEFAULT_CHAIN); repeated here so this module imports nothing that
# imports runner.
BLANK_CHAIN = "A"

# PDBQT hydrogen atom types; hydrogens do not count towards residue distances.
_PDBQT_HYDROGEN_TYPES = ("H", "HD", "HS")


def detect_flexible_residues(pdbqt_file: str, binding_site_coords: Tuple[float, float, float],
                             radius: float = defaults.AUTO_FLEX_RADIUS,
                             max_residues: int = defaults.AUTO_FLEX_MAX_RESIDUES) -> List[str]:
    """Auto-detect residues near the binding site that can be made flexible.

    Residues are ranked by their closest heavy atom to the site (ties: by
    residue key); residues without rotatable side-chain bonds (GLY, ALA,
    PRO) are skipped. Raises if the receptor cannot be read.
    """
    cx, cy, cz = binding_site_coords
    residues = {}
    with open(pdbqt_file, 'r', encoding="utf-8") as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            if line[17:20].strip() not in SIDE_CHAIN_TORSIONS:
                continue
            if line[77:79].strip() in _PDBQT_HYDROGEN_TYPES:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                logger.warning(f"[!] Skipping malformed receptor line: {line.rstrip()}")
                continue
            dist = math.sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2)
            if dist <= radius:
                res_key = format_residue_key(line[21], line[22:26], line[26:27])
                residues[res_key] = min(dist, residues.get(res_key, dist))

    sorted_residues = sorted(residues.items(), key=lambda x: (x[1], x[0]))
    return [res[0] for res in sorted_residues[:max_residues]]


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
# Bond if closer than the sum of covalent radii plus this tolerance (A), as
# Open Babel's ConnectTheDots does.
_BOND_TOLERANCE = 0.45
# Covalent radius used for elements missing from _ATOMIC_RADII (that of carbon).
_DEFAULT_ATOMIC_RADIUS = 0.76

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


# Residue keys: [chain]resnum[icode] with a letter chain ('A45', 'a45',
# 'A100B'), chain:resnum[icode] for any chain character ('1:45', 'A:45'), or
# resnum[icode] alone for the blank chain ('45'). Chains are case-sensitive.
_RESIDUE_KEY = re.compile(
    r"^(?:(?P<sep_chain>[^\s:,]):(?P<sep_num>-?\d+)(?P<sep_icode>[A-Za-z]?)"
    r"|(?P<chain>[A-Za-z]?)(?P<num>-?\d+)(?P<icode>[A-Za-z]?))$")


def _parse_residue_key(key: str) -> Optional[Tuple[str, str]]:
    """Parse a residue key into (chain, resnum-with-icode), or None.

    'A45' -> ('A', '45'); 'a45' -> ('a', '45') (chain IDs are
    case-sensitive, as for --chain); 'A100B' -> ('A', '100B'); '1:45' ->
    ('1', '45') (the separator form addresses any chain, including digit
    chains); '45' -> ('', '45'), the blank chain. A blank chain is read as
    chain A everywhere (see ``chain_matches``), so 'A45' and '45' address
    the same residue of a blank-chain receptor.
    """
    match = _RESIDUE_KEY.match(key.strip())
    if not match:
        return None
    if match.group("sep_num") is not None:
        return (match.group("sep_chain"),
                match.group("sep_num") + match.group("sep_icode").upper())
    return match.group("chain"), match.group("num") + match.group("icode").upper()


def _record_chain(chain_column: str) -> str:
    """A PDB/PDBQT chain column as AutoDocker reads it: as written, a blank
    column as ``BLANK_CHAIN``."""
    return chain_column.strip() or BLANK_CHAIN


def chain_matches(chain_column: str, key_chain: str) -> bool:
    """True if a record's chain column (PDB column 22) is the chain of a
    parsed residue key. Exact and case-sensitive; a blank column and an
    empty key chain both mean ``BLANK_CHAIN``, as in receptor preparation."""
    return _record_chain(chain_column) == (key_chain or BLANK_CHAIN)


def residue_identity(key: str) -> Optional[Tuple[str, str]]:
    """(chain as read, resnum-with-icode) of a residue key, so that '45' and
    'A45' (the blank chain read as A) are the same residue; None if invalid."""
    parsed = _parse_residue_key(key)
    if parsed is None:
        return None
    chain, resnum = parsed
    return chain or BLANK_CHAIN, resnum


def format_residue_key(chain_column: str, resnum: str, icode: str = "") -> str:
    """Residue key of a record that ``_parse_residue_key`` reads back: chain
    + resnum (+ icode) for a letter chain, resnum alone for a blank chain,
    chain:resnum for any other chain character (e.g. a digit)."""
    chain = chain_column.strip()
    rest = f"{resnum.strip()}{icode.strip()}"
    if not chain or chain.isalpha():
        return f"{chain}{rest}"
    return f"{chain}:{rest}"


def missing_residues(lines: Iterable[str], selected_chains: Sequence[str],
                     keys: Iterable[str]) -> List[str]:
    """The residue keys among ``keys`` with no ATOM/HETATM record in
    ``lines`` within ``selected_chains`` (chain IDs as reported by
    ``get_chains``, a blank chain as ``BLANK_CHAIN``). Uses the same matching
    as flexible-residue building; invalid keys count as missing."""
    present = set()
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) > 26:
            chain = _record_chain(line[21])
            if chain in selected_chains:
                present.add((chain, line[22:27].strip()))
    missing = []
    for key in keys:
        identity = residue_identity(key)
        if identity is None or not any(
                c == identity[0] and _resnum_matches(r, identity[1]) for c, r in present):
            missing.append(key)
    return missing


def _residue_atoms_from_pdbqt(pdbqt_file: str,
                              key: str) -> List[Dict]:
    """Collect PDBQT atoms belonging to a residue identified by a key such
    as 'A45' (see ``_parse_residue_key``; chains match exactly)."""
    parsed = _parse_residue_key(key)
    if parsed is None:
        logger.warning(f"[!] Invalid residue key '{key}' (expected e.g. A45)")
        return []
    chain, resnum = parsed
    atoms = []
    with open(pdbqt_file, "r", encoding="utf-8") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_resnum = line[22:27].strip()
            if not chain_matches(line[21], chain):
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
    radius = (_ATOMIC_RADII.get(a["element"], _DEFAULT_ATOMIC_RADIUS)
              + _ATOMIC_RADII.get(b["element"], _DEFAULT_ATOMIC_RADIUS) + _BOND_TOLERANCE)
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
        logger.warning("[!] Cannot build flex residue: no CA atom found")
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
        logger.warning("[!] Residue has no flexible side-chain atoms; skipping")
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
        logger.warning("[!] No side-chain atoms bonded to CA; skipping residue")
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

    Residues are given as e.g. ``['A45', 'B102', '1:45']`` (see
    ``_parse_residue_key``); a residue given twice (also as '45' and 'A45'
    for a blank chain) is built once and both keys count as built. Returns the flex file path, or None (and no file) if no
    residue could be converted. The file's REMARK lists the residues
    actually written. When ``built`` is given, the keys of those residues
    are appended to it, so only they are removed from the rigid receptor.
    """
    if not residues:
        return None

    flex_file = os.path.join(outdir, outputs.FLEX_PDBQT)
    blocks = []
    written_keys = []
    aliases = []  # keys naming a residue given earlier under another key
    seen = {}
    for key in residues:
        key = key.strip()
        identity = residue_identity(key)
        if identity is not None and identity in seen:
            logger.warning(f"[!] Residue '{key}' given more than once (as "
                           f"'{seen[identity]}'); using it once")
            if key != seen[identity]:
                aliases.append((key, seen[identity]))
            continue
        if identity is not None:
            seen[identity] = key
        atoms = _residue_atoms_from_pdbqt(receptor_pdbqt, key)
        if not atoms:
            logger.warning(
                f"[!] Residue '{key}' not found in receptor; skipping")
            continue
        block = _build_flex_residue_block(atoms)
        if not block:
            logger.warning(f"[!] Could not build flex block for '{key}'")
            continue
        blocks.extend(block)
        written_keys.append(key)

    written = len(written_keys)
    if written:
        with open(flex_file, "w", encoding="utf-8") as dst:
            dst.write("REMARK flexible residues: " + ", ".join(written_keys) + "\n")
            dst.writelines(blocks)
        if built is not None:
            built.extend(written_keys)
            built.extend(alias for alias, first in aliases if first in written_keys)

    if written == 0:
        logger.warning("[!] No flexible residues could be built; flex docking disabled")
        return None

    logger.info(
        f"[*] Flexible receptor PDBQT written: {flex_file} "
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
        identity = residue_identity(key)
        if identity is not None:
            keys.append(identity)

    flex_atoms = None
    if flex_file and os.path.exists(flex_file):
        flex_atoms = set()
        with open(flex_file, encoding="utf-8") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    flex_atoms.add((line[6:11].strip(), _record_chain(line[21]),
                                    line[22:27].strip()))

    removed = 0
    kept = []
    with open(receptor_pdbqt, "r", encoding="utf-8") as src:
        for line in src:
            if line.startswith(("ATOM", "HETATM")):
                chain = _record_chain(line[21])
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

    with open(outpath, "w", encoding="utf-8") as dst:
        dst.writelines(kept)

    if removed:
        logger.info(
            f"[*] Rigid receptor: removed {removed} flexible-residue atoms "
            f"({outpath})")
    return outpath
