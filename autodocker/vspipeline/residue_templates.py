#!/usr/bin/env python3
"""Ring hydrogens of Trp and His from residue templates.

Why
---
Open Babel (3.1.1 and 3.2.1) often fails to kekulize the indole and
imidazole rings of a protein PDB ("Failed to kekulize aromatic bonds"), and
its pH model then adds no ring hydrogen: on the 39 CASF-2016 + 10 DUD-E
receptors every Trp NE1 (280/280) was written as acceptor ``NA`` without a
hydrogen, and about half of the His (196/396 named HIS) had no ring hydrogen
at all, with ND1 typed ``NA`` (or ``N``, neither donor nor acceptor for
Vina, e.g. 1HSG His 69). AutoDock Vina makes an N a donor only when an
``HD`` atom lies within its bond distance (``model::assign_bonds``,
1.1 x covalent radii; ``bonded_to_HD``), so these residues lost their donors.

What
----
``apply_residue_templates`` edits the receptor PDBQT that Open Babel wrote
(after ``obabel -p``, the hydrogen stabilisation and the residue-label
restoration, before the frame translation):

* every Open Babel hydrogen bonded to Trp NE1 or His ND1/NE2 is removed (the
  template decides them all, so the result does not depend on whether Open
  Babel perceived the ring);
* **Trp**: NE1 always gets HE1;
* **His**: the tautomer is chosen by the first rule that applies
  (recorded per residue as ``decided_by``):

  1. ``residue_name``: HID/HSD -> HID (H on ND1), HIE/HSE -> HIE (H on NE2),
     HIP/HSP -> HIP (both); the residue is written as HIS (atom names are
     already the standard ones);
  2. ``input_hydrogens``: the input PDB (first model, any altloc) has a
     hydrogen within Vina's N-H bond distance (1.232 A) of ND1 and/or NE2:
     ND1 -> HID, NE2 -> HIE, both -> HIP. The pipeline removes input
     hydrogens before protonation; they are read from the input file here;
  3. ``metal_coordination``: exactly one ring N lies within
     ``METAL_COORDINATION_DISTANCE`` of a transition-metal/zinc HETATM of the
     input (whether or not ``--keep-hetero`` keeps it): that N stays
     unprotonated and the other one gets the hydrogen;
  4. ``default``: HIE (H on NE2), the more common neutral tautomer in
     solution and the default of Amber/CHARMM preparation tools.

  The choice does not depend on ``--ph``: no His is charged (HIP) unless the
  input says so by name or hydrogens. No hydrogen-bond-network optimisation
  (Reduce/PDB2PQR style flips) is done; that would need a search whose
  outcome depends on many neighbours.
* each hydrogen is placed in the ring plane on the external bisector of the
  two ring bonds of its N (CD1-NE1-CE2, CG-ND1-CE1, CD2-NE2-CE1) at
  ``N_H_LENGTH`` (0.9845 A: Open Babel's sp2 N-H length, 0.95 x 0.71 + 0.31,
  the length of every ND1-H Open Babel itself adds), named HE1/HD1/HE2, typed
  ``HD`` with charge +0.000, and written right after its N. Open Babel's
  Gasteiger charges are computed before ``-p`` adds hydrogens, so every polar
  hydrogen it writes has charge +0.000 and its N the charge of the
  unprotonated atom; the template hydrogens follow that convention (Vina does
  not read charges);
* parent N types: protonated ring N -> ``N`` (Vina N_D), unprotonated His N
  -> ``NA`` (N_A). No other atom's type, charge or coordinates change.
* atom serials are renumbered consecutively within each block of ATOM/HETATM
  records (Open Babel starts every block at 1).

Everything is a function of the coordinates and residue labels, so the
output is byte-identical run to run. The returned summary goes into
``run_manifest.json`` (``receptor_residue_templates``).
"""
import logging
import math
from typing import Dict, List, Optional, Set, Tuple

from vspipeline.hydrogens import PDBQT_HYDROGEN_TYPES, SpatialHash, pdbqt_element, vina_bond_cutoff

logger = logging.getLogger(__name__)

# Residue names of histidine and the tautomer they fix (None: not fixed).
HIS_RESIDUE_NAMES = {
    "HIS": None,
    "HID": "HID", "HSD": "HID",
    "HIE": "HIE", "HSE": "HIE",
    "HIP": "HIP", "HSP": "HIP",
}
TRP_RESIDUE_NAME = "TRP"
DEFAULT_HIS_STATE = "HIE"

# Ring nitrogens that carry the hydrogen in each His state.
HIS_PROTONATED = {"HID": ("ND1",), "HIE": ("NE2",), "HIP": ("ND1", "NE2")}
# N -> (its two ring neighbours, hydrogen name)
RING_N = {
    ("TRP", "NE1"): (("CD1", "CE2"), "HE1"),
    ("HIS", "ND1"): (("CG", "CE1"), "HD1"),
    ("HIS", "NE2"): (("CD2", "CE1"), "HE2"),
}
# Ring atoms a residue needs for its template.
TEMPLATE_ATOMS = {
    template: frozenset(n for (t, n) in RING_N if t == template)
    | frozenset(a for (t, _), (nbrs, _) in RING_N.items() if t == template for a in nbrs)
    for template in ("TRP", "HIS")
}

# Open Babel's N-H bond length for an sp2 nitrogen (OBAtom CorrectedBondRad:
# 0.95 x 0.71 + 0.31); the length of the ND1-H Open Babel adds itself.
N_H_LENGTH = 0.9845
# Charge written on the template hydrogens: Open Babel's Gasteiger charges are
# computed before -p adds hydrogens, so all its polar hydrogens carry 0.000.
TEMPLATE_H_CHARGE = 0.0

# Rule 3: metals that coordinate His ring nitrogens, and the N-metal distance
# (A) taken as coordination (Zn/Fe/Cu/Ni/Co-N(His) is 1.9-2.4 A).
COORDINATING_METALS = frozenset({"ZN", "FE", "CU", "NI", "CO", "MN", "CD", "HG"})
METAL_COORDINATION_DISTANCE = 2.8

# Chain ID of records with a blank chain column (protein.BLANK_CHAIN_ID).
_BLANK_CHAIN = "A"

POLICY = ("Trp NE1 always protonated (HE1). His tautomer by precedence: residue name "
          "(HID/HSD, HIE/HSE, HIP/HSP) > input hydrogens on ND1/NE2 > metal "
          "coordination of one ring N (that N unprotonated) > default HIE; "
          "pH-independent. Ring hydrogens in the ring plane on the external "
          f"bisector, N-H {N_H_LENGTH} A, type HD, charge +0.000.")

Xyz = Tuple[float, float, float]


def _xyz(line: str) -> Xyz:
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def _element(line: str) -> str:
    """Element of a PDB record: columns 77-78, else the first letter of the
    atom name (hydrogens and the ring atoms read here)."""
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element:
        return element
    name = line[12:16].strip().lstrip("0123456789")
    return name[:1].upper()


def _residue_label(key: str) -> str:
    """'A69', 'B12A' from columns 22-27 (chain, resnum, insertion code)."""
    chain = key[0].strip() or _BLANK_CHAIN
    return f"{chain}{key[1:5].strip()}{key[5:].strip()}"


def _unit(v):
    length = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    return (v[0] / length, v[1] / length, v[2] / length)


def ring_hydrogen_position(n: Xyz, a: Xyz, b: Xyz, length: float = N_H_LENGTH) -> Xyz:
    """Hydrogen on ring atom ``n`` bonded to ``a`` and ``b``: in their plane,
    on the external bisector of a-n-b, at ``length``."""
    ua = _unit((n[0] - a[0], n[1] - a[1], n[2] - a[2]))
    ub = _unit((n[0] - b[0], n[1] - b[1], n[2] - b[2]))
    d = _unit((ua[0] + ub[0], ua[1] + ub[1], ua[2] + ub[2]))
    return (n[0] + length * d[0], n[1] + length * d[1], n[2] + length * d[2])


# ---------------------------------------------------------------------------
# Input evidence (rules 2 and 3)
# ---------------------------------------------------------------------------

def read_input_evidence(input_pdb: str) -> Tuple[Dict[str, Set[str]], List[Xyz]]:
    """From the first model of ``input_pdb``: {residue key (columns 22-27) of
    a His: ring N names with an input hydrogen bonded} and the coordinates of
    the coordinating-metal HETATM atoms."""
    ring_n: Dict[Tuple[str, str], List[Xyz]] = {}
    hydrogens: Dict[Tuple[str, str], List[Xyz]] = {}
    metals: List[Xyz] = []
    with open(input_pdb, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
                continue
            resname = line[17:20].strip()
            element = _element(line)
            if line.startswith("HETATM"):
                # Single-atom ions: element column, else the residue name (ZN).
                if (line[76:78].strip().upper() or resname.upper()) in COORDINATING_METALS:
                    metals.append(_xyz(line))
                continue
            if resname not in HIS_RESIDUE_NAMES:
                continue
            key = (line[21:27], line[16])  # residue, altloc
            if element in ("H", "D"):
                hydrogens.setdefault(key, []).append(_xyz(line))
            elif line[12:16].strip() in ("ND1", "NE2"):
                ring_n.setdefault(key, []).append((line[12:16].strip(), _xyz(line)))
    cutoff = vina_bond_cutoff("N")
    found: Dict[str, Set[str]] = {}
    for (residue, altloc), atoms in ring_n.items():
        # An H with a blank altloc applies to every altloc of the N, and vice versa.
        hs = list(hydrogens.get((residue, altloc), ()))
        if altloc != " ":
            hs += hydrogens.get((residue, " "), [])
        else:
            hs += [h for (r, a), xs in hydrogens.items() if r == residue and a != " " for h in xs]
        for name, xyz in atoms:
            if any(math.dist(xyz, h) < cutoff for h in hs):
                found.setdefault(residue, set()).add(name)
    return found, metals


def choose_his_state(resname: str, input_h: Optional[Set[str]],
                     coordinated: Set[str]) -> Tuple[str, str]:
    """(state, decided_by) for one His; see the module docstring."""
    named = HIS_RESIDUE_NAMES.get(resname)
    if named:
        return named, "residue_name"
    if input_h:
        if input_h >= {"ND1", "NE2"}:
            return "HIP", "input_hydrogens"
        return ("HID" if "ND1" in input_h else "HIE"), "input_hydrogens"
    if len(coordinated) == 1:
        return ("HIE" if "ND1" in coordinated else "HID"), "metal_coordination"
    return DEFAULT_HIS_STATE, "default"


# ---------------------------------------------------------------------------
# PDBQT edit
# ---------------------------------------------------------------------------

def _is_atom(line: str) -> bool:
    return line.startswith(("ATOM", "HETATM"))


def _atom_type(line: str) -> str:
    return line[77:79].strip()


def _retype(line: str, atom_type: str) -> str:
    return line[:77] + f"{atom_type:<2s}" + line[79:]


def _hydrogen_line(parent: str, name: str, xyz: Xyz) -> str:
    """PDBQT record of a template hydrogen, labelled like its parent N."""
    return (f"{parent[:12]}{name:>4s}{parent[16:30]}"
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{parent[54:70]}"
            f"{TEMPLATE_H_CHARGE:+6.3f} HD\n")


def _renumber(lines: List[str]) -> None:
    """Consecutive serials within each run of ATOM/HETATM records, starting
    from the run's first serial (Open Babel starts every block at 1)."""
    serial = None
    for i, line in enumerate(lines):
        if not _is_atom(line):
            serial = None
            continue
        if serial is None:
            try:
                serial = int(line[6:11])
            except ValueError:
                serial = 1
        lines[i] = f"{line[:6]}{serial:5d}{line[11:]}"
        serial += 1


def apply_residue_templates(pdbqt_file: str, input_pdb: str) -> Dict:
    """Give every Trp and His of ``pdbqt_file`` its template ring hydrogens
    (see the module docstring); ``input_pdb`` is the receptor as given by the
    user (rules 2 and 3). Rewrites the file and returns the summary recorded
    in run_manifest.json."""
    with open(pdbqt_file, encoding="utf-8") as f:
        lines = f.readlines()
    input_h, metals = read_input_evidence(input_pdb)

    # Residues: (resname, columns 22-27) -> {atom name: line index}, heavy atoms only.
    residues: Dict[Tuple[str, str], Dict[str, int]] = {}
    order: List[Tuple[str, str]] = []
    # Cell edge >= the largest Vina H bond distance (H-metal, 2.33 A).
    heavy = SpatialHash(2.5)
    hydrogen_idx = []
    for i, line in enumerate(lines):
        if not _is_atom(line):
            continue
        if _atom_type(line) in PDBQT_HYDROGEN_TYPES:
            hydrogen_idx.append(i)
            continue
        heavy.add(_xyz(line), i)
        resname = line[17:20].strip()
        if resname in HIS_RESIDUE_NAMES or resname == TRP_RESIDUE_NAME:
            key = (resname, line[21:27])
            if key not in residues:
                residues[key] = {}
                order.append(key)
            residues[key].setdefault(line[12:16].strip(), i)

    metal_hash = SpatialHash(METAL_COORDINATION_DISTANCE)
    for m in metals:
        metal_hash.add(m, None)

    protonate: Dict[int, str] = {}      # N line index -> hydrogen name
    retype: Dict[int, str] = {}         # N line index -> new type
    ring_n_idx: Set[int] = set()
    his_records = []
    skipped = []
    trp_added = 0
    for key in order:
        resname, label = key
        atoms = residues[key]
        template = "TRP" if resname == TRP_RESIDUE_NAME else "HIS"
        if not TEMPLATE_ATOMS[template] <= atoms.keys():
            skipped.append(f"{resname} {_residue_label(label)}")
            continue
        if template == "TRP":
            protonate[atoms["NE1"]] = "HE1"
            retype[atoms["NE1"]] = "N"
            ring_n_idx.add(atoms["NE1"])
            trp_added += 1
            continue
        coordinated = {n for n in ("ND1", "NE2")
                       if any(math.dist(_xyz(lines[atoms[n]]), m) < METAL_COORDINATION_DISTANCE
                              for m, _ in metal_hash.near(_xyz(lines[atoms[n]])))}
        state, rule = choose_his_state(resname, input_h.get(label), coordinated)
        for n in ("ND1", "NE2"):
            idx = atoms[n]
            ring_n_idx.add(idx)
            if n in HIS_PROTONATED[state]:
                protonate[idx] = RING_N[("HIS", n)][1]
                retype[idx] = "N"
            else:
                retype[idx] = "NA"
        his_records.append({"residue": _residue_label(label), "input_name": resname,
                            "state": state, "decided_by": rule})

    # Open Babel hydrogens on the template nitrogens: those whose nearest heavy
    # atom is a template ring N within Vina's bond distance.
    cutoff = vina_bond_cutoff("N")
    removed: Set[int] = set()
    for i in hydrogen_idx:
        xyz = _xyz(lines[i])
        near = [(math.dist(xyz, p), j) for p, j in heavy.near(xyz)]
        if not near:
            continue
        d, j = min(near)
        if j in ring_n_idx and d < cutoff:
            removed.add(i)

    # New hydrogens, and a check that Vina bonds each to its N only.
    added: Dict[int, str] = {}
    close_contacts = []
    for idx, h_name in protonate.items():
        line = lines[idx]
        resname = line[17:20].strip()
        template = "TRP" if resname == TRP_RESIDUE_NAME else "HIS"
        atoms = residues[(resname, line[21:27])]
        nbrs = RING_N[(template, line[12:16].strip())][0]
        xyz = ring_hydrogen_position(_xyz(line), _xyz(lines[atoms[nbrs[0]]]),
                                     _xyz(lines[atoms[nbrs[1]]]))
        h_line = _hydrogen_line(line, h_name, xyz)
        bonded = [j for p, j in heavy.near(_xyz(h_line))
                  if math.dist(_xyz(h_line), p)
                  < vina_bond_cutoff(pdbqt_element(_atom_type(lines[j])))]
        if bonded != [idx]:
            close_contacts.append(f"{resname} {_residue_label(line[21:27])} {h_name}")
        added[idx] = h_line

    renamed = sum(1 for resname, _ in order if resname in HIS_RESIDUE_NAMES and resname != "HIS")
    out = []
    for i, line in enumerate(lines):
        if i in removed:
            continue
        if i in added:
            out.append(_retype(line, retype[i]))
            out.append(added[i])
            continue
        if i in retype:
            line = _retype(line, retype[i])
        out.append(line)
    for i, line in enumerate(out):
        if _is_atom(line) and line[17:20] != "HIS" and line[17:20].strip() in HIS_RESIDUE_NAMES:
            out[i] = line[:17] + "HIS" + line[20:]
    _renumber(out)
    with open(pdbqt_file, "w", encoding="utf-8") as f:
        f.writelines(out)

    states = {s: sum(r["state"] == s for r in his_records) for s in ("HID", "HIE", "HIP")}
    decided = {r: sum(h["decided_by"] == r for h in his_records)
               for r in ("residue_name", "input_hydrogens", "metal_coordination", "default")}
    summary = {
        "policy": POLICY,
        "trp_he1_added": trp_added,
        "his_states": states,
        "his_decided_by": decided,
        "his_residues": his_records,
        "his_residue_names_written_as_HIS": renamed,
        "openbabel_ring_hydrogens_replaced": len(removed),
        "skipped_incomplete_residues": skipped,
        "hydrogens_within_bond_distance_of_another_atom": close_contacts,
    }
    his_text = "no His"
    if his_records:
        his_text = ("His " + ", ".join(f"{s} {n}" for s, n in states.items())
                    + " (decided by " + ", ".join(f"{r.replace('_', ' ')} {n}"
                                                  for r, n in decided.items() if n) + ")")
    logger.info(f"[*] Residue templates: HE1 on {trp_added} Trp; {his_text}; "
                f"replaced {len(removed)} Open Babel ring hydrogen(s)")
    if skipped:
        logger.warning(f"[!] Residue templates: {len(skipped)} Trp/His residue(s) lack "
                       "ring atoms and keep Open Babel's hydrogens: " + ", ".join(skipped))
    if close_contacts:
        logger.warning(
            f"[!] Residue templates: {len(close_contacts)} ring hydrogen(s) within Vina's "
            "bond distance of a second atom (clashing input geometry): "
            + ", ".join(close_contacts))
    return summary
