#!/usr/bin/env python3
"""PDB records: the reader shared by receptor cleaning and HETATM handling
(first model, one alternate location per residue) and element perception.
"""
import re
from collections import defaultdict
from typing import List, Sequence, Tuple

from vspipeline.hydrogens import PDBQT_HYDROGEN_TYPES

# Chain ID given to ATOM/HETATM records whose chain column is blank.
BLANK_CHAIN_ID = "A"

HYDROGEN_ELEMENTS = frozenset({"H", "D"})


def _element_from_atom_name(field: str) -> str:
    """Element from a PDB atom-name field (columns 13-16) when the element
    column is blank. PDB right-justifies the element symbol in columns 13-14,
    so C-alpha is ' CA ' (carbon) while calcium is 'CA  '; four-character
    hydrogen names such as 'HG21' start in column 13."""
    field = field.ljust(4)[:4]
    if field[0].isdigit():
        return field[1].upper()
    if field[0] == "H" and field.strip() == field and len(field.strip()) == 4:
        return "H"
    return re.sub(r"[^A-Za-z]", "", field[:2]).upper()


def pdb_record_element(line: str) -> str:
    """Element symbol (upper case) of a PDB ATOM/HETATM record.

    Uses the element column (77-78) and falls back to the atom name. In ATOM
    records (amino acids) a name starting with H or D, after any leading
    digit, is always a hydrogen ('HA', 'HG1', '1HD2', 'HNZ1'); the
    calcium/mercury ambiguity of the name only exists for HETATM records.
    """
    element = re.sub(r"[^A-Za-z]", "", line[76:78]).upper() if len(line) > 76 else ""
    if element:
        return element
    if line.startswith("ATOM"):
        name = line[12:16].strip().lstrip("0123456789")
        if name[:1].upper() in HYDROGEN_ELEMENTS:
            return name[0].upper()
    return _element_from_atom_name(line[12:16])


def is_hydrogen_record(line: str) -> bool:
    """True for a PDB ATOM/HETATM record of a hydrogen (or deuterium)."""
    return pdb_record_element(line) in HYDROGEN_ELEMENTS


def _is_pdbqt_hydrogen(line: str) -> bool:
    return line[77:79].strip() in PDBQT_HYDROGEN_TYPES


def _first_model_lines(pdb_file: str) -> List[str]:
    """Lines of a PDB file up to the end of its first model.

    NMR ensembles: keeping every model would superimpose copies of the
    receptor.
    """
    lines = []
    with open(pdb_file, encoding="utf-8") as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            lines.append(line)
    return lines


def _residue_position(line: str) -> Tuple[str, str]:
    """(chain, resnum + insertion code) of an ATOM/HETATM record. The residue
    name is left out on purpose: with microheterogeneity two residue types
    share one position and must not both be kept."""
    return line[21:22], line[22:27]


def _occupancy(line: str) -> float:
    try:
        return float(line[54:60])
    except ValueError:
        return 0.0


def select_altlocs(lines: Sequence[str]) -> Tuple[List[str], int]:
    """Keep one alternate location per residue.

    For each residue position the altloc label with the highest total
    occupancy wins (ties: the label seen first), and all of that residue's
    atoms use it, so a residue is never assembled from atoms of different
    conformers. Records without an altloc are always kept. The altloc column
    of the kept records is blanked. Returns (lines, dropped record count).
    """
    occupancy = defaultdict(lambda: defaultdict(float))
    first_seen = defaultdict(dict)
    for i, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")) and len(line) > 16 and line[16] != " ":
            pos = _residue_position(line)
            occupancy[pos][line[16]] += _occupancy(line)
            first_seen[pos].setdefault(line[16], i)
    chosen = {
        pos: max(labels, key=lambda lab: (labels[lab], -first_seen[pos][lab]))
        for pos, labels in occupancy.items()
    }
    kept, dropped = [], 0
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) > 16 and line[16] != " ":
            if line[16] != chosen[_residue_position(line)]:
                dropped += 1
                continue
            line = line[:16] + " " + line[17:]
        kept.append(line)
    return kept, dropped


def read_pdb_records(pdb_file: str) -> Tuple[List[str], int]:
    """First-model lines of ``pdb_file`` with one altloc per residue.

    The single reader for receptor cleaning and for HETATM detection and
    appending, so all of them see the same atoms. Returns (lines, dropped
    altloc record count).
    """
    return select_altlocs(_first_model_lines(pdb_file))


def record_chain(line: str) -> str:
    """Chain ID of a record, with a blank chain read as ``BLANK_CHAIN_ID``."""
    return (line[21:22].strip() or BLANK_CHAIN_ID)


def _has_blank_chain(pdb_file: str) -> bool:
    """True if an ATOM record of ``pdb_file`` has a blank chain ID."""
    with open(pdb_file, errors="replace", encoding="utf-8") as f:
        return any(line.startswith("ATOM") and len(line) > 21 and not line[21].strip()
                   for line in f)


def _heavy_atom_coords(path: str) -> Tuple[List[float], List[float], List[float]]:
    """Coordinates of the non-hydrogen ATOM/HETATM records of a PDB/PDBQT."""
    pdbqt = path.lower().endswith(".pdbqt")
    xs, ys, zs = [], [], []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if _is_pdbqt_hydrogen(line) if pdbqt else is_hydrogen_record(line):
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError as e:
                raise ValueError(f"Malformed coordinates in {path}: {line.rstrip()}") from e
            xs.append(x)
            ys.append(y)
            zs.append(z)
    if not xs:
        raise ValueError(f"No heavy-atom ATOM/HETATM records in {path}")
    return xs, ys, zs
