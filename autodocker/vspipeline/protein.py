#!/usr/bin/env python3
"""Receptor preparation: PDB cleaning, protonation, pocket ranking, grid box.

Cross-module references go through the ``runner`` namespace at call time so
the public entry point and its monkeypatch surface stay unchanged.

The module is organised as small helpers grouped by concern, used by the
``ProteinPreparation`` class:

* PDB record reading (first model, one alternate location per residue,
  element perception),
* protonation (Open Babel) and deterministic hydrogen placement,
* fpocket pockets (parsing and ranking),
* box geometry and the grid/PyMOL writers.
"""
import csv
import logging
import math
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# The runner module owns the tool-path globals (OBABEL, COMMAND_TIMEOUT) and
# every public name; during the initial import of runner this binds the
# partially-initialized module object, which is fine because attributes are
# only accessed at call time.
from vspipeline import defaults  # noqa: E402
from vspipeline import outputs  # noqa: E402
import runner  # noqa: E402
from vspipeline.hydrogens import (  # noqa: E402,F401  (some re-exported)
    H_PARENT_MAX_DISTANCE,
    OB_BOND_TOLERANCE,
    OB_COVALENT_RADII,
    OB_MIN_BOND_DISTANCE,
    OB_UNKNOWN_COVALENT_RADIUS,
    PDBQT_HYDROGEN_TYPES,
    VINA_BOND_ALLOWANCE,
    VINA_COVALENT_RADII,
    VINA_METAL_COVALENT_RADIUS,
    SpatialHash as _SpatialHash,
    pdbqt_element as _pdbqt_element,
    same_bytes,
    stabilize_hydrogens,
    vina_bond_cutoff as _vina_bond_cutoff,
)
from vspipeline.utils import elements_without_gasteiger  # noqa: E402

logger = logging.getLogger(__name__)

# Modified amino acids deposited as HETATM; they are part of the protein and
# dropping them would leave holes in the backbone or side chains.
MODIFIED_RESIDUES = {
    "MSE", "SEP", "TPO", "PTR", "CSO", "CSD", "CME", "MLY", "HYP", "KCX",
    "PCA", "CSS", "OCS", "SCH", "NEP", "ALY", "M3L", "TYS",
}

# Selenomethionine is written as methionine (Se -> S): AutoDock Vina has no
# selenium atom type and Open Babel no Gasteiger parameters for Se, so a SeMet
# structure could not be converted at all. The side chains are isosteric
# (C-Se 1.95 A vs C-S 1.81 A); coordinates are kept as deposited.
SELENOMETHIONINE = "MSE"
_SEMET_TO_MET_ATOM = {"SE": " SD "}

# A receptor with fewer ATOM records than this is rejected as truncated.
MIN_PROTEIN_ATOM_LINES = 100
# Default padding (A) added on each side of the pocket atoms' bounding box.
DEFAULT_PADDING = defaults.PADDING
# Minimum edge (A) of the whole-receptor box used with --no-fpocket.
MIN_RECEPTOR_BOX_SIZE = 24.0
# Chain ID given to ATOM/HETATM records whose chain column is blank.
BLANK_CHAIN_ID = "A"
# fpocket's documented druggability threshold: "A druggability score at 0.5
# (the threshold) indicates that binding of prodrugs or druglike molecules can
# be possible" (fpocket doc/GETTINGSTARTED.md, "Druggability Score").
FPOCKET_DRUGGABILITY_THRESHOLD = 0.5

HYDROGEN_ELEMENTS = frozenset({"H", "D"})


# ---------------------------------------------------------------------------
# PDB records
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# PDBQT post-processing
# ---------------------------------------------------------------------------

def _translate_pdbqt(path: str, offset: Tuple[float, float, float]) -> None:
    """Subtract ``offset`` from every ATOM/HETATM coordinate of a PDBQT file."""
    dx, dy, dz = offset
    out = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                x = float(line[30:38]) - dx
                y = float(line[38:46]) - dy
                z = float(line[46:54]) - dz
                line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            out.append(line)
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(out)


def _strip_remark_name_paths(path: str) -> None:
    """Reduce Open Babel's ``REMARK  Name = <input path>`` to the file name,
    so the receptor PDBQT does not depend on the output directory."""
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    prefix = "REMARK  Name = "
    lines = [prefix + os.path.basename(l[len(prefix):].rstrip("\n")) + "\n"
             if l.startswith(prefix) else l for l in lines]
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(lines)


def _semet_to_met(line: str) -> Tuple[str, Tuple[str, str]]:
    """A selenomethionine (MSE) record as methionine: residue name MET, the
    selenium as ' SD ' with element S. Returns (line, residue position)."""
    name = line[12:16]
    element = line[76:78] if len(line) > 76 else ""
    if name.strip() in _SEMET_TO_MET_ATOM:
        name = _SEMET_TO_MET_ATOM[name.strip()]
        element = " S"
    line = line[:12] + name + line[16:17] + "MET" + line[20:]
    if len(line) > 76:
        line = line[:76] + element + line[78:]
    return line, _residue_position(line)


def _write_without_residues(src: str, dst: str, residues: set) -> None:
    """Copy a PDB without the HETATM residues whose ``pocket._residue_key``
    is in ``residues``."""
    from vspipeline import pocket
    with open(src, encoding="utf-8") as f:
        lines = f.readlines()
    kept = []
    for line in lines:
        parsed = pocket._parse_hetatm(line)
        if parsed is not None and pocket._residue_key(parsed) in residues:
            continue
        kept.append(line)
    with open(dst, "w", encoding="utf-8") as f:
        f.writelines(kept)


def _unsupported_elements_hint(pdb_file: str) -> str:
    """' (elements without Open Babel Gasteiger parameters: ...)' for the
    records of ``pdb_file``, or '' when every element is supported."""
    found = defaultdict(set)
    try:
        with open(pdb_file, encoding="utf-8") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    found[pdb_record_element(line)].add(
                        f"{line[17:20].strip()} {record_chain(line)}{line[22:27].strip()}")
    except OSError:
        return ""
    bad = elements_without_gasteiger(found)
    if not bad:
        return ""
    where = "; ".join(f"{e} in " + ", ".join(sorted(found[e])[:3])
                      + (" ..." if len(found[e]) > 3 else "") for e in bad)
    return (f" (element(s) without Open Babel Gasteiger parameters: {where}; remove "
            "those residues or run without --keep-hetero)")


def _xyz_key(line: str) -> Tuple[float, float, float]:
    return (round(float(line[30:38]), 3), round(float(line[38:46]), 3),
            round(float(line[46:54]), 3))


def _restore_residue_labels(pdbqt_file: str, pdb_file: str) -> int:
    """Copy atom names and residue labels from ``pdb_file`` back onto the
    PDBQT that Open Babel wrote from it.

    ``obabel -p`` re-perceives residues, and on a receptor with a chain break
    it renames, renumbers and re-chains most atoms. Heavy-atom coordinates
    are untouched, so each heavy atom is matched to its input record by
    coordinates; each added hydrogen takes the residue labels of the nearest
    heavy atom within ``H_PARENT_MAX_DISTANCE`` (the atom it is bonded to).
    Raises if the heavy-atom counts differ, a heavy atom has no match or a
    hydrogen has no heavy atom nearby. Returns the number of relabelled atoms.
    """
    labels = {}
    duplicates = 0
    n_input_heavy = 0
    with open(pdb_file, encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                if is_hydrogen_record(line):
                    continue
                n_input_heavy += 1
                # columns 13-27: atom name, altloc, resname, chain, resnum, icode
                label = (line[:6], line[12:27])
                previous = labels.setdefault(_xyz_key(line), label)
                duplicates += previous != label
    if duplicates:
        logger.warning(
            f"[!] {duplicates} atom(s) in {pdb_file} share coordinates with an "
            "earlier atom; their PDBQT records take the earlier atom's labels")

    with open(pdbqt_file, encoding="utf-8") as f:
        lines = f.readlines()
    atoms = [(i, l) for i, l in enumerate(lines) if l.startswith(("ATOM", "HETATM"))]
    n_output_heavy = sum(1 for _, l in atoms if not _is_pdbqt_hydrogen(l))
    if n_output_heavy != n_input_heavy:
        raise RuntimeError(
            f"{pdbqt_file} has {n_output_heavy} heavy atoms but {pdb_file} has "
            f"{n_input_heavy}; Open Babel added or dropped atoms")

    heavy = _SpatialHash(H_PARENT_MAX_DISTANCE)
    changed = 0
    for i, line in atoms:
        if _is_pdbqt_hydrogen(line):
            continue
        found = labels.get(_xyz_key(line))
        if found is None:
            raise RuntimeError(
                f"Receptor atom at ({line[30:54].strip()}) in {pdbqt_file} has no "
                f"matching record in {pdb_file}; cannot restore residue labels")
        record, label = found
        new = record + line[6:12] + label + line[27:]
        changed += new != line
        lines[i] = new
        heavy.add(_xyz_key(line), found)
    for i, line in atoms:
        if not _is_pdbqt_hydrogen(line):
            continue
        xyz = _xyz_key(line)
        nearby = [(math.dist(p, xyz), parent) for p, parent in heavy.near(xyz)]
        nearby = [entry for entry in nearby if entry[0] <= H_PARENT_MAX_DISTANCE]
        if not nearby:
            raise RuntimeError(
                f"Hydrogen at ({line[30:54].strip()}) has no heavy atom within "
                f"{H_PARENT_MAX_DISTANCE} A in {pdbqt_file}")
        record, label = min(nearby, key=lambda entry: entry[0])[1]
        # Hydrogen keeps its own name ("H"); record type and residue fields
        # come from its heavy atom.
        new = record + line[6:17] + label[5:] + line[27:]
        changed += new != line
        lines[i] = new
    with open(pdbqt_file, "w", encoding="utf-8") as f:
        f.writelines(lines)
    return changed


def _count_polar_atoms(pdbqt_file: str) -> Tuple[int, int]:
    """(polar hydrogens, N/O heavy atoms) in a PDBQT file."""
    hd = no = 0
    with open(pdbqt_file, encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                t = line[77:79].strip()
                hd += t == "HD"
                no += t in ("N", "NA", "NS", "OA", "OS", "O")
    return hd, no


# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------

def _log_tool_stderr(tool: str, stderr: str) -> None:
    """Log a tool's stderr at WARNING, collapsing repeated messages.

    Open Babel prints one block per atom for some PDB warnings (e.g. a blank
    element column), separated by '=====' lines; identical blocks (ignoring
    the atom number) are logged once with a count.
    """
    if not stderr or not stderr.strip():
        return
    blocks = [b.strip() for b in re.split(r"^=+\s*$", stderr, flags=re.M) if b.strip()]
    counts = {}
    for block in blocks:
        key = re.sub(r"\(atom \d+\)", "(atom N)", " ".join(block.split()))
        counts[key] = counts.get(key, 0) + 1
    for message, count in counts.items():
        if re.fullmatch(r"\d+ molecules? converted", message):
            logger.debug(f"{tool}: {message}")
            continue
        suffix = f" [x{count}]" if count > 1 else ""
        logger.warning(f"[!] {tool}: {message}{suffix}")


def _run_tool(cmd: List[str], tool: str, log_stderr: bool = True) -> str:
    """Run an external tool through ``runner.run``; its stderr goes to the
    log (unless ``log_stderr`` is False). Returns stdout. Raises RuntimeError
    with the tool's stderr on failure."""
    try:
        result = runner.run(cmd, capture=True)
    except subprocess.CalledProcessError as e:
        stderr = (e.stderr or "").strip()
        raise RuntimeError(
            f"{tool} failed (exit code {e.returncode}): {' '.join(cmd)}"
            + (f"\n{stderr[-2000:]}" if stderr else "")) from e
    stdout, stderr = (result or ("", ""))[:2]
    if log_stderr:
        _log_tool_stderr(tool, stderr)
    return stdout or ""


def _run_obabel(args: List[str], outputs: Sequence[str] = (), log_stderr: bool = True) -> str:
    """Open Babel through ``runner.run_obabel`` (single-threaded, partial
    outputs removed on failure); its stderr goes to the log unless
    ``log_stderr`` is False. Returns stdout."""
    stdout, stderr = runner.run_obabel(list(args), outputs=outputs, return_stderr=True)
    if log_stderr:
        _log_tool_stderr("obabel", stderr)
    return stdout


def _fpocket_version(fpocket: str) -> Optional[str]:
    """fpocket's version from its usage text (it has no --version flag)."""
    try:
        proc = subprocess.run([fpocket, "-h"], capture_output=True, text=True,
                              timeout=runner.COMMAND_TIMEOUT)
    except (OSError, subprocess.SubprocessError):
        return None
    match = re.search(r"fpocket\s+(\d[\w.]*)", re.sub(r"\x1b\[[0-9;]*m", "",
                                                     proc.stdout + proc.stderr))
    return match.group(1) if match else None


# ---------------------------------------------------------------------------
# Pockets and box geometry
# ---------------------------------------------------------------------------

def bounding_box(xs: Sequence[float], ys: Sequence[float], zs: Sequence[float],
                 padding: float, min_size: float = 0.0
                 ) -> Tuple[float, float, float, float, float, float]:
    """Box centred on the midpoint of the coordinates' bounding box, with
    ``padding`` added on each side (so the padding is symmetric) and every
    edge at least ``min_size``."""
    cx = (min(xs) + max(xs)) / 2.0
    cy = (min(ys) + max(ys)) / 2.0
    cz = (min(zs) + max(zs)) / 2.0
    sx = max(max(xs) - min(xs) + 2 * padding, min_size)
    sy = max(max(ys) - min(ys) + 2 * padding, min_size)
    sz = max(max(zs) - min(zs) + 2 * padding, min_size)
    return cx, cy, cz, sx, sy, sz


def is_druggable(pocket: Dict) -> bool:
    """fpocket druggability score at or above the documented 0.5 threshold."""
    return pocket.get("druggability_score", 0) >= FPOCKET_DRUGGABILITY_THRESHOLD


def rank_pockets(pockets: List[Dict]) -> List[Dict]:
    """Pockets best first: fpocket score, then druggability score, then the
    lower pocket number. fpocket's volume is a Monte Carlo estimate seeded
    from the clock, so it is not used."""
    return sorted(pockets, key=lambda p: (-p.get("score", -1),
                                          -p.get("druggability_score", -1),
                                          p["number"]))


class ProteinPreparation:
    """Receptor preparation: chain selection, cleaning, protonation, pocket
    detection and the docking box."""

    def __init__(self, pdb_file: str, outdir: str):
        self.pdb_file = pdb_file
        self.outdir = outdir
        self.pdb_clean = os.path.join(outdir, outputs.PROTEIN_CLEAN_PDB)
        self.receptor_pdbqt = os.path.join(outdir, outputs.RECEPTOR_PDBQT)
        self.grid_conf = os.path.join(outdir, outputs.GRID_CONF)
        self.grid_box_script = os.path.join(outdir, outputs.GRID_BOX_PML)
        self.pocket_summary_file = os.path.join(outdir, outputs.POCKET_SUMMARY_CSV)

    def validate(self):
        """Validate protein PDB file."""
        if not os.path.exists(self.pdb_file):
            raise FileNotFoundError(f"Protein file not found: {self.pdb_file}")

        if not os.access(self.pdb_file, os.R_OK):
            raise PermissionError(f"Cannot read protein file: {self.pdb_file}")

        try:
            with open(self.pdb_file, encoding="utf-8") as f:
                count = sum(1 for line in f if line.startswith("ATOM"))
        except IOError as e:
            raise IOError(f"Error reading protein file: {e}")

        if count < MIN_PROTEIN_ATOM_LINES:
            raise ValueError(
                f"Invalid protein: only {count} ATOM lines "
                f"(need >= {MIN_PROTEIN_ATOM_LINES})")

        logger.info(f"[*] Protein validated: {count} ATOM lines")

    def select_chain(self) -> str:
        """Choose a chain for docking when --chain is not given.

        Non-interactive: with several protein chains the first one (in
        sorted order) is used and a warning lists the alternatives.
        """
        chains = runner.get_chains(self.pdb_file)

        if len(chains) == 1:
            note = ""
            if chains[0] == BLANK_CHAIN_ID and _has_blank_chain(self.pdb_file):
                note = (f" (blank chain ID in the PDB, read as {BLANK_CHAIN_ID}; "
                        f"--chain {BLANK_CHAIN_ID} and residue keys such as "
                        f"{BLANK_CHAIN_ID}45 or 45 address it)")
            logger.info(f"[*] Single chain detected: {chains[0]}{note}")
            return chains[0]

        logger.warning(
            f"[!] WARNING: {len(chains)} protein chains present "
            f"({', '.join(chains)}); docking against chain '{chains[0]}' only. "
            "Pass --chain (e.g. --chain all or --chain A,B) to choose explicitly.")
        return chains[0]

    def _write_clean_pdb(self, selected_chains: Sequence[str], keep_hetero: bool) -> None:
        """Write ``pdb_clean``: first model, one altloc per residue, the
        selected chains, no hydrogens, HETATM records only with
        ``keep_hetero`` (modified amino acids are always kept as ATOM;
        selenomethionine is written as methionine, see SELENOMETHIONINE)."""
        try:
            lines, skipped_altloc = read_pdb_records(self.pdb_file)
        except IOError as e:
            raise IOError(f"Failed to read {self.pdb_file}: {e}")

        kept_atoms = skipped_hetero = removed_h = 0
        semet = set()
        other_chains = defaultdict(lambda: [0, 0])  # chain -> [ATOM, HETATM] records
        last_was_atom = False
        out = []
        for line in lines:
            is_atom = line.startswith("ATOM")
            is_het = line.startswith("HETATM")
            if is_het and line[17:20].strip() in MODIFIED_RESIDUES:
                line = "ATOM  " + line[6:]
                is_atom, is_het = True, False
            if is_atom or (keep_hetero and is_het):
                if len(line) <= 21 or record_chain(line) not in selected_chains:
                    if len(line) > 21:
                        other_chains[record_chain(line)][0 if is_atom else 1] += 1
                    continue
                # Input hydrogens are dropped: Open Babel's pH model would
                # protonate the already protonated groups once more and put
                # the extra hydrogen in a random direction.
                if is_hydrogen_record(line):
                    removed_h += 1
                    continue
                if is_atom and line[17:20] == SELENOMETHIONINE:
                    line, residue = _semet_to_met(line)
                    semet.add(residue)
                out.append(line)
                kept_atoms += 1
                last_was_atom = True
            elif is_het:
                skipped_hetero += 1
            elif line.startswith("TER") and last_was_atom:
                out.append(line)
                last_was_atom = False
        out.append("END\n")

        if kept_atoms == 0:
            raise ValueError(
                f"Selected chain(s) '{', '.join(selected_chains)}' have no atoms "
                f"in {self.pdb_file}")
        try:
            with open(self.pdb_clean, "w", encoding="utf-8") as dst:
                dst.writelines(out)
        except IOError as e:
            raise IOError(
                f"Failed to filter chain(s) {', '.join(selected_chains)}: {e}")

        if skipped_altloc:
            logger.info(f"[*] Kept one alternate location per residue (highest "
                        f"occupancy); dropped {skipped_altloc} alternate-location "
                        "atom records")
        if removed_h:
            logger.warning(
                f"[!] Removed {removed_h} hydrogen records from the input; the receptor "
                "is re-protonated by Open Babel at the requested pH (input "
                "protonation states are not kept)")
        if semet:
            logger.warning(
                f"[!] Wrote {len(semet)} selenomethionine (MSE) residue(s) as "
                "methionine (MET, Se -> S): AutoDock Vina has no selenium atom type "
                "and Open Babel no Gasteiger parameters for Se; coordinates are kept")
        if skipped_hetero and not keep_hetero:
            logger.info(
                f"[*] Removed {skipped_hetero} HETATM records during receptor cleaning")
        het_only = sorted(c for c, (n_atom, _) in other_chains.items() if n_atom == 0)
        if keep_hetero and het_only:
            logger.warning(
                "[!] --keep-hetero: not keeping HETATM records of chain(s) without "
                "protein atoms, which --chain never selects: " + ", ".join(
                    f"{c} ({other_chains[c][1]} records)" for c in het_only))
        dropped_het = sum(n_het for c, (n_atom, n_het) in other_chains.items() if n_atom)
        if keep_hetero and dropped_het:
            logger.info(f"[*] --keep-hetero: {dropped_het} HETATM records of chains "
                        "outside the selection are not kept")

    def _metal_ions(self) -> List[Tuple[str, str, str]]:
        """Residue keys (see ``pocket._residue_key``) of the single-atom metal
        ions in ``pdb_clean`` that AutoDock Vina can type."""
        from vspipeline import pocket
        ions = []
        for key, atoms in pocket._hetatm_residues(self.pdb_clean).items():
            if len(atoms) == 1 and atoms[0]["element"] in pocket.VINA_METAL_TYPE_MAP:
                ions.append(key)
        return ions

    def _protonate(self, ph: float, source: Optional[str] = None) -> None:
        """Convert ``source`` (default ``pdb_clean``) to ``receptor_pdbqt``
        with Open Babel, adding hydrogens for ``ph``, with deterministic
        hydrogen positions."""
        source = source or self.pdb_clean
        probe = self.receptor_pdbqt + ".probe"
        # -p adds hydrogens for the given pH. Vina's atom typing needs them:
        # without polar hydrogens every N is typed as an acceptor and the
        # receptor has no hydrogen-bond donors. The PDBQT writer (-xr) keeps
        # the polar hydrogens and merges the non-polar ones.
        base = ["-ipdb", source, "-opdbqt", "-O"]
        options = ["-xr", "-p", f"{ph:g}", "--partialcharge", "gasteiger"]
        try:
            try:
                _run_obabel(base + [self.receptor_pdbqt] + options,
                            outputs=[self.receptor_pdbqt])
            except RuntimeError as e:
                raise RuntimeError(f"{e}{_unsupported_elements_hint(source)}") from e
            # A second conversion exposes every hydrogen Open Babel placed at
            # random (see stabilize_hydrogens). Open Babel seeds its generator
            # from the clock (usec % 70092), so two runs can draw the same
            # seed; when the two files are identical a third run is compared.
            _run_obabel(base + [probe] + options, outputs=[probe], log_stderr=False)
            if same_bytes(self.receptor_pdbqt, probe):
                _run_obabel(base + [probe] + options, outputs=[probe], log_stderr=False)
            moved = stabilize_hydrogens(self.receptor_pdbqt, probe)
        finally:
            if os.path.exists(probe):
                os.remove(probe)
        if moved:
            logger.info(f"[*] Placed {moved} hydrogens deterministically where Open "
                        "Babel's placement is random (isolated atoms, two-atom "
                        "fragments, over-coordinated atoms)")
        _strip_remark_name_paths(self.receptor_pdbqt)

    def prepare_receptor(self, chain: str = BLANK_CHAIN_ID, keep_hetero: bool = False,
                         ph: float = defaults.PH) -> List[str]:
        """Clean protein and convert to PDBQT."""
        available_chains = runner.get_chains(self.pdb_file)
        selected_chains = runner._parse_chain_selection(chain, available_chains)
        logger.info(
            f"[*] Preparing receptor (chain(s): {', '.join(selected_chains)})...")

        if os.path.realpath(self.pdb_file) == os.path.realpath(self.pdb_clean):
            raise ValueError(
                f"Receptor input {self.pdb_file} is the file the cleaned receptor "
                "is written to; move it out of the output directory")

        self._write_clean_pdb(selected_chains, keep_hetero)
        self.selected_chains = selected_chains

        try:
            # Metal ions kept by --keep-hetero stay in pdb_clean (and in the
            # frame offset) but are not given to Open Babel, which has no
            # Gasteiger parameters for them and would skip the whole receptor;
            # they are appended with their Vina metal types afterwards.
            ions = self._metal_ions() if keep_hetero else []
            with tempfile.TemporaryDirectory(prefix="autodocker_receptor_") as tmp:
                source = self.pdb_clean
                if ions:
                    # Same file name, so the PDBQT's REMARK Name is unchanged.
                    source = os.path.join(tmp, os.path.basename(self.pdb_clean))
                    _write_without_residues(self.pdb_clean, source, set(ions))
                self._protonate(ph, source)
                # obabel -p re-perceives residues; put the input's labels back.
                relabelled = _restore_residue_labels(self.receptor_pdbqt, source)
            logger.info(f"[*] Restored input residue labels on {relabelled} receptor atoms")
            # Centre on the heavy atoms of pdb_clean (the frame offset written
            # to frame_offset.txt).
            _translate_pdbqt(self.receptor_pdbqt, self._receptor_centroid())
            if ions:
                appended = runner._append_hetatm_to_receptor(
                    self.receptor_pdbqt, self.pdb_clean, ions, selected_chains,
                    offset=self.frame_offset())
                logger.info(
                    f"[*] --keep-hetero: added {appended} metal ion(s) with their Vina "
                    "metal types (not charged by Open Babel): "
                    + ", ".join(f"{n} {c}{r}" for c, r, n in ions))

            # Receptor PDBQT must not contain ligand torsion tags (ROOT/BRANCH/TORSDOF).
            runner._sanitize_receptor_pdbqt(self.receptor_pdbqt)

            if not runner._ensure_pdbqt_has_charges(self.receptor_pdbqt):
                raise RuntimeError(
                    f"Open Babel assigned no partial charges to {self.receptor_pdbqt}")
            n_hd, n_polar = _count_polar_atoms(self.receptor_pdbqt)
            if n_hd == 0 and n_polar > 0:
                raise RuntimeError(
                    f"Receptor {self.receptor_pdbqt} has no polar hydrogens after "
                    f"protonation at pH {ph:g}; Vina would see no H-bond donors")

            logger.info(f"[*] Receptor prepared: {self.receptor_pdbqt} "
                        f"(protonated at pH {ph:g}, {n_hd} polar hydrogens)")
        except Exception as e:
            logger.error(f"Failed to prepare receptor: {e}")
            raise

        return selected_chains

    def rebuild_receptor_keep_hetatm(self, keep_residues) -> str:
        """Append the selected HETATM residues to the receptor PDBQT.

        ``keep_residues`` is a sequence of ``(chain, resnum, resname)``
        triples identifying the HETATM residues to retain (waters, metal
        ions, cofactors). Appending (rather than re-converting the whole
        receptor) preserves the protein charges already assigned by
        OpenBabel. Returns the receptor PDBQT path.
        """
        if not os.path.exists(self.receptor_pdbqt):
            raise FileNotFoundError(
                f"Receptor PDBQT not found: {self.receptor_pdbqt}")

        selected_chains = getattr(self, "selected_chains", None) or [BLANK_CHAIN_ID]
        appended = runner._append_hetatm_to_receptor(
            self.receptor_pdbqt, self.pdb_file, keep_residues, selected_chains,
            offset=self.frame_offset())

        if appended == 0:
            logger.warning(
                "[!] No HETATM atoms were added to the receptor (none found, or "
                "all already present)")
            return self.receptor_pdbqt

        runner._sanitize_receptor_pdbqt(self.receptor_pdbqt)
        logger.info(
            f"[*] Receptor updated keeping {appended} HETATM atoms: {self.receptor_pdbqt}")
        return self.receptor_pdbqt

    def _box_source(self) -> str:
        if self.receptor_pdbqt and os.path.exists(self.receptor_pdbqt):
            return self.receptor_pdbqt
        return self.pdb_clean if os.path.exists(self.pdb_clean) else self.pdb_file

    def centroid_grid(self, padding: float = DEFAULT_PADDING
                      ) -> Tuple[float, float, float, float, float, float]:
        """Box around the whole receptor (used with --no-fpocket).

        Centred on the midpoint of the receptor's heavy-atom bounding box,
        ``padding`` on each side and every edge at least
        ``MIN_RECEPTOR_BOX_SIZE``. Raises if the receptor cannot be read.

        Public API (formerly the private ``_protein_centroid_grid``).
        """
        source = self._box_source()
        try:
            xs, ys, zs = _heavy_atom_coords(source)
        except (ValueError, IOError) as e:
            raise RuntimeError(f"Cannot build the receptor box from {source}: {e}") from e
        cx, cy, cz, sx, sy, sz = bounding_box(xs, ys, zs, padding, MIN_RECEPTOR_BOX_SIZE)
        logger.info(
            f"[*] Whole-receptor grid centred on the receptor bounding box "
            f"({cx:.2f}, {cy:.2f}, {cz:.2f}), size ({sx:.1f}, {sy:.1f}, {sz:.1f})")
        return cx, cy, cz, sx, sy, sz

    def _run_fpocket(self, target: str) -> str:
        """Run fpocket on ``target``; returns its output directory. Raises
        when fpocket is missing or fails."""
        fpocket = runner.find_tool("fpocket")
        if fpocket is None:
            raise RuntimeError(
                "fpocket not found on PATH; install fpocket, or pass --no-fpocket "
                "to dock in a box around the whole receptor, or give the box with "
                "--grid-center/--grid-size")
        version = _fpocket_version(fpocket)
        logger.info(f"[*] fpocket version: {version or '(unknown)'} ({fpocket})")
        pocket_root = Path(target).with_suffix("").as_posix() + "_out"
        shutil.rmtree(pocket_root, ignore_errors=True)
        _run_tool([fpocket, "-f", target], "fpocket")
        return pocket_root

    def _collect_pockets(self, pocket_root: str, target: str) -> List[Dict]:
        pocket_dir = os.path.join(pocket_root, "pockets")
        if not os.path.isdir(pocket_dir):
            raise RuntimeError(f"fpocket wrote no pocket directory ({pocket_dir})")
        info_file = os.path.join(pocket_root, f"{Path(target).stem}_info.txt")
        if not os.path.exists(info_file):
            raise RuntimeError(f"fpocket wrote no pocket scores ({info_file})")
        pocket_info = self._parse_fpocket_info(info_file)

        pockets = []
        for f in sorted(os.listdir(pocket_dir)):
            match = re.match(r"pocket(\d+)_atm\.pdb$", f)
            if not match:
                continue
            number = int(match.group(1))
            if number not in pocket_info:
                raise RuntimeError(f"Pocket {number} has no entry in {info_file}")
            p = os.path.join(pocket_dir, f)
            xs, ys, zs = runner.parse_pdb_coords(p)
            data = dict(pocket_info[number])
            data.update({
                "number": number,
                "file": f,
                "path": p,
                "box_volume": (max(xs) - min(xs)) * (max(ys) - min(ys)) * (max(zs) - min(zs)),
            })
            pockets.append(data)
        if not pockets:
            raise RuntimeError(
                f"fpocket found no pockets in {target}; pass --grid-center and "
                "--grid-size, or --no-fpocket for a box around the whole receptor")
        return pockets

    def detect_pocket(self, pocket_spec: Optional[str] = None,
                      padding: float = DEFAULT_PADDING
                      ) -> Tuple[float, float, float, float, float, float]:
        """Detect binding pockets with fpocket and build a grid from one or
        more of them. Raises when fpocket is unavailable, fails or finds no
        pocket (use ``centroid_grid`` for a deliberate whole-receptor box)."""
        logger.info("[*] Running fpocket...")
        fpocket_target = self.pdb_clean if os.path.exists(
            self.pdb_clean) else self.pdb_file
        pocket_root = self._run_fpocket(fpocket_target)
        pockets = rank_pockets(self._collect_pockets(pocket_root, fpocket_target))
        self._write_pocket_summary(pockets)

        logger.info(f"[*] Found {len(pockets)} pockets. Top candidates:")
        for pocket in pockets[:5]:
            logger.info(
                "  pocket%-2s score=%6.3f druggability=%5.3f volume=%7.1f %s"
                % (
                    pocket["number"],
                    pocket.get("score", 0),
                    pocket.get("druggability_score", 0),
                    pocket.get("volume", 0),
                    "druggable" if is_druggable(pocket) else "low-druggability",
                )
            )

        selected_numbers = self._parse_pocket_selection(pocket_spec, pockets)
        selected = [p for p in pockets if p["number"] in selected_numbers]
        selected.sort(key=lambda p: selected_numbers.index(p["number"]))

        selected_names = ", ".join(f"pocket{p['number']}" for p in selected)
        logger.info(
            f"[*] Pocket selection: {selected_names} (padding={padding:g} A)")
        cx, cy, cz, sx, sy, sz = self._get_pocket_info(
            [p["path"] for p in selected], padding=padding)
        dx, dy, dz = self._receptor_centroid()
        return cx - dx, cy - dy, cz - dz, sx, sy, sz

    def frame_offset(self) -> Tuple[float, float, float]:
        """Translation from the receptor PDBQT/docking frame to the input PDB frame.

        Add it to docked-pose coordinates to overlay them on the input PDB;
        subtract it to move input-PDB coordinates into the docking frame.
        """
        return self._receptor_centroid()

    def write_frame_offset(self) -> Optional[str]:
        """Record the frame offset in ``frame_offset.txt`` next to grid.conf."""
        path = os.path.join(self.outdir, outputs.FRAME_OFFSET_TXT)
        try:
            dx, dy, dz = self.frame_offset()
            with open(path, "w", encoding="utf-8") as f:
                f.write(
                    "# receptor.pdbqt, grid.conf and docked poses are translated so\n"
                    "# the mean of the receptor heavy atoms is at the origin. Add this\n"
                    "# offset to their coordinates to overlay them on the input PDB, e.g. with\n"
                    "#   python3 scripts/restore_frame.py <output_dir> docked/*_out.pdbqt\n"
                    f"offset_x = {dx}\noffset_y = {dy}\noffset_z = {dz}\n")
        except (ValueError, IOError) as e:
            logger.warning(f"[!] Could not write frame offset: {e}")
            return None
        return path

    def _receptor_centroid(self) -> Tuple[float, float, float]:
        """Unweighted mean of the cleaned receptor's heavy atoms.

        ``prepare_receptor`` translates the receptor PDBQT by this vector so
        it is centred at the origin, while ``pdb_clean`` and fpocket's pocket
        files keep the original coordinates; the fpocket-derived grid centre
        is shifted by it to land in the receptor PDBQT frame.
        """
        source = self.pdb_clean if os.path.exists(
            self.pdb_clean) else self.pdb_file
        xs, ys, zs = _heavy_atom_coords(source)
        return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)

    def _parse_fpocket_info(self, info_file: str) -> Dict[int, Dict]:
        """Parse fpocket's *_info.txt file into a pocket metadata map.

        Raises if the file cannot be read.
        """
        info = {}
        current = None
        try:
            with open(info_file, encoding="utf-8") as f:
                for line in f:
                    pocket_match = re.match(
                        r"Pocket\s+(\d+)\s*:", line.strip())
                    if pocket_match:
                        current = int(pocket_match.group(1))
                        info[current] = {}
                        continue

                    if current is None or ":" not in line:
                        continue

                    key, value = line.split(":", 1)
                    key = key.strip().lower().replace(" ", "_").replace(".", "")
                    key = key.replace("-", "_")
                    try:
                        info[current][key] = float(value.strip().split()[0])
                    except (ValueError, IndexError):
                        continue
        except IOError as e:
            raise RuntimeError(f"Could not read fpocket info file {info_file}: {e}") from e

        return info

    def _write_pocket_summary(self, pockets: List[Dict]):
        try:
            with open(self.pocket_summary_file, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow([
                    "Pocket", "Score", "Druggability_Score", "Volume",
                    "Box_Volume", "Could_Be_Druggable", "File"
                ])
                for pocket in sorted(pockets, key=lambda p: p["number"]):
                    writer.writerow([
                        pocket["number"],
                        pocket.get("score", ""),
                        pocket.get("druggability_score", ""),
                        pocket.get("volume", ""),
                        f"{pocket.get('box_volume', 0):.3f}",
                        "YES" if is_druggable(pocket) else "NO",
                        pocket["file"],
                    ])
            logger.info(
                f"[*] Pocket summary saved: {self.pocket_summary_file}")
        except IOError as e:
            logger.warning(f"Could not save pocket summary: {e}")

    def _parse_pocket_selection(self, pocket_spec: Optional[str], pockets: List[Dict]) -> List[int]:
        available = {p["number"] for p in pockets}
        if not pocket_spec or pocket_spec.lower() == "auto":
            return [pockets[0]["number"]]

        selected = []
        for token in re.split(r"[,\s]+", pocket_spec.strip()):
            if not token:
                continue
            try:
                pocket_number = int(token)
            except ValueError as e:
                raise ValueError(
                    f"Invalid pocket selection '{token}'. Use numbers like 1 or 1,3.") from e
            if pocket_number not in available:
                raise ValueError(
                    f"Pocket {pocket_number} was not found. Available: {sorted(available)}")
            selected.append(pocket_number)

        if not selected:
            raise ValueError("No pockets selected")
        return selected

    def _get_pocket_info(self, pocket_files: List[str], padding: float = DEFAULT_PADDING
                         ) -> Tuple[float, float, float, float, float, float]:
        """Grid centre and size from one or more pocket PDB files: the
        midpoint of the pocket atoms' bounding box, ``padding`` on each side."""
        xs, ys, zs = [], [], []
        for pocket_file in pocket_files:
            px, py, pz = runner.parse_pdb_coords(pocket_file)
            xs.extend(px)
            ys.extend(py)
            zs.extend(pz)

        cx, cy, cz, sx, sy, sz = bounding_box(xs, ys, zs, padding)
        logger.debug(f"Grid center: ({cx:.2f}, {cy:.2f}, {cz:.2f}), "
                     f"size: ({sx:.2f}, {sy:.2f}, {sz:.2f})")
        return cx, cy, cz, sx, sy, sz

    def write_grid(self, cx: float, cy: float, cz: float,
                   sx: float, sy: float, sz: float):
        """Write AutoDock Vina grid configuration."""
        try:
            with open(self.grid_conf, "w", encoding="utf-8") as f:
                f.write(f"""center_x = {cx}
center_y = {cy}
center_z = {cz}
size_x = {sx}
size_y = {sy}
size_z = {sz}
""")
            logger.info(f"[*] Grid ready: {self.grid_conf}")
            self.write_grid_box_script(cx, cy, cz, sx, sy, sz)
            self.write_frame_offset()
        except IOError as e:
            logger.error(f"Failed to write grid file: {e}")
            raise

    def write_grid_box_script(self, cx: float, cy: float, cz: float,
                              sx: float, sy: float, sz: float):
        """Write a PyMOL script that visualizes the final docking box."""
        _dir = os.path.dirname(os.path.abspath(self.grid_box_script))
        pdb_relpath = os.path.relpath(self.receptor_pdbqt, _dir)

        x0, x1 = cx - sx / 2, cx + sx / 2
        y0, y1 = cy - sy / 2, cy + sy / 2
        z0, z1 = cz - sz / 2, cz + sz / 2
        vertices = [
            (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
            (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
        ]
        edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ]
        # CGO: COLOR r, g, b (green), LINEWIDTH w, BEGIN LINES, VERTEX x, y, z ..., END
        cgo = [("COLOR", "0.0", "1.0", "0.0"), ("LINEWIDTH", "3.0"), ("BEGIN", "LINES")]
        for a, b in edges:
            for x, y, z in (vertices[a], vertices[b]):
                cgo.append(("VERTEX", f"{x:7.3f}", f"{y:7.3f}", f"{z:7.3f}"))
        cgo.append(("END",))
        box = ",\n       ".join(", ".join(item) for item in cgo)

        pml_lines = (
            f'load {pdb_relpath}, receptor',
            r'show cartoon, receptor',
            r'',
            r'python',
            r'from pymol import cmd',
            r'from pymol.cgo import *',
            r'',
            f'box = [{box}]',
            r'cmd.load_cgo(box, "docking_grid")',
            r'python end',
            r'',
            r'zoom all',
            r''
        )

        try:
            with open(self.grid_box_script, "w", encoding="utf-8") as f:
                f.write('\n'.join(pml_lines))
                logger.info(
                    f"[*] Grid visualization script saved: {self.grid_box_script}")
        except IOError as e:
            logger.warning(f"Could not write grid visualization script: {e}")

    # Backward-compatible alias (deprecated): the former private name.
    _protein_centroid_grid = centroid_grid
