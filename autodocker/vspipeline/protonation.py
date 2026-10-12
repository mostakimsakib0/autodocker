#!/usr/bin/env python3
"""Receptor protonation with Open Babel and the repairs of its output.

``obabel -p`` adds the hydrogens for the requested pH and writes the PDBQT;
a probe conversion with another OB_RANDOM_SEED exposes the hydrogens it
placed at random, which are re-placed deterministically (hydrogens.py), and
the residue labels it re-perceived are restored from the input. Open Babel
runs through ``utils.run_obabel``, called through its module so tests can
replace it.
"""
import logging
import math
import os
from collections import defaultdict
from typing import List, Optional, Sequence, Tuple

from vspipeline import utils
from vspipeline.hydrogens import (
    H_PARENT_MAX_DISTANCE,
    SpatialHash,
    probe_seed,
    same_bytes,
    stabilize_hydrogens,
)
from vspipeline.pdb_records import (
    _is_pdbqt_hydrogen,
    is_hydrogen_record,
    pdb_record_element,
    record_chain,
)

logger = logging.getLogger(__name__)


def _run_obabel(args: List[str], outputs: Sequence[str] = (), log_stderr: bool = True,
                seed: Optional[int] = None) -> str:
    """Open Babel through ``utils.run_obabel`` (single-threaded,
    OB_RANDOM_SEED = ``seed``, partial outputs removed on failure); its
    stderr goes to the log unless ``log_stderr`` is False. Returns stdout."""
    stdout, stderr = utils.run_obabel(list(args), outputs=outputs, return_stderr=True,
                                      seed=seed)
    if log_stderr:
        utils.log_tool_stderr("obabel", stderr)
    return stdout


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
    bad = utils.elements_without_gasteiger(found)
    if not bad:
        return ""
    where = "; ".join(f"{e} in " + ", ".join(sorted(found[e])[:3])
                      + (" ..." if len(found[e]) > 3 else "") for e in bad)
    return (f" (element(s) without Open Babel Gasteiger parameters: {where}; remove "
            "those residues or run without --keep-hetero)")


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


def protonate(source: str, receptor_pdbqt: str, ph: float, seed: int) -> None:
    """Convert ``source`` (a cleaned receptor PDB) to ``receptor_pdbqt``
    with Open Babel, adding hydrogens for ``ph``, with deterministic
    hydrogen positions; OB_RANDOM_SEED = ``seed``."""
    probe = receptor_pdbqt + ".probe"
    # -p adds hydrogens for the given pH. Vina's atom typing needs them:
    # without polar hydrogens every N is typed as an acceptor and the
    # receptor has no hydrogen-bond donors. The PDBQT writer (-xr) keeps
    # the polar hydrogens and merges the non-polar ones.
    base = ["-ipdb", source, "-opdbqt", "-O"]
    options = ["-xr", "-p", f"{ph:g}", "--partialcharge", "gasteiger"]
    try:
        try:
            _run_obabel(base + [receptor_pdbqt] + options,
                        outputs=[receptor_pdbqt], seed=seed)
        except RuntimeError as e:
            raise RuntimeError(f"{e}{_unsupported_elements_hint(source)}") from e
        # A second conversion with a different OB_RANDOM_SEED exposes
        # every hydrogen Open Babel placed at random (see
        # stabilize_hydrogens and run_obabel_deterministic). Open Babel
        # 3.1 ignores the seed and seeds from the clock (usec % 70092),
        # so two runs can draw the same seed; when the two files are
        # identical a third run (another seed) is compared.
        _run_obabel(base + [probe] + options, outputs=[probe], log_stderr=False,
                    seed=probe_seed(seed, 1))
        if same_bytes(receptor_pdbqt, probe):
            _run_obabel(base + [probe] + options, outputs=[probe], log_stderr=False,
                        seed=probe_seed(seed, 2))
        moved = stabilize_hydrogens(receptor_pdbqt, probe)
    finally:
        if os.path.exists(probe):
            os.remove(probe)
    if moved:
        logger.info(f"[*] Placed {moved} hydrogens deterministically where Open "
                    "Babel's placement is random (isolated atoms, two-atom "
                    "fragments, over-coordinated atoms)")
    _strip_remark_name_paths(receptor_pdbqt)


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

    heavy = SpatialHash(H_PARENT_MAX_DISTANCE)
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
