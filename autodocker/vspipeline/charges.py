#!/usr/bin/env python3
"""PDBQT charge inspection and Gasteiger charge assignment.

Split out of runner.py. Cross-module references (run_obabel) go through the
``runner`` namespace at call time so the monkeypatch surface stays unchanged.

Note: the Vina, Vinardo and default SMINA scoring functions do not use
partial charges, so all-zero charges are reported as a warning by the docking
step, not treated as a failure.
"""
import logging
import os
from typing import List, Optional, Tuple

# The runner module owns the tool-path globals; attributes are only accessed
# at call time, so binding the partially-initialized module is fine.
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


def _pdbqt_has_nonzero_charges(pdbqt_file: str) -> bool:
    """True if the PDBQT has at least one ATOM/HETATM with a nonzero charge.

    Only inspects the file. False for all-zero charges, no atoms, or an
    unreadable file.
    """
    try:
        with open(pdbqt_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    charge = _extract_pdbqt_charge(line)
                    if charge is not None and charge != 0.0:
                        return True
        return False
    except Exception as e:
        logger.warning(f"Charge check failed: {e}")
        return False


# Backwards-compatible name (re-exported by runner.py, monkeypatched by tests).
_ensure_pdbqt_has_charges = _pdbqt_has_nonzero_charges


def _pdbqt_has_atoms(pdbqt_file: str) -> bool:
    """Return True when a PDBQT contains at least one ATOM/HETATM record."""
    try:
        with open(pdbqt_file, 'r', encoding='utf-8') as f:
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
        with open(pdbqt_file, "r", encoding='utf-8') as f:
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

        with open(pdbqt_file, "w", encoding='utf-8') as f:
            f.writelines(kept)

        if removed:
            logger.info(
                f"[*] Sanitized receptor PDBQT: removed {removed} ligand-only tags")
    except Exception as e:
        raise RuntimeError(
            f"Failed to sanitize receptor PDBQT {pdbqt_file}: {e}")


_OBABEL_INPUT_FORMATS = {"sdf": "sdf", "sd": "sdf", "mol": "mol", "mol2": "mol2",
                         "pdbqt": "pdbqt"}

# Coordinates whose z values all lie within this of 0 (or whose atoms all
# share one position) are treated as 2D / missing geometry.
_FLAT_TOLERANCE = 1e-4


def _read_coordinates(path: str) -> List[Tuple[float, float, float]]:
    """Atom coordinates of a PDB/PDBQT, SDF/MOL (V2000) or MOL2 file
    (first molecule only); [] when none can be read."""
    ext = os.path.splitext(path)[1].lower().lstrip(".")
    coords: List[Tuple[float, float, float]] = []
    try:
        with open(path, encoding='utf-8') as f:
            lines = f.read().splitlines()
    except OSError:
        return coords
    if ext in ("sdf", "sd", "mol"):
        if len(lines) > 3:
            try:
                n_atoms = int(lines[3][0:3])
            except ValueError:
                return coords
            for line in lines[4:4 + n_atoms]:
                try:
                    coords.append((float(line[0:10]), float(line[10:20]), float(line[20:30])))
                except ValueError:
                    continue
        return coords
    if ext == "mol2":
        in_atoms = False
        for line in lines:
            if line.startswith("@<TRIPOS>"):
                if in_atoms:
                    break
                in_atoms = line.startswith("@<TRIPOS>ATOM")
                continue
            parts = line.split()
            if in_atoms and len(parts) >= 5:
                try:
                    coords.append((float(parts[2]), float(parts[3]), float(parts[4])))
                except ValueError:
                    continue
        return coords
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
        elif line.startswith("ENDMDL"):
            break
    return coords


def _has_3d_coordinates(path: str) -> bool:
    """True when the file holds atoms with a non-planar (z != 0) geometry."""
    coords = _read_coordinates(path)
    if not coords:
        return False
    if all(abs(c[2]) < _FLAT_TOLERANCE for c in coords):
        return False
    return len(set(coords)) > 1


def _fix_pdbqt_charges(pdbqt_file: str, source_file: str,
                       gen3d: Optional[bool] = None) -> None:
    """Rebuild a ligand PDBQT with Open Babel Gasteiger charges, never
    discarding 3D coordinates.

    The geometry comes from, in order: the existing ``pdbqt_file`` when it
    already has 3D coordinates (rebuilt from itself); else ``source_file``.
    ``--gen3d`` is added when ``gen3d`` is True, or, when ``gen3d`` is None,
    when the chosen input has no 3D coordinates (2D SDF, flat PDB), so a 2D
    source is never docked flat.

    All-zero charges after the rebuild (pure hydrocarbons once nonpolar
    hydrogens are merged, e.g. benzene) are accepted with a warning: the
    Vina, Vinardo and default SMINA scoring functions do not use partial
    charges. Raises RuntimeError only when Open Babel fails or writes no atoms.
    """
    logger.warning(f"[!] Recomputing Gasteiger charges for: {pdbqt_file}")

    rebuild_in_place = os.path.exists(pdbqt_file) and _has_3d_coordinates(pdbqt_file)
    input_file = pdbqt_file if rebuild_in_place else source_file
    ext = os.path.splitext(input_file)[1].lower().lstrip(".")
    in_fmt = "pdbqt" if rebuild_in_place else _OBABEL_INPUT_FORMATS.get(ext, "pdb")
    if gen3d is None:
        gen3d = not rebuild_in_place and not _has_3d_coordinates(input_file)

    target = pdbqt_file + ".charges.tmp.pdbqt" if rebuild_in_place else pdbqt_file
    # -h: hydrogens are added again (the PDBQT keeps polar H only, and a
    # source rebuild would otherwise lose the H added by the first conversion).
    args = [f"-i{in_fmt}", input_file, "-opdbqt", "-O", target,
            "--partialcharge", "gasteiger", "-h"]
    if gen3d:
        logger.warning(f"[!] {input_file} has no 3D coordinates; generating them (--gen3d)")
        args.append("--gen3d")

    if not os.path.exists(input_file):
        raise RuntimeError(
            f"Cannot assign valid charges for {pdbqt_file}: source {input_file} not found")
    if not rebuild_in_place and os.path.exists(pdbqt_file):
        os.remove(pdbqt_file)  # rebuilt from the source; never keep the old file
    try:
        # Single-threaded; a partial output is deleted when Open Babel fails.
        # Without --gen3d the conversion is run through the deterministic
        # hydrogen wrapper (-h places some hydrogens at random); --gen3d
        # output differs in its heavy atoms on every run, so it runs once.
        if gen3d:
            runner.run_obabel(args, outputs=[target])
        else:
            from vspipeline.hydrogens import run_obabel_deterministic
            run_obabel_deterministic(args, target, kind="ligand")
        if rebuild_in_place:
            os.replace(target, pdbqt_file)
    except Exception as e:
        if rebuild_in_place and os.path.exists(target):
            os.remove(target)
        raise RuntimeError(
            f"Cannot assign valid charges for {pdbqt_file}: Open Babel failed ({e})")

    if not os.path.exists(pdbqt_file) or not _pdbqt_has_atoms(pdbqt_file):
        raise RuntimeError(
            f"Cannot assign valid charges for {pdbqt_file}: Open Babel wrote no atoms")
    if _pdbqt_has_nonzero_charges(pdbqt_file):
        logger.info("[*] Charges assigned via OpenBabel Gasteiger")
    else:
        logger.warning(
            f"[!] All Gasteiger charges are zero in {pdbqt_file} (e.g. a pure "
            "hydrocarbon); accepted, since Vina-family scoring ignores charges")
