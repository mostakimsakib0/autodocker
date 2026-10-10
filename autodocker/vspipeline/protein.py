#!/usr/bin/env python3
"""PROTEIN PREPARATION.

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

# Modified amino acids deposited as HETATM; they are part of the protein and
# dropping them would leave holes in the backbone or side chains.
MODIFIED_RESIDUES = {
    "MSE", "SEP", "TPO", "PTR", "CSO", "CSD", "CME", "MLY", "HYP", "KCX",
    "PCA", "CSS", "OCS", "SCH", "NEP", "ALY", "M3L", "TYS",
}


def _translate_pdbqt(path: str, offset: Tuple[float, float, float]) -> None:
    """Subtract ``offset`` from every ATOM/HETATM coordinate of a PDBQT file."""
    dx, dy, dz = offset
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                x = float(line[30:38]) - dx
                y = float(line[38:46]) - dy
                z = float(line[46:54]) - dz
                line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
            out.append(line)
    with open(path, "w") as f:
        f.writelines(out)


def _restore_residue_labels(pdbqt_file: str, pdb_file: str) -> int:
    """Copy atom names and residue labels from ``pdb_file`` back onto the
    PDBQT that Open Babel wrote from it.

    ``obabel -p`` re-perceives residues, and on a receptor with a chain break
    it renames, renumbers and re-chains most atoms. Heavy-atom coordinates
    are untouched, so each heavy atom is matched to its input record by
    coordinates; each added hydrogen takes the labels of the heavy atom it is
    bonded to. Returns the number of relabelled atoms; raises if a heavy atom
    has no match.
    """
    def key(line):
        return (round(float(line[30:38]), 3), round(float(line[38:46]), 3),
                round(float(line[46:54]), 3))

    labels = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                # columns 13-27: atom name, altloc, resname, chain, resnum, icode
                labels.setdefault(key(line), (line[:6], line[12:27]))

    with open(pdbqt_file) as f:
        lines = f.readlines()
    atoms = [(i, l) for i, l in enumerate(lines) if l.startswith(("ATOM", "HETATM"))]
    heavy_xyz = []
    changed = 0
    for i, line in atoms:
        if line[77:79].strip() in ("HD", "H"):
            continue
        found = labels.get(key(line))
        if found is None:
            raise RuntimeError(
                f"Receptor atom at ({line[30:54].strip()}) in {pdbqt_file} has no "
                f"matching record in {pdb_file}; cannot restore residue labels")
        record, label = found
        new = record + line[6:12] + label + line[27:]
        changed += new != line
        lines[i] = new
        heavy_xyz.append((key(line), label))
    # Spatial hash (2 A cells) so each hydrogen only checks nearby heavy atoms.
    cell = 2.0
    grid = defaultdict(list)
    for xyz, label in heavy_xyz:
        grid[tuple(int(math.floor(c / cell)) for c in xyz)].append((xyz, label))
    for i, line in atoms:
        if line[77:79].strip() not in ("HD", "H"):
            continue
        x, y, z = key(line)
        cx, cy, cz = (int(math.floor(c / cell)) for c in (x, y, z))
        nearby = [hk for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
                  for hk in grid.get((cx + dx, cy + dy, cz + dz), ())]
        if not nearby:
            raise RuntimeError(f"Hydrogen at ({line[30:54].strip()}) has no heavy atom "
                               f"within {cell} A in {pdbqt_file}")
        near = min(nearby, key=lambda hk: (hk[0][0] - x) ** 2 + (hk[0][1] - y) ** 2
                   + (hk[0][2] - z) ** 2)
        # Hydrogen keeps its own name ("H"); residue fields come from its heavy atom.
        new = line[:17] + near[1][5:] + line[27:]
        changed += new != line
        lines[i] = new
    with open(pdbqt_file, "w") as f:
        f.writelines(lines)
    return changed


def _count_polar_atoms(pdbqt_file: str) -> Tuple[int, int]:
    """(polar hydrogens, N/O heavy atoms) in a PDBQT file."""
    hd = no = 0
    with open(pdbqt_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                t = line[77:79].strip()
                hd += t == "HD"
                no += t in ("N", "NA", "NS", "OA", "OS", "O")
    return hd, no


class ProteinPreparation:
    """Advanced protein preparation with chain selection."""

    def __init__(self, pdb_file: str, outdir: str):
        self.pdb_file = pdb_file
        self.outdir = outdir
        self.pdb_clean = os.path.join(outdir, "protein_clean.pdb")
        self.receptor_pdbqt = os.path.join(outdir, "receptor.pdbqt")
        self.grid_conf = os.path.join(outdir, "grid.conf")
        self.grid_box_script = os.path.join(outdir, "grid_box.pml")
        self.pocket_summary_file = os.path.join(outdir, "pocket_summary.csv")

    def validate(self):
        """Validate protein PDB file."""
        if not os.path.exists(self.pdb_file):
            raise FileNotFoundError(f"Protein file not found: {self.pdb_file}")

        if not os.access(self.pdb_file, os.R_OK):
            raise PermissionError(f"Cannot read protein file: {self.pdb_file}")

        try:
            count = sum(1 for l in open(self.pdb_file) if l.startswith("ATOM"))
        except IOError as e:
            raise IOError(f"Error reading protein file: {e}")

        if count < 100:
            raise ValueError(
                f"Invalid protein: only {count} ATOM lines (need >= 100)")

        logger.info(f"[✔] Protein validated: {count} ATOM lines")

    def select_chain(self) -> str:
        """Choose a chain for docking. Non-interactive: never prompts."""
        chains = runner.get_chains(self.pdb_file)

        if len(chains) == 1:
            logger.info(f"[*] Single chain detected: {chains[0]}")
            return chains[0]

        logger.info(f"[*] Multiple chains detected: {', '.join(chains)}")
        logger.warning(
            f"[!] Multiple chains present; using first chain '{chains[0]}'. "
            "Specify --chain to choose explicitly.")
        return chains[0]

    def prepare_receptor(self, chain: str = "A", keep_hetero: bool = False,
                         ph: float = 7.4) -> List[str]:
        """Clean protein and convert to PDBQT."""
        available_chains = runner.get_chains(self.pdb_file)
        selected_chains = runner._parse_chain_selection(chain, available_chains)
        logger.info(
            f"[*] Preparing receptor (chain(s): {', '.join(selected_chains)})...")

        if os.path.realpath(self.pdb_file) == os.path.realpath(self.pdb_clean):
            raise ValueError(
                f"Receptor input {self.pdb_file} is the file the cleaned receptor "
                "is written to; move it out of the output directory")

        kept_atoms = 0
        skipped_hetero = 0
        skipped_altloc = 0
        seen_atoms = set()
        try:
            with open(self.pdb_file, "r") as src, open(self.pdb_clean, "w") as dst:
                for line in src:
                    # NMR ensembles: keep the first model only; stacking every
                    # model superimposes copies of the receptor.
                    if line.startswith("ENDMDL"):
                        break
                    is_atom = line.startswith("ATOM")
                    is_het = line.startswith("HETATM")
                    if is_het and line[17:20].strip() in MODIFIED_RESIDUES:
                        line = "ATOM  " + line[6:]
                        is_atom, is_het = True, False
                    if is_atom or (keep_hetero and is_het):
                        if len(line) <= 21:
                            continue
                        chain = line[21].strip() or "A"
                        if chain not in selected_chains:
                            continue
                        # Alternate locations: keep the first conformer of
                        # each atom so side chains are not duplicated.
                        atom_key = (line[21], line[22:27], line[12:16])
                        if line[16] != " ":
                            if atom_key in seen_atoms:
                                skipped_altloc += 1
                                continue
                            line = line[:16] + " " + line[17:]
                        seen_atoms.add(atom_key)
                        dst.write(line)
                        kept_atoms += 1
                    elif is_het:
                        skipped_hetero += 1
                    elif line.startswith("TER"):
                        dst.write(line)
                dst.write("END\n")
        except IOError as e:
            raise IOError(
                f"Failed to filter chain(s) {', '.join(selected_chains)}: {e}")

        if kept_atoms == 0:
            raise ValueError(
                f"Selected chain(s) '{', '.join(selected_chains)}' have no atoms in {self.pdb_file}")

        if skipped_altloc:
            logger.info(f"[*] Kept the first alternate location; dropped "
                        f"{skipped_altloc} alternate-location atom records")

        try:
            # -p adds hydrogens for the given pH. Vina's atom typing needs
            # them: without polar hydrogens every N is typed as an acceptor
            # and the receptor has no hydrogen-bond donors. The PDBQT writer
            # (-xr) keeps the polar hydrogens and merges the non-polar ones.
            runner.run([
                runner.OBABEL,
                "-ipdb", self.pdb_clean,
                "-opdbqt", "-O", self.receptor_pdbqt,
                "-xr", "-p", f"{ph:g}", "--partialcharge", "gasteiger",
            ])
            # obabel -p re-perceives residues; put the input's labels back.
            _restore_residue_labels(self.receptor_pdbqt, self.pdb_clean)
            # Center on the heavy atoms of pdb_clean (the frame offset written
            # to frame_offset.txt), not on all atoms as obabel -c would now
            # that hydrogens are present.
            _translate_pdbqt(self.receptor_pdbqt, self._receptor_centroid())

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

            logger.info(f"[✔] Receptor prepared: {self.receptor_pdbqt} "
                        f"(protonated at pH {ph:g}, {n_hd} polar hydrogens)")
            if skipped_hetero and not keep_hetero:
                logger.info(
                    f"[*] Removed {skipped_hetero} HETATM records during receptor cleaning")
        except Exception as e:
            logger.error(f"Failed to prepare receptor: {e}")
            raise

        self.selected_chains = selected_chains

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

        selected_chains = getattr(self, "selected_chains", None) or ["A"]
        appended = runner._append_hetatm_to_receptor(
            self.receptor_pdbqt, self.pdb_file, keep_residues, selected_chains,
            offset=self.frame_offset())

        if appended == 0:
            logger.warning(
                "[!] None of the requested HETATM residues could be added to the "
                "receptor; docking without them")
            return self.receptor_pdbqt

        runner._sanitize_receptor_pdbqt(self.receptor_pdbqt)
        logger.info(
            f"[✔] Receptor updated keeping {appended} HETATM atoms: {self.receptor_pdbqt}")
        return self.receptor_pdbqt

    def centroid_grid(self, padding: float = 6.0) -> Tuple[float, float, float, float, float, float]:
        """Build a grid centered on the receptor bounding box when fpocket is unavailable.

        Public API (formerly the private ``_protein_centroid_grid``).
        """
        source = self.receptor_pdbqt
        if not source or not os.path.exists(source):
            source = self.pdb_clean if os.path.exists(
                self.pdb_clean) else self.pdb_file
        try:
            xs, ys, zs = runner.parse_pdb_coords(source)
        except (ValueError, IOError):
            return 0, 0, 0, 24, 24, 24
        cx = (min(xs) + max(xs)) / 2.0
        cy = (min(ys) + max(ys)) / 2.0
        cz = (min(zs) + max(zs)) / 2.0
        sx = max(max(xs) - min(xs) + 2 * padding, 24.0)
        sy = max(max(ys) - min(ys) + 2 * padding, 24.0)
        sz = max(max(zs) - min(zs) + 2 * padding, 24.0)
        logger.info(
            f"[✔] Fallback grid centered on receptor centroid "
            f"({cx:.2f}, {cy:.2f}, {cz:.2f}), size ({sx:.1f}, {sy:.1f}, {sz:.1f})")
        return cx, cy, cz, sx, sy, sz

    def detect_pocket(self, pocket_spec: Optional[str] = None,
                      padding: float = 6.0) -> Tuple[float, float, float, float, float, float]:
        """Detect binding pockets using fpocket and build a grid from one or more pockets."""
        logger.info("[*] Running fpocket...")
        fpocket_target = self.pdb_clean if os.path.exists(
            self.pdb_clean) else self.pdb_file
        try:
            shutil.rmtree(Path(fpocket_target).with_suffix(
                "").as_posix() + "_out", ignore_errors=True)
            runner.run(["fpocket", "-f", fpocket_target])
        except Exception as e:
            logger.warning(
                f"fpocket failed: {e}. Using centroid fallback grid.")
            return self.centroid_grid(padding)

        pocket_root = Path(fpocket_target).with_suffix("").as_posix() + "_out"
        pocket_dir = os.path.join(pocket_root, "pockets")

        if not os.path.exists(pocket_dir):
            logger.warning(
                f"Pocket directory not found: {pocket_dir}. Using centroid fallback grid.")
            return self.centroid_grid(padding)

        info_file = os.path.join(
            pocket_root, f"{Path(fpocket_target).stem}_info.txt"
        )
        pocket_info = self._parse_fpocket_info(info_file)
        pockets = []

        for f in os.listdir(pocket_dir):
            match = re.match(r"pocket(\d+)_atm\.pdb$", f)
            if not match:
                continue

            number = int(match.group(1))
            p = os.path.join(pocket_dir, f)

            try:
                xs, ys, zs = runner.parse_pdb_coords(p)
                size = (max(xs) - min(xs)) * \
                    (max(ys) - min(ys)) * (max(zs) - min(zs))
                data = dict(pocket_info.get(number, {}))
                data.update({
                    "number": number,
                    "file": f,
                    "path": p,
                    "box_volume": size,
                })
                pockets.append(data)
            except ValueError:
                continue

        if not pockets:
            logger.warning("No valid pockets found. Using fallback grid.")
            return 0, 0, 0, 24, 24, 24

        # FIXED: Sort by SCORE first (binding pocket quality), then druggability
        # Previously: sorted by druggability_score first, which selected worst-scoring
        # pocket (Pocket 4) over best-scoring pocket (Pocket 1) when druggability was slightly higher.
        # New behavior: Best-scoring pocket wins, with druggability as tiebreaker.
        pockets.sort(key=lambda p: (
            p.get("score", -1),
            p.get("druggability_score", -1),
            p.get("volume", -1),
            p.get("box_volume", -1),
        ), reverse=True)
        self._write_pocket_summary(pockets)

        logger.info(f"[*] Found {len(pockets)} pockets. Top candidates:")
        for pocket in pockets[:5]:
            druggable = "druggable" if pocket.get(
                "druggability_score", 0) > 0.15 else "low-druggability"
            logger.info(
                "  pocket%-2s score=%6.3f druggability=%5.3f volume=%7.1f %s"
                % (
                    pocket["number"],
                    pocket.get("score", 0),
                    pocket.get("druggability_score", 0),
                    pocket.get("volume", 0),
                    druggable,
                )
            )

        selected_numbers = self._parse_pocket_selection(pocket_spec, pockets)
        selected = [p for p in pockets if p["number"] in selected_numbers]
        selected.sort(key=lambda p: selected_numbers.index(p["number"]))

        selected_names = ", ".join(f"pocket{p['number']}" for p in selected)
        logger.info(
            f"[✔] Pocket selection: {selected_names} (padding={padding:g} Å)")
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
        path = os.path.join(self.outdir, "frame_offset.txt")
        try:
            dx, dy, dz = self.frame_offset()
            with open(path, "w") as f:
                f.write(
                    "# receptor.pdbqt, grid.conf and docked poses are centered at the\n"
                    "# origin (obabel -c). Add this offset to their coordinates to\n"
                    "# overlay them on the input PDB, e.g. with\n"
                    "#   python3 scripts/restore_frame.py <output_dir> docked/*_out.pdbqt\n"
                    f"offset_x = {dx}\noffset_y = {dy}\noffset_z = {dz}\n")
        except (ValueError, IOError, ZeroDivisionError) as e:
            logger.warning(f"[!] Could not write frame offset: {e}")
            return None
        return path

    def _receptor_centroid(self) -> Tuple[float, float, float]:
        """Center of mass of the prepared receptor.

        The receptor PDBQT is centered at the origin by ``obabel -c`` during
        ``prepare_receptor``, while ``pdb_clean`` and fpocket's pocket files
        keep the original coordinates. The fpocket-derived grid center must be
        shifted by this centroid to land in the centered receptor frame.
        """
        source = self.pdb_clean if os.path.exists(
            self.pdb_clean) else self.pdb_file
        xs, ys, zs = runner.parse_pdb_coords(source)
        return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)

    def _parse_fpocket_info(self, info_file: str) -> Dict[int, Dict]:
        """Parse fpocket's *_info.txt file into a pocket metadata map."""
        info = {}
        if not os.path.exists(info_file):
            return info

        current = None
        try:
            with open(info_file) as f:
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
            logger.warning(
                f"Could not parse fpocket info file {info_file}: {e}")

        return info

    def _write_pocket_summary(self, pockets: List[Dict]):
        try:
            with open(self.pocket_summary_file, "w", newline="") as f:
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
                        "YES" if pocket.get(
                            "druggability_score", 0) > 0.15 else "NO",
                        pocket["file"],
                    ])
            logger.info(
                f"[✔] Pocket summary saved: {self.pocket_summary_file}")
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

    def _get_pocket_info(self, pocket_files: List[str], padding: float = 6.0) -> Tuple[float, float, float, float, float, float]:
        """Calculate grid center and size from one or more pocket PDB files."""
        xs, ys, zs = [], [], []
        for pocket_file in pocket_files:
            px, py, pz = runner.parse_pdb_coords(pocket_file)
            xs.extend(px)
            ys.extend(py)
            zs.extend(pz)

        cx = sum(xs) / len(xs)
        cy = sum(ys) / len(ys)
        cz = sum(zs) / len(zs)

        sx = max(xs) - min(xs) + (2 * padding)
        sy = max(ys) - min(ys) + (2 * padding)
        sz = max(zs) - min(zs) + (2 * padding)

        logger.debug(f"Grid center: ({cx:.2f}, {cy:.2f}, {cz:.2f}), "
                     f"size: ({sx:.2f}, {sy:.2f}, {sz:.2f})")

        return cx, cy, cz, sx, sy, sz

    def write_grid(self, cx: float, cy: float, cz: float,
                   sx: float, sy: float, sz: float):
        """Write AutoDock Vina grid configuration."""
        try:
            with open(self.grid_conf, "w") as f:
                f.write(f"""center_x = {cx}
center_y = {cy}
center_z = {cz}
size_x = {sx}
size_y = {sy}
size_z = {sz}
""")
            logger.info(f"[✔] Grid ready: {self.grid_conf}")
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

        def generate_box(): return ',\n       '.join(
            map(lambda x: ', '.join(x), (
                ('COLOR', '0.0', '1.0'),
                ('LINEWIDTH', '3.0'),
                ('BEGIN', 'LINES'),
                *(
                    vertex
                    for edge in edges
                    for vertex in map(
                        lambda points: (
                            'VERTEX',
                            f'{points[0]:7.3f}',
                            f'{points[1]:7.3f}',
                            f'{points[2]:7.3f}',
                        ),
                        (vertices[edge[0]], vertices[edge[1]])
                    )
                ),
                ('END',),
            ))
        )

        pml_lines = (
            f'load {pdb_relpath}, receptor',
            r'show cartoon, receptor',
            r'',
            r'python',
            r'from pymol import cmd',
            r'from pymol.cgo import *',
            r'',
            f'box = [{generate_box()}]',
            r'cmd.load_cgo(box, "docking_grid")',
            r'python end',
            r'',
            r'zoom all',
            r''
        )

        try:
            with open(self.grid_box_script, "w") as f:
                f.write('\n'.join(pml_lines))
                logger.info(
                    f"[✔] Grid visualization script saved: {self.grid_box_script}")
        except IOError as e:
            logger.warning(f"Could not write grid visualization script: {e}")

    # Backward-compatible alias (deprecated): the former private name.
    _protein_centroid_grid = centroid_grid
