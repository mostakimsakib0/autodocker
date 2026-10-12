#!/usr/bin/env python3
"""Docking box geometry and the files that describe the box: grid.conf
(AutoDock Vina), grid_box.pml (PyMOL) and frame_offset.txt."""
import logging
import os
from typing import Sequence, Tuple

from vspipeline import utils

logger = logging.getLogger(__name__)

# (center x, y, z, size x, y, z), Angstrom.
Box = Tuple[float, float, float, float, float, float]


def bounding_box(xs: Sequence[float], ys: Sequence[float], zs: Sequence[float],
                 padding: float, min_size: float = 0.0) -> Box:
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


def pocket_box(pocket_files: Sequence[str], padding: float) -> Box:
    """Grid centre and size from one or more pocket PDB files: the
    midpoint of the pocket atoms' bounding box, ``padding`` on each side."""
    xs, ys, zs = [], [], []
    for pocket_file in pocket_files:
        px, py, pz = utils.parse_pdb_coords(pocket_file)
        xs.extend(px)
        ys.extend(py)
        zs.extend(pz)

    cx, cy, cz, sx, sy, sz = bounding_box(xs, ys, zs, padding)
    logger.debug(f"Grid center: ({cx:.2f}, {cy:.2f}, {cz:.2f}), "
                 f"size: ({sx:.2f}, {sy:.2f}, {sz:.2f})")
    return cx, cy, cz, sx, sy, sz


def write_grid_conf(path: str, cx: float, cy: float, cz: float,
                    sx: float, sy: float, sz: float) -> None:
    """Write the AutoDock Vina grid configuration (raises on a write error)."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"""center_x = {cx}
center_y = {cy}
center_z = {cz}
size_x = {sx}
size_y = {sy}
size_z = {sz}
""")


def write_frame_offset(path: str, dx: float, dy: float, dz: float) -> None:
    """Write frame_offset.txt (raises on a write error)."""
    with open(path, "w", encoding="utf-8") as f:
        f.write(
            "# receptor.pdbqt, grid.conf and docked poses are translated so\n"
            "# the mean of the receptor heavy atoms is at the origin. Add this\n"
            "# offset to their coordinates to overlay them on the input PDB, e.g. with\n"
            "#   python3 scripts/restore_frame.py <output_dir> docked/*_out.pdbqt\n"
            f"offset_x = {dx}\noffset_y = {dy}\noffset_z = {dz}\n")


def write_grid_box_script(path: str, receptor_pdbqt: str, cx: float, cy: float, cz: float,
                          sx: float, sy: float, sz: float) -> None:
    """Write a PyMOL script that shows ``receptor_pdbqt`` and the docking
    box (a write error is only logged)."""
    _dir = os.path.dirname(os.path.abspath(path))
    pdb_relpath = os.path.relpath(receptor_pdbqt, _dir)

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
        with open(path, "w", encoding="utf-8") as f:
            f.write('\n'.join(pml_lines))
            logger.info(
                f"[*] Grid visualization script saved: {path}")
    except IOError as e:
        logger.warning(f"Could not write grid visualization script: {e}")
