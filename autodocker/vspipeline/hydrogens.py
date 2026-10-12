#!/usr/bin/env python3
"""Deterministic hydrogen positions for Open Babel output.

Open Babel adds each hydrogen along ``OBBuilder::GetNewBondVector``. When
the local geometry does not fix that direction (no reference atom two bonds
away, or a hybridisation/degree combination the builder has no rule for) it
calls ``vector3::randomUnitVector()``. Open Babel 3.1 seeds that generator
from the clock (``OBRandom::TimeSeed``: ``srand(usec % 70092)``), so the
hydrogen lands elsewhere on every run although heavy atoms, atom order,
types and charges are identical. Open Babel >= 3.2 seeds it from the
environment variable ``OB_RANDOM_SEED`` (random when unset; upstream PR
#2241): with a fixed seed the output is repeatable, but such a hydrogen is
still placed by a random draw, not by the geometry, and moves with the seed.
Seen with ``obabel -p`` on receptors (waters, atoms of broken residues,
two-atom fragments) and with ``obabel -h`` on ligands read from MOL2 (e.g.
the amino N of the CASF-2016 1nc3 ligand).

``stabilize_hydrogens`` compares two conversions of the same input, finds the
hydrogens that differ (plus the structurally undefined groups) and re-places
every hydrogen of those heavy atoms from heavy-atom geometry only, keeping
each one within AutoDock Vina's bond distance of its own heavy atom and of no
other. Only coordinates change: atom order, serials, names, types, charges
and PDBQT torsion-tree records are untouched, so a ligand hydrogen stays in
the branch of its heavy atom. ``run_obabel_deterministic`` wraps an Open
Babel call with this check. The main conversion runs with the run's seed S
and the probe conversions with *different* seeds (S+1, then S+2 when the
first two agree; see ``probe_seed``): a hydrogen that moves when only the
seed changes was placed at random. (Two conversions with the same seed
would agree on Open Babel >= 3.2 and hide it.) With Open Babel 3.1, which
ignores the variable, the probes differ by their clock seeds instead.

Supported files: PDB/PDBQT (ATOM/HETATM records) and V2000 SDF/MOL. Other
formats are accepted only if the conversions are byte-identical.
"""
import logging
import math
import os
from collections import defaultdict
from typing import Callable, List, Optional, Sequence, Tuple

from vspipeline import defaults, utils

logger = logging.getLogger(__name__)

# Hybridisation from the mean bond angle at the parent atom: roughly the
# midpoints between the ideal angles of sp3 (109.47), sp2 (120) and sp (180)
# centres. Kept at the exact values the hydrogen placement was validated with.
SP2_SP3_ANGLE_CUTOFF = 114.7  # degrees
SP_SP2_ANGLE_CUTOFF = 150.0   # degrees

PDBQT_HYDROGEN_TYPES = frozenset({"H", "HD", "HS"})

# Open Babel covalent radii (src/elementtable.h) and the bond tolerance of
# OBMol::ConnectTheDots, used to reproduce Open Babel's bond perception.
OB_COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07,
    "S": 1.05, "CL": 1.02, "BR": 1.20, "I": 1.39, "SE": 1.20, "NA": 1.66,
    "K": 2.03, "MG": 1.41, "CA": 1.76, "MN": 1.39, "FE": 1.32, "CO": 1.26,
    "NI": 1.24, "CU": 1.32, "ZN": 1.22, "CD": 1.44, "HG": 1.32,
}
OB_UNKNOWN_COVALENT_RADIUS = 1.6
OB_BOND_TOLERANCE = 0.45
OB_MIN_BOND_DISTANCE = 0.4

# AutoDock Vina covalent radii (src/lib/atom_constants.h) and the 1.1
# allowance of model::assign_bonds; Vina calls an N/O a donor when it is
# bonded, by this distance rule, to an HD atom.
VINA_COVALENT_RADII = {
    "H": 0.37, "C": 0.77, "N": 0.75, "O": 0.73, "P": 1.06, "S": 1.02,
    "F": 0.71, "I": 1.33, "CL": 0.99, "BR": 1.14, "SI": 1.11, "MG": 1.30,
    "MN": 1.39, "ZN": 1.31, "CA": 1.74, "FE": 1.25,
}
VINA_METAL_COVALENT_RADIUS = 1.75
VINA_BOND_ALLOWANCE = 1.1

# A hydrogen's parent is the nearest heavy atom within this distance (A).
H_PARENT_MAX_DISTANCE = 2.0

_PDBQT_TYPE_ELEMENTS = {
    "A": "C", "NA": "N", "NS": "N", "OA": "O", "OS": "O", "SA": "S",
    "HD": "H", "HS": "H",
}



class SpatialHash:
    """Points bucketed in cubic cells; ``near`` returns every point within
    one cell edge of the query (and some farther ones)."""

    def __init__(self, cell: float):
        self.cell = cell
        self.cells = defaultdict(list)

    def _cell(self, xyz):
        return tuple(int(math.floor(c / self.cell)) for c in xyz)

    def add(self, xyz, item):
        self.cells[self._cell(xyz)].append((xyz, item))

    def near(self, xyz):
        cx, cy, cz = self._cell(xyz)
        return [entry for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
                for entry in self.cells.get((cx + dx, cy + dy, cz + dz), ())]


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

_TETRAHEDRAL = math.radians(109.4712206)
_SPIN_STEP_DEGREES = 30


def pdbqt_element(atom_type: str) -> str:
    return _PDBQT_TYPE_ELEMENTS.get(atom_type, atom_type.upper())


def _vsub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vadd(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _vscale(a, s):
    return (a[0] * s, a[1] * s, a[2] * s)


def _vdot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _vcross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _vnorm(a):
    length = math.sqrt(_vdot(a, a))
    return _vscale(a, 1.0 / length) if length > 1e-9 else None


_LAB_AXES = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))


def _perpendicular(axis):
    """Unit vector perpendicular to ``axis`` built from the lab axis least
    parallel to it (ties: x before y before z)."""
    ref = min(_LAB_AXES, key=lambda e: abs(_vdot(e, axis)))
    return _vnorm(_vsub(ref, _vscale(axis, _vdot(ref, axis))))


def _cone(axis, half_angle, n, spin_degrees):
    """``n`` unit vectors at ``half_angle`` from ``axis``, 120 deg apart."""
    b = _perpendicular(axis)
    c = _vcross(axis, b)
    out = []
    for k in range(n):
        phi = math.radians(spin_degrees + 120.0 * k)
        side = _vadd(_vscale(b, math.cos(phi)), _vscale(c, math.sin(phi)))
        out.append(_vnorm(_vadd(_vscale(axis, math.cos(half_angle)),
                                _vscale(side, math.sin(half_angle)))))
    return out


# Directions sampled on the unit sphere when searching free space around an
# isolated atom for its hydrogens.
_SPHERE_POINTS = 256


def _fibonacci_sphere(n=_SPHERE_POINTS):
    golden = math.pi * (3.0 - math.sqrt(5.0))
    points = []
    for i in range(n):
        z = 1.0 - 2.0 * (i + 0.5) / n
        r = math.sqrt(max(0.0, 1.0 - z * z))
        points.append((r * math.cos(golden * i), r * math.sin(golden * i), z))
    return points


_SPHERE = _fibonacci_sphere()


def _spread_directions(bonds, n):
    """Greedily add ``n`` directions that keep the largest angle to the
    existing bond directions (fixed sampling of the sphere)."""
    bonds = list(bonds)
    out = []
    for _ in range(n):
        best = min(_SPHERE, key=lambda p: max((_vdot(p, b) for b in bonds), default=-1.0))
        out.append(best)
        bonds.append(best)
    return out


def _candidate_h_directions(bonds: List[Tuple[float, float, float]], n_h: int):
    """Candidate sets of ``n_h`` hydrogen directions, in a fixed order, for an
    atom whose heavy-atom bonds point along the unit vectors ``bonds``.

    Hydrogens take the free tetrahedral positions; where those leave a free
    rotation (no or one heavy neighbour) every 30 deg spin is a candidate.
    Hydrogens beyond the free tetrahedral positions (over-coordinated atoms)
    are spread away from the other bonds.
    """
    deg = len(bonds)
    candidates = []
    if deg == 0:
        for axis in _LAB_AXES + tuple(_vscale(a, -1.0) for a in _LAB_AXES):
            for spin in range(0, 360, _SPIN_STEP_DEGREES):
                slots = [axis] + _cone(axis, _TETRAHEDRAL, 3, spin)
                candidates.append(slots)
    elif deg == 1:
        away = _vscale(bonds[0], -1.0)
        for spin in range(0, 360, _SPIN_STEP_DEGREES):
            candidates.append(_cone(away, math.pi - _TETRAHEDRAL, 3, spin))
    elif deg == 2:
        bisector = _vnorm(_vscale(_vadd(bonds[0], bonds[1]), -1.0)) or _perpendicular(bonds[0])
        normal = _vnorm(_vcross(bonds[0], bonds[1])) or _vcross(bonds[0], bisector)
        half = _TETRAHEDRAL / 2.0
        up = _vnorm(_vadd(_vscale(bisector, math.cos(half)), _vscale(normal, math.sin(half))))
        down = _vnorm(_vsub(_vscale(bisector, math.cos(half)), _vscale(normal, math.sin(half))))
        candidates = [[up, down], [down, up]]
    elif deg == 3:
        candidates = [[_vnorm(_vscale(_vadd(_vadd(bonds[0], bonds[1]), bonds[2]), -1.0))
                       or _perpendicular(bonds[0])]]
    else:
        candidates = [[]]
    result = []
    for slots in candidates:
        chosen = list(slots[:n_h])
        if len(chosen) < n_h:
            chosen += _spread_directions(list(bonds) + chosen, n_h - len(chosen))
        result.append(chosen)
    return result



def vina_bond_cutoff(element: str) -> float:
    """Distance below which Vina bonds a hydrogen to an atom of ``element``."""
    radius = VINA_COVALENT_RADII.get(element, VINA_METAL_COVALENT_RADIUS)
    return VINA_BOND_ALLOWANCE * (VINA_COVALENT_RADII["H"] + radius)


def _ob_bonded(el_a, xyz_a, el_b, xyz_b) -> bool:
    cutoff = (OB_COVALENT_RADII.get(el_a, OB_UNKNOWN_COVALENT_RADIUS)
              + OB_COVALENT_RADII.get(el_b, OB_UNKNOWN_COVALENT_RADIUS) + OB_BOND_TOLERANCE)
    return OB_MIN_BOND_DISTANCE < math.dist(xyz_a, xyz_b) < cutoff


# ---------------------------------------------------------------------------
# Coordinate files
# ---------------------------------------------------------------------------

def file_format(path: str) -> Optional[str]:
    """'pdbqt', 'pdb' or 'sdf' from the file extension, else None."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".pdbqt":
        return "pdbqt"
    if ext in (".pdb", ".ent"):
        return "pdb"
    if ext in (".sdf", ".sd", ".mol"):
        return "sdf"
    return None


def _pdb_element(line: str) -> str:
    element = "".join(c for c in line[76:78] if c.isalpha()).upper()
    if element:
        return element
    name = line[12:16]
    if line.startswith("ATOM") or (name[0] == "H" and len(name.strip()) == 4):
        return name.strip().lstrip("0123456789")[:1].upper()
    return "".join(c for c in name[:2] if c.isalpha()).upper()


class _Atoms:
    """Atom records of a coordinate file: (line index, element, is_h, xyz,
    rest of the record, molecule index); writes coordinates back in place."""

    def __init__(self, path: str, fmt: str):
        self.fmt = fmt
        with open(path, encoding="utf-8") as f:
            self.lines = f.readlines()
        self.records = []
        if fmt in ("pdbqt", "pdb"):
            self._parse_pdb()
        elif fmt == "sdf":
            self._parse_sdf(path)
        else:
            raise ValueError(f"Unsupported coordinate format: {fmt}")

    def _parse_pdb(self):
        for i, line in enumerate(self.lines):
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if self.fmt == "pdbqt":
                atom_type = line[77:79].strip()
                is_h = atom_type in PDBQT_HYDROGEN_TYPES
                element = pdbqt_element(atom_type)
            else:
                element = _pdb_element(line)
                is_h = element in ("H", "D")
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            self.records.append((i, element, is_h, xyz, line[:30] + line[54:], 0))

    def _parse_sdf(self, path):
        lines = self.lines
        i = molecule = 0
        while i + 3 < len(lines):
            counts = lines[i + 3]
            if "V3000" in counts:
                raise ValueError(f"V3000 molfiles are not supported: {path}")
            n_atoms = int(counts[0:3])
            for k in range(i + 4, i + 4 + n_atoms):
                line = lines[k]
                element = line[31:34].strip().upper()
                xyz = (float(line[0:10]), float(line[10:20]), float(line[20:30]))
                self.records.append((k, element, element in ("H", "D"), xyz,
                                     line[30:], molecule))
            j = i + 4 + n_atoms
            while j < len(lines) and not lines[j].startswith("$$$$"):
                j += 1
            i, molecule = j + 1, molecule + 1

    def set_xyz(self, line_index: int, xyz) -> None:
        line = self.lines[line_index]
        x, y, z = xyz
        if self.fmt == "sdf":
            self.lines[line_index] = f"{x:10.4f}{y:10.4f}{z:10.4f}{line[30:]}"
        else:
            self.lines[line_index] = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"

    def write(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as f:
            f.writelines(self.lines)


# ---------------------------------------------------------------------------
# Ligand placement: keep Open Babel's hybridisation, fix the rotation
# ---------------------------------------------------------------------------

# OBAtom::CorrectedBondRad scales the covalent radius of sp2 and sp atoms.
_BOND_SCALE = {1: 0.90, 2: 0.95, 3: 1.0}


def _hybridisation(xyz, nbr_xyz, h_xyz) -> int:
    """1, 2 or 3 from the mean H-X-neighbour angle Open Babel used (the angle
    is deterministic even when the rotation is not)."""
    angles = []
    for h in h_xyz:
        u = _vnorm(_vsub(h, xyz))
        for n in nbr_xyz:
            v = _vnorm(_vsub(n, xyz))
            if u and v:
                angles.append(math.degrees(math.acos(max(-1.0, min(1.0, _vdot(u, v))))))
    if not angles:
        return 3
    mean = sum(angles) / len(angles)
    if mean > SP_SP2_ANGLE_CUTOFF:
        return 1
    return 2 if mean > SP2_SP3_ANGLE_CUTOFF else 3


def _ligand_candidates(xyz, q_xyz, ref_xyz, hyb, n_h):
    """Candidate direction sets for the hydrogens of an atom at ``xyz`` with
    one heavy neighbour at ``q_xyz``: sp2 in the plane of the neighbour's own
    neighbour ``ref_xyz`` (conjugation plane), sp3 staggered with the first
    hydrogen anti to it, sp linear; then every 30 deg rotation."""
    away = _vnorm(_vsub(xyz, q_xyz))
    side = None
    if ref_xyz is not None:
        r = _vsub(ref_xyz, q_xyz)
        side = _vnorm(_vsub(r, _vscale(away, _vdot(r, away))))
    side = side or _perpendicular(away)
    third = _vcross(away, side)
    if hyb == 1:
        return [[away] * n_h]
    if hyb == 2:
        half, slots, start = math.radians(60.0), 2, 0
    else:
        half, slots, start = math.pi - _TETRAHEDRAL, 3, 180
    spacing = 360.0 / slots
    out = []
    for step in range(0, 360, _SPIN_STEP_DEGREES):
        dirs = []
        for k in range(slots):
            phi = math.radians(start + step + spacing * k)
            s = _vadd(_vscale(side, math.cos(phi)), _vscale(third, math.sin(phi)))
            dirs.append(_vnorm(_vadd(_vscale(away, math.cos(half)), _vscale(s, math.sin(half)))))
        chosen = dirs[:n_h]
        if len(chosen) < n_h:
            chosen += _spread_directions([_vscale(away, -1.0)] + chosen, n_h - len(chosen))
        out.append(chosen)
    return out


# ---------------------------------------------------------------------------
# Detection and re-placement
# ---------------------------------------------------------------------------

def stabilize_hydrogens(path: str, probe_path: Optional[str] = None,
                        kind: str = "receptor", fmt: Optional[str] = None) -> int:
    """Give Open Babel's randomly placed hydrogens a deterministic position.

    ``probe_path`` is a second conversion of the same input. A hydrogen group
    (all hydrogens of one heavy atom) is re-placed when any of its hydrogens
    differs between the two files, when its heavy atom has no heavy
    neighbour and two or more hydrogens, or when it belongs to a two-atom
    fragment. Everything except hydrogen coordinates must be identical in
    both files, else RuntimeError. ``fmt`` defaults to the format of
    ``path``'s extension (see ``file_format``).

    ``kind="receptor"``: hydrogens take free tetrahedral positions at the
    sum of Open Babel covalent radii; where those can rotate, the 30 deg
    rotation keeping them farthest from other heavy atoms wins.
    ``kind="ligand"``: an atom with one heavy neighbour keeps the
    hybridisation Open Babel gave it (from its H-X-Y angle) and Open Babel's
    X-H length for it; sp2 hydrogens lie in the plane of the neighbour's
    lowest-numbered other neighbour, sp3 ones are staggered with the first
    anti to it, and the first such rotation that keeps every hydrogen out of
    Vina's bond distance of other heavy atoms is used. Other ligand groups
    are placed as for receptors.

    Only coordinate columns change. Returns the number of hydrogens moved.
    """
    if kind not in ("receptor", "ligand"):
        raise ValueError(f"kind must be 'receptor' or 'ligand', not {kind!r}")
    fmt = fmt or file_format(path)
    if fmt is None:
        raise ValueError(f"Cannot tell the coordinate format of {path}")
    atoms = _Atoms(path, fmt)
    probe = _Atoms(probe_path, fmt) if probe_path is not None else None
    if probe is not None and len(probe.records) != len(atoms.records):
        raise RuntimeError(
            f"Open Babel wrote {len(atoms.records)} and {len(probe.records)} atoms in "
            "two conversions of the same input")

    heavy, hydrogens = [], []
    for k, (i, element, is_h, xyz, sig, mol) in enumerate(atoms.records):
        moved, alt = False, xyz
        if probe is not None:
            other_i, _, _, other_xyz, other_sig, _ = probe.records[k]
            moved = other_xyz != xyz
            if other_sig != sig or (moved and not is_h):
                raise RuntimeError(
                    "Open Babel wrote different heavy atoms in two conversions of "
                    f"the same input:\n  {atoms.lines[i].rstrip()}\n  "
                    f"{probe.lines[other_i].rstrip()}")
            alt = other_xyz
        if is_h:
            hydrogens.append((i, xyz, alt, moved, mol))
        else:
            heavy.append((i, element, xyz, mol))

    if not hydrogens:
        return 0
    max_bond = 2 * max(OB_COVALENT_RADII.values()) + OB_BOND_TOLERANCE
    grid = SpatialHash(max_bond)
    for h in range(len(heavy)):
        grid.add(heavy[h][2], h)

    # Parent: the heavy atom closest to the hydrogen in both conversions (a
    # randomly placed hydrogen may pass close to another atom in one of them).
    groups = defaultdict(list)
    for i, xyz, alt, moved, mol in hydrogens:
        best = None
        for p, h in grid.near(xyz):
            if heavy[h][3] != mol:
                continue
            d = max(math.dist(p, xyz), math.dist(p, alt))
            if d <= H_PARENT_MAX_DISTANCE and (best is None or d < best[0]):
                best = (d, h)
        if best is None:
            raise RuntimeError(
                f"Hydrogen at {xyz} has no heavy atom within "
                f"{H_PARENT_MAX_DISTANCE} A in {path}")
        groups[best[1]].append((i, moved, xyz))

    def heavy_neighbours(h):
        _, el, xyz, mol = heavy[h]
        return [n for p, n in grid.near(xyz)
                if n != h and heavy[n][3] == mol and _ob_bonded(el, xyz, heavy[n][1], p)]

    moved_total = 0
    for h in sorted(groups):
        members = groups[h]
        nbrs = heavy_neighbours(h)
        undefined = ((not nbrs and len(members) >= 2)
                     or (len(nbrs) == 1 and heavy_neighbours(nbrs[0]) == [h]))
        if not (undefined or any(moved for _, moved, _ in members)):
            continue
        _, el, xyz, mol = heavy[h]
        others = [(heavy[n][1], p) for p, n in grid.near(xyz)
                  if n != h and heavy[n][3] == mol]

        def clearance_of(positions):
            value = min((math.dist(pos, p) - vina_bond_cutoff(e)
                         for pos in positions for e, p in others), default=math.inf)
            return round(value, 4)

        radius = OB_COVALENT_RADII.get(el, OB_UNKNOWN_COVALENT_RADIUS)
        best = None
        if kind == "ligand" and len(nbrs) == 1 and not undefined:
            q = nbrs[0]
            hyb = _hybridisation(xyz, [heavy[q][2]], [m[2] for m in members])
            length = OB_COVALENT_RADII["H"] + radius * _BOND_SCALE[hyb]
            refs = [n for n in heavy_neighbours(q) if n != h]
            ref_xyz = heavy[min(refs)][2] if refs else None
            for directions in _ligand_candidates(xyz, heavy[q][2], ref_xyz, hyb, len(members)):
                positions = [_vadd(xyz, _vscale(d, length)) for d in directions]
                clearance = clearance_of(positions)
                if best is None or clearance > best[0]:
                    best = (clearance, positions)
                if clearance > 0:
                    break
        else:
            length = radius + OB_COVALENT_RADII["H"]
            bonds = [_vnorm(_vsub(heavy[n][2], xyz)) for n in nbrs]
            bonds = [b for b in bonds if b is not None]
            for directions in _candidate_h_directions(bonds, len(members)):
                positions = [_vadd(xyz, _vscale(d, length)) for d in directions]
                clearance = clearance_of(positions)
                if best is None or clearance > best[0]:
                    best = (clearance, positions)
        clearance, positions = best
        if clearance <= 0:
            logger.warning(
                f"[!] A hydrogen of the {el} atom at {xyz} is within Vina's bond "
                "distance of another heavy atom in every orientation; check the "
                "input geometry there")
        for (i, _, _), pos in zip(members, positions):
            atoms.set_xyz(i, pos)
            moved_total += 1

    atoms.write(path)
    return moved_total


# ---------------------------------------------------------------------------
# Open Babel wrapper
# ---------------------------------------------------------------------------

# Seeds are positive 32-bit integers (see runner.SEED_RANGE); the probe seeds
# wrap around inside that range.
_SEED_MAX = 2 ** 31 - 1


def probe_seed(seed: Optional[int], k: int) -> int:
    """OB_RANDOM_SEED of the k-th probe conversion (k >= 1): ``seed + k``,
    wrapped into 1..2**31-1, so never equal to ``seed`` for k < 2**31-1."""
    base = defaults.SEED if seed is None else int(seed)
    return (base - 1 + k) % _SEED_MAX + 1


def same_bytes(a: str, b: str) -> bool:
    """Byte-for-byte comparison of two files. (filecmp.cmp caches results by
    size and mtime, so a probe rewritten within one timestamp tick with the
    same size could be reported equal to its earlier content.)"""
    with open(a, "rb") as fa, open(b, "rb") as fb:
        return fa.read() == fb.read()


def _probe_path(output: str) -> str:
    root, ext = os.path.splitext(output)
    return f"{root}.h-probe{ext}"


def _default_obabel(args: Sequence[str], outputs: Sequence[str] = (),
                    seed: Optional[int] = None) -> None:
    """One Open Babel conversion through ``utils.run_obabel``."""
    utils.run_obabel(list(args), outputs=list(outputs), seed=seed)


def run_obabel_deterministic(args: Sequence[str], output: str, kind: str = "ligand",
                             run_obabel: Optional[Callable] = None,
                             seed: Optional[int] = None) -> int:
    """Run an Open Babel conversion two or three times and make ``output``'s
    hydrogen positions deterministic.

    ``args`` are Open Babel arguments (without the executable) that write
    ``output``; the path must appear in ``args`` exactly once (after -O).
    ``run_obabel(args, outputs=[...], seed=...)`` performs one conversion
    with ``OB_RANDOM_SEED`` = ``seed`` (default: ``utils.run_obabel``, which
    runs obabel with OMP_NUM_THREADS=1). ``output`` is written
    with ``seed`` (default ``defaults.SEED``); the second (and third)
    conversion writes a temporary sibling file (the probe) with the seed
    ``probe_seed(seed, 1)`` (and ``probe_seed(seed, 2)``).

    Open Babel >= 3.2 draws a different random direction for each seed, so
    a randomly placed hydrogen differs between the output and the first
    probe. Open Babel 3.1 ignores OB_RANDOM_SEED and seeds from the clock
    (``usec % 70092``), so two runs can draw the same seed and place a
    random hydrogen identically. When the first two files are identical a
    third conversion is therefore compared as well; only if it is identical
    too is the output accepted as it is. With Open Babel 3.1 a randomly
    placed hydrogen can then survive only if all three runs drew the same
    clock seed: with run-to-run timing jitter far above one microsecond,
    about 1 in 70092**2 (~2e-10) per conversion, against ~1.4e-5 with two
    runs. If the files differ, ``stabilize_hydrogens``
    re-places the hydrogens that differ (and the structurally undefined
    groups); if anything else differs (heavy atoms, e.g. from ``--gen3d``)
    ``output`` is deleted and RuntimeError raised. Returns the number of
    hydrogens moved.

    The structurally undefined groups are deliberately not re-placed when
    all runs agree: Open Babel places some of them deterministically and
    better (e.g. the linear H of HCN), and an input without random hydrogens
    must keep the output of the first conversion unchanged.
    """
    args = list(args)
    if args.count(output) != 1:
        raise ValueError(f"Output path {output!r} must appear exactly once in {args}")
    run = run_obabel or _default_obabel
    seed = defaults.SEED if seed is None else int(seed)
    probe = _probe_path(output)
    probe_args = [probe if a == output else a for a in args]
    run(args, outputs=[output], seed=seed)
    try:
        run(probe_args, outputs=[probe], seed=probe_seed(seed, 1))
        if same_bytes(output, probe):
            # Open Babel 3.1 can draw the same clock seed twice (see above);
            # compare a third run.
            run(probe_args, outputs=[probe], seed=probe_seed(seed, 2))
            if same_bytes(output, probe):
                return 0
        fmt = file_format(output)
        if fmt is None:
            if os.path.exists(output):
                os.remove(output)
            raise RuntimeError(
                f"Open Babel wrote a different {output} in two runs with different "
                "random seeds; "
                "its format is not supported for hydrogen re-placement")
        try:
            moved = stabilize_hydrogens(output, probe, kind=kind, fmt=fmt)
        except (RuntimeError, ValueError):
            if os.path.exists(output):
                os.remove(output)
            raise
        if moved:
            logger.debug(f"Placed {moved} hydrogens of {output} deterministically")
        return moved
    finally:
        if os.path.exists(probe):
            os.remove(probe)
