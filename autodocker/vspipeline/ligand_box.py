"""Docking box from a reference ligand (--box-ligand, --box-from-het).

The box is built around a ligand known to bind the site, instead of around
an fpocket pocket: smina's --autobox_ligand, with AutoDocker's padding rule.

* Center: midpoint of the bounding box of the reference ligand's heavy atoms
  (as for fpocket pockets and the whole-receptor box, gridbox.bounding_box).
* Size, per axis: max(extent + 2 x padding, MIN_LIGAND_BOX_EDGE), with
  padding = --padding (defaults.PADDING, 6 A): max(22.5, extent + 12) A by
  default. See docs/features_2_0.md for the choice.
* Frame: the reference coordinates are in the input receptor PDB frame. The
  docking frame (receptor.pdbqt, grid.conf, docked poses) is that frame
  minus ProteinPreparation.frame_offset(), exactly as for
  --grid-center-frame input.

The reference comes from a file (SDF/MOL/MOL2/PDB/PDBQT, first molecule or
model) or from a HETATM group of the input receptor PDB
(RESNAME[:CHAIN[:RESNUM]]); the group is then removed from the receptor
(write_receptor_without_group), also with --keep-hetero.

Pure Python, no Open Babel: the box depends only on the coordinates in the
file. This module imports nothing that imports it, so docking, protein and
pipeline can all use it.
"""
import hashlib
import math
import os
import re
from collections import OrderedDict
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from vspipeline import defaults
from vspipeline.fpocket import rank_pockets
from vspipeline.gridbox import bounding_box
from vspipeline.pdb_records import is_hydrogen_record, read_pdb_records, record_chain

Coord = Tuple[float, float, float]

# Smallest box edge (Angstrom) around a reference ligand; the benchmarks'
# rule (scripts/_benchlib.py). A small reference (a fragment, a metal ion)
# still gives a box in which a drug-like library ligand (up to ~20 A long)
# can turn.
MIN_LIGAND_BOX_EDGE = defaults.LIGAND_BOX_MIN_SIZE
# Padding on each side: --padding (defaults.PADDING, the fpocket-pocket padding).
LIGAND_BOX_PADDING = defaults.PADDING
# Warn when no receptor heavy atom is this close (Angstrom) to the reference
# ligand: a bound ligand touches the receptor at ~3-4 A, so a larger gap
# usually means the ligand is in another frame (another PDB entry, not
# superposed) or near a chain that is not docked.
CONTACT_WARN_DISTANCE = 5.0
# A receptor HETATM group "is" the reference ligand when this fraction of the
# reference heavy atoms lies within SAME_ATOM_DISTANCE of one of its atoms.
SAME_ATOM_DISTANCE = 1.0
SAME_GROUP_FRACTION = 0.5

LIGAND_EXTENSIONS = {".sdf": "sdf", ".mol": "sdf", ".mol2": "mol2",
                     ".pdb": "pdb", ".ent": "pdb", ".pdbqt": "pdbqt"}
WATER_RESNAMES = frozenset({"HOH", "WAT", "H2O", "DOD", "SOL", "TIP", "T3P", "T4P", "SPC"})
_PDBQT_HYDROGEN_TYPES = frozenset({"H", "HD", "HS"})

SOURCE_FILE = "ligand_file"   # --box-ligand
SOURCE_HET = "het_group"      # --box-from-het


class LigandBoxError(ValueError):
    """A reference ligand that cannot give a docking box (the runner reports
    it as a usage error)."""


# ---------------------------------------------------------------------------
# Reading reference ligands
# ---------------------------------------------------------------------------

def ligand_format(path: str) -> str:
    """'sdf', 'mol2', 'pdb' or 'pdbqt' from the file extension."""
    ext = os.path.splitext(path)[1].lower()
    if ext not in LIGAND_EXTENSIONS:
        raise LigandBoxError(
            f"--box-ligand {path}: unsupported format {ext or '(no extension)'}; use "
            + ", ".join(sorted(LIGAND_EXTENSIONS)))
    return LIGAND_EXTENSIONS[ext]


def _is_hydrogen_symbol(symbol: str) -> bool:
    return symbol.strip().capitalize() in ("H", "D")


def _sdf_heavy_atoms(lines: List[str], path: str) -> Tuple[List[Coord], int]:
    """Heavy atoms of the first record of an SDF/MOL file (V2000 or V3000);
    returns (coordinates, number of records)."""
    n_records = sum(1 for line in lines if line.startswith("$$$$")) or 1
    if len(lines) < 4:
        raise LigandBoxError(f"--box-ligand {path}: not a MOL/SDF file (no counts line)")
    counts = lines[3]
    coords: List[Coord] = []
    if "V3000" in counts:
        in_atoms = False
        for line in lines[4:]:
            if line.startswith("$$$$"):
                break
            text = line.strip()
            if text == "M  V30 BEGIN ATOM":
                in_atoms = True
            elif text == "M  V30 END ATOM":
                break
            elif in_atoms and text.startswith("M  V30"):
                parts = text.split()
                # M V30 index symbol x y z aamap ...
                try:
                    symbol, x, y, z = parts[3], float(parts[4]), float(parts[5]), float(parts[6])
                except (IndexError, ValueError):
                    raise LigandBoxError(f"--box-ligand {path}: malformed V3000 atom "
                                         f"line: {text}") from None
                if not _is_hydrogen_symbol(symbol):
                    coords.append((x, y, z))
        return coords, n_records
    try:
        n_atoms = int(counts[0:3])
    except ValueError:
        raise LigandBoxError(f"--box-ligand {path}: not a MOL/SDF file (bad counts "
                             f"line {counts.rstrip()!r})") from None
    for line in lines[4:4 + n_atoms]:
        try:
            x, y, z = float(line[0:10]), float(line[10:20]), float(line[20:30])
        except ValueError:
            raise LigandBoxError(f"--box-ligand {path}: malformed atom line: "
                                 f"{line.rstrip()}") from None
        if not _is_hydrogen_symbol(line[31:34]):
            coords.append((x, y, z))
    return coords, n_records


def _mol2_heavy_atoms(lines: List[str], path: str) -> Tuple[List[Coord], int]:
    """Heavy atoms of the first molecule of a MOL2 file."""
    n_records = sum(1 for line in lines if line.startswith("@<TRIPOS>MOLECULE")) or 1
    coords: List[Coord] = []
    in_atoms = seen_atoms = False
    for line in lines:
        if line.startswith("@<TRIPOS>"):
            if seen_atoms:
                break
            in_atoms = line.startswith("@<TRIPOS>ATOM")
            seen_atoms = in_atoms
            continue
        if not in_atoms or not line.strip():
            continue
        parts = line.split()
        try:
            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
            atom_type = parts[5]
        except (IndexError, ValueError):
            raise LigandBoxError(f"--box-ligand {path}: malformed MOL2 atom line: "
                                 f"{line.rstrip()}") from None
        if not _is_hydrogen_symbol(atom_type.split(".")[0]):
            coords.append((x, y, z))
    return coords, n_records


def _pdb_heavy_atoms(lines: List[str], pdbqt: bool) -> Tuple[List[Coord], int]:
    """Heavy ATOM/HETATM atoms of the first model of a PDB/PDBQT file."""
    n_models = sum(1 for line in lines if line.startswith("MODEL")) or 1
    coords: List[Coord] = []
    for line in lines:
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if pdbqt:
            if line[77:79].strip() in _PDBQT_HYDROGEN_TYPES:
                continue
        elif is_hydrogen_record(line):
            continue
        try:
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError:
            raise LigandBoxError(f"malformed coordinates: {line.rstrip()}") from None
    return coords, n_models


def read_ligand_heavy_atoms(path: str) -> Tuple[List[Coord], int]:
    """Heavy-atom coordinates of a reference ligand file and the number of
    molecules/models in it (only the first is used). Raises LigandBoxError
    for unreadable files, no heavy atoms or 2D coordinates."""
    fmt = ligand_format(path)
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
    except OSError as e:
        raise LigandBoxError(f"--box-ligand: cannot read {path}: {e}") from None
    try:
        if fmt == "sdf":
            coords, n = _sdf_heavy_atoms(lines, path)
        elif fmt == "mol2":
            coords, n = _mol2_heavy_atoms(lines, path)
        else:
            coords, n = _pdb_heavy_atoms(lines, pdbqt=fmt == "pdbqt")
    except LigandBoxError as e:
        message = str(e)
        if not message.startswith("--box-ligand"):
            message = f"--box-ligand {path}: {message}"
        raise LigandBoxError(message) from None
    if not coords:
        raise LigandBoxError(f"--box-ligand {path}: no heavy atoms in the first "
                             f"{'molecule' if fmt in ('sdf', 'mol2') else 'model'}")
    if not all(math.isfinite(v) for c in coords for v in c):
        raise LigandBoxError(f"--box-ligand {path}: non-finite coordinates")
    if len(coords) > 1 and all(c[2] == 0.0 for c in coords):
        raise LigandBoxError(
            f"--box-ligand {path}: the ligand has 2D coordinates (every z is 0); the box "
            "needs the ligand's 3D position in the receptor's frame")
    return coords, n


# ---------------------------------------------------------------------------
# HETATM groups of the receptor PDB
# ---------------------------------------------------------------------------

@dataclass(frozen=True, order=True)
class HetGroup:
    """One HETATM residue: chain (blank read as A), residue number with
    insertion code, residue name."""
    chain: str
    resid: str
    resname: str

    @property
    def spec(self) -> str:
        """The --box-from-het value that selects exactly this group."""
        return f"{self.resname}:{self.chain}:{self.resid}"

    @property
    def is_water(self) -> bool:
        return self.resname in WATER_RESNAMES


@dataclass(frozen=True)
class HetSpec:
    """A parsed --box-from-het value; None fields match anything."""
    resname: str
    chain: Optional[str] = None
    resid: Optional[str] = None

    def matches(self, group: HetGroup) -> bool:
        return (group.resname == self.resname
                and (self.chain is None or group.chain == self.chain)
                and (self.resid is None or group.resid.upper() == self.resid.upper()))

    def __str__(self) -> str:
        parts = [self.resname, self.chain or "", self.resid or ""]
        while len(parts) > 1 and not parts[-1]:
            parts.pop()
        return ":".join(parts)


_RESNAME_RE = re.compile(r"[A-Za-z0-9]{1,3}")
_RESID_RE = re.compile(r"-?[0-9]{1,4}[A-Za-z]?")


def parse_het_spec(text: str) -> HetSpec:
    """Parse RESNAME[:CHAIN[:RESNUM]] (e.g. MK1, MK1:B, MK1:B:902, MK1:B:902A
    with an insertion code, MK1::902 for any chain). The residue name is
    matched case-insensitively, the chain ID exactly (as --chain)."""
    raw = (text or "").strip()
    parts = raw.split(":")
    if not raw or len(parts) > 3:
        raise LigandBoxError(
            f"--box-from-het must be RESNAME[:CHAIN[:RESNUM]] (e.g. MK1, MK1:B or "
            f"MK1:B:902), got {text!r}")
    resname = parts[0].strip().upper()
    chain = parts[1].strip() if len(parts) > 1 else ""
    resid = parts[2].strip() if len(parts) > 2 else ""
    if not _RESNAME_RE.fullmatch(resname):
        raise LigandBoxError(
            f"--box-from-het: residue name {parts[0]!r} must be 1-3 letters or digits "
            "(the PDB residue name, columns 18-20)")
    if len(chain) > 1:
        raise LigandBoxError(f"--box-from-het: chain {chain!r} must be one character")
    if resid and not _RESID_RE.fullmatch(resid):
        raise LigandBoxError(
            f"--box-from-het: residue number {resid!r} must be a number with an optional "
            "insertion code (e.g. 902 or 902A)")
    if resname in WATER_RESNAMES:
        raise LigandBoxError(f"--box-from-het {raw}: {resname} is water, not a ligand")
    return HetSpec(resname, chain or None, resid or None)


def _record_group(line: str) -> HetGroup:
    return HetGroup(record_chain(line), line[22:27].strip(), line[17:20].strip().upper())


def het_groups(pdb_file: str) -> "OrderedDict[HetGroup, List[Coord]]":
    """HETATM groups of ``pdb_file`` in file order, with their heavy-atom
    coordinates: first model, one alternate location per residue (the
    receptor-cleaning rule, pdb_records.read_pdb_records). A group of hydrogens
    only has no coordinates."""
    lines, _ = read_pdb_records(pdb_file)
    groups: "OrderedDict[HetGroup, List[Coord]]" = OrderedDict()
    for line in lines:
        if not line.startswith("HETATM") or len(line) < 54:
            continue
        coords = groups.setdefault(_record_group(line), [])
        if is_hydrogen_record(line):
            continue
        try:
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except ValueError:
            raise LigandBoxError(f"{pdb_file}: malformed coordinates: {line.rstrip()}") from None
    return groups


def describe_het_groups(groups: Dict[HetGroup, List[Coord]], limit: int = 20) -> str:
    """'MK1:B:902 (45 heavy atoms), ...; and 127 water residues' for errors."""
    ligands = [g for g in groups if not g.is_water]
    waters = len(groups) - len(ligands)
    if not ligands:
        text = "no HETATM group other than water"
    else:
        shown = [f"{g.spec} ({len(groups[g])} heavy atoms)" for g in ligands[:limit]]
        if len(ligands) > limit:
            shown.append(f"... {len(ligands) - limit} more")
        text = ", ".join(shown)
    if waters:
        text += f"; and {waters} water residue(s)"
    return text


def select_het_group(pdb_file: str, spec: str) -> Tuple[HetGroup, List[Coord]]:
    """The one HETATM group of ``pdb_file`` that ``spec`` selects, with its
    heavy-atom coordinates. Raises LigandBoxError listing the available
    groups when none or several match, or when the group has no heavy atom."""
    het = parse_het_spec(spec)
    try:
        groups = het_groups(pdb_file)
    except (OSError, UnicodeDecodeError) as e:
        raise LigandBoxError(f"--box-from-het: cannot read {pdb_file}: {e}") from None
    matches = [g for g in groups if het.matches(g)]
    if not matches:
        hint = ""
        same_name = [g for g in groups if g.resname == het.resname]
        if same_name:
            hint = " (" + het.resname + " is present as " + ", ".join(
                g.spec for g in same_name) + ")"
        raise LigandBoxError(
            f"--box-from-het {het}: no such HETATM group in {pdb_file}{hint}. "
            f"Available: {describe_het_groups(groups)}")
    if len(matches) > 1:
        raise LigandBoxError(
            f"--box-from-het {het} matches {len(matches)} HETATM groups in {pdb_file}: "
            + ", ".join(g.spec for g in matches)
            + f"; add the chain and residue number, e.g. --box-from-het {matches[0].spec}")
    group = matches[0]
    if not groups[group]:
        raise LigandBoxError(f"--box-from-het {group.spec}: the group has no heavy atoms")
    return group, groups[group]


def write_receptor_without_group(src: str, dst: str, group: HetGroup) -> int:
    """Copy ``src`` to ``dst`` without the HETATM (and ANISOU) records of
    ``group``, in every model and alternate location; every other line is
    kept byte for byte. Receptor preparation then reads ``dst``, so the
    reference ligand is never part of the receptor (also with --keep-hetero,
    --detect-cofactors). Returns the number of records removed."""
    removed = 0
    out = []
    with open(src, encoding="utf-8", newline="") as f:
        for line in f:
            if (line.startswith(("HETATM", "ANISOU")) and len(line) >= 27
                    and _record_group(line) == group):
                removed += 1
                continue
            out.append(line)
    if not removed:
        raise LigandBoxError(f"{group.spec}: no HETATM records to remove from {src}")
    with open(dst, "w", encoding="utf-8", newline="") as f:
        f.writelines(out)
    return removed


def write_reference_pdb(path: str, src: str, group: HetGroup) -> None:
    """The group's HETATM records (first model, chosen altloc) as a PDB file,
    in the input frame, for viewing next to the input receptor."""
    lines, _ = read_pdb_records(src)
    records = [line if line.endswith("\n") else line + "\n" for line in lines
               if line.startswith("HETATM") and len(line) >= 27
               and _record_group(line) == group]
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"REMARK   AutoDocker --box-from-het {group.spec}\n")
        f.writelines(records)
        f.write("END\n")


# ---------------------------------------------------------------------------
# Box geometry and frames
# ---------------------------------------------------------------------------

def to_docking_frame(point: Sequence[float], offset: Sequence[float]) -> Coord:
    """Input-PDB coordinates -> docking frame (receptor.pdbqt, grid.conf):
    subtract ProteinPreparation.frame_offset(), as --grid-center-frame input."""
    return (point[0] - offset[0], point[1] - offset[1], point[2] - offset[2])


def to_input_frame(point: Sequence[float], offset: Sequence[float]) -> Coord:
    """Docking frame -> input-PDB coordinates (inverse of to_docking_frame)."""
    return (point[0] + offset[0], point[1] + offset[1], point[2] + offset[2])


def ligand_box(coords: Sequence[Coord], padding: float = LIGAND_BOX_PADDING,
               min_edge: float = MIN_LIGAND_BOX_EDGE) -> Tuple[Coord, Coord]:
    """(center, size) of the box around ``coords``: bounding-box midpoint,
    per-axis edge max(extent + 2 x padding, min_edge)."""
    if not coords:
        raise LigandBoxError("no reference coordinates for the docking box")
    if padding < 0 or min_edge < 0:
        raise LigandBoxError("padding and minimum edge must be >= 0")
    xs, ys, zs = zip(*coords)
    cx, cy, cz, sx, sy, sz = bounding_box(xs, ys, zs, padding, min_edge)
    return (cx, cy, cz), (sx, sy, sz)


def extent(coords: Sequence[Coord]) -> Coord:
    xs, ys, zs = zip(*coords)
    return (max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))


def points_in_box(coords: Sequence[Coord], center: Sequence[float],
                  size: Sequence[float]) -> int:
    """How many of ``coords`` lie inside the box (edges included)."""
    half = [s / 2.0 for s in size]
    return sum(1 for c in coords
               if all(abs(c[i] - center[i]) <= half[i] + 1e-9 for i in range(3)))


def min_distance(a: Sequence[Coord], b: Sequence[Coord]) -> float:
    """Smallest distance between a point of ``a`` and a point of ``b``
    (inf when either is empty). Brute force: a ligand against a receptor is
    ~50 x 10^4 pairs."""
    best = math.inf
    for p in a:
        for q in b:
            d = (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2
            if d < best:
                best = d
    return math.sqrt(best)


def matching_het_groups(coords: Sequence[Coord], groups: Dict[HetGroup, List[Coord]],
                        distance: float = SAME_ATOM_DISTANCE,
                        fraction: float = SAME_GROUP_FRACTION) -> List[HetGroup]:
    """Non-water HETATM groups that occupy the reference ligand's place: at
    least ``fraction`` of its heavy atoms within ``distance`` of a group atom
    (the reference was taken from the receptor's own entry). Used to warn
    that --keep-hetero keeps that ligand in the receptor."""
    found = []
    d2 = distance * distance
    for group, atoms in groups.items():
        if group.is_water or not atoms:
            continue
        close = sum(1 for p in coords if any(
            (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2 <= d2 for q in atoms))
        if close >= fraction * len(coords):
            found.append(group)
    return found


# ---------------------------------------------------------------------------
# fpocket transparency: which pocket holds the reference ligand
# ---------------------------------------------------------------------------

def read_alpha_spheres(vert_pqr: str) -> List[Tuple[Coord, float]]:
    """(center, radius) of the alpha spheres in fpocket's pocketN_vert.pqr."""
    spheres = []
    with open(vert_pqr, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            parts = line.split()
            try:
                x, y, z, radius = (float(parts[-5]), float(parts[-4]), float(parts[-3]),
                                   float(parts[-1]))
            except (IndexError, ValueError):
                continue
            spheres.append(((x, y, z), radius))
    return spheres


def pocket_overlap(coords: Sequence[Coord], spheres: Sequence[Tuple[Coord, float]]
                   ) -> float:
    """Fraction of the reference heavy atoms inside at least one alpha sphere."""
    if not coords:
        return 0.0
    inside = sum(1 for p in coords if any(
        (p[0] - c[0]) ** 2 + (p[1] - c[1]) ** 2 + (p[2] - c[2]) ** 2 <= r * r
        for c, r in spheres))
    return inside / len(coords)


def pocket_overlaps(coords: Sequence[Coord], pockets: List[Dict],
                    center: Sequence[float]) -> List[Dict]:
    """Overlap of the reference ligand (input frame) with every fpocket pocket
    (dicts as ProteinPreparation._collect_pockets builds them; their files are
    in the input frame too). Returns one entry per pocket in fpocket rank
    order (fpocket.rank_pockets: score, then druggability): number, rank,
    score, druggability score, fraction of reference heavy atoms inside the
    pocket's alpha spheres, and the distance from the pocket's sphere
    centroid to the box center."""
    rows = []
    for rank, p in enumerate(rank_pockets(pockets), 1):
        vert = p["path"][:-len("_atm.pdb")] + "_vert.pqr" if p["path"].endswith(
            "_atm.pdb") else None
        spheres = read_alpha_spheres(vert) if vert and os.path.exists(vert) else []
        row = {"pocket": p["number"], "rank": rank, "score": p.get("score"),
               "druggability_score": p.get("druggability_score"),
               "overlap_fraction": None, "center_distance": None}
        if spheres:
            n = float(len(spheres))
            centroid = tuple(sum(s[0][i] for s in spheres) / n for i in range(3))
            row["overlap_fraction"] = round(pocket_overlap(coords, spheres), 4)
            row["center_distance"] = round(math.dist(centroid, center), 3)
        rows.append(row)
    return rows


def best_overlapping_pocket(rows: List[Dict]) -> Optional[Dict]:
    """The pocket holding most of the reference ligand (ties: better fpocket
    rank), or None when no pocket overlaps it."""
    hits = [r for r in rows if r["overlap_fraction"]]
    if not hits:
        return None
    return min(hits, key=lambda r: (-r["overlap_fraction"], r["rank"]))


def fpocket_summary(rows: List[Dict], limit: int = 5) -> Dict:
    """run_manifest.json entry: the top-ranked pocket, the best-overlapping
    one, and whether fpocket's default choice (rank 1) would have been it."""
    best = best_overlapping_pocket(rows)
    top = rows[0] if rows else None
    return {
        "pockets_found": len(rows),
        "top_ranked_pocket": top,
        "best_overlapping_pocket": best,
        "top_pocket_contains_reference": bool(best and top and best["pocket"] == top["pocket"]),
        "overlapping_pockets": [r for r in rows if r["overlap_fraction"]][:limit],
        "overlap_definition": ("fraction of the reference ligand's heavy atoms inside "
                               "the pocket's fpocket alpha spheres (pocketN_vert.pqr)"),
    }


# ---------------------------------------------------------------------------
# The reference box of a run
# ---------------------------------------------------------------------------

def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _rounded(values: Sequence[float], digits: int = 3) -> List[float]:
    return [round(float(v), digits) for v in values]


@dataclass(frozen=True)
class ReferenceBox:
    """A docking box built around a reference ligand, in the input frame."""
    source: str                       # SOURCE_FILE or SOURCE_HET
    reference: str                    # file path, or the group's spec
    coords: Tuple[Coord, ...]         # reference heavy atoms, input frame
    center: Coord                     # input frame
    size: Coord
    padding: float
    min_edge: float
    sha256: Optional[str] = None      # of the --box-ligand file
    group: Optional[HetGroup] = None  # --box-from-het
    n_records: int = 1                # molecules/models in the file (first used)

    @property
    def extent(self) -> Coord:
        return extent(self.coords)

    def docking_center(self, offset: Sequence[float]) -> Coord:
        return to_docking_frame(self.center, offset)

    def describe(self) -> str:
        what = (f"--box-ligand {self.reference}" if self.source == SOURCE_FILE
                else f"--box-from-het {self.reference} (HETATM group of the receptor)")
        return f"{what}, {len(self.coords)} heavy atoms"

    def log_lines(self, offset: Sequence[float]) -> List[str]:
        """The INFO lines the runner logs for the box."""
        ex = self.extent
        dc = self.docking_center(offset)
        lines = [
            f"[*] Docking box: reference ligand {self.describe()}",
            f"[*] Reference ligand extent: {ex[0]:.2f}, {ex[1]:.2f}, {ex[2]:.2f} A; box = "
            f"max(extent + 2 x {self.padding:g}, {self.min_edge:g}) A per axis",
            f"[*] Grid center (input-PDB frame): {self.center[0]:.2f}, {self.center[1]:.2f}, "
            f"{self.center[2]:.2f}",
            f"[*] Grid center (docking frame): {dc[0]:.2f}, {dc[1]:.2f}, {dc[2]:.2f}",
            f"[*] Grid size: {self.size[0]:.2f}, {self.size[1]:.2f}, {self.size[2]:.2f}",
        ]
        if self.n_records > 1:
            lines.append(f"[*] {self.reference} holds {self.n_records} molecules/models; "
                         "the box uses the first")
        return lines

    def manifest(self, offset: Sequence[float], grid_size_override: Optional[Coord] = None,
                 fpocket: Optional[Dict] = None,
                 receptor_contact: Optional[float] = None) -> Dict:
        """The run_manifest.json "box" entry."""
        size = grid_size_override or self.size
        entry = {
            "source": self.source,
            "reference": self.reference,
            "reference_sha256": self.sha256,
            "het_group": (None if self.group is None else
                          {"resname": self.group.resname, "chain": self.group.chain,
                           "resid": self.group.resid, "spec": self.group.spec}),
            "reference_heavy_atoms": len(self.coords),
            "reference_records_in_file": self.n_records,
            "reference_extent": _rounded(self.extent),
            "rule": "center = midpoint of the reference heavy atoms' bounding box; "
                    "size = max(extent + 2 x padding, min_edge) per axis",
            "padding": self.padding,
            "min_edge": self.min_edge,
            "center_input_frame": _rounded(self.center),
            "center_docking_frame": _rounded(self.docking_center(offset)),
            "frame_offset": _rounded(offset),
            "size_computed": _rounded(self.size),
            "size": _rounded(size),
            "size_source": "--grid-size" if grid_size_override else "reference ligand",
            "reference_atoms_in_box": points_in_box(self.coords, self.center, size),
            "closest_receptor_heavy_atom": (None if receptor_contact is None
                                            else round(receptor_contact, 3)),
            "fpocket": fpocket,
        }
        return entry


def box_from_ligand_file(path: str, padding: float = LIGAND_BOX_PADDING,
                         min_edge: float = MIN_LIGAND_BOX_EDGE) -> ReferenceBox:
    """--box-ligand: box around the first molecule/model of ``path``."""
    if not os.path.isfile(path):
        raise LigandBoxError(f"--box-ligand: no such file: {path}")
    coords, n = read_ligand_heavy_atoms(path)
    center, size = ligand_box(coords, padding, min_edge)
    return ReferenceBox(SOURCE_FILE, path, tuple(coords), center, size, padding,
                        min_edge, sha256=_sha256(path), n_records=n)


def box_from_het(pdb_file: str, spec: str, padding: float = LIGAND_BOX_PADDING,
                 min_edge: float = MIN_LIGAND_BOX_EDGE) -> ReferenceBox:
    """--box-from-het: box around a HETATM group of the input receptor PDB."""
    group, coords = select_het_group(pdb_file, spec)
    center, size = ligand_box(coords, padding, min_edge)
    return ReferenceBox(SOURCE_HET, group.spec, tuple(coords), center, size, padding,
                        min_edge, group=group)


def check_options(box_ligand: Optional[str], box_from_het: Optional[str],
                  pockets: Optional[str], grid_center: Optional[str]) -> None:
    """Combination rules of the reference-ligand box options (raises
    LigandBoxError). With a reference box:

    * --box-ligand and --box-from-het exclude each other;
    * --grid-center is refused (the reference gives the center);
      --grid-center-frame then is too (runner: it needs --grid-center);
    * --pockets is refused (no pocket selects the box);
    * --grid-size is accepted and replaces the computed size (center from
      the reference, edge fixed, e.g. a benchmark's 22.5 A cube);
    * --no-fpocket is accepted: it skips the fpocket run that otherwise only
      reports which pocket holds the reference (pocket_summary.csv, manifest);
    * --padding sets the padding of the rule.
    """
    if box_ligand is None and box_from_het is None:
        return
    if box_ligand is not None and box_from_het is not None:
        raise LigandBoxError("--box-ligand and --box-from-het cannot be combined "
                             "(give one reference ligand)")
    option = "--box-ligand" if box_ligand is not None else "--box-from-het"
    if grid_center is not None:
        raise LigandBoxError(f"{option} cannot be combined with --grid-center (the box "
                             "center is the reference ligand's; --grid-size may set the size)")
    if pockets is not None:
        raise LigandBoxError(f"{option} cannot be combined with --pockets (the box comes "
                             "from the reference ligand, not from an fpocket pocket)")
    if box_ligand is not None:
        if not box_ligand.strip():
            raise LigandBoxError("--box-ligand must not be empty")
        ligand_format(box_ligand)
        if not os.path.isfile(box_ligand):
            raise LigandBoxError(f"--box-ligand: no such file: {box_ligand}")
    else:
        parse_het_spec(box_from_het)
