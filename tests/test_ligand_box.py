"""Docking box from a reference ligand (vspipeline.ligand_box).

Fixtures (tests/data/ligand_box, see make_fixtures.py): a binding-site
fragment of PDB 1HSG with its crystal ligand MK1 (indinavir, B 902), and MK1
alone as PDB, SDF, MOL2 and PDBQT.
"""
import math
import os
import shutil

import pytest

from vspipeline import ligand_box as lb
from vspipeline import protein

DATA = os.path.join(os.path.dirname(__file__), "data", "ligand_box")
SITE = os.path.join(DATA, "1hsg_site.pdb")
MK1_FILES = [os.path.join(DATA, f"mk1.{ext}") for ext in ("pdb", "sdf", "mol2", "pdbqt")]
# The paper's 1HSG receptor (MK1 and waters removed), when available locally.
PAPER_RECEPTOR = os.path.expanduser(
    "~/Downloads/Research/AutoDocker_Paper/rerun/v115_image/protein.pdb")
OBABEL = shutil.which("obabel")


def _mk1_coords():
    coords = []
    with open(os.path.join(DATA, "mk1.pdb"), encoding="utf-8") as f:
        for line in f:
            if line.startswith("HETATM"):
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return coords


def _het(serial, name, resname, chain, resnum, x, y, z, element="C", altloc=" ",
         icode=" ", occ=1.0):
    return (f"HETATM{serial:5d} {name:<4s}{altloc}{resname:>3s} {chain}{resnum:>4s}{icode}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f} 20.00          {element:>2s}\n")


def _atom(serial, x, y, z, chain="A", resnum=1):
    return (f"ATOM  {serial:5d}  CA  ALA {chain}{resnum:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C\n")


# ---------------------------------------------------------------------------
# Reading reference ligands
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", MK1_FILES, ids=lambda p: os.path.basename(p))
def test_every_format_gives_the_crystal_heavy_atoms(path):
    coords, n = lb.read_ligand_heavy_atoms(path)
    assert n == 1
    expected = _mk1_coords()
    assert len(coords) == len(expected) == 45
    assert sorted(coords) == pytest.approx(sorted(expected), abs=1e-3)


def test_hydrogens_are_ignored(tmp_path):
    sdf = tmp_path / "h.sdf"
    sdf.write_text(
        "h\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    1.0000 C   0  0\n"
        "    5.0000    0.0000    1.0000 O   0  0\n"
        "   50.0000   50.0000   50.0000 H   0  0\n"
        "  1  2  1  0\n  1  3  1  0\nM  END\n$$$$\n", encoding="utf-8")
    coords, _ = lb.read_ligand_heavy_atoms(str(sdf))
    assert coords == [(0.0, 0.0, 1.0), (5.0, 0.0, 1.0)]


def test_v3000_sdf(tmp_path):
    sdf = tmp_path / "v3.sdf"
    sdf.write_text(
        "v3\n  test\n\n  0  0  0     0  0            999 V3000\n"
        "M  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 0 0 0\nM  V30 BEGIN ATOM\n"
        "M  V30 1 C 1.0 2.0 3.0 0\nM  V30 2 N 4.0 5.0 6.0 0\nM  V30 3 H 9.0 9.0 9.0 0\n"
        "M  V30 END ATOM\nM  V30 END CTAB\nM  END\n$$$$\n", encoding="utf-8")
    coords, _ = lb.read_ligand_heavy_atoms(str(sdf))
    assert coords == [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)]


def test_multi_record_files_use_the_first(tmp_path):
    text = open(os.path.join(DATA, "mk1.sdf"), encoding="utf-8").read()
    shifted = text.replace("    9.2800   23.7630", "  109.2800   23.7630")
    path = tmp_path / "two.sdf"
    path.write_text(text + shifted, encoding="utf-8")
    coords, n = lb.read_ligand_heavy_atoms(str(path))
    assert n == 2 and max(c[0] for c in coords) < 100

    pdbqt = tmp_path / "poses.pdbqt"
    body = open(os.path.join(DATA, "mk1.pdbqt"), encoding="utf-8").read()
    pdbqt.write_text("MODEL 1\n" + body + "ENDMDL\nMODEL 2\n"
                     "ATOM      1  C   UNL     1     500.000 500.000 500.000  0.00  0.00"
                     "     0.000 C \nENDMDL\n", encoding="utf-8")
    coords, n = lb.read_ligand_heavy_atoms(str(pdbqt))
    assert n == 2 and len(coords) == 45


@pytest.mark.parametrize("name,content,message", [
    ("lig.xyz", "1\n\nC 0 0 0\n", "unsupported format"),
    ("flat.sdf", "f\n  t\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n"
                 "    0.0000    0.0000    0.0000 C   0  0\n"
                 "    1.5000    0.0000    0.0000 C   0  0\n  1  2  1  0\nM  END\n",
     "2D coordinates"),
    ("onlyh.pdb", "HETATM    1  H1  LIG A   1       1.000   1.000   1.000  1.00  0.00"
                  "           H\n", "no heavy atoms"),
    ("bad.sdf", "x\n", "not a MOL/SDF"),
])
def test_unusable_reference_files(tmp_path, name, content, message):
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")
    with pytest.raises(lb.LigandBoxError, match=message):
        lb.read_ligand_heavy_atoms(str(path))


# ---------------------------------------------------------------------------
# Box rule and frames
# ---------------------------------------------------------------------------

def test_box_rule_padding_and_minimum_edge():
    coords = [(0.0, 0.0, 0.0), (20.0, 4.0, 1.0)]
    center, size = lb.ligand_box(coords, padding=6.0, min_edge=22.5)
    assert center == (10.0, 2.0, 0.5)
    assert size == (32.0, 22.5, 22.5)          # 20 + 12; max(4 + 12, 22.5); ...
    assert lb.ligand_box(coords, padding=4.0, min_edge=0.0)[1] == (28.0, 12.0, 9.0)
    assert lb.ligand_box([(1.0, 2.0, 3.0)])[1] == (22.5, 22.5, 22.5)   # a single ion
    with pytest.raises(lb.LigandBoxError):
        lb.ligand_box([])


def test_default_rule_is_the_benchmark_rule():
    """Defaults: padding = defaults.PADDING (6 A), minimum edge 22.5 A, the
    rule of scripts/_benchlib.ligand_box: max(22.5, extent + 12)."""
    from vspipeline import defaults
    assert lb.LIGAND_BOX_PADDING == defaults.PADDING == 6.0
    assert lb.MIN_LIGAND_BOX_EDGE == 22.5


@pytest.mark.parametrize("path", MK1_FILES, ids=lambda p: os.path.basename(p))
def test_box_contains_the_crystal_ligand_with_padding(path):
    box = lb.box_from_ligand_file(path)
    mk1 = _mk1_coords()
    assert lb.points_in_box(mk1, box.center, box.size) == len(mk1)
    for i in range(3):
        lo = min(c[i] for c in mk1) - (box.center[i] - box.size[i] / 2)
        hi = (box.center[i] + box.size[i] / 2) - max(c[i] for c in mk1)
        assert lo >= 6.0 - 1e-6 and hi >= 6.0 - 1e-6
        assert lo == pytest.approx(hi)        # bounding-box midpoint: symmetric margins
    assert box.sha256 and len(box.sha256) == 64
    assert box.source == lb.SOURCE_FILE


def test_frame_round_trip():
    offset = (12.5, -3.25, 7.0)
    point = (1.0, 2.0, 3.0)
    docking = lb.to_docking_frame(point, offset)
    assert docking == (1.0 - 12.5, 2.0 + 3.25, 3.0 - 7.0)
    assert lb.to_input_frame(docking, offset) == pytest.approx(point)


def test_docking_frame_matches_grid_center_frame_input(tmp_path):
    """The reference box center goes through the same conversion as
    --grid-center X --grid-center-frame input: X - frame_offset()."""
    prep = protein.ProteinPreparation(SITE, str(tmp_path))
    offset = prep.frame_offset()               # mean of the receptor heavy atoms
    box = lb.box_from_het(SITE, "MK1")
    cx, cy, cz = box.center
    dx, dy, dz = offset
    assert box.docking_center(offset) == (cx - dx, cy - dy, cz - dz)
    entry = box.manifest(offset)
    assert entry["center_input_frame"] == [round(v, 3) for v in box.center]
    assert entry["center_docking_frame"] == [round(c - d, 3) for c, d in zip(box.center, offset)]
    assert entry["reference_atoms_in_box"] == 45
    assert entry["het_group"] == {"resname": "MK1", "chain": "B", "resid": "902",
                                  "spec": "MK1:B:902"}


# ---------------------------------------------------------------------------
# HETATM groups
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("spec", ["MK1", "mk1", "MK1:B", "MK1:B:902", "MK1::902"])
def test_box_from_het_1hsg(spec):
    box = lb.box_from_het(SITE, spec)
    assert box.group == lb.HetGroup("B", "902", "MK1")
    assert box.reference == "MK1:B:902"
    file_box = lb.box_from_ligand_file(os.path.join(DATA, "mk1.pdb"))
    assert box.center == pytest.approx(file_box.center)
    assert box.size == pytest.approx(file_box.size)


@pytest.mark.parametrize("spec,message", [
    ("XYZ", r"no such HETATM group .*Available: MK1:B:902 \(45 heavy atoms\); and "
            r"\d+ water residue"),
    ("MK1:A", r"MK1 is present as MK1:B:902"),
    ("MK1:B:901", r"no such HETATM group"),
    ("HOH", r"is water"),
    ("MK1:B:902:X", r"RESNAME\[:CHAIN\[:RESNUM\]\]"),
    ("TOOLONG", r"1-3 letters"),
    ("MK1:BB", r"one character"),
    ("MK1:B:x9", r"residue number"),
    ("", r"RESNAME"),
])
def test_het_selection_errors(spec, message):
    with pytest.raises(lb.LigandBoxError, match=message):
        lb.box_from_het(SITE, spec)


def test_ambiguous_het_lists_the_candidates(tmp_path):
    pdb = tmp_path / "two.pdb"
    pdb.write_text(_atom(1, 0, 0, 0) + _het(2, "C1", "LIG", "A", "301", 1, 1, 1)
                   + _het(3, "C1", "LIG", "B", "301", 9, 9, 9) + "END\n", encoding="utf-8")
    with pytest.raises(lb.LigandBoxError, match=r"matches 2 HETATM groups.*LIG:A:301, "
                                               r"LIG:B:301.*--box-from-het LIG:A:301"):
        lb.box_from_het(str(pdb), "LIG")
    assert lb.box_from_het(str(pdb), "LIG:B").center == (9.0, 9.0, 9.0)


def test_insertion_codes_and_altlocs(tmp_path):
    """Insertion codes tell groups apart; of an altloc group only the
    higher-occupancy conformer is used (the receptor-cleaning rule)."""
    pdb = tmp_path / "alt.pdb"
    pdb.write_text(
        _atom(1, 0, 0, 0)
        + _het(2, "C1", "LIG", "A", "301", 1, 1, 1)
        + _het(3, "C1", "LIG", "A", "301", 30, 30, 30, icode="A")
        + _het(4, "C1", "LIG", "A", "302", 2, 0, 0, altloc="A", occ=0.3)
        + _het(5, "C2", "LIG", "A", "302", 4, 0, 0, altloc="A", occ=0.3)
        + _het(6, "C1", "LIG", "A", "302", 2, 10, 0, altloc="B", occ=0.7)
        + _het(7, "C2", "LIG", "A", "302", 4, 10, 0, altloc="B", occ=0.7)
        + "END\n", encoding="utf-8")
    assert lb.box_from_het(str(pdb), "LIG:A:301").center == (1.0, 1.0, 1.0)
    assert lb.box_from_het(str(pdb), "LIG:A:301A").center == (30.0, 30.0, 30.0)
    assert lb.box_from_het(str(pdb), "lig:A:301a").center == (30.0, 30.0, 30.0)
    box = lb.box_from_het(str(pdb), "LIG:A:302")
    assert box.coords == ((2.0, 10.0, 0.0), (4.0, 10.0, 0.0))


def test_blank_chain_is_chain_a(tmp_path):
    pdb = tmp_path / "blank.pdb"
    pdb.write_text(_atom(1, 0, 0, 0, chain=" ") + _het(2, "C1", "LIG", " ", "1", 3, 3, 3)
                   + "END\n", encoding="utf-8")
    assert lb.box_from_het(str(pdb), "LIG:A:1").group.chain == "A"


def test_receptor_without_group_keeps_everything_else(tmp_path):
    dst = tmp_path / "receptor_no_ref.pdb"
    removed = lb.write_receptor_without_group(SITE, str(dst), lb.HetGroup("B", "902", "MK1"))
    assert removed == 45
    src_lines = open(SITE, encoding="utf-8").read().splitlines(keepends=True)
    kept = [line for line in src_lines if not (line.startswith("HETATM") and "MK1" in line)]
    assert open(dst, encoding="utf-8").read() == "".join(kept)
    assert "MK1" not in lb.describe_het_groups(lb.het_groups(str(dst)))
    with pytest.raises(lb.LigandBoxError, match="no HETATM records"):
        lb.write_receptor_without_group(str(dst), str(tmp_path / "x.pdb"),
                                        lb.HetGroup("B", "902", "MK1"))


def test_receptor_without_group_removes_every_altloc_and_model(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(
        "MODEL        1\n" + _atom(1, 0, 0, 0)
        + _het(2, "C1", "LIG", "A", "5", 1, 1, 1, altloc="A")
        + "ANISOU    2  C1 ALIG A   5     100    100    100      0      0      0       C\n"
        + _het(3, "C1", "LIG", "A", "5", 1, 2, 1, altloc="B")
        + _het(4, "O", "HOH", "A", "6", 5, 5, 5, element="O")
        + "ENDMDL\nMODEL        2\n" + _atom(1, 0, 0, 0)
        + _het(2, "C1", "LIG", "A", "5", 1, 1, 1) + "ENDMDL\nEND\n", encoding="utf-8")
    dst = tmp_path / "out.pdb"
    assert lb.write_receptor_without_group(str(src), str(dst),
                                           lb.HetGroup("A", "5", "LIG")) == 4
    text = dst.read_text(encoding="utf-8")
    assert "LIG" not in text and "HOH" in text and text.count("ATOM  ") == 2


def test_reference_pdb_has_the_group_only(tmp_path):
    path = tmp_path / "reference_ligand.pdb"
    lb.write_reference_pdb(str(path), SITE, lb.HetGroup("B", "902", "MK1"))
    coords, _ = lb.read_ligand_heavy_atoms(str(path))
    assert sorted(coords) == sorted(_mk1_coords())


def test_reference_from_file_matches_receptor_het_group():
    """--box-ligand with the receptor's own ligand: the group is found (the
    runner warns that --keep-hetero would keep it in the receptor)."""
    coords, _ = lb.read_ligand_heavy_atoms(os.path.join(DATA, "mk1.sdf"))
    groups = lb.het_groups(SITE)
    assert lb.matching_het_groups(coords, groups) == [lb.HetGroup("B", "902", "MK1")]
    moved = [(x + 5.0, y, z) for x, y, z in coords]
    assert lb.matching_het_groups(moved, groups) == []


def test_receptor_contact_distance():
    mk1 = _mk1_coords()
    xs, ys, zs = protein._heavy_atom_coords(SITE)
    receptor = [c for c in zip(xs, ys, zs) if c not in set(mk1)]
    contact = lb.min_distance(mk1, receptor)
    assert 2.0 < contact < lb.CONTACT_WARN_DISTANCE
    far = [(x + 100.0, y, z) for x, y, z in mk1]
    assert lb.min_distance(far, receptor) > lb.CONTACT_WARN_DISTANCE
    assert lb.min_distance([], receptor) == math.inf


# ---------------------------------------------------------------------------
# Option rules
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kwargs,message", [
    (dict(box_ligand="mk1.sdf", box_from_het="MK1"), "cannot be combined"),
    (dict(box_ligand="mk1.sdf", grid_center="1,2,3"), "--grid-center"),
    (dict(box_from_het="MK1", grid_center="1,2,3"), "--grid-center"),
    (dict(box_from_het="MK1", pockets="1"), "--pockets"),
    (dict(box_ligand="missing.sdf"), "no such file"),
    (dict(box_ligand="mk1.cif"), "unsupported format"),
    (dict(box_from_het="MK1:B:1:2"), "RESNAME"),
])
def test_option_conflicts(kwargs, message):
    args = dict(box_ligand=None, box_from_het=None, pockets=None, grid_center=None)
    args.update(kwargs)
    if args["box_ligand"] == "mk1.sdf":
        args["box_ligand"] = os.path.join(DATA, "mk1.sdf")
    with pytest.raises(lb.LigandBoxError, match=message):
        lb.check_options(**args)


def test_accepted_option_combinations():
    lb.check_options(None, None, "1", "1,2,3")           # no reference box: not our rules
    lb.check_options(os.path.join(DATA, "mk1.sdf"), None, None, None)
    lb.check_options(None, "MK1:B:902", None, None)


def test_grid_size_override_in_manifest():
    box = lb.box_from_het(SITE, "MK1")
    entry = box.manifest((0.0, 0.0, 0.0), grid_size_override=(20.0, 20.0, 20.0))
    assert entry["size"] == [20.0, 20.0, 20.0] and entry["size_source"] == "--grid-size"
    assert entry["size_computed"] == [round(v, 3) for v in box.size]
    assert entry["reference_atoms_in_box"] <= 45


def test_log_lines_name_both_frames():
    box = lb.box_from_het(SITE, "MK1")
    text = "\n".join(box.log_lines((1.0, 2.0, 3.0)))
    assert "MK1:B:902" in text and "input-PDB frame" in text and "docking frame" in text
    assert "max(extent + 2 x 6, 22.5)" in text


# ---------------------------------------------------------------------------
# fpocket transparency
# ---------------------------------------------------------------------------

def _write_pocket(root, number, spheres):
    pockets = root / "pockets"
    pockets.mkdir(exist_ok=True)
    atm = pockets / f"pocket{number}_atm.pdb"
    atm.write_text(_atom(1, *spheres[0][0]), encoding="utf-8")
    lines = [f"ATOM  {i + 1:5d}    C STP  {number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
             f"    0.00   {r:5.2f}\n" for i, ((x, y, z), r) in enumerate(spheres)]
    (pockets / f"pocket{number}_vert.pqr").write_text("".join(lines), encoding="utf-8")
    return str(atm)


def test_pocket_overlap_names_the_pocket_holding_the_reference(tmp_path):
    box = lb.box_from_het(SITE, "MK1")
    mk1 = list(box.coords)
    centroid = tuple(sum(c[i] for c in mk1) / len(mk1) for i in range(3))
    pockets = [
        # fpocket's top pocket, 30 A away from the ligand
        {"number": 1, "score": 0.9, "druggability_score": 0.8,
         "path": _write_pocket(tmp_path, 1, [((centroid[0] + 30, centroid[1], centroid[2]),
                                              4.0)])},
        # pocket 2 has spheres over every MK1 atom
        {"number": 2, "score": 0.5, "druggability_score": 0.9,
         "path": _write_pocket(tmp_path, 2, [(c, 3.5) for c in mk1])},
        # pocket 3 covers part of the ligand
        {"number": 3, "score": 0.1, "druggability_score": 0.1,
         "path": _write_pocket(tmp_path, 3, [(mk1[0], 0.5)])},
    ]
    rows = lb.pocket_overlaps(mk1, pockets, box.center)
    assert [r["pocket"] for r in rows] == [1, 2, 3]
    assert rows[0]["overlap_fraction"] == 0.0 and rows[1]["overlap_fraction"] == 1.0
    assert 0 < rows[2]["overlap_fraction"] < 1
    summary = lb.fpocket_summary(rows)
    assert summary["best_overlapping_pocket"]["pocket"] == 2
    assert summary["top_ranked_pocket"]["pocket"] == 1
    assert summary["top_pocket_contains_reference"] is False
    assert [r["pocket"] for r in summary["overlapping_pockets"]] == [2, 3]
    assert lb.best_overlapping_pocket([dict(rows[0])]) is None


# ---------------------------------------------------------------------------
# Real receptor preparation (Open Babel)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(OBABEL is None, reason="obabel not available on PATH")
def test_prepared_receptor_has_no_reference_and_box_holds_it(tmp_path):
    """--box-from-het MK1 with --keep-hetero: the receptor is prepared from
    the copy without MK1, so receptor.pdbqt has no MK1 atom, and the box,
    moved into the docking frame, still holds every MK1 atom (moved the same
    way)."""
    out = tmp_path / "out"
    out.mkdir()
    box = lb.box_from_het(SITE, "MK1:B:902")
    stripped = str(out / "receptor_no_reference.pdb")
    lb.write_receptor_without_group(SITE, stripped, box.group)
    prep = protein.ProteinPreparation(stripped, str(out))
    prep.prepare_receptor("all", keep_hetero=True)
    with open(prep.receptor_pdbqt, encoding="utf-8") as f:
        receptor = f.read()
    assert "MK1" not in receptor and "HOH" in receptor     # other HETATMs kept
    offset = prep.frame_offset()
    center = box.docking_center(offset)
    mk1_docking = [lb.to_docking_frame(c, offset) for c in _mk1_coords()]
    assert lb.points_in_box(mk1_docking, center, box.size) == 45
    rec_xyz = [(float(line[30:38]), float(line[38:46]), float(line[46:54]))
               for line in receptor.splitlines() if line.startswith(("ATOM", "HETATM"))]
    assert lb.min_distance(mk1_docking, rec_xyz) < lb.CONTACT_WARN_DISTANCE


@pytest.mark.skipif(not os.path.exists(PAPER_RECEPTOR),
                    reason="the paper's 1HSG receptor is not available locally")
def test_paper_receptor_box_from_crystal_mk1():
    """The case-study receptor (1HSG without MK1): the box built from the
    crystal MK1 contains it and sits on the receptor (same frame)."""
    box = lb.box_from_ligand_file(os.path.join(DATA, "mk1.sdf"))
    mk1 = _mk1_coords()
    assert lb.points_in_box(mk1, box.center, box.size) == len(mk1)
    xs, ys, zs = protein._heavy_atom_coords(PAPER_RECEPTOR)
    receptor = list(zip(xs, ys, zs))
    assert lb.min_distance(mk1, receptor) < lb.CONTACT_WARN_DISTANCE
    assert lb.points_in_box(receptor, box.center, box.size) > 100
    assert lb.het_groups(PAPER_RECEPTOR) == {} or "MK1" not in lb.describe_het_groups(
        lb.het_groups(PAPER_RECEPTOR))
