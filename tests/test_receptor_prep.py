"""Receptor preparation: deterministic protonation, PDB record filtering,
fail-loud pocket detection and box geometry."""
import logging
import math
import os
import re

import pytest

import vspipeline.protein as prot
from vspipeline import tools, utils

DATA = os.path.join(os.path.dirname(__file__), "data")
# Residues 1-12 of the DUD-E AKT1 receptor (input hydrogens, blank chain and
# element columns) plus a broken lysine (CE-NZ), an isolated arginine NH1 and
# a water. Before the fix Open Babel placed ~10 of its hydrogens at random.
AKT1_FRAGMENT = os.path.join(DATA, "akt1_fragment_protonated.pdb")

needs_obabel = pytest.mark.skipif(tools.settings().obabel is None, reason="obabel not available")


def _pdb_atom(serial, name, resname, chain, resnum, x, y, z, element="",
              altloc=" ", occupancy=1.0, record="ATOM"):
    return (f"{record:<6s}{serial:5d} {name:<4s}{altloc}{resname:>3s} {chain:1s}"
            f"{resnum:>4s}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f} 20.00"
            f"          {element:>2s}\n")


def _protein(n=120, chain="A"):
    return [_pdb_atom(i + 1, "C", "ALA", chain, str(i // 5 + 1), float(i), float(i % 7),
                      float(i % 3), "C") for i in range(n)]


def _atoms(path):
    out = []
    for line in open(path, encoding="utf-8"):
        if line.startswith(("ATOM", "HETATM")):
            out.append((line[77:79].strip(),
                        (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                        line))
    return out


def _vina_bonded_heavy(h_xyz, heavy):
    """Heavy atoms within Vina's covalent bond distance of a hydrogen."""
    out = []
    for t, xyz, line in heavy:
        el = prot._pdbqt_element(t)
        if math.dist(h_xyz, xyz) < prot._vina_bond_cutoff(el):
            out.append(line)
    return out


# --- deterministic protonation ---------------------------------------------

@needs_obabel
def test_protonation_is_byte_identical_across_runs(tmp_path):
    outputs = []
    for run in ("run1", "run2", "run3"):
        out = tmp_path / run
        out.mkdir()
        pp = prot.ProteinPreparation(AKT1_FRAGMENT, str(out))
        pp.prepare_receptor("all", keep_hetero=True)
        outputs.append(open(pp.receptor_pdbqt, "rb").read())
    assert outputs[0] == outputs[1] == outputs[2]
    # No absolute (output-directory) path in the file.
    assert str(tmp_path).encode() not in outputs[0]
    assert b"REMARK  Name = protein_clean.pdb" in outputs[0]


@needs_obabel
def test_every_polar_hydrogen_bonded_to_exactly_one_heavy_atom(tmp_path):
    pp = prot.ProteinPreparation(AKT1_FRAGMENT, str(tmp_path))
    pp.prepare_receptor("all", keep_hetero=True)
    atoms = _atoms(pp.receptor_pdbqt)
    heavy = [a for a in atoms if a[0] not in prot.PDBQT_HYDROGEN_TYPES]
    for t, xyz, line in atoms:
        if t == "HD":
            assert len(_vina_bonded_heavy(xyz, heavy)) == 1, line


@needs_obabel
def test_input_hydrogens_are_replaced_not_stacked(tmp_path, caplog):
    """A pre-protonated lysine NZ must end with three hydrogens, not four."""
    pp = prot.ProteinPreparation(AKT1_FRAGMENT, str(tmp_path))
    with caplog.at_level(logging.WARNING):
        pp.prepare_receptor("all")
    assert "Removed" in caplog.text and "hydrogen records" in caplog.text
    assert not any(prot.is_hydrogen_record(l) for l in open(pp.pdb_clean, encoding="utf-8")
                   if l.startswith(("ATOM", "HETATM")))
    atoms = _atoms(pp.receptor_pdbqt)
    nz = [a for a in atoms if a[2][12:16] == " NZ " and a[2][22:26] == "  11"][0]
    n_h = sum(1 for t, xyz, _ in atoms if t == "HD" and math.dist(xyz, nz[1]) < 1.2)
    assert n_h == 3


def _pdbqt_line(serial, name, x, y, z, atom_type, resnum=1):
    return (f"ATOM  {serial:5d} {name:<4s} UNL A{resnum:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  0.00  0.00    +0.000 {atom_type:<2s}\n")


def test_stabilize_hydrogens_replaces_random_groups(tmp_path):
    # An isolated water and an N-H that differs between two conversions.
    common = [
        _pdbqt_line(1, "O", 0.0, 0.0, 0.0, "OA"),
        _pdbqt_line(4, "C", 10.0, 0.0, 0.0, "C", 2),
        _pdbqt_line(5, "C", 11.5, 0.0, 0.0, "C", 2),
        _pdbqt_line(6, "N", 12.2, 1.2, 0.0, "N", 2),
    ]
    first = common + [_pdbqt_line(2, "H", 0.97, 0.0, 0.0, "HD"),
                      _pdbqt_line(3, "H", -0.2, 0.9, 0.2, "HD"),
                      _pdbqt_line(7, "H", 12.9, 1.6, 0.4, "HD", 2)]
    second = common + [_pdbqt_line(2, "H", 0.97, 0.0, 0.0, "HD"),
                       _pdbqt_line(3, "H", -0.3, 0.1, 0.9, "HD"),
                       _pdbqt_line(7, "H", 12.0, 2.2, 0.1, "HD", 2)]
    a, b = tmp_path / "a.pdbqt", tmp_path / "b.pdbqt"
    a.write_text("".join(first))
    b.write_text("".join(second))
    a2, b2 = tmp_path / "a2.pdbqt", tmp_path / "b2.pdbqt"
    a2.write_text("".join(second))
    b2.write_text("".join(first))
    assert prot.stabilize_hydrogens(str(a), str(b)) == 3
    assert prot.stabilize_hydrogens(str(a2), str(b2)) == 3
    # Same result whichever random placement Open Babel produced first.
    assert a.read_text() == a2.read_text()
    atoms = _atoms(str(a))
    by_serial = {int(l[6:11]): xyz for _, xyz, l in atoms}
    assert math.dist(by_serial[2], by_serial[1]) == pytest.approx(0.97, abs=2e-3)
    assert math.dist(by_serial[7], by_serial[6]) == pytest.approx(1.02, abs=2e-3)


def test_stabilize_hydrogens_leaves_deterministic_groups(tmp_path):
    lines = [_pdbqt_line(1, "C", 0.0, 0.0, 0.0, "C"),
             _pdbqt_line(2, "C", 1.5, 0.0, 0.0, "C"),
             _pdbqt_line(3, "O", 2.2, 1.2, 0.0, "OA"),
             _pdbqt_line(4, "H", 3.1, 1.0, 0.3, "HD")]
    a, b = tmp_path / "a.pdbqt", tmp_path / "b.pdbqt"
    a.write_text("".join(lines))
    b.write_text("".join(lines))
    assert prot.stabilize_hydrogens(str(a), str(b)) == 0
    assert a.read_text() == "".join(lines)


def test_stabilize_hydrogens_rejects_different_heavy_atoms(tmp_path):
    a, b = tmp_path / "a.pdbqt", tmp_path / "b.pdbqt"
    a.write_text(_pdbqt_line(1, "C", 0.0, 0.0, 0.0, "C"))
    b.write_text(_pdbqt_line(1, "C", 0.1, 0.0, 0.0, "C"))
    with pytest.raises(RuntimeError, match="different heavy atoms"):
        prot.stabilize_hydrogens(str(a), str(b))


# --- PDB records -----------------------------------------------------------

@pytest.mark.parametrize("line,expected", [
    (_pdb_atom(1, "HA", "ALA", "A", "1", 0, 0, 0), True),          # blank element
    (_pdb_atom(1, "HNZ1", "LYS", "A", "1", 0, 0, 0), True),
    (_pdb_atom(1, "1HD2", "LEU", "A", "1", 0, 0, 0), True),
    (_pdb_atom(1, "CA", "ALA", "A", "1", 0, 0, 0), False),
    (_pdb_atom(1, "HG", "HG", "A", "1", 0, 0, 0, record="HETATM"), False),  # mercury
    (_pdb_atom(1, "D1", "ALA", "A", "1", 0, 0, 0, "D"), True),
    (_pdb_atom(1, "N", "ALA", "A", "1", 0, 0, 0, "N"), False),
])
def test_is_hydrogen_record(line, expected):
    assert prot.is_hydrogen_record(line) is expected


def test_altloc_one_label_per_residue_highest_occupancy():
    lines = [
        _pdb_atom(1, "CA", "SER", "A", "5", 0, 0, 0, "C", "A", 0.30),
        _pdb_atom(2, "CA", "SER", "A", "5", 0, 0, 1, "C", "B", 0.70),
        _pdb_atom(3, "CB", "SER", "A", "5", 1, 0, 0, "C", "A", 0.30),
        _pdb_atom(4, "CB", "SER", "A", "5", 1, 0, 1, "C", "B", 0.70),
        # An atom present only in conformer A must not be mixed in.
        _pdb_atom(5, "OG", "SER", "A", "5", 2, 0, 0, "O", "A", 0.30),
        _pdb_atom(6, "N", "SER", "A", "5", 3, 0, 0, "N"),
    ]
    kept, dropped = prot.select_altlocs(lines)
    assert dropped == 3
    assert [l[12:16].strip() for l in kept] == ["CA", "CB", "N"]
    assert all(l[16] == " " for l in kept)
    assert float(kept[0][46:54]) == 1.0  # conformer B


def test_altloc_tie_takes_first_label_and_resolves_microheterogeneity():
    lines = [
        _pdb_atom(1, "CB", "SER", "A", "7", 0, 0, 0, "C", "A", 0.50),
        _pdb_atom(2, "OG", "SER", "A", "7", 0, 0, 1, "O", "A", 0.50),
        _pdb_atom(3, "CB", "THR", "A", "7", 1, 0, 0, "C", "B", 0.50),
        _pdb_atom(4, "OG1", "THR", "A", "7", 1, 0, 1, "O", "B", 0.50),
    ]
    kept, dropped = prot.select_altlocs(lines)
    assert dropped == 2
    assert {l[17:20] for l in kept} == {"SER"}


def test_clean_pdb_first_model_selected_chain_ter(tmp_path):
    lines = (["MODEL        1\n"] + _protein(60, "A") + ["TER\n"] + _protein(60, "B")
             + ["TER\n", "ENDMDL\n", "MODEL        2\n"] + _protein(60, "A") + ["ENDMDL\n"])
    p = tmp_path / "in.pdb"
    p.write_text("".join(lines))
    pp = prot.ProteinPreparation(str(p), str(tmp_path / "out"))
    os.makedirs(pp.outdir)
    pp._write_clean_pdb(["A"], keep_hetero=False)
    clean = open(pp.pdb_clean, encoding="utf-8").read().splitlines()
    assert sum(1 for l in clean if l.startswith("ATOM")) == 60
    assert sum(1 for l in clean if l.startswith("TER")) == 1


def test_select_chain_warns_with_chain_list(tmp_path, caplog):
    p = tmp_path / "in.pdb"
    p.write_text("".join(_protein(60, "B") + _protein(60, "A")))
    pp = prot.ProteinPreparation(str(p), str(tmp_path))
    with caplog.at_level(logging.WARNING):
        assert pp.select_chain() == "A"
    assert "A, B" in caplog.text and "--chain" in caplog.text


def test_restore_labels_rejects_atom_count_mismatch(tmp_path):
    pdb = tmp_path / "c.pdb"
    pdb.write_text(_pdb_atom(1, "CA", "ALA", "A", "1", 0, 0, 0, "C")
                   + _pdb_atom(2, "CB", "ALA", "A", "1", 1.5, 0, 0, "C"))
    q = tmp_path / "r.pdbqt"
    q.write_text(_pdbqt_line(1, "C", 0.0, 0.0, 0.0, "C"))
    with pytest.raises(RuntimeError, match="heavy atoms"):
        prot._restore_residue_labels(str(q), str(pdb))


def test_restore_labels_hydrogen_search_radius(tmp_path):
    pdb = tmp_path / "c.pdb"
    pdb.write_text(_pdb_atom(1, "N", "ALA", "A", "1", 0, 0, 0, "N"))
    q = tmp_path / "r.pdbqt"
    q.write_text(_pdbqt_line(1, "N", 0.0, 0.0, 0.0, "N")
                 + _pdbqt_line(2, "H", 2.5, 0.0, 0.0, "HD"))
    with pytest.raises(RuntimeError, match="within 2.0 A"):
        prot._restore_residue_labels(str(q), str(pdb))
    q.write_text(_pdbqt_line(1, "N", 0.0, 0.0, 0.0, "N")
                 + _pdbqt_line(2, "H", 1.0, 0.0, 0.0, "HD"))
    assert prot._restore_residue_labels(str(q), str(pdb)) == 2
    h = [l for l in open(q, encoding="utf-8") if l[77:79] == "HD"][0]
    assert h[17:26] == "ALA A   1"


def test_restore_labels_warns_on_duplicate_coordinates(tmp_path, caplog):
    pdb = tmp_path / "c.pdb"
    pdb.write_text(_pdb_atom(1, "CA", "ALA", "A", "1", 0, 0, 0, "C")
                   + _pdb_atom(2, "CA", "GLY", "A", "2", 0, 0, 0, "C"))
    q = tmp_path / "r.pdbqt"
    q.write_text(_pdbqt_line(1, "C", 0.0, 0.0, 0.0, "C")
                 + _pdbqt_line(2, "C", 0.0, 0.0, 0.0, "C"))
    with caplog.at_level(logging.WARNING):
        prot._restore_residue_labels(str(q), str(pdb))
    assert "share coordinates" in caplog.text


# --- pockets, boxes ----------------------------------------------------------

def test_rank_pockets_ties_broken_by_pocket_number():
    pockets = [
        {"number": 3, "score": 0.5, "druggability_score": 0.2, "volume": 900.0},
        {"number": 1, "score": 0.5, "druggability_score": 0.2, "volume": 100.0},
        {"number": 2, "score": 0.5, "druggability_score": 0.3, "volume": 50.0},
        {"number": 4, "score": 0.7, "druggability_score": 0.0, "volume": 10.0},
    ]
    assert [p["number"] for p in prot.rank_pockets(pockets)] == [4, 2, 1, 3]
    assert [p["number"] for p in prot.rank_pockets(pockets[::-1])] == [4, 2, 1, 3]


def test_druggability_threshold_is_fpocket_documented_value():
    assert prot.FPOCKET_DRUGGABILITY_THRESHOLD == 0.5
    assert prot.is_druggable({"druggability_score": 0.5})
    assert not prot.is_druggable({"druggability_score": 0.3})


def test_bounding_box_padding_is_symmetric():
    xs, ys, zs = [0.0, 1.0, 1.0, 10.0], [0.0, 0.0, 0.0, 2.0], [5.0, 5.0, 5.0, 5.0]
    cx, cy, cz, sx, sy, sz = prot.bounding_box(xs, ys, zs, padding=4.0)
    assert (cx, cy, cz) == (5.0, 1.0, 5.0)
    for c, s, vals in ((cx, sx, xs), (cy, sy, ys), (cz, sz, zs)):
        assert min(vals) - (c - s / 2) == pytest.approx(4.0)
        assert (c + s / 2) - max(vals) == pytest.approx(4.0)


def test_pocket_box_centred_on_bounding_box(tmp_path):
    p = tmp_path / "pocket1_atm.pdb"
    p.write_text("".join(_pdb_atom(i, "C", "ALA", "A", "1", x, 0.0, 0.0, "C")
                         for i, x in enumerate((0.0, 1.0, 2.0, 12.0), 1)))
    pp = prot.ProteinPreparation("x.pdb", str(tmp_path))
    cx, _, _, sx, _, _ = pp._get_pocket_info([str(p)], padding=6.0)
    assert cx == 6.0 and sx == 24.0


def test_centroid_grid_uses_heavy_atoms_only(tmp_path):
    rec = tmp_path / "receptor.pdbqt"
    rec.write_text(_pdbqt_line(1, "C", 0.0, 0.0, 0.0, "C")
                   + _pdbqt_line(2, "C", 30.0, 30.0, 30.0, "C")
                   + _pdbqt_line(3, "H", 31.0, 31.0, 31.0, "HD"))
    pp = prot.ProteinPreparation("x.pdb", str(tmp_path))
    pp.receptor_pdbqt = str(rec)
    cx, cy, cz, sx, sy, sz = pp.centroid_grid(padding=1.0)
    assert (cx, cy, cz) == (15.0, 15.0, 15.0)
    assert sx == 32.0


def test_centroid_grid_raises_when_receptor_unreadable(tmp_path):
    pp = prot.ProteinPreparation(str(tmp_path / "missing.pdb"), str(tmp_path))
    with pytest.raises(RuntimeError, match="Cannot build the receptor box"):
        pp.centroid_grid()


def _fake_fpocket(pockets, info=True):
    """utils.run replacement writing an fpocket-like output directory."""
    def fake_run(cmd, capture=False, **kwargs):
        target = cmd[-1]
        root = os.path.splitext(target)[0] + "_out"
        os.makedirs(os.path.join(root, "pockets"), exist_ok=True)
        text = []
        for number, score, drug, coords in pockets:
            with open(os.path.join(root, "pockets", f"pocket{number}_atm.pdb"), "w", encoding="utf-8") as f:
                for i, (x, y, z) in enumerate(coords, 1):
                    f.write(_pdb_atom(i, "C", "ALA", "A", "1", x, y, z, "C"))
            text.append(f"Pocket {number} :\n\tScore : \t{score}\n"
                        f"\tDruggability Score : \t{drug}\n\tVolume : \t{100 * number}\n\n")
        if info:
            with open(os.path.join(root, os.path.basename(root[:-4]) + "_info.txt"), "w", encoding="utf-8") as f:
                f.write("".join(text))
        return ("", "", 0) if capture else None
    return fake_run


def _prep_with_receptor(tmp_path):
    p = tmp_path / "prot.pdb"
    p.write_text("".join(_protein(120)))
    out = tmp_path / "out"
    out.mkdir()
    pp = prot.ProteinPreparation(str(p), str(out))
    return pp


def test_detect_pocket_tie_break_and_symmetric_box(tmp_path, monkeypatch):
    pp = _prep_with_receptor(tmp_path)
    monkeypatch.setattr(tools, "find_tool", lambda *names: "/usr/bin/fpocket")
    monkeypatch.setattr(prot, "_fpocket_version", lambda path: "4.0")
    box3 = [(0.0, 0.0, 0.0), (10.0, 2.0, 4.0), (1.0, 1.0, 1.0)]
    monkeypatch.setattr(utils, "run", _fake_fpocket([
        (3, 0.4, 0.1, box3), (2, 0.4, 0.1, box3), (1, 0.1, 0.9, box3)]))
    cx, cy, cz, sx, sy, sz = pp.detect_pocket(padding=2.0)
    dx, dy, dz = pp._receptor_centroid()
    assert (cx + dx, cy + dy, cz + dz) == pytest.approx((5.0, 1.0, 2.0))
    assert (sx, sy, sz) == pytest.approx((14.0, 6.0, 8.0))
    rows = open(pp.pocket_summary_file, encoding="utf-8").read().splitlines()
    assert rows[1].startswith("1,") and rows[1].split(",")[5] == "YES"


def test_detect_pocket_selects_lowest_number_on_full_tie(tmp_path, monkeypatch):
    pp = _prep_with_receptor(tmp_path)
    monkeypatch.setattr(tools, "find_tool", lambda *names: "/usr/bin/fpocket")
    monkeypatch.setattr(prot, "_fpocket_version", lambda path: "4.0")
    far = [(50.0, 50.0, 50.0), (52.0, 52.0, 52.0)]
    near = [(0.0, 0.0, 0.0), (2.0, 2.0, 2.0)]
    monkeypatch.setattr(utils, "run", _fake_fpocket([
        (12, 0.3, 0.3, far), (2, 0.3, 0.3, near), (10, 0.3, 0.3, far)]))
    cx, cy, cz, _, _, _ = pp.detect_pocket(padding=1.0)
    dx, dy, dz = pp._receptor_centroid()
    assert (cx + dx, cy + dy, cz + dz) == pytest.approx((1.0, 1.0, 1.0))


def test_detect_pocket_fails_without_fpocket(tmp_path, monkeypatch):
    pp = _prep_with_receptor(tmp_path)
    monkeypatch.setattr(tools, "find_tool", lambda *names: None)
    with pytest.raises(RuntimeError, match="--no-fpocket"):
        pp.detect_pocket()


def test_detect_pocket_fails_when_no_pockets(tmp_path, monkeypatch):
    pp = _prep_with_receptor(tmp_path)
    monkeypatch.setattr(tools, "find_tool", lambda *names: "/usr/bin/fpocket")
    monkeypatch.setattr(prot, "_fpocket_version", lambda path: None)
    monkeypatch.setattr(utils, "run", _fake_fpocket([]))
    with pytest.raises(RuntimeError, match="no pockets"):
        pp.detect_pocket()


def test_detect_pocket_fails_without_info_file(tmp_path, monkeypatch):
    pp = _prep_with_receptor(tmp_path)
    monkeypatch.setattr(tools, "find_tool", lambda *names: "/usr/bin/fpocket")
    monkeypatch.setattr(prot, "_fpocket_version", lambda path: None)
    monkeypatch.setattr(utils, "run", _fake_fpocket(
        [(1, 0.3, 0.3, [(0.0, 0.0, 0.0)])], info=False))
    with pytest.raises(RuntimeError, match="pocket scores"):
        pp.detect_pocket()


def _parse_cgo(pml_text):
    """Tokens of the CGO list in a grid_box.pml script."""
    body = re.search(r"box = \[(.*?)\]", pml_text, re.S).group(1)
    return [t.strip() for t in body.replace("\n", " ").split(",") if t.strip()]


def test_grid_box_cgo_is_well_formed(tmp_path):
    pp = prot.ProteinPreparation("x.pdb", str(tmp_path))
    pp.write_grid_box_script(1.0, 2.0, 3.0, 20.0, 22.0, 24.0)
    tokens = _parse_cgo(open(pp.grid_box_script, encoding="utf-8").read())
    # Each CGO opcode takes a fixed number of floats.
    arity = {"COLOR": 3, "LINEWIDTH": 1, "BEGIN": 1, "VERTEX": 3, "END": 0}
    ops, i, vertices = [], 0, []
    while i < len(tokens):
        op = tokens[i]
        assert op in arity, op
        args = tokens[i + 1:i + 1 + arity[op]]
        if op == "BEGIN":
            assert args == ["LINES"]
        else:
            values = [float(a) for a in args]
            if op == "VERTEX":
                vertices.append(values)
        ops.append(op)
        i += 1 + arity[op]
    assert ops[:3] == ["COLOR", "LINEWIDTH", "BEGIN"] and ops[-1] == "END"
    assert len(vertices) == 24
    xs = sorted({v[0] for v in vertices})
    assert xs == pytest.approx([-9.0, 11.0])
