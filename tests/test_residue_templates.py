"""Trp/His ring hydrogens from residue templates (vspipeline.residue_templates)."""
import math
import os
import shutil

import pytest

import vspipeline.protein as prot
from vspipeline import residue_templates as rt
from vspipeline import tools
from vspipeline.hydrogens import pdbqt_element, vina_bond_cutoff

DATA = os.path.join(os.path.dirname(__file__), "data", "templates")
# 1HSG chain A residues 1-10 (Trp 6) and 60-75 (His 69), no hydrogens.
FRAGMENT = os.path.join(DATA, "1hsg_trp_his_fragment.pdb")

needs_obabel = pytest.mark.skipif(tools.settings().obabel is None, reason="obabel not available")


def _xyz(line):
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def _atoms(path):
    return [line for line in open(path, encoding="utf-8")
            if line.startswith(("ATOM", "HETATM"))]


def _residue(lines, resnum, chain="A"):
    return {line[12:16].strip(): line for line in lines
            if line[21] == chain and line[22:26].strip() == str(resnum)}


def _vina_bonded(h_line, lines):
    """Heavy atoms Vina bonds to the hydrogen (1.1 x covalent radii)."""
    return [line for line in lines if line[77:79].strip() not in ("H", "HD", "HS")
            and math.dist(_xyz(h_line), _xyz(line))
            < vina_bond_cutoff(pdbqt_element(line[77:79].strip()))]


def _ring_h(lines, n_line):
    return [line for line in lines if line[77:79].strip() == "HD"
            and math.dist(_xyz(line), _xyz(n_line)) < vina_bond_cutoff("N")]


def _pdb_line(name, resname, resnum, xyz, element, record="ATOM", chain="A", serial=1):
    return (f"{record:<6s}{serial:5d} {name:<4s} {resname:>3s} {chain}{resnum:>4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00          {element:>2s}\n")


def _pdbqt_from_fragment(path, rename=None, ob_his_h=None):
    """A PDBQT like Open Babel's: heavy atoms, ring N typed NA, no Trp HE1;
    ``ob_his_h`` adds an Open Babel-style hydrogen on that His N (at the end
    of the block, named H, as Open Babel writes it)."""
    out = ["REMARK  Name = protein_clean.pdb\n"]
    his = {}
    for i, line in enumerate(_atoms(FRAGMENT)):
        name = line[12:16].strip()
        element = line[76:78].strip()
        t = {"O": "OA", "S": "SA"}.get(element, element)
        if name in ("NE1", "ND1", "NE2"):
            t = "NA"
        if rename and line[17:20] == "HIS":
            line = line[:17] + rename + line[20:]
        if line[17:20] in rt.HIS_RESIDUE_NAMES:
            his[name] = _xyz(line)
        out.append(f"{line[:6]}{i + 1:5d}{line[11:54]}  1.00  0.00    -0.100 {t:<2s}\n")
    if ob_his_h:
        nbrs = rt.RING_N[("HIS", ob_his_h)][0]
        h = rt.ring_hydrogen_position(his[ob_his_h], his[nbrs[0]], his[nbrs[1]], 0.98)
        parent = [line for line in out if line[12:16].strip() == ob_his_h][0]
        out.append(f"ATOM  {len(out):5d}  H  {parent[16:30]}{h[0]:8.3f}{h[1]:8.3f}{h[2]:8.3f}"
                   "  1.00  0.00    +0.000 HD\n")
    out.append("TER\n")
    return "".join(out)


# --- geometry and rules --------------------------------------------------

def test_ring_hydrogen_is_on_the_external_bisector_in_the_ring_plane():
    n, a, b = (0.0, 0.0, 0.0), (-1.0, 1.0, 0.0), (-1.0, -1.0, 0.0)
    h = rt.ring_hydrogen_position(n, a, b)
    assert h == pytest.approx((rt.N_H_LENGTH, 0.0, 0.0))
    assert math.dist(h, n) < vina_bond_cutoff("N")


@pytest.mark.parametrize("resname, input_h, coordinated, expected", [
    ("HIS", None, set(), ("HIE", "default")),
    ("HID", None, set(), ("HID", "residue_name")),
    ("HSD", {"NE2"}, set(), ("HID", "residue_name")),     # name beats input H
    ("HIE", None, {"NE2"}, ("HIE", "residue_name")),      # name beats metal
    ("HSE", None, set(), ("HIE", "residue_name")),
    ("HIP", None, set(), ("HIP", "residue_name")),
    ("HSP", None, set(), ("HIP", "residue_name")),
    ("HIS", {"ND1"}, set(), ("HID", "input_hydrogens")),
    ("HIS", {"NE2"}, {"NE2"}, ("HIE", "input_hydrogens")),  # input H beats metal
    ("HIS", {"ND1", "NE2"}, set(), ("HIP", "input_hydrogens")),
    ("HIS", None, {"NE2"}, ("HID", "metal_coordination")),
    ("HIS", None, {"ND1"}, ("HIE", "metal_coordination")),
    ("HIS", None, {"ND1", "NE2"}, ("HIE", "default")),
])
def test_his_state_precedence(resname, input_h, coordinated, expected):
    assert rt.choose_his_state(resname, input_h, coordinated) == expected


def test_input_evidence_reads_ring_hydrogens_altlocs_and_metals(tmp_path):
    src = _atoms(FRAGMENT)
    his = _residue(src, 69)
    nd1 = _xyz(his["ND1"])
    h = rt.ring_hydrogen_position(nd1, _xyz(his["CG"]), _xyz(his["CE1"]), 1.01)
    hb = rt.ring_hydrogen_position(nd1, _xyz(his["CG"]), _xyz(his["CE1"]), 2.5)
    lines = src + [
        _pdb_line("HD1", "HIS", 69, h, "H"),
        _pdb_line("HE1", "TRP", 6, h, "H"),            # not a His: ignored
        _pdb_line("ZN", "ZN", 201, hb, "ZN", record="HETATM"),
        "ENDMDL\n",
        _pdb_line("HE2", "HIS", 69, _xyz(his["NE2"]), "H"),  # second model: ignored
    ]
    path = tmp_path / "in.pdb"
    path.write_text("".join(lines))
    found, metals = rt.read_input_evidence(str(path))
    assert found == {"A  69 ": {"ND1"}}
    assert metals == [pytest.approx(hb, abs=1e-3)]


# --- the PDBQT edit (no Open Babel) --------------------------------------

def test_templates_on_openbabel_like_pdbqt(tmp_path):
    q = tmp_path / "receptor.pdbqt"
    q.write_text(_pdbqt_from_fragment(FRAGMENT, ob_his_h="ND1"))
    before = _atoms(str(q))
    summary = rt.apply_residue_templates(str(q), FRAGMENT)
    after = _atoms(str(q))

    assert summary["trp_he1_added"] == 1
    assert summary["his_states"] == {"HID": 0, "HIE": 1, "HIP": 0}
    assert summary["his_residues"] == [{"residue": "A69", "input_name": "HIS",
                                        "state": "HIE", "decided_by": "default"}]
    assert summary["openbabel_ring_hydrogens_replaced"] == 1
    assert summary["skipped_incomplete_residues"] == []
    assert summary["hydrogens_within_bond_distance_of_another_atom"] == []

    trp, his = _residue(after, 6), _residue(after, 69)
    assert trp["NE1"][77:79] == "N " and his["NE2"][77:79] == "N "
    assert his["ND1"][77:79] == "NA"
    assert _ring_h(after, his["ND1"]) == []
    for n_line, h_name in ((trp["NE1"], "HE1"), (his["NE2"], "HE2")):
        (h_line,) = _ring_h(after, n_line)
        assert h_line[12:16] == f" {h_name}" and h_line[70:79] == "+0.000 HD"
        assert h_line[16:30] == n_line[16:30]
        assert math.dist(_xyz(h_line), _xyz(n_line)) == pytest.approx(rt.N_H_LENGTH, abs=2e-3)
        assert _vina_bonded(h_line, after) == [n_line]
        # Written right after its nitrogen.
        assert after[after.index(n_line) + 1] == h_line
    # Serials consecutive from 1; every other record unchanged except serial.
    assert [int(line[6:11]) for line in after] == list(range(1, len(after) + 1))
    changed = {line[11:] for line in after} ^ {line[11:] for line in before}
    assert {line[1:5].strip() for line in changed} == {"NE1", "NE2", "HE1", "HE2", "H"}


def test_templates_are_idempotent(tmp_path):
    q = tmp_path / "receptor.pdbqt"
    q.write_text(_pdbqt_from_fragment(FRAGMENT, ob_his_h="ND1"))
    rt.apply_residue_templates(str(q), FRAGMENT)
    once = q.read_bytes()
    summary = rt.apply_residue_templates(str(q), FRAGMENT)
    assert q.read_bytes() == once
    assert summary["openbabel_ring_hydrogens_replaced"] == 2


@pytest.mark.parametrize("name, state, protonated", [
    ("HID", "HID", ("ND1",)), ("HSE", "HIE", ("NE2",)), ("HIP", "HIP", ("ND1", "NE2"))])
def test_his_variant_names_are_used_and_written_as_his(tmp_path, name, state, protonated):
    q = tmp_path / "receptor.pdbqt"
    q.write_text(_pdbqt_from_fragment(FRAGMENT, rename=name))
    pdb = tmp_path / "in.pdb"
    pdb.write_text(open(FRAGMENT).read().replace(" HIS A  69", f" {name} A  69"))
    summary = rt.apply_residue_templates(str(q), str(pdb))
    assert summary["his_residues"][0]["state"] == state
    assert summary["his_residues"][0]["input_name"] == name
    assert summary["his_residue_names_written_as_HIS"] == 1
    after = _atoms(str(q))
    assert not any(line[17:20] == name for line in after)
    his = _residue(after, 69)
    for n in ("ND1", "NE2"):
        assert len(_ring_h(after, his[n])) == (n in protonated)
        assert his[n][77:79] == ("N " if n in protonated else "NA")


def test_incomplete_ring_is_skipped_and_reported(tmp_path):
    q = tmp_path / "receptor.pdbqt"
    text = "".join(line for line in _pdbqt_from_fragment(FRAGMENT).splitlines(True)
                   if not (line[17:20] == "TRP" and line[12:16] == " CE2"))
    q.write_text(text)
    summary = rt.apply_residue_templates(str(q), FRAGMENT)
    assert summary["trp_he1_added"] == 0
    assert summary["skipped_incomplete_residues"] == ["TRP A6"]
    assert _residue(_atoms(str(q)), 6)["NE1"][77:79] == "NA"


# --- receptor preparation with Open Babel --------------------------------

def _prepare(src, out, keep_hetero=False):
    out.mkdir()
    pp = prot.ProteinPreparation(str(src), str(out))
    pp.prepare_receptor("all", keep_hetero=keep_hetero)
    return pp


@needs_obabel
def test_receptor_preparation_applies_templates_deterministically(tmp_path):
    runs = [_prepare(FRAGMENT, tmp_path / f"r{i}") for i in range(3)]
    data = [open(pp.receptor_pdbqt, "rb").read() for pp in runs]
    assert data[0] == data[1] == data[2]
    pp = runs[0]
    assert pp.residue_templates["trp_he1_added"] == 1
    assert pp.residue_templates["his_states"]["HIE"] == 1
    lines = _atoms(pp.receptor_pdbqt)
    trp, his = _residue(lines, 6), _residue(lines, 69)
    assert trp["NE1"][77:79] == "N " and len(_ring_h(lines, trp["NE1"])) == 1
    assert his["NE2"][77:79] == "N " and len(_ring_h(lines, his["NE2"])) == 1
    assert his["ND1"][77:79] == "NA" and _ring_h(lines, his["ND1"]) == []
    for line in lines:
        if line[77:79] == "HD":
            assert len(_vina_bonded(line, lines)) == 1, line
    # The frame offset is still the centroid of the heavy atoms of pdb_clean.
    heavy = [_xyz(line) for line in lines if line[77:79] != "HD"]
    centroid = [sum(c) / len(heavy) for c in zip(*heavy)]
    assert centroid == pytest.approx([0.0, 0.0, 0.0], abs=2e-3)


@needs_obabel
def test_input_hydrogen_and_metal_rules_in_preparation(tmp_path):
    src = _atoms(FRAGMENT)
    his = _residue(src, 69)
    nd1, ne2 = _xyz(his["ND1"]), _xyz(his["NE2"])
    h_nd1 = rt.ring_hydrogen_position(nd1, _xyz(his["CG"]), _xyz(his["CE1"]), 1.01)
    zn_ne2 = rt.ring_hydrogen_position(ne2, _xyz(his["CD2"]), _xyz(his["CE1"]), 2.1)

    protonated = tmp_path / "protonated.pdb"
    protonated.write_text("".join(src + [_pdb_line("HD1", "HIS", 69, h_nd1, "H")]))
    pp = _prepare(protonated, tmp_path / "inputH")
    assert pp.residue_templates["his_residues"][0]["state"] == "HID"
    assert pp.residue_templates["his_residues"][0]["decided_by"] == "input_hydrogens"

    zinc = tmp_path / "zinc.pdb"
    zinc.write_text("".join(src + [_pdb_line("ZN", "ZN", 201, zn_ne2, "ZN", record="HETATM")]))
    for keep_hetero in (False, True):
        pp = _prepare(zinc, tmp_path / f"zn{keep_hetero}", keep_hetero=keep_hetero)
        rec = pp.residue_templates["his_residues"][0]
        assert (rec["state"], rec["decided_by"]) == ("HID", "metal_coordination")
        lines = _atoms(pp.receptor_pdbqt)
        h = _residue(lines, 69)
        assert h["NE2"][77:79] == "NA" and _ring_h(lines, h["NE2"]) == []
        assert h["ND1"][77:79] == "N " and len(_ring_h(lines, h["ND1"])) == 1


@needs_obabel
def test_template_receptor_is_accepted_by_vina(tmp_path):
    vina = shutil.which("vina")
    if vina is None:
        pytest.skip("vina not available")
    import subprocess
    pp = _prepare(FRAGMENT, tmp_path / "r")
    lig = tmp_path / "lig.pdbqt"
    lig.write_text("ROOT\nATOM      1  C   LIG     1       0.000   0.000   0.000  1.00  0.00"
                   "     0.000 C \nENDROOT\nTORSDOF 0\n")
    out = subprocess.run([vina, "--receptor", pp.receptor_pdbqt, "--ligand", str(lig),
                          "--score_only", "--center_x", "0", "--center_y", "0",
                          "--center_z", "0", "--size_x", "30", "--size_y", "30",
                          "--size_z", "30"], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr + out.stdout
