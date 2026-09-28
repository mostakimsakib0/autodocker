import os

import pytest

import runner
import vspipeline.protein as prot
from tests.test_flexible import _atom  # reuse PDB ATOM line builder
from tests.test_pocket import _het  # reuse HETATM line builder


def _pdb(path, lines):
    path.write_text("\n".join(lines) + "\n")


def _protein_pdb(n=120, chain="A"):
    lines = []
    for i in range(n):
        lines.append(_atom(i + 1, "C", "ALA", chain, str((i % 10) + 1),
                           float(i), float(i), float(i), "C"))
    return lines


def test_validate_ok(tmp_path):
    p = tmp_path / "prot.pdb"
    _pdb(p, _protein_pdb(120))
    pp = prot.ProteinPreparation(str(p), str(tmp_path))
    pp.validate()  # should not raise


def test_validate_too_few(tmp_path):
    p = tmp_path / "prot.pdb"
    _pdb(p, _protein_pdb(10))
    pp = prot.ProteinPreparation(str(p), str(tmp_path))
    with pytest.raises(ValueError):
        pp.validate()


def test_validate_missing(tmp_path):
    pp = prot.ProteinPreparation(str(tmp_path / "nope.pdb"), str(tmp_path))
    with pytest.raises(FileNotFoundError):
        pp.validate()


def test_select_chain(tmp_path):
    p = tmp_path / "prot.pdb"
    _pdb(p, _protein_pdb(120, chain="A"))
    pp = prot.ProteinPreparation(str(p), str(tmp_path))
    assert pp.select_chain() == "A"


def test_select_chain_multiple(tmp_path):
    lines = _protein_pdb(60, chain="A") + _protein_pdb(60, chain="B")
    p = tmp_path / "prot.pdb"
    _pdb(p, lines)
    pp = prot.ProteinPreparation(str(p), str(tmp_path))
    # non-interactive: picks first chain
    assert pp.select_chain() == "A"


@pytest.mark.skipif(runner.OBABEL is None, reason="obabel not available")
def test_prepare_receptor_real(tmp_path):
    p = tmp_path / "prot.pdb"
    _pdb(p, _protein_pdb(120))
    pp = prot.ProteinPreparation(str(p), str(tmp_path))
    chains = pp.prepare_receptor(chain="A")
    assert chains == ["A"]
    assert os.path.exists(pp.receptor_pdbqt)
    assert runner._ensure_pdbqt_has_charges(pp.receptor_pdbqt)


def test_protein_centroid_grid(tmp_path):
    p = tmp_path / "prot.pdb"
    _pdb(p, _protein_pdb(120))
    pp = prot.ProteinPreparation(str(p), str(tmp_path))
    cx, cy, cz, sx, sy, sz = pp.centroid_grid()
    assert sx >= 24.0 and sy >= 24.0 and sz >= 24.0
    assert (cx, cy, cz) == (59.5, 59.5, 59.5)  # avg of 0..119


def test_receptor_centroid(tmp_path):
    p = tmp_path / "prot.pdb"
    _pdb(p, _protein_pdb(120))
    pp = prot.ProteinPreparation(str(p), str(tmp_path))
    c = pp._receptor_centroid()
    assert c == (59.5, 59.5, 59.5)


def test_detect_pocket_fallback(tmp_path, monkeypatch):
    p = tmp_path / "prot.pdb"
    _pdb(p, _protein_pdb(120))

    def fake_run(*a, **k):
        raise RuntimeError("fpocket missing")

    monkeypatch.setattr(runner, "run", fake_run)
    pp = prot.ProteinPreparation(str(p), str(tmp_path))
    grid = pp.detect_pocket(padding=6.0)
    # fallback grid: non-zero size, centroid on receptor
    assert grid[3] >= 24.0


def test_parse_fpocket_info(tmp_path):
    info = tmp_path / "prot_info.txt"
    info.write_text(
        "Pocket 1 :\n"
        "  Score : 0.85\n"
        "  Druggability Score : 0.52\n"
        "  Volume : 1234.5\n"
        "Pocket 2 :\n"
        "  Score : 0.60\n"
        "  Druggability Score : 0.20\n"
        "  Volume : 800.0\n"
    )
    pp = prot.ProteinPreparation("x.pdb", str(tmp_path))
    parsed = pp._parse_fpocket_info(str(info))
    assert parsed[1]["score"] == 0.85
    assert parsed[1]["druggability_score"] == 0.52
    assert parsed[1]["volume"] == 1234.5


def test_parse_pocket_selection():
    pockets = [{"number": 1}, {"number": 2}, {"number": 3}]
    pp = prot.ProteinPreparation("x.pdb", ".")
    assert pp._parse_pocket_selection("auto", pockets) == [1]
    assert pp._parse_pocket_selection("1,3", pockets) == [1, 3]
    with pytest.raises(ValueError):
        pp._parse_pocket_selection("9", pockets)


def test_write_pocket_summary(tmp_path):
    pp = prot.ProteinPreparation("x.pdb", str(tmp_path))
    pockets = [{"number": 1, "score": 0.8, "druggability_score": 0.5,
                "volume": 100.0, "box_volume": 200.0, "file": "p1.pdb"}]
    pp._write_pocket_summary(pockets)
    assert os.path.exists(pp.pocket_summary_file)
    txt = open(pp.pocket_summary_file).read()
    assert "Score" in txt and "0.8" in txt


def test_get_pocket_info(tmp_path):
    p = tmp_path / "pocket.pdb"
    _pdb(p, _protein_pdb(20))
    pp = prot.ProteinPreparation("x.pdb", str(tmp_path))
    cx, cy, cz, sx, sy, sz = pp._get_pocket_info([str(p)], padding=6.0)
    assert sx >= 6.0


def test_write_grid(tmp_path):
    pp = prot.ProteinPreparation("x.pdb", str(tmp_path))
    pp.write_grid(1.0, 2.0, 3.0, 25.0, 26.0, 27.0)
    assert os.path.exists(pp.grid_conf)
    txt = open(pp.grid_conf).read()
    assert "center_x = 1.0" in txt and "size_z = 27.0" in txt
    assert os.path.exists(pp.grid_box_script)


def test_write_grid_box_script(tmp_path):
    pp = prot.ProteinPreparation("x.pdb", str(tmp_path))
    pp.write_grid_box_script(0, 0, 0, 20, 20, 20)
    txt = open(pp.grid_box_script).read()
    assert "VERTEX" in txt


def test_rebuild_receptor_keep_hetatm(tmp_path):
    # receptor pdbqt (no ROOT tags) and a source pdb with a water HETATM
    rec = tmp_path / "receptor.pdbqt"
    rec.write_text(_atom(1, "C", "ALA", "A", "1", 0.0, 0.0, 0.0, "C"))
    pdb = tmp_path / "prot.pdb"
    pdb.write_text(_het("O", "HOH", "A", "1", 0.0, 0.0, 0.0, element="O"))
    pp = prot.ProteinPreparation(str(pdb), str(tmp_path))
    pp.receptor_pdbqt = str(rec)
    pp.selected_chains = ["A"]
    out = pp.rebuild_receptor_keep_hetatm([("A", "1", "HOH")])
    assert os.path.exists(out)
    assert "HETATM" in open(str(rec)).read()


def test_kept_hetatm_lands_in_receptor_frame(tmp_path):
    # Receptor PDBQT is centered by obabel -c; a kept water must move with it.
    p = tmp_path / "prot.pdb"
    water = _het("O", "HOH", "A", "500", 60.0, 59.0, 59.0, element="O")
    _pdb(p, _protein_pdb(120) + [water])
    pp = prot.ProteinPreparation(str(p), str(tmp_path))
    pp.prepare_receptor(chain="A")
    pp.rebuild_receptor_keep_hetatm([("A", "500", "HOH")])
    lines = open(pp.receptor_pdbqt).read().splitlines()
    het = [l for l in lines if l.startswith("HETATM") and "HOH" in l][0]
    xyz = (float(het[30:38]), float(het[38:46]), float(het[46:54]))
    # centroid of the protein is (59.5, 59.5, 59.5)
    assert xyz == pytest.approx((0.5, -0.5, -0.5), abs=0.01)


def test_frame_offset_file_restores_input_frame(tmp_path):
    import subprocess
    import sys
    p = tmp_path / "prot.pdb"
    _pdb(p, _protein_pdb(120))
    pp = prot.ProteinPreparation(str(p), str(tmp_path))
    pp.prepare_receptor(chain="A")
    pp.write_grid(0, 0, 0, 20, 20, 20)
    assert os.path.exists(tmp_path / "frame_offset.txt")
    script = os.path.join(os.path.dirname(__file__), "..", "scripts", "restore_frame.py")
    subprocess.run([sys.executable, script, str(tmp_path), pp.receptor_pdbqt], check=True)
    restored = str(pp.receptor_pdbqt).replace(".pdbqt", "_inputframe.pdbqt")
    first = [l for l in open(restored) if l.startswith("ATOM")][0]
    assert (float(first[30:38]), float(first[38:46]), float(first[46:54])) == \
        pytest.approx((0.0, 0.0, 0.0), abs=0.01)
