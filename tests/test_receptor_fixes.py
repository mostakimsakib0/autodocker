"""Regression tests for the receptor preparation review findings (R1-R5)
and the clock-seed collision of the hydrogen check."""
import logging
import os

import pytest

import runner
import vspipeline.flexible as fl
import vspipeline.hydrogens as hyd
import vspipeline.pocket as pocket
import vspipeline.protein as prot
import vspipeline.utils as utils

DATA = os.path.join(os.path.dirname(__file__), "data", "prep")
SEMET = os.path.join(DATA, "1hsg_semet_fragment.pdb")
ZN = os.path.join(DATA, "1hsg_zn_fragment.pdb")

needs_obabel = pytest.mark.skipif(runner.OBABEL is None, reason="obabel not available")


def _atoms(path):
    return [line for line in open(path, encoding="utf-8")
            if line.startswith(("ATOM", "HETATM"))]


def _pdb_atom(serial, name, resname, chain, resnum, x, y, z, element="C",
              record="ATOM", icode=" "):
    return (f"{record:<6s}{serial:5d} {name:<4s} {resname:>3s} {chain:1s}"
            f"{resnum:>4s}{icode}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00"
            f"          {element:>2s}\n")


def _pdbqt_atom(serial, name, resname, chain, resnum, x, y, z, atom_type="C"):
    return (f"ATOM  {serial:5d} {name:<4s} {resname:>3s} {chain:1s}{resnum:>4s}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    +0.000 {atom_type}\n")


# ---------------------------------------------------------------------------
# R1: selenomethionine, metal ions, Open Babel converting nothing
# ---------------------------------------------------------------------------

def _fake_run(stderr, write=None):
    def run(cmd, capture=False, extra_env=None):
        if write is not None and "-O" in cmd:
            with open(cmd[cmd.index("-O") + 1], "w", encoding="utf-8") as f:
                f.write(write)
        return ("", stderr, 0)
    return run


def test_run_obabel_fails_when_zero_molecules_converted(tmp_path, monkeypatch):
    out = tmp_path / "x.pdbqt"
    monkeypatch.setattr(runner, "run", _fake_run(
        "*** Open Babel Warning\n0 molecules converted\n", write=""))
    with pytest.raises(RuntimeError, match="converted 0 molecules") as exc:
        utils.run_obabel(["-ipdb", "a.pdb", "-opdbqt", "-O", str(out)], outputs=[str(out)])
    assert isinstance(exc.value, utils.ObabelNoOutput)
    assert not out.exists()


def test_run_obabel_fails_on_missing_or_empty_output(tmp_path, monkeypatch):
    out = tmp_path / "x.pdbqt"
    monkeypatch.setattr(runner, "run", _fake_run("1 molecule converted\n"))
    with pytest.raises(RuntimeError, match="wrote no output file x.pdbqt"):
        utils.run_obabel(["-O", str(out)], outputs=[str(out)])
    monkeypatch.setattr(runner, "run", _fake_run("1 molecule converted\n", write=""))
    with pytest.raises(RuntimeError, match="empty x.pdbqt"):
        utils.run_obabel(["-O", str(out)], outputs=[str(out)])
    assert not out.exists()
    monkeypatch.setattr(runner, "run", _fake_run("1 molecule converted\n", write="ATOM\n"))
    utils.run_obabel(["-O", str(out)], outputs=[str(out)])
    assert out.read_text() == "ATOM\n"


def test_semet_records_are_written_as_met(tmp_path):
    pp = prot.ProteinPreparation(SEMET, str(tmp_path))
    pp._write_clean_pdb(["A"], keep_hetero=False)
    clean = _atoms(pp.pdb_clean)
    assert not any("MSE" in line or line[76:78] == "SE" for line in clean)
    sd = [line for line in clean if line[12:16] == " SD " and line[22:26] == "  36"]
    assert len(sd) == 1 and sd[0].startswith("ATOM  ") and sd[0][17:20] == "MET"
    assert sd[0][76:78] == " S" and sd[0][30:54] == "  15.123  31.051  22.667"


@needs_obabel
def test_semet_receptor_is_prepared_as_met(tmp_path):
    """R1: Open Babel has no Gasteiger parameters for Se; a SeMet receptor
    gave an empty receptor.pdbqt. It is now prepared exactly like the same
    structure deposited with methionine."""
    met_pdb = tmp_path / "met.pdb"
    lines = []
    for line in open(SEMET, encoding="utf-8"):
        if line[17:20] == "MSE":
            line, _ = prot._semet_to_met("ATOM  " + line[6:])
        lines.append(line)
    met_pdb.write_text("".join(lines))
    outputs = []
    for name, src in (("semet", SEMET), ("met", str(met_pdb))):
        out = tmp_path / name
        out.mkdir()
        pp = prot.ProteinPreparation(src, str(out))
        pp.prepare_receptor("A")
        outputs.append(open(pp.receptor_pdbqt, "rb").read())
    assert outputs[0] == outputs[1]
    assert b" SD  MET A  36" in outputs[0]


@needs_obabel
def test_keep_hetero_metal_ions_are_appended_with_vina_types(tmp_path):
    """R1: with --keep-hetero a zinc ion made Open Babel convert nothing. Ions
    now bypass Open Babel and are appended with their Vina metal type; the
    two ions and two waters differing only by insertion code are all kept
    (R4)."""
    pp = prot.ProteinPreparation(ZN, str(tmp_path))
    pp.prepare_receptor("all", keep_hetero=True)
    atoms = _atoms(pp.receptor_pdbqt)
    zn = [line for line in atoms if line[77:79].strip() == "Zn"]
    assert [line[22:27] for line in zn] == [" 201 ", " 201A"]
    waters = [line for line in atoms if line[17:20] == "HOH" and line[77:79].strip() == "OA"]
    assert [line[22:27] for line in waters] == [" 301 ", " 301A"]
    dx, dy, dz = pp.frame_offset()
    assert abs(float(zn[0][30:38]) - (14.0 - dx)) < 1e-3
    # The ions are in pdb_clean (and in the frame offset) as kept residues.
    assert sum(1 for line in _atoms(pp.pdb_clean) if line[17:20] == " ZN") == 2
    # Keeping them again (--detect-metals) adds nothing.
    assert pp.rebuild_receptor_keep_hetatm([("A", "201", "ZN"), ("A", "201A", "ZN")]) \
        == pp.receptor_pdbqt
    assert len([line for line in _atoms(pp.receptor_pdbqt)
                if line[77:79].strip() == "Zn"]) == 2


@needs_obabel
def test_unsupported_hetero_element_fails_loudly_and_is_named(tmp_path):
    src = tmp_path / "xse.pdb"
    src.write_text(open(SEMET, encoding="utf-8").read().replace("MSE", "XSE"))
    pp = prot.ProteinPreparation(str(src), str(tmp_path / "out"))
    os.makedirs(pp.outdir)
    with pytest.raises(RuntimeError, match=r"converted 0 molecules.*SE in XSE A36"):
        pp.prepare_receptor("A", keep_hetero=True)
    assert not os.path.exists(pp.receptor_pdbqt)


# ---------------------------------------------------------------------------
# R2: flexible-residue keys
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("key,expected", [
    ("A45", ("A", "45")), ("a45", ("a", "45")), ("A100b", ("A", "100B")),
    ("45", ("", "45")), ("1:45", ("1", "45")), ("A:45", ("A", "45")),
    ("1:-3A", ("1", "-3A")), ("A-3", ("A", "-3")), ("145", ("", "145")),
    (":45", None), ("AB45", None), ("1:", None), ("A 45", None),
])
def test_parse_residue_key_forms(key, expected):
    assert fl._parse_residue_key(key) == expected


def test_blank_chain_constants_agree():
    assert fl.BLANK_CHAIN == prot.BLANK_CHAIN_ID == utils.DEFAULT_CHAIN


def _two_chain_receptor(tmp_path, chains=("A", "a")):
    lines = []
    serial = 1
    for k, chain in enumerate(chains):
        for name, xyz in (("N", (0, 0, 0)), ("CA", (1.45, 0, 0)), ("C", (2.0, 1.4, 0)),
                          ("O", (3.2, 1.6, 0)), ("CB", (1.9, -0.8, 1.2)),
                          ("OG", (3.3, -0.9, 1.3))):
            x, y, z = xyz
            atom_type = {"N": "NA", "O": "OA", "OG": "OA"}.get(name, "C")
            lines.append(_pdbqt_atom(serial, name, "SER", chain, "45", x + 10 * k, y, z,
                                     atom_type))
            serial += 1
    path = tmp_path / "rec.pdbqt"
    path.write_text("".join(lines))
    return str(path)


def test_flexible_residue_chains_are_case_sensitive(tmp_path):
    """R2: 'A45' collected the atoms of chains A and a."""
    rec = _two_chain_receptor(tmp_path)
    upper = fl._residue_atoms_from_pdbqt(rec, "A45")
    lower = fl._residue_atoms_from_pdbqt(rec, "a45")
    assert len(upper) == len(lower) == 6
    assert {a["serial"] for a in upper}.isdisjoint({a["serial"] for a in lower})
    out = tmp_path / "out"
    out.mkdir()
    built = []
    flex = fl.build_flexible_residue_pdbqt(rec, ["a45"], str(out), built=built)
    assert built == ["a45"]
    rigid = fl._remove_residues_from_pdbqt(rec, built, str(tmp_path / "rigid.pdbqt"),
                                           flex_file=flex)
    kept = _atoms(rigid)
    assert all(line[21] == "A" for line in kept if line[12:16].strip() in ("CA", "CB", "OG"))
    assert sum(1 for line in kept if line[21] == "a") == 3  # backbone N, C, O


def test_digit_chain_needs_separator_form(tmp_path):
    rec = _two_chain_receptor(tmp_path, chains=("1", "B"))
    assert len(fl._residue_atoms_from_pdbqt(rec, "1:45")) == 6
    assert fl._residue_atoms_from_pdbqt(rec, "145") == []  # blank chain, residue 145
    assert fl.format_residue_key("1", "  45", " ") == "1:45"
    assert fl._parse_residue_key(fl.format_residue_key("1", "  45", "B")) == ("1", "45B")
    keys = fl.detect_flexible_residues(rec, (1.5, 0.0, 0.0), radius=3.0, max_residues=5)
    assert keys == ["1:45"]


def test_blank_chain_answers_to_A45_and_45(tmp_path):
    """R2: a blank-chain receptor is reported as chain A; 'A46' and '46' both
    address its residues, and naming both builds the residue once."""
    rec = _two_chain_receptor(tmp_path, chains=(" ",))
    assert len(fl._residue_atoms_from_pdbqt(rec, "A45")) == 6
    assert len(fl._residue_atoms_from_pdbqt(rec, "45")) == 6
    out = tmp_path / "out"
    out.mkdir()
    built = []
    flex = fl.build_flexible_residue_pdbqt(rec, ["45", "A45"], str(out), built=built)
    assert built == ["45", "A45"]
    assert open(flex, encoding="utf-8").read().count("BEGIN_RES") == 1
    # The runner's check that every requested residue was built.
    assert {fl._parse_residue_key(k) for k in built} >= {
        fl._parse_residue_key(k) for k in ["45", "A45"]}


def test_missing_residues_uses_flexible_matching(tmp_path):
    lines = [_pdb_atom(1, "CA", "SER", "a", "45", 0, 0, 0),
             _pdb_atom(2, "CA", "SER", "1", "7", 0, 0, 0),
             _pdb_atom(3, "CA", "SER", "B", "100", 0, 0, 0, icode="B")]
    assert fl.missing_residues(lines, ["a", "1", "B"], ["a45", "A45", "1:7", "100", "B100B"]) \
        == ["A45", "100"]
    blank = [_pdb_atom(1, "CA", "SER", " ", "46", 0, 0, 0)]
    assert fl.missing_residues(blank, ["A"], ["46", "A46", "B46", "bad"]) == ["B46", "bad"]


# ---------------------------------------------------------------------------
# R3: blank chain IDs next to named chains
# ---------------------------------------------------------------------------

def _write(path, lines):
    path.write_text("".join(lines) + "END\n")
    return str(path)


def test_get_chains_reports_blank_chain_as_A(tmp_path):
    lines = ([_pdb_atom(i + 1, "CA", "ALA", " ", str(i + 1), i, 0, 0) for i in range(3)]
             + [_pdb_atom(i + 10, "CA", "ALA", "B", str(i + 1), i, 5, 0) for i in range(3)])
    path = _write(tmp_path / "blank_b.pdb", lines)
    assert utils.get_chains(path) == ["A", "B"]
    # --chain all keeps both chains (it kept only chain B).
    pp = prot.ProteinPreparation(path, str(tmp_path))
    pp._write_clean_pdb(utils._parse_chain_selection("all", utils.get_chains(path)), False)
    assert len(_atoms(pp.pdb_clean)) == 6


def test_get_chains_refuses_blank_next_to_chain_A(tmp_path):
    lines = [_pdb_atom(1, "CA", "ALA", " ", "1", 0, 0, 0),
             _pdb_atom(2, "CA", "ALA", "A", "2", 1, 0, 0)]
    path = _write(tmp_path / "blank_a.pdb", lines)
    with pytest.raises(utils.ChainAmbiguityError, match="blank chain ID"):
        utils.get_chains(path)
    with pytest.raises(ValueError, match="cannot be told apart"):
        prot.ProteinPreparation(path, str(tmp_path)).prepare_receptor("all")


def test_select_chain_logs_blank_chain(tmp_path, caplog):
    path = _write(tmp_path / "blank.pdb",
                  [_pdb_atom(1, "CA", "ALA", " ", "1", 0, 0, 0)])
    with caplog.at_level(logging.INFO):
        assert prot.ProteinPreparation(path, str(tmp_path)).select_chain() == "A"
    assert "blank chain ID in the PDB, read as A" in caplog.text


# ---------------------------------------------------------------------------
# R4: insertion codes of kept HETATM residues; R5: HETATM-only chains
# ---------------------------------------------------------------------------

def test_hetatm_residue_keys_include_insertion_code():
    keys = list(pocket._hetatm_residues(ZN))
    assert keys == [("A", "201", "ZN"), ("A", "201A", "ZN"),
                    ("A", "301", "HOH"), ("A", "301A", "HOH")]
    assert [(m["resnum"]) for m in pocket.detect_metal_ions(ZN)] == ["201", "201A"]
    waters = pocket.detect_water_molecules(ZN, (14.0, 22.0, 25.0), distance_threshold=10)
    assert [runner._keep_residue_key(w) for w in waters] == [
        ("A", "301", "HOH"), ("A", "301A", "HOH")]


def test_append_keeps_only_the_requested_insertion_code(tmp_path):
    rec = tmp_path / "rec.pdbqt"
    rec.write_text(_pdbqt_atom(1, "CA", "ALA", "A", "1", 0, 0, 0))
    n = pocket._append_hetatm_to_receptor(str(rec), ZN, [("A", "301A", "HOH")], ["A"])
    assert n == 1
    assert _atoms(str(rec))[-1][17:27] == "HOH A 301A"


def test_hetatm_only_chain_dropped_with_a_message(tmp_path, caplog):
    lines = ([_pdb_atom(i + 1, "CA", "ALA", "A", str(i + 1), i, 0, 0) for i in range(3)]
             + [_pdb_atom(20, "O", "HOH", "W", "1", 0, 9, 0, "O", record="HETATM")])
    path = _write(tmp_path / "hetchain.pdb", lines)
    pp = prot.ProteinPreparation(path, str(tmp_path))
    with caplog.at_level(logging.WARNING):
        pp._write_clean_pdb(utils._parse_chain_selection("all", utils.get_chains(path)), True)
    assert "chain(s) without protein atoms" in caplog.text and "W (1 records)" in caplog.text


# ---------------------------------------------------------------------------
# Hydrogen check: probes with other seeds; Open Babel 3.1 clock-seed collision
# ---------------------------------------------------------------------------

def _h_pdb(x):
    return (_pdb_atom(1, "O", "HOH", " ", "1", 0.0, 0.0, 0.0, "O", record="HETATM")
            + _pdb_atom(2, "H", "HOH", " ", "1", 0.957, 0.0, 0.0, "H", record="HETATM")
            + _pdb_atom(3, "H", "HOH", " ", "1", x, 0.927, 0.0, "H", record="HETATM"))


def test_third_conversion_when_first_two_agree(tmp_path):
    """Two runs with the same (Open Babel 3.1: clock) seed give identical
    random hydrogens; a third run that differs exposes them, and they are
    re-placed. The output is written with the run seed, the probes with the
    next two seeds."""
    xs = iter([-0.240, -0.240, 0.300])
    calls = []

    def fake(args, outputs=(), seed=None):
        calls.append((list(args), seed))
        with open(outputs[0], "w", encoding="utf-8") as f:
            f.write(_h_pdb(next(xs)))

    out = str(tmp_path / "w.pdb")
    moved = hyd.run_obabel_deterministic(["-O", out], out, run_obabel=fake, seed=7)
    assert len(calls) == 3 and moved == 2
    assert [seed for _, seed in calls] == [7, 8, 9]
    assert not os.path.exists(hyd._probe_path(out))


def test_probe_seeds_differ_from_the_run_seed():
    assert [hyd.probe_seed(42, k) for k in (1, 2)] == [43, 44]
    top = 2 ** 31 - 1  # largest --seed: the probes wrap to 1 and 2
    assert [hyd.probe_seed(top, k) for k in (1, 2)] == [1, 2]
    assert hyd.probe_seed(None, 1) == runner.DEFAULT_SEED + 1


def test_three_identical_conversions_are_kept_unchanged(tmp_path):
    calls = []

    def fake(args, outputs=(), seed=None):
        calls.append(list(args))
        with open(outputs[0], "w", encoding="utf-8") as f:
            f.write(_h_pdb(-0.240))

    out = str(tmp_path / "w.pdb")
    assert hyd.run_obabel_deterministic(["-O", out], out, run_obabel=fake) == 0
    assert len(calls) == 3
    assert open(out, encoding="utf-8").read() == _h_pdb(-0.240)


def test_receptor_protonation_compares_a_third_run(tmp_path, monkeypatch):
    """The receptor check also runs Open Babel a third time when the first
    two outputs agree, and re-places hydrogens that differ only then."""
    xs = iter([-0.240, -0.240, 0.300])
    calls = []

    def fake(args, outputs=(), log_stderr=True, seed=None):
        calls.append((list(args), seed))
        with open(outputs[0], "w", encoding="utf-8") as f:
            x = next(xs)
            f.write(f"ATOM      1  O   HOH A   1       0.000   0.000   0.000  0.00  0.00"
                    f"    -0.400 OA\n"
                    f"ATOM      2  H   HOH A   1       0.957   0.000   0.000  0.00  0.00"
                    f"    +0.200 HD\n"
                    f"ATOM      3  H   HOH A   1    {x:8.3f}   0.927   0.000  0.00  0.00"
                    f"    +0.200 HD\n")
        return ""

    monkeypatch.setattr(prot, "_run_obabel", fake)
    pp = prot.ProteinPreparation(str(tmp_path / "in.pdb"), str(tmp_path), seed=11)
    pp._protonate(7.4, source=str(tmp_path / "in.pdb"))
    assert [seed for _, seed in calls] == [11, 12, 13]
    lines = _atoms(pp.receptor_pdbqt)
    assert lines[2][30:38] != "  -0.240"  # re-placed deterministically
    assert not os.path.exists(pp.receptor_pdbqt + ".probe")
