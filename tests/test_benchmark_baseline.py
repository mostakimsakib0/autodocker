"""Tests of scripts/benchmark_baseline.py (standard Vina workflow baseline).

Fixtures (tests/data/baseline):
  casf/3lka  : the fetch_casf_rcsb.py files of CASF-2016 3lka with the
               receptor cut to the residues within 9 A of the ligand
  dude/akt1  : DUD-E akt1 crystal ligand, receptor cut to 9 A, 2 actives
               and 3 decoys (one MOL2 entry each)

The preparation/docking smoke tests need reduce, Meeko (mk_prepare_*.py),
Vina, Open Babel (obrms) and RDKit. Meeko is taken from
$AUTODOCKER_MEEKO_BIN (a virtual environment's bin/, as in
docs/reproducibility.md) or PATH; the tests are skipped when it does not run.
The EasyDock smoke test runs when $AUTODOCKER_EASYDOCK_BIN names the
easydock executable.
"""
import ast
import csv
import gzip
import importlib.util
import json
import os
import shutil
import subprocess
import sys

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
DATA = os.path.join(os.path.dirname(__file__), "data", "baseline")
CASF = os.path.join(DATA, "casf")
DUDE = os.path.join(DATA, "dude")


def _load(name):
    spec = importlib.util.spec_from_file_location(
        f"bench_{name}", os.path.join(ROOT, "scripts", f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    # Registered so multiprocessing can pickle the script's worker functions.
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


bb = _load("benchmark_baseline")
br = bb.br
be = bb.be


def _meeko_bin():
    """Directory holding working Meeko scripts, or None."""
    env = os.environ.get("AUTODOCKER_MEEKO_BIN")
    script = (os.path.join(env, "mk_prepare_ligand.py") if env
              else shutil.which("mk_prepare_ligand.py"))
    if not script or not os.path.exists(script):
        return None
    try:
        proc = subprocess.run([script, "--help"], capture_output=True, timeout=120)
    except (OSError, subprocess.SubprocessError):
        return None
    return os.path.dirname(script) if proc.returncode == 0 else None


MEEKO_BIN = _meeko_bin()
needs_chain = pytest.mark.skipif(
    not (MEEKO_BIN and shutil.which("vina") and shutil.which("obabel")
         and shutil.which("obrms")),
    reason="Meeko (AUTODOCKER_MEEKO_BIN), vina, obabel or obrms not available")


# ---------------------------------------------------------------------------
# Schema parity with the AutoDocker benchmarks
# ---------------------------------------------------------------------------

def _redock_header():
    """The results.csv header written by benchmark_redock.main()."""
    with open(os.path.join(ROOT, "scripts", "benchmark_redock.py"), encoding="utf-8") as f:
        tree = ast.parse(f.read())
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "writerow"
                and node.args and isinstance(node.args[0], ast.List)):
            values = [e.value for e in node.args[0].elts if isinstance(e, ast.Constant)]
            if values and values[0] == "pdbid":
                return values
    raise AssertionError("results.csv header not found in benchmark_redock.py")


def test_redock_columns_match_benchmark_redock():
    assert bb.REDOCK_COLUMNS == _redock_header()


def test_redock_results_and_summary_schema(tmp_path):
    recs = [{"pdbid": "1abc", "status": "ok", "affinity": -7.1, "affinity_mode1": -7.1,
             "rmsd_top1": 1.2345, "rmsd_best": 0.9, "rmsd_modes": [1.2345, None, 0.9],
             "n_modes_rmsd_failed": 1, "reused": False, "error": None},
            {"pdbid": "2xyz", "status": "error", "affinity": None, "affinity_mode1": None,
             "rmsd_top1": None, "rmsd_best": None, "rmsd_modes": [], "reused": False,
             "error": "boom"}]
    path = bb.write_redock_results(recs, str(tmp_path))
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[0] == _redock_header()
    assert rows[1] == ["1abc", "ok", "-7.1", "-7.1", "1.234", "0.900", "1.234 NA 0.900",
                       "1", "False", ""]
    assert rows[2][-1] == "boom"


def test_enrichment_ligand_selection_matches_benchmark_enrichment(tmp_path, monkeypatch):
    """Without --ligands-from the baseline samples the same compounds, with
    the same entries and file names, as benchmark_enrichment.py."""
    def fake_convert(mol2_path, pdbqt_path):
        with open(pdbqt_path, "w", encoding="utf-8") as f:
            f.write("x")
    monkeypatch.setattr(be, "mol2_to_pdbqt", fake_convert)
    seed = bb.target_seed("akt1")
    assert seed == int(__import__("hashlib").md5(b"akt1").hexdigest()[:8], 16)
    for kind, prefix, s in (("actives", "A_", seed), ("decoys", "D_", seed ^ be.DECOY_SEED_SALT)):
        gz = os.path.join(DUDE, "akt1", f"{kind}_final.mol2.gz")
        _, _, info = be.sample_ligands(gz, str(tmp_path / kind), prefix, 2, s,
                                       str(tmp_path / f"{kind}_scratch"))
        sel, failures, counts = bb.select_entries(gz, prefix, 2, s,
                                                  lambda block, name: name + ".sdf")
        assert failures == []
        assert [(r["ligand"], r["compound"], r["entry_index"], r["n_entries"]) for r in sel] \
            == [(r["ligand"], r["compound"], r["entry_index"], r["n_entries"])
                for r in info["selected"]]
        assert counts == {"file_entries": info["file_entries"],
                          "file_compounds": info["file_compounds"]}


def test_select_entries_from_ligand_table_keeps_failures():
    gz = os.path.join(DUDE, "akt1", "decoys_final.mol2.gz")
    blocks = be.read_molecules(gz)
    table = [{"ligand": f"D_{be.molecule_name(b)}", "compound": be.molecule_name(b),
              "entry_index": i, "n_entries": 1} for i, b in enumerate(blocks)]

    def convert(block, name):
        if name == table[1]["ligand"]:
            raise RuntimeError("unreadable")
        return name + ".sdf"
    sel, failures, _ = bb.select_entries(gz, "D_", 99, 0, convert, table)
    assert [r["ligand"] for r in sel] == [r["ligand"] for r in table]
    assert sel[1]["sdf"] is None and failures[0]["ligand"] == table[1]["ligand"]
    bad = [dict(table[0], compound="OTHER")]
    with pytest.raises(RuntimeError, match="is not OTHER"):
        bb.select_entries(gz, "D_", 99, 0, convert, bad)


def _pose(path, score):
    path.write_text(f"MODEL 1\nREMARK VINA RESULT:    {score}      0.000      0.000\n"
                    "ATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  0.00"
                    "    +0.000 C\nENDMDL\n", encoding="utf-8")
    return str(path)


def test_ranking_uses_the_pipeline_format_and_statuses(tmp_path):
    poses = {"A_good": _pose(tmp_path / "a.pdbqt", "-8.250"),
             "D_weak": _pose(tmp_path / "b.pdbqt", "-3.100"),
             "D_out": _pose(tmp_path / "c.pdbqt", "0.500")}
    ranking = bb.write_ranking(str(tmp_path / "ranking.csv"), poses, {"A_fail": "x"}, "Vina")
    with open(ranking, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        assert tuple(reader.fieldnames) == bb.vr.RANKING_COLUMNS
        rows = [(r["Ligand"], r["Binding_Affinity"], r["Status"]) for r in reader]
    assert rows == [("A_good", "-8.250", "OK"), ("D_weak", "-3.100", "WEAK"),
                    ("D_out", "0.500", "OUT_OF_WINDOW"), ("A_fail", "FAILED", "FAILED")]
    names, scores, info = be.collect_scores(ranking, ["A_good", "A_fail", "D_weak", "D_out"])
    assert dict(zip(names, scores)) == {"A_fail": None, "A_good": -8.25, "D_out": 0.5,
                                        "D_weak": -3.1}
    _, scores, _ = be.collect_scores(ranking, ["A_good", "D_out"], out_of_window="fail")
    assert scores == [-8.25, None]


# ---------------------------------------------------------------------------
# Preparation helpers
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name, element", [
    (" N  ", "N"), (" CA ", "C"), ("CA  ", "CA"), ("ZN  ", "ZN"), (" OG1", "O"),
    ("HD21", "H"), ("1HB ", "H"), (" SD ", "S"), ("FE  ", "FE")])
def test_element_from_name(name, element):
    assert bb.element_from_name(name) == element


def test_fill_pdb_elements_keeps_existing_ones(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(
        "ATOM      1  N   ARG     1      -2.419  -7.160  -1.307\n"
        "ATOM      2  CA  ARG     1      -1.648  -7.209  -2.588  1.00  0.00           C\n"
        "HETATM    3 ZN    ZN A 301       1.000   2.000   3.000  1.00  0.00\n"
        "END\n", encoding="utf-8")
    dst = tmp_path / "out.pdb"
    assert bb.fill_pdb_elements(str(src), str(dst)) == 2
    lines = dst.read_text(encoding="utf-8").splitlines()
    assert lines[0][76:78] == " N" and lines[0][:54] == \
        "ATOM      1  N   ARG     1      -2.419  -7.160  -1.307"
    assert lines[1][76:78] == " C"
    assert lines[2][76:78] == "ZN"
    assert lines[3] == "END"


def test_fill_pdb_elements_restores_dude_residue_names(tmp_path):
    """DUD-E's LEV/GLZ/ASQ/MEU become LEU/GLY/ASP/MET (Meeko has no template
    for them); other names, including the real TPO, are left alone."""
    src = tmp_path / "in.pdb"
    src.write_text(
        "ATOM    724  N   LEV    83       3.745  34.040   5.822\n"
        "ATOM    725  CA  GLZ A  27       5.140  33.982   6.264\n"
        "ATOM    726  CA  ASQ A  29       5.140  33.982   6.264\n"
        "ATOM    727  SD  MEU   793       5.140  33.982   6.264\n"
        "ATOM    728  OG1 TPO   308       5.140  33.982   6.264\n", encoding="utf-8")
    dst = tmp_path / "out.pdb"
    renamed = {}
    bb.fill_pdb_elements(str(src), str(dst), renamed)
    names = [line[17:20] for line in dst.read_text(encoding="utf-8").splitlines()]
    assert names == ["LEU", "GLY", "ASP", "MET", "TPO"]
    assert renamed == {"LEV": 1, "GLZ": 1, "ASQ": 1, "MEU": 1}


def test_script_interpreter_forms(tmp_path):
    plain = tmp_path / "plain" / "tool.py"
    plain.parent.mkdir()
    plain.write_text("#!/opt/py/bin/python3 -u\n", encoding="utf-8")
    assert bb.script_interpreter(str(plain)) == "/opt/py/bin/python3"
    pip = tmp_path / "pip" / "tool.py"
    pip.parent.mkdir()
    pip.write_text("#!/bin/sh\n'''exec' /venv/bin/python3 \"$0\" \"$@\"\n' '''\n",
                   encoding="utf-8")
    assert bb.script_interpreter(str(pip)) == "/venv/bin/python3"
    venv = tmp_path / "venv"
    venv.mkdir()
    (venv / "python3").write_text("", encoding="utf-8")
    (venv / "tool.py").write_text("#!/usr/bin/python3\n", encoding="utf-8")
    assert bb.script_interpreter(str(venv / "tool.py")) == str(venv / "python3")


def test_vina_command_uses_the_protocol(tmp_path):
    tools = {"vina": "/opt/vina"}
    opts = {"exhaustiveness": 8, "binding_modes": 9, "energy_range": 3.0, "seed": 42,
            "threads": 1}
    cmd = bb.vina_command(tools, "r.pdbqt", "l.pdbqt", "o.pdbqt", (1.0, 2.5, -3.25), 22.5,
                          opts)
    assert cmd == ["/opt/vina", "--receptor", "r.pdbqt", "--ligand", "l.pdbqt",
                   "--out", "o.pdbqt", "--center_x", "1.000", "--center_y", "2.500",
                   "--center_z", "-3.250", "--size_x", "22.500", "--size_y", "22.500",
                   "--size_z", "22.500", "--exhaustiveness", "8", "--num_modes", "9",
                   "--energy_range", "3.0", "--seed", "42", "--cpu", "1"]


def test_baseline_box_equals_the_autodocker_benchmark_box():
    crystal = os.path.join(CASF, "3lka", "3lka_ligand.mol2")
    cx, cy, cz, size = br.ligand_grid(crystal)
    center, edge = _load("_benchlib").ligand_box(_load("_benchlib").mol2_heavy_coords(crystal))
    assert (cx, cy, cz) == pytest.approx(center) and size == edge == 22.5


@pytest.mark.skipif(not shutil.which("obabel"), reason="obabel not available")
def test_redock_start_ligand_is_the_autodocker_start_conformer(tmp_path):
    pytest.importorskip("rdkit")
    opts = {"start_conformer": "rdkit", "conformer_seed": 42}
    sdf, meta = bb.redock_start_sdf("3lka", os.path.join(CASF, "3lka"), str(tmp_path), opts)
    assert meta["seed"] == 42 and meta["method"].startswith("ETKDGv3")
    # Same coordinates as the MOL2 benchmark_redock.py hands to AutoDocker.
    mol2 = os.path.join(tmp_path, "input", "3lka_ligand.mol2")
    lib = _load("_benchlib")
    mol2_xyz = [(a["x"], a["y"], a["z"]) for a in lib.parse_mol2_atoms(mol2)]
    from rdkit import Chem
    mol = Chem.MolFromMolFile(sdf, removeHs=False)
    sdf_xyz = [tuple(p) for p in mol.GetConformer().GetPositions()]
    assert len(sdf_xyz) == len(mol2_xyz) == 21
    assert all(a == pytest.approx(b, abs=1e-4) for a, b in zip(sdf_xyz, mol2_xyz))
    assert any(a.GetAtomicNum() == 1 for a in mol.GetAtoms())


# ---------------------------------------------------------------------------
# Engine availability
# ---------------------------------------------------------------------------

def test_easydock_unavailable_is_reported(tmp_path):
    out = tmp_path / "out"
    code = bb.main(["redock", "--data-dir", CASF, "--pdbids", "3lka", "--outdir", str(out),
                    "--engine", "easydock", "--easydock-bin", str(tmp_path / "nope")])
    assert code == bb.EXIT_ENGINE_UNAVAILABLE
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert summary["available"] is False and "not found" in summary["reason"]
    broken = tmp_path / "easydock"
    broken.write_text("#!/bin/sh\necho 'ModuleNotFoundError: No module named rdkit' >&2\n"
                      "exit 1\n", encoding="utf-8")
    broken.chmod(0o755)
    code = bb.main(["redock", "--data-dir", CASF, "--pdbids", "3lka", "--outdir", str(out),
                    "--engine", "easydock", "--easydock-bin", str(broken)])
    assert code == bb.EXIT_ENGINE_UNAVAILABLE
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert "No module named rdkit" in summary["reason"]


def test_unknown_residue_names_and_ids(tmp_path):
    text = ("Error: Creation of data structure for receptor failed.\n"
            "Template generation failed for unknown residues: {'CD', 'NI'}, which appear "
            "to be linking fragments.\n")
    assert bb.unknown_residue_names(text) == {"CD", "NI"}
    assert bb.unknown_residue_names("Template matching failed for: ['A:136']") == set()
    pdb = tmp_path / "r.pdb"
    pdb.write_text(
        "ATOM      1  N   ALA A   1      0.000   0.000   0.000  1.00  0.00           N\n"
        "HETATM    2 CD    CD A 201      1.000   0.000   0.000  1.00  0.00          CD\n"
        "HETATM    3 CD    CD A 201A     2.000   0.000   0.000  1.00  0.00          CD\n"
        "HETATM    4 NI    NI B1003      3.000   0.000   0.000  1.00  0.00          NI\n",
        encoding="utf-8")
    # The insertion code is part of the ID: Meeko does not delete A:201A for A:201.
    assert bb.residue_ids(str(pdb), {"CD", "NI"}) == ["A:201", "A:201A", "B:1003"]
    assert bb.residue_ids(str(pdb), set()) == []


def test_macrocycle_pseudo_atoms_are_detected(tmp_path):
    mac = os.path.join(DATA, "macrocycle", "2x00_ligand_out.pdbqt")
    assert bb.has_macrocycle_pseudo_atoms(mac)
    plain = tmp_path / "plain.pdbqt"
    plain.write_text("MODEL 1\n"
                     "ATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  0.00"
                     "     0.000 C \n"
                     "ENDMDL\n", encoding="utf-8")
    assert not bb.has_macrocycle_pseudo_atoms(str(plain))


@pytest.mark.skipif(not (MEEKO_BIN and shutil.which("obrms")
                         and os.path.exists(os.path.join(MEEKO_BIN or "", "mk_export.py"))),
                    reason="Meeko (AUTODOCKER_MEEKO_BIN) with mk_export.py, or obrms, not available")
def test_macrocycle_poses_get_an_rmsd(tmp_path):
    # Meeko opens the 15-membered ring of 2x00's ligand (G0/CG0 pseudo atoms);
    # obrms alone returns inf for such poses.
    mac = os.path.join(DATA, "macrocycle")
    out = str(tmp_path / "2x00_ligand_out.pdbqt")
    shutil.copy(os.path.join(mac, "2x00_ligand_out.pdbqt"), out)
    tools = {"mk_export": os.path.join(MEEKO_BIN, "mk_export.py")}
    rmsds, errors = bb.pose_rmsds(os.path.join(mac, "2x00_ligand.mol2"), out, str(tmp_path), tools)
    assert errors == []
    assert len(rmsds) == 2
    assert all(r is not None and 0.0 <= r < 20.0 for r in rmsds)


@needs_chain
def test_receptor_with_untemplated_ions_is_prepared(tmp_path):
    # Cd has no Meeko template: plain and -x preparation both fail, so the ion
    # residues (both insertion codes) are deleted with -d and recorded.
    pdb = tmp_path / "cd.pdb"
    lines = []
    with open(os.path.join(os.path.dirname(__file__), "data", "prep", "1hsg_zn_fragment.pdb"),
              encoding="utf-8") as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20] == " ZN":
                line = line[:12] + " CD " + line[16] + " CD" + line[20:76].ljust(56) + "CD\n"
            lines.append(line)
    pdb.write_text("".join(lines), encoding="utf-8")
    tools = bb.resolve_tools(bb.build_parser().parse_args(
        ["enrichment", "--workdir", str(tmp_path), "--outdir", str(tmp_path / "out"),
         "--meeko-bin", MEEKO_BIN]))
    out, rec = bb.prepare_receptor(str(pdb), str(tmp_path / "prep"), tools)
    assert os.path.getsize(out) > 0
    assert rec["deleted_untemplated_residues"] == ["A:201", "A:201A"]
    assert " -d A:201,A:201A -x" in rec["mk_prepare_receptor"]


def test_missing_preparation_tools_are_reported(tmp_path):
    out = tmp_path / "out"
    code = bb.main(["enrichment", "--workdir", DUDE, "--targets", "akt1", "--outdir",
                    str(out), "--meeko-bin", str(tmp_path)])
    assert code == bb.EXIT_ENGINE_UNAVAILABLE
    reason = json.loads((out / "summary.json").read_text(encoding="utf-8"))["reason"]
    assert "mk_prepare_receptor" in reason and "mk_prepare_ligand" in reason


# ---------------------------------------------------------------------------
# Preparation + docking smoke tests (tiny fixtures, exhaustiveness 1)
# ---------------------------------------------------------------------------

@pytest.mark.slow
@needs_chain
def test_redock_smoke(tmp_path):
    out = tmp_path / "redock"
    code = bb.main(["redock", "--data-dir", CASF, "--pdbids", "3lka", "--outdir", str(out),
                    "--meeko-bin", MEEKO_BIN, "--exhaustiveness", "1", "--binding-modes", "3",
                    "--processes", "1"])
    assert code == 0
    with open(out / "results.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert list(rows[0]) == _redock_header()
    (row,) = rows
    assert row["status"] == "ok", row["error"]
    assert float(row["affinity"]) < -3 and len(row["rmsd_modes"].split()) <= 3
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    expected = set(br.summarize([], 2.0)) | {"protocol"}
    assert set(summary) == expected
    proto = summary["protocol"]
    assert proto["seed"] == 42 and proto["exhaustiveness"] == 1
    assert proto["tools"]["meeko_python"]["packages"]["meeko"]
    work = out / "work" / "3lka"
    log = (work / "3lka_ligand_vina.log").read_text(encoding="utf-8")
    assert log.startswith("# Command: ") and "--seed 42" in log
    params = json.loads((work / "params.json").read_text(encoding="utf-8"))
    assert params["box"]["edge"] == 22.5
    assert params["receptor_preparation"]["receptor_pdbqt_sha256"]
    rec = (work / "receptor" / "receptor.pdbqt").read_text(encoding="utf-8")
    assert " HD" in rec  # polar hydrogens from reduce, typed by Meeko
    # A second run reuses the docking (same parameters and inputs).
    code = bb.main(["redock", "--data-dir", CASF, "--pdbids", "3lka", "--outdir", str(out),
                    "--meeko-bin", MEEKO_BIN, "--exhaustiveness", "1", "--binding-modes", "3",
                    "--processes", "1"])
    with open(out / "results.csv", newline="", encoding="utf-8") as f:
        again = list(csv.DictReader(f))[0]
    assert again["reused"] == "True" and again["rmsd_modes"] == row["rmsd_modes"]


@pytest.mark.slow
@needs_chain
def test_enrichment_smoke(tmp_path):
    work = tmp_path / "dude"
    shutil.copytree(DUDE, work)
    out = tmp_path / "enrichment"
    code = bb.main(["enrichment", "--workdir", str(work), "--targets", "akt1",
                    "--outdir", str(out), "--meeko-bin", MEEKO_BIN, "--exhaustiveness", "1",
                    "--max-actives", "1", "--max-decoys", "2", "--processes", "2",
                    "--bootstrap", "10"])
    assert code == 0
    with open(out / "enrichment.csv", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        assert reader.fieldnames == be.SUMMARY_COLS
        (row,) = list(reader)
    assert row["status"] == "ok", row["error"]
    assert (row["n_actives"], row["n_decoys"], row["total"]) == ("1", "2", "3")
    assert row["dude_actives_compounds"] == "2" and row["dude_decoys_entries"] == "3"
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert set(summary) == {"protocol", "primary_metrics", "secondary_metrics_caveats",
                            "targets"}
    assert summary["protocol"]["binding_modes"] == 1
    with open(out / "akt1" / "metrics.csv", newline="", encoding="utf-8") as f:
        metrics = list(csv.DictReader(f))
    assert [m["status"] for m in metrics] == ["OK", "OK", "OK"]
    with open(out / "akt1" / "ligands.csv", newline="", encoding="utf-8") as f:
        ligands = list(csv.DictReader(f))
    assert [r["ligand"][:2] for r in ligands] == ["A_", "D_", "D_"]
    # The baseline's ranking can be re-read exactly like AutoDocker's.
    names, scores, _ = be.collect_scores(str(out / "akt1" / "ranking.csv"))
    assert sorted(names) == sorted(r["ligand"] for r in ligands) and None not in scores
    # --ligands-from: the same compounds again, now from the ligands.csv.
    out2 = tmp_path / "enrichment2"
    code = bb.main(["enrichment", "--workdir", str(work), "--targets", "akt1",
                    "--outdir", str(out2), "--meeko-bin", MEEKO_BIN, "--exhaustiveness", "1",
                    "--ligands-from", str(out), "--processes", "2", "--bootstrap", "10"])
    assert code == 0
    assert (out2 / "akt1" / "ligands.csv").read_text(encoding="utf-8") == \
        (out / "akt1" / "ligands.csv").read_text(encoding="utf-8")


def test_dude_fixture_is_readable():
    for kind, n in (("actives", 2), ("decoys", 3)):
        with gzip.open(os.path.join(DUDE, "akt1", f"{kind}_final.mol2.gz"), "rt") as f:
            assert f.read().count("@<TRIPOS>MOLECULE") == n


@pytest.mark.skipif(not os.environ.get("AUTODOCKER_EASYDOCK_BIN"),
                    reason="set AUTODOCKER_EASYDOCK_BIN to the easydock executable")
@pytest.mark.slow
@needs_chain
def test_easydock_smoke(tmp_path):
    out = tmp_path / "easydock"
    code = bb.main(["redock", "--data-dir", CASF, "--pdbids", "3lka", "--outdir", str(out),
                    "--meeko-bin", MEEKO_BIN, "--exhaustiveness", "1", "--binding-modes", "2",
                    "--processes", "1", "--engine", "easydock",
                    "--easydock-bin", os.environ["AUTODOCKER_EASYDOCK_BIN"]])
    assert code == 0
    with open(out / "results.csv", newline="", encoding="utf-8") as f:
        (row,) = list(csv.DictReader(f))
    assert row["status"] == "ok", row["error"]
    assert len(row["rmsd_modes"].split()) == 2
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert summary["protocol"]["tools"]["easydock"]["packages"]["easydock"]
