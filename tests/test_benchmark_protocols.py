"""Regression tests for the 2.0.0 benchmark protocol fixes.

Fixtures (tests/data/bench/, no network):
  * dude/akt1_actives_head.mol2.gz: the first 11 entries of DUD-E's akt1
    actives_final.mol2.gz = 7 CHEMBL compounds, 4 of them with two
    protonation/tautomer entries.
  * casf/<pdbid>_modelserver.sdf: the RCSB ModelServer SDF of the ligand
    (H-less, CCD bond orders); casf/<pdbid>.pdb: a cut-down entry (24
    protein atoms + that ligand as HETATM records). 1y6r: the CCD bond
    orders lack a ring double bond; 1vso: phosphate written
    charge-separated in the core-set SMILES; 1nc3: an N-H the old
    sdf -> mol2 conversion lost.
"""
import argparse
import csv
import gzip
import importlib.util
import os
import shutil
from collections import Counter

import pytest

HERE = os.path.dirname(__file__)
SCRIPTS = os.path.join(HERE, "..", "scripts")
DATA = os.path.join(HERE, "data", "bench")
CORE_SET = os.path.join(HERE, "..", "benchmarks", "data", "casf2016_core_set.csv")


def _load(name):
    spec = importlib.util.spec_from_file_location(
        f"benchproto_{name}", os.path.join(SCRIPTS, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


be = _load("benchmark_enrichment")
br = _load("benchmark_redock")
fc = _load("fetch_casf_rcsb")
bl = _load("_benchlib")

needs_obabel = pytest.mark.skipif(not shutil.which("obabel"), reason="obabel not available")


def _rdkit():
    return pytest.importorskip("rdkit")


# ---------------------------------------------------------------------------
# Bug 1: DUD-E sampling by compound, not by MOL2 entry
# ---------------------------------------------------------------------------

AKT1_HEAD = os.path.join(DATA, "dude", "akt1_actives_head.mol2.gz")
AKT1_HEAD_ENTRIES = ["CHEMBL1081852", "CHEMBL469561", "CHEMBL469561", "CHEMBL576386",
                     "CHEMBL576386", "CHEMBL577134", "CHEMBL575824", "CHEMBL577379",
                     "CHEMBL577379", "CHEMBL573096", "CHEMBL573096"]


def _fake_convert(mol2, pdbqt):
    with open(mol2, encoding="utf-8") as src, open(pdbqt, "w", encoding="utf-8") as out:
        out.write(src.read())


def _sample(tmp_path, gz, n, seed=7, prefix="A_"):
    return be.sample_ligands(str(gz), str(tmp_path / "lig"), prefix, n, seed,
                             str(tmp_path / "scratch"))


class TestCompoundSampling:
    def test_fixture_has_several_entries_per_compound(self):
        blocks = be.read_molecules(AKT1_HEAD)
        assert [be.molecule_name(b) for b in blocks] == AKT1_HEAD_ENTRIES
        groups = be.group_compounds(blocks)
        assert len(groups) == 7 and groups["CHEMBL469561"] == [1, 2]

    def test_sample_holds_unique_compounds(self, tmp_path, monkeypatch):
        monkeypatch.setattr(be, "mol2_to_pdbqt", _fake_convert)
        # asking for more than there are compounds: every compound once
        paths, failures, info = _sample(tmp_path, AKT1_HEAD, 50)
        assert not failures and len(paths) == 7
        assert (info["file_entries"], info["file_compounds"]) == (11, 7)
        assert info["file_compounds_with_several_entries"] == 4
        compounds = [r["compound"] for r in info["selected"]]
        assert sorted(compounds) == sorted(set(AKT1_HEAD_ENTRIES))
        # one file per compound, named after it
        assert sorted(os.path.basename(p) for p in paths) == sorted(
            f"A_{c}.pdbqt" for c in set(AKT1_HEAD_ENTRIES))

    def test_first_entry_of_each_compound_is_used(self, tmp_path, monkeypatch):
        monkeypatch.setattr(be, "mol2_to_pdbqt", _fake_convert)
        blocks = be.read_molecules(AKT1_HEAD)
        paths, _, info = _sample(tmp_path, AKT1_HEAD, 50)
        for path, rec in zip(paths, info["selected"]):
            assert rec["entry_index"] == AKT1_HEAD_ENTRIES.index(rec["compound"])
            assert rec["n_entries"] == AKT1_HEAD_ENTRIES.count(rec["compound"])
            assert open(path, encoding="utf-8").read() == blocks[rec["entry_index"]]

    def test_sample_is_seeded_and_order_independent(self, tmp_path, monkeypatch):
        monkeypatch.setattr(be, "mol2_to_pdbqt", _fake_convert)
        blocks = be.read_molecules(AKT1_HEAD)
        rev = tmp_path / "reversed.mol2.gz"
        with gzip.open(rev, "wt") as f:
            f.write("".join(reversed(blocks)))
        _, _, a = _sample(tmp_path / "a", AKT1_HEAD, 4, seed=11)
        _, _, b = _sample(tmp_path / "b", AKT1_HEAD, 4, seed=11)
        _, _, c = _sample(tmp_path / "c", rev, 4, seed=11)
        pick = lambda info: [r["compound"] for r in info["selected"]]  # noqa: E731
        assert pick(a) == pick(b) == pick(c)
        assert len(set(pick(a))) == 4

    def test_scattered_decoy_duplicates_are_merged(self, tmp_path, monkeypatch):
        # DUD-E decoy files also repeat a ZINC ID, not always contiguously
        monkeypatch.setattr(be, "mol2_to_pdbqt", _fake_convert)
        gz = tmp_path / "decoys.mol2.gz"
        with gzip.open(gz, "wt") as f:
            f.write("".join(f"@<TRIPOS>MOLECULE\n{name}\n 1 0 0\nSMALL\n# {i}\n"
                            for i, name in enumerate(["ZINC1", "ZINC2", "ZINC1", "ZINC3"])))
        paths, _, info = _sample(tmp_path, gz, 10, prefix="D_")
        assert len(paths) == 3 and info["file_compounds"] == 3
        z1 = next(r for r in info["selected"] if r["compound"] == "ZINC1")
        assert (z1["entry_index"], z1["n_entries"]) == (0, 2)
        assert "# 0" in open(os.path.join(tmp_path, "lig", "D_ZINC1.pdbqt"),
                             encoding="utf-8").read()

    @needs_obabel
    def test_real_conversion(self, tmp_path):
        paths, failures, _ = _sample(tmp_path, AKT1_HEAD, 2)
        assert not failures and len(paths) == 2
        assert all(any(ln.startswith(("ATOM", "HETATM")) for ln in open(p, encoding="utf-8"))
                   for p in paths)

    def test_counts_and_rule_recorded_per_target(self, tmp_path, monkeypatch):
        monkeypatch.setattr(be, "mol2_to_pdbqt", _fake_convert)
        tdir = tmp_path / "work" / "akt1"
        tdir.mkdir(parents=True)
        shutil.copy(AKT1_HEAD, tdir / be.DUD_E_ACTIVES)
        with gzip.open(tdir / be.DUD_E_DECOYS, "wt") as f:
            f.write("".join(f"@<TRIPOS>MOLECULE\nZINC{i % 30}\n 1 0 0\nSMALL\n"
                            for i in range(40)))
        (tdir / be.DUD_E_RECEPTOR).write_text("ATOM\n")
        outdir = tmp_path / "out"
        rec = be.process_target("akt1", str(tmp_path / "work"), {
            "max_actives": 5, "max_decoys": 50, "box_source": "fpocket",
            "prepare_only": True, "outdir": str(outdir)})
        assert rec["status"] == "prepared", rec["error"]
        assert (rec["n_actives"], rec["n_decoys"]) == (5, 30)
        assert (rec["dude_actives_entries"], rec["dude_actives_compounds"]) == (11, 7)
        assert (rec["dude_decoys_entries"], rec["dude_decoys_compounds"]) == (40, 30)
        assert any("only 30 unique decoys" in w for w in rec["warnings"])
        with open(outdir / "akt1" / "ligands.csv", newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 35 and len({r["compound"] for r in rows}) == 35


# ---------------------------------------------------------------------------
# Bug 5: option-B defaults and primary metrics
# ---------------------------------------------------------------------------

class TestEnrichmentProtocol:
    def test_main_uses_option_b(self, tmp_path, monkeypatch):
        seen = {}

        def fake_process(target, workdir, opts):
            seen.update(opts)
            return {"target": target, "status": "prepared", "n_actives": 0,
                    "n_decoys": 0, "error": None}

        monkeypatch.setattr(be, "process_target", fake_process)
        assert be.main(["--targets", "akt1", "--prepare-only", "--ingest", str(tmp_path),
                        "--outdir", str(tmp_path / "out")]) == 0
        assert be.DEFAULT_MAX_ACTIVES == 50 and be.DEFAULT_MAX_DECOYS == 500
        assert (seen["max_actives"], seen["max_decoys"]) == (50, 500)
        import json
        with open(tmp_path / "out" / "summary.json", encoding="utf-8") as f:
            summary = json.load(f)
        assert summary["primary_metrics"] == ["auc", "ef5"]
        assert set(summary["secondary_metrics_caveats"]) == {"ef1", "bedroc"}
        assert summary["protocol"]["sampling_unit"].startswith("compound")
        assert "first" in summary["protocol"]["entry_per_compound"]

    def test_secondary_metrics_carry_warnings_at_50_of_550(self):
        names = [f"A_{i}" for i in range(50)] + [f"D_{i}" for i in range(500)]
        scores = [float(i) for i in range(550)]
        row = be.compute_stats(names, scores, 1, 0)
        assert row["ef1_n_top"] == 6 and row["ef5_n_top"] == 28
        text = " ".join(row["warnings"])
        assert "EF1% uses only n_top=6" in text
        assert "alpha*Ra = 7.3" in text
        assert not any(w.startswith("EF5%") for w in row["warnings"])

    def test_primary_metrics_lead_the_table(self):
        cols = be.SUMMARY_COLS
        assert cols.index("auc") < cols.index("ef5") < cols.index("ef1") < cols.index("bedroc")


# ---------------------------------------------------------------------------
# Bugs 2/3: one protonated ligand for both start modes; template bond orders
# ---------------------------------------------------------------------------

def _core_smiles(pdbid):
    with open(CORE_SET, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["pdbid"].lower() == pdbid:
                return row["smiles"].split("</s>")[-1].replace(" ", "")
    raise KeyError(pdbid)


def _fake_fetch(pdbid):
    def fetch(url, retries=4):
        name = f"{pdbid}.pdb" if "files.rcsb.org" in url else f"{pdbid}_modelserver.sdf"
        with open(os.path.join(DATA, "casf", name), "rb") as f:
            return f.read()
    return fetch


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """fetch_casf_rcsb.build() on the three fixture complexes (no network)."""
    pytest.importorskip("rdkit")
    if not shutil.which("obabel"):
        pytest.skip("obabel not available")
    out = tmp_path_factory.mktemp("casf")
    results = {}
    original = fc.fetch
    try:
        for pdbid in ("1y6r", "1vso", "1nc3"):
            fc.fetch = _fake_fetch(pdbid)
            results[pdbid] = fc.build(pdbid, _core_smiles(pdbid), str(out))
    finally:
        fc.fetch = original
    return out, results


def _mol2_elements(path):
    return Counter(a["element"].upper() for a in bl.parse_mol2_atoms(path))


class TestCasfLigand:
    @pytest.mark.parametrize("pdbid", ["1y6r", "1vso", "1nc3"])
    def test_both_start_modes_are_the_core_set_molecule(self, built, pdbid):
        out, results = built
        info, reason = results[pdbid]
        assert reason is None and "conformer_error" not in info, (reason, info)
        cdir = out / pdbid
        crystal = str(cdir / f"{pdbid}_ligand.mol2")
        rdkit = str(cdir / f"{pdbid}_ligand_rdkit.mol2")
        sdf = str(cdir / f"{pdbid}_ligand.sdf")
        key = bl.inchikey(crystal)
        assert key == info["inchikey_ligand"] == bl.inchikey(sdf) == bl.inchikey(rdkit)
        assert bl.same_compound_and_protonation(key, info["inchikey_core_set"])
        assert info["full_key_match"]
        # identical topology incl. every hydrogen in both start files
        assert _mol2_elements(crystal) == _mol2_elements(rdkit)
        assert _mol2_elements(crystal)["H"] > 0
        # crystal-start coordinates are the crystal's
        from rdkit import Chem
        ms = Chem.MolFromMolFile(str(cdir / f"{pdbid}_ligand_modelserver.sdf"),
                                 sanitize=False)
        ref = sorted(tuple(round(v, 3) for v in ms.GetConformer().GetAtomPosition(i))
                     for i in range(ms.GetNumAtoms()))
        heavy = sorted((round(x, 3), round(y, 3), round(z, 3))
                       for x, y, z in bl.mol2_heavy_coords(crystal))
        assert heavy == ref

    def test_ccd_bond_orders_replaced_by_template(self, built):
        # 1y6r: the CCD bond orders give a different molecule (a ring double
        # bond missing); the old graph-only fallback accepted it as it was.
        _, results = built
        info, _ = results["1y6r"]
        assert not bl.same_connectivity(info["inchikey_pdb"], info["inchikey_core_set"])
        assert info["identity_check"] == "heavy-atom-graph+smiles-bond-orders"
        assert info["inchikey_ligand"] == info["inchikey_core_set"]

    def test_lost_hydrogens_are_kept(self, built):
        # 1nc3: obabel sdf -> mol2 on the H-less SDF dropped an N-H; the
        # protonated SDF and the MOL2 now agree (checked by InChIKey above)
        # and the MOL2 has as many hydrogens as the RDKit molecule.
        out, results = built
        from rdkit import Chem
        mol = Chem.MolFromMolFile(str(out / "1nc3" / "1nc3_ligand.sdf"), removeHs=False)
        n_h = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
        assert _mol2_elements(str(out / "1nc3" / "1nc3_ligand.mol2"))["H"] == n_h

    def test_phosphate_written_as_dianion(self, built):
        _, results = built
        info, _ = results["1vso"]
        assert info["core_set_smiles_normalized"]
        assert "[P+]" not in info["core_set_smiles"]

    def test_pipeline_conversion_keeps_explicit_hydrogens(self, built, tmp_path):
        """The runner converts a MOL2 ligand with obabel -h (library.py);
        with every hydrogen explicit it must add none."""
        from rdkit import Chem
        from vspipeline.hydrogens import run_obabel_deterministic
        out, _ = built
        for pdbid in ("1y6r", "1vso", "1nc3"):
            mol2 = str(out / pdbid / f"{pdbid}_ligand.mol2")
            pdbqt = str(tmp_path / f"{pdbid}.pdbqt")
            run_obabel_deterministic(["-imol2", mol2, "-opdbqt", "-O", pdbqt, "-h"], pdbqt)
            types = [ln.split()[-1] for ln in open(pdbqt, encoding="utf-8")
                     if ln.startswith(("ATOM", "HETATM"))]
            mol = Chem.MolFromMolFile(str(out / pdbid / f"{pdbid}_ligand.sdf"),
                                      removeHs=False)
            polar_h = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1
                          and a.GetNeighbors()[0].GetAtomicNum() != 6)
            assert types.count("HD") == polar_h, pdbid
            assert bl.inchikey(mol2, extra_args=("-h",)) == bl.inchikey(mol2)

    def test_wrong_isomer_is_rejected(self):
        _rdkit()
        with open(os.path.join(DATA, "casf", "1nc3_modelserver.sdf"), encoding="utf-8") as f:
            sdf = f.read()
        # same formula, CH2OH and OH swapped on the sugar ring
        isomer = "Nc1ncnc2c([C@@H]3O[C@H](O)[C@@H](CO)[C@H]3O)n[nH]c12"
        with pytest.raises(fc.LigandMismatch):
            fc.protonated_ligand(sdf, isomer)
        fc.protonated_ligand(sdf, _core_smiles("1nc3"))  # the right one maps

    @needs_obabel
    def test_wrong_isomer_skips_the_complex(self, tmp_path, monkeypatch):
        _rdkit()
        monkeypatch.setattr(fc, "fetch", _fake_fetch("1nc3"))
        isomer = "Nc1ncnc2c([C@@H]3O[C@H](O)[C@@H](CO)[C@H]3O)n[nH]c12"
        info, reason = fc.build("1nc3", isomer, str(tmp_path))
        assert info is None and "none is the core-set molecule" in reason
        assert not (tmp_path / "1nc3" / "1nc3_ligand.mol2").exists()

    def test_normalize_core_set_smiles(self):
        _rdkit()
        smi, changed = fc.normalize_core_set_smiles("CCO[P+](=O)(=O)[O-]")
        assert changed
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smi)
        p = next(a for a in mol.GetAtoms() if a.GetSymbol() == "P")
        assert p.GetFormalCharge() == 0 and Chem.GetFormalCharge(mol) == -2
        assert sum(a.GetTotalNumHs() for a in mol.GetAtoms()) == 5   # unchanged
        assert fc.normalize_core_set_smiles("CC(=O)[O-]") == ("CC(=O)[O-]", False)
        assert fc.normalize_core_set_smiles("C[N+](=O)[O-]")[1] is False

    def test_redock_refuses_different_start_molecule(self):
        keys = {"crystal.mol2": "AAAAAAAAAAAAAA-BBBBBBBBSA-N",
                "start.mol2": "AAAAAAAAAAAAAA-BBBBBBBBSA-M"}
        with pytest.raises(RuntimeError, match="not the crystal ligand"):
            br.check_same_molecule("crystal.mol2", "start.mol2", key_fn=keys.get)
        br.check_same_molecule("crystal.mol2", "crystal.mol2", key_fn=keys.get)


# ---------------------------------------------------------------------------
# Bug 4: receptor keeps modified residues; metals only as ions
# ---------------------------------------------------------------------------

def _het(serial, name, resname, resnum, element, chain="A"):
    return (f"HETATM{serial:5d} {name:<4s} {resname:>3s} {chain}{resnum:4d}    "
            f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}")


class TestReceptor:
    RECORDS = [
        "ATOM      1  CA  ALA A  10       0.000   0.000   0.000  1.00 20.00           C",
        _het(2, "SE", "MSE", 11, "SE"), _het(3, "CA", "MSE", 11, "C"),
        _het(4, "P", "TPO", 12, "P"),
        "TER",
        _het(5, "FE", "HEM", 500, "FE"), _het(6, "NA", "HEM", 500, "N"),
        _het(7, "ZN", "ZN", 600, "ZN"), _het(8, "CA", "CA", 601, "CA"),
        _het(9, "O", "HOH", 700, "O"), _het(10, "C1", "LIG", 800, "C"),
    ]

    def _resnames(self, lines):
        return [ln[17:20].strip() for ln in lines if ln.startswith("HETATM")]

    def test_pipeline_modified_residues(self):
        mods = bl.pipeline_modified_residues()
        from vspipeline import protein
        assert mods == frozenset(protein.MODIFIED_RESIDUES) and "MSE" in mods

    def test_metals_mode(self):
        mods = bl.pipeline_modified_residues()
        out = fc.receptor_records(self.RECORDS, "metals", ("LIG", "A", "800"), mods)
        # MSE/TPO kept as protein; Zn and Ca ions kept; the heme's Fe is not
        # an ion and goes with the heme; water and the ligand dropped
        assert self._resnames(out) == ["MSE", "MSE", "TPO", "ZN", "CA"]

    def test_none_and_all_modes(self):
        mods = bl.pipeline_modified_residues()
        none = fc.receptor_records(self.RECORDS, "none", ("LIG", "A", "800"), mods)
        assert self._resnames(none) == ["MSE", "MSE", "TPO"]
        every = fc.receptor_records(self.RECORDS, "all", ("LIG", "A", "800"), mods)
        assert self._resnames(every) == ["MSE", "MSE", "TPO", "HEM", "HEM", "ZN", "CA"]

    def test_modified_residue_is_never_the_ligand(self):
        mods = bl.pipeline_modified_residues()
        groups = [k[0] for k, _ in fc.het_groups(self.RECORDS, mods)]
        assert "MSE" not in groups and "TPO" not in groups and "LIG" in groups


# ---------------------------------------------------------------------------
# Bug 7: fetch and redock select the same complexes
# ---------------------------------------------------------------------------

class TestSelection:
    def _fetch_ids(self, tmp_path, monkeypatch, *extra):
        _rdkit()
        seen = []
        monkeypatch.setattr(fc, "build", lambda pdbid, *a, **k: (seen.append(pdbid)
                                                                  or (None, "test")))
        fc.main(["--outdir", str(tmp_path / "fetch"), *extra])
        return seen

    def _redock_ids(self, *extra):
        ap = argparse.Namespace(pdbids=None, first=None, sample=bl.BENCH_SAMPLE_SIZE,
                                sample_seed=bl.BENCH_SEED, ids_csv=br.DEFAULT_IDS)
        for flag, value in zip(extra[::2], extra[1::2]):
            setattr(ap, flag.lstrip("-").replace("-", "_"), int(value))
        return br.select_ids(ap)

    def test_default_sample_is_identical(self, tmp_path, monkeypatch):
        fetched = self._fetch_ids(tmp_path, monkeypatch)
        assert fetched == self._redock_ids() and len(fetched) == bl.BENCH_SAMPLE_SIZE

    def test_first_uses_the_csv_in_both(self, tmp_path, monkeypatch):
        fetched = self._fetch_ids(tmp_path, monkeypatch, "--first", "5")
        assert fetched == self._redock_ids("--first", "5")
        assert fetched == sorted(bl.read_pdbid_rows(br.DEFAULT_IDS))[:5]

    def test_pdbids_are_lower_cased(self):
        ap = argparse.Namespace(pdbids=" 1Y6R,1nc3 ,1Y6R,", first=None)
        assert br.select_ids(ap) == ["1nc3", "1y6r"]
