#!/usr/bin/env python3
"""Build a re-docking benchmark set for CASF-2016 core-set IDs from the RCSB PDB.

The PDBbind/CASF-2016 files need a registered download. This script builds
an openly reproducible substitute from the deposited RCSB entries of the
same PDB IDs. The structures are NOT the PDBbind-processed files: PDBbind
selects, completes and protonates its own receptor/ligand files, while this
script uses the deposited coordinates as they are (first model). Results on
this set are therefore not directly comparable with published CASF numbers.

  * IDs and ligand SMILES : ``benchmarks/data/casf2016_core_set.csv``;
                            a seeded random sample (``--sample 40
                            --sample-seed 42``) or ``--first N`` (first N
                            IDs of the CSV in sorted order); both select
                            exactly what ``benchmark_redock.py`` selects
                            with the same options
  * Alternate locations   : per residue, the blank/'A' conformer is kept
                            (the highest-occupancy one when there is no 'A')
  * Receptor              : ``{pdbid}_protein.pdb`` = ATOM/TER records plus
                            the modified amino acids the pipeline keeps as
                            protein (``vspipeline.protein.MODIFIED_RESIDUES``:
                            MSE, TPO, SEP, ...) plus, with ``--hetatm metals``
                            (default), metal ions = HETATM residues made of a
                            single metal atom (Zn, Mg, Ca, Mn, Fe, ...; the
                            rule of the runner's --detect-metals, so the Fe of
                            a heme is dropped with the heme);
                            ``--hetatm all`` keeps every non-water HETATM
                            except the chosen ligand instance; ``--hetatm
                            none`` keeps the protein only
  * Ligand                : one crystal instance of the HET group whose
                            heavy-atom formula matches the core-set SMILES
                            (coordinates from the RCSB ModelServer, kept as
                            ``{pdbid}_ligand_modelserver.sdf``). One RDKit
                            molecule is built from it: crystal coordinates
                            and heavy-atom graph; bond orders, formal charges
                            and hydrogen counts (the protonation state) from
                            the core-set SMILES via AssignBondOrdersFromTemplate
                            (the crystal graph must be isomorphic to the
                            SMILES graph; the SMILES' charge-separated
                            ``[P+](=O)(=O)[O-]`` phosphates, a valence-6
                            phosphorus, are first rewritten as the
                            ``P(=O)([O-])[O-]`` dianion with the same atoms
                            and hydrogens); stereo from the crystal
                            coordinates; explicit hydrogens added once. It is
                            written as ``{pdbid}_ligand.sdf`` and, via Open
                            Babel, ``{pdbid}_ligand.mol2`` (crystal-start
                            ligand and RMSD reference). The complex is skipped
                            unless the MOL2's InChIKey equals the SDF's (also
                            after ``obabel -h``, as the pipeline converts it)
                            and matches the core-set InChIKey in connectivity
                            block and protonation flag. The full-key match
                            (incl. stereo) and whether the CCD bond orders
                            alone already gave the right connectivity
                            (``identity_check``) are reported
  * Start conformer       : ``{pdbid}_ligand_rdkit.mol2`` (+ ``.json``) =
                            seeded ETKDGv3 + MMFF conformer of the SAME
                            protonated molecule, for ``benchmark_redock.py
                            --start-conformer rdkit`` (lets the re-docking run
                            where RDKit is missing). Its InChIKey must equal
                            the crystal MOL2's, else it is not written

RDKit is required (bond orders, protonation, conformers).

Ligands that are polymers (peptides) or span several residues have no single
matching HET group and are skipped. ``manifest.csv`` records the chosen
component, chain and residue, or the reason for skipping.

Usage:
  python3 scripts/fetch_casf_rcsb.py --outdir DATA_DIR [--sample 40 --sample-seed 42]
then
  python3 scripts/benchmark_redock.py --data-dir DATA_DIR ...
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _benchlib  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_IDS = os.path.join(HERE, "..", "benchmarks", "data", "casf2016_core_set.csv")
PDB_URL = "https://files.rcsb.org/download/{pdbid}.pdb"
LIGAND_URL = ("https://models.rcsb.org/v1/{pdbid}/ligand?label_comp_id={comp}"
              "&auth_asym_id={chain}&auth_seq_id={resnum}&encoding=sdf")
WATERS = {"HOH", "WAT", "DOD"}
# Metal ions kept in the receptor with --hetatm metals (single-atom residues
# only; see receptor_records).
METAL_ELEMENTS = {"ZN", "MG", "CA", "MN", "FE", "CU", "NI", "CO", "CD", "HG"}
# Never chosen as "the ligand" (solvent, buffer, ions).
SKIP_RESNAMES = WATERS | {"SO4", "PO4", "GOL", "EDO", "ACT", "DMS", "PEG",
                          "CL", "NA", "K", "IOD", "BR"} | METAL_ELEMENTS


def fetch(url, retries=4):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "autodocker-bench"})
            with urllib.request.urlopen(req, timeout=_benchlib.DOWNLOAD_TIMEOUT) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise
        except (urllib.error.URLError, ConnectionError, TimeoutError):
            pass
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"could not download {url}")


# ---------------------------------------------------------------------------
# PDB records
# ---------------------------------------------------------------------------

def _element(line):
    return (line[76:78].strip() or line[12:14].strip()).upper()


def _occupancy(line):
    try:
        return float(line[54:60])
    except ValueError:
        return 0.0


def first_model_records(pdb_text):
    """ATOM/HETATM/TER lines of the first model."""
    out = []
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM", "TER")):
            out.append(line)
    return out


def filter_altlocs(lines):
    """Keep one alternate conformer per residue.

    Atoms with a blank altloc are always kept. For a residue (chain, number,
    insertion code) that has alternate locations, the 'A' conformer is kept
    when present, otherwise the altloc with the highest mean occupancy (ties
    broken alphabetically). Kept lines get a blank altloc column. Non-atom
    lines (TER) pass through.
    """
    occ = defaultdict(lambda: defaultdict(list))
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and line[16] != " ":
            occ[(line[21], line[22:27])][line[16]].append(_occupancy(line))
    chosen = {}
    for res, by_alt in occ.items():
        if "A" in by_alt:
            chosen[res] = "A"
        else:
            chosen[res] = min(by_alt, key=lambda a: (-sum(by_alt[a]) / len(by_alt[a]), a))
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            alt = line[16]
            if alt != " ":
                if chosen.get((line[21], line[22:27])) != alt:
                    continue
                line = line[:16] + " " + line[17:]
        out.append(line)
    return out


def het_groups(records, modified_residues=frozenset()):
    """Heavy-atom element counts of every candidate ligand HET residue
    (solvent, buffer, ions and modified amino acids are never candidates)."""
    groups = defaultdict(Counter)
    order = []
    for line in records:
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        if resname in SKIP_RESNAMES or resname in modified_residues:
            continue
        element = _element(line)
        if element in ("H", "D"):
            continue
        key = (resname, line[21], line[22:26].strip())
        if key not in groups:
            order.append(key)
        groups[key][element] += 1
    return [(k, groups[k]) for k in order]


def _het_residue_key(line):
    return (line[17:20].strip(), line[21], line[22:26].strip(), line[26])


def receptor_records(records, hetatm="metals", ligand_key=None,
                     modified_residues=frozenset()):
    """Receptor lines: ATOM/TER, the modified amino acids of the chain, and
    the HETATM records selected by ``hetatm``:

    * ``none``: no other HETATM records;
    * ``metals``: metal ions, i.e. HETATM residues made of one metal atom
      (the rule of the pipeline's --detect-metals); a metal inside a larger
      HET group (the Fe of HEM, the Fe/S of SF4) is not an ion and is
      dropped with its group;
    * ``all``: every non-water HETATM except the residue ``ligand_key`` =
      (resname, chain, resnum).

    ``modified_residues`` (the pipeline's MODIFIED_RESIDUES: MSE, TPO, ...)
    are part of the protein and are kept in every mode; the pipeline turns
    them into ATOM records.
    """
    heavy_atoms = Counter(_het_residue_key(line) for line in records
                          if line.startswith("HETATM")
                          and _element(line) not in ("H", "D"))
    out = []
    for line in records:
        if line.startswith(("ATOM", "TER")):
            out.append(line)
        elif line.startswith("HETATM"):
            resname = line[17:20].strip()
            key = (resname, line[21], line[22:26].strip())
            if resname in WATERS or key == ligand_key:
                continue
            if resname in modified_residues or hetatm == "all":
                out.append(line)
            elif (hetatm == "metals" and _element(line) in METAL_ELEMENTS
                  and heavy_atoms[_het_residue_key(line)] == 1):
                out.append(line)
    return out


# ---------------------------------------------------------------------------
# Ligand identity
# ---------------------------------------------------------------------------

def formula_counts(formula):
    counts, elem, num = Counter(), "", ""
    for ch in formula + "Z":
        if ch.isupper() or ch == "Z":
            if elem:
                counts[elem.upper()] += int(num) if num else 1
            elem, num = ch, ""
        elif ch.islower():
            elem += ch
        elif ch.isdigit():
            num += ch
    counts.pop("H", None)
    return counts


def smiles_formula(smiles):
    """Heavy-atom element counts of a SMILES string."""
    out = subprocess.run(["obabel", f"-:{smiles}", "-d", "-otxt", "--append", "formula"],
                         capture_output=True, text=True).stdout.strip()
    return formula_counts(out)


class LigandMismatch(ValueError):
    """The crystal ligand cannot be made into the core-set molecule."""


def normalize_core_set_smiles(smiles):
    """Repair the one non-physical charge pattern of the core-set SMILES.

    The core-set CSV writes deprotonated phosphate groups charge-separated
    as ``[P+](=O)(=O)[O-]``: a phosphorus with valence 6 (9 molecules of
    the 285). It is rewritten as ``P(=O)([O-])[O-]``, the dianion with the
    same atoms and hydrogens (no hydrogen is added or removed, only the
    formal charges and one P=O bond order change). Other molecules are
    returned unchanged. Returns ``(smiles, changed)``.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles, False
    rw, changed = Chem.RWMol(mol), False
    Chem.Kekulize(rw, clearAromaticFlags=True)
    for atom in rw.GetAtoms():
        if not (atom.GetSymbol() == "P" and atom.GetFormalCharge() == 1
                and atom.GetTotalValence() == 6):
            continue
        for bond in atom.GetBonds():
            other = bond.GetOtherAtom(atom)
            if (bond.GetBondType() == Chem.BondType.DOUBLE and other.GetSymbol() == "O"
                    and other.GetDegree() == 1):
                bond.SetBondType(Chem.BondType.SINGLE)
                other.SetFormalCharge(-1)
                atom.SetFormalCharge(0)
                changed = True
                break
    if not changed:
        return smiles, False
    Chem.SanitizeMol(rw)
    return Chem.MolToSmiles(rw), True


def protonated_ligand(sdf_text, smiles):
    """The crystal ligand as one RDKit molecule with explicit hydrogens.

    Coordinates and the heavy-atom graph come from the crystal (ModelServer
    SDF); bond orders, formal charges and hydrogen counts (i.e. the
    protonation state) come from the core-set SMILES, assigned with
    ``AllChem.AssignBondOrdersFromTemplate`` after the crystal's own bond
    orders and charges are reset, so the result never depends on how the
    CCD (or a reader) drew them. Stereo is taken from the crystal
    coordinates. Hydrogens are added once, with coordinates.

    Raises LigandMismatch when the heavy-atom graphs differ (atom or bond
    count, or no isomorphism with matching elements).
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    crystal = Chem.MolFromMolBlock(sdf_text, sanitize=False, removeHs=False)
    if crystal is None or crystal.GetNumConformers() == 0:
        raise LigandMismatch("RDKit cannot read the crystal ligand SDF")
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise LigandMismatch(f"RDKit cannot read the core-set SMILES {smiles}")
    heavy = Chem.RWMol(crystal)
    for idx in sorted((a.GetIdx() for a in heavy.GetAtoms() if a.GetAtomicNum() == 1),
                      reverse=True):
        heavy.RemoveAtom(idx)
    for bond in heavy.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in heavy.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetIsAromatic(False)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
        atom.SetNumRadicalElectrons(0)
    heavy = heavy.GetMol()
    heavy.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(heavy)
    if (heavy.GetNumAtoms(), heavy.GetNumBonds()) != (template.GetNumAtoms(),
                                                      template.GetNumBonds()):
        raise LigandMismatch(
            f"heavy-atom graph differs from the core-set SMILES ({heavy.GetNumAtoms()} "
            f"atoms/{heavy.GetNumBonds()} bonds vs {template.GetNumAtoms()}/"
            f"{template.GetNumBonds()})")
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, heavy)
    except (ValueError, RuntimeError) as e:  # no match, or sanitization failed
        raise LigandMismatch(f"core-set SMILES does not map onto the crystal "
                             f"heavy-atom graph: {e}") from None
    Chem.AssignStereochemistryFrom3D(mol)
    mol = Chem.AddHs(mol, addCoords=True)
    return mol


def write_ligand_files(mol, sdf_path, mol2_path, name):
    """Write ``mol`` (explicit H, crystal coordinates) as SDF (RDKit) and
    MOL2 (Open Babel ``sdf -> mol2``; the MOL2 keeps every hydrogen)."""
    from rdkit import Chem
    mol = Chem.Mol(mol)
    mol.SetProp("_Name", name)
    with open(sdf_path, "w", encoding="utf-8") as f:
        f.write(Chem.MolToMolBlock(mol, kekulize=True))
        f.write("$$$$\n")
    subprocess.run(["obabel", "-isdf", sdf_path, "-omol2", "-O", mol2_path],
                   check=True, capture_output=True, timeout=_benchlib.TOOL_TIMEOUT)


def check_ligand_files(sdf_path, mol2_path, ref_key):
    """InChIKeys of the written ligand; returns ``(keys, problem)``.

    The MOL2 must be the same molecule as the SDF it was written from (full
    InChIKey), also after ``obabel -h`` (the pipeline's MOL2 -> PDBQT
    conversion adds hydrogens with -h; with every hydrogen explicit it must
    add none), and must be the core-set compound in the core-set
    protonation state (see _benchlib.same_compound_and_protonation).
    """
    keys = {"sdf": _benchlib.inchikey(sdf_path),
            "mol2": _benchlib.inchikey(mol2_path),
            "mol2_obabel_h": _benchlib.inchikey(mol2_path, extra_args=("-h",)),
            "core_set": ref_key}
    if not keys["sdf"] or keys["mol2"] != keys["sdf"]:
        return keys, (f"MOL2 written by Open Babel is not the built molecule "
                      f"(SDF {keys['sdf']}, MOL2 {keys['mol2']})")
    if keys["mol2_obabel_h"] != keys["mol2"]:
        return keys, (f"obabel -h changes the MOL2 ligand ({keys['mol2']} -> "
                      f"{keys['mol2_obabel_h']})")
    if not _benchlib.same_compound_and_protonation(keys["mol2"], ref_key):
        return keys, (f"ligand InChIKey {keys['mol2']} is not the core-set compound/"
                      f"protonation state {ref_key}")
    return keys, None


def normalize_sdf_elements(sdf):
    """Write two-letter element symbols as Cl/Br, not CL/BR.

    The ModelServer writes them in capitals, which Open Babel reads as an
    unknown (dummy) atom, so the ligand loses its halogens.
    """
    lines = sdf.split("\n")
    try:
        n_atoms = int(lines[3][:3])
    except (IndexError, ValueError):
        return sdf
    for i in range(4, 4 + n_atoms):
        line = lines[i]
        elem = line[31:34]
        fixed = elem.strip().capitalize().ljust(3)
        lines[i] = line[:31] + fixed + line[34:]
    return "\n".join(lines)


def build(pdbid, smiles, outdir, hetatm="metals", conformer_seed=_benchlib.BENCH_SEED,
          modified_residues=None):
    """Write one complex; returns ``(info, None)`` or ``(None, reason)``."""
    if modified_residues is None:
        modified_residues = _benchlib.pipeline_modified_residues()
    pdb_text = fetch(PDB_URL.format(pdbid=pdbid.upper())).decode()
    records = filter_altlocs(first_model_records(pdb_text))
    target = smiles_formula(smiles)
    matches = [key for key, counts in het_groups(records, modified_residues)
               if counts == target]
    if not matches:
        return None, "no HET group matches the ligand formula (peptide or multi-residue ligand)"
    smiles, smiles_normalized = normalize_core_set_smiles(smiles)
    ref_key = _benchlib.inchikey(smiles, is_smiles=True)
    if not ref_key:
        return None, f"Open Babel computes no InChIKey for the core-set SMILES {smiles}"
    cdir = os.path.join(outdir, pdbid)
    os.makedirs(cdir, exist_ok=True)
    raw_sdf = os.path.join(cdir, f"{pdbid}_ligand_modelserver.sdf")
    sdf_path = os.path.join(cdir, f"{pdbid}_ligand.sdf")
    mol2 = os.path.join(cdir, f"{pdbid}_ligand.mol2")
    chosen, keys, ccd_key, problems = None, None, None, []
    for comp, chain, resnum in matches:
        sdf = fetch(LIGAND_URL.format(pdbid=pdbid.lower(), comp=comp, chain=chain,
                                      resnum=resnum)).decode()
        sdf = normalize_sdf_elements(sdf.split("$$$$")[0] + "$$$$\n")
        with open(raw_sdf, "w", encoding="utf-8") as f:
            f.write(sdf)
        ccd_key = _benchlib.inchikey(raw_sdf)
        try:
            mol = protonated_ligand(sdf, smiles)
            write_ligand_files(mol, sdf_path, mol2, pdbid)
            keys, problem = check_ligand_files(sdf_path, mol2, ref_key)
        except (LigandMismatch, subprocess.SubprocessError) as e:
            problem = str(e)
        if problem is None:
            chosen = (comp, chain, resnum)
            break
        problems.append(f"{comp} {chain}{resnum}: {problem}")
    if chosen is None:
        for path in (raw_sdf, sdf_path, mol2):
            if os.path.exists(path):
                os.remove(path)
        return None, (f"formula matches {len(matches)} HET group(s) but none is the "
                      f"core-set molecule ({'; '.join(problems)[:600]})")
    with open(os.path.join(cdir, f"{pdbid}_protein.pdb"), "w", encoding="utf-8") as f:
        for line in receptor_records(records, hetatm, ligand_key=chosen,
                                     modified_residues=modified_residues):
            f.write(line + "\n")
        f.write("END\n")
    info = {"comp": chosen[0], "chain": chosen[1], "resnum": chosen[2],
            "instances": len(matches),
            "inchikey_pdb": ccd_key,
            "inchikey_ligand": keys["mol2"],
            "inchikey_core_set": ref_key,
            "core_set_smiles": smiles,
            "core_set_smiles_normalized": smiles_normalized,
            # How the crystal heavy atoms were identified as the core-set
            # compound; bond orders/protonation always come from the SMILES.
            "identity_check": ("ccd-inchikey-connectivity"
                               if _benchlib.same_connectivity(ccd_key, ref_key)
                               else "heavy-atom-graph+smiles-bond-orders"),
            "full_key_match": keys["mol2"] == ref_key}
    rd_mol2 = os.path.join(cdir, f"{pdbid}_ligand_rdkit.mol2")
    rd_json = rd_mol2[:-len(".mol2")] + ".json"
    try:
        meta = _benchlib.rdkit_start_conformer(sdf_path, rd_mol2, conformer_seed)
        meta["inchikey"] = _benchlib.inchikey(rd_mol2)
        if meta["inchikey"] != keys["mol2"]:
            raise RuntimeError(f"RDKit start conformer {meta['inchikey']} differs from "
                               f"the crystal ligand {keys['mol2']}")
        with open(rd_json, "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)
    except Exception as e:  # the crystal-start benchmark still works
        info["conformer_error"] = str(e)[:200]
        for path in (rd_mol2, rd_json):
            if os.path.exists(path):
                os.remove(path)
    return info, None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ids", default=DEFAULT_IDS, help="CSV with pdbid,smiles columns")
    ap.add_argument("--sample", type=int, default=_benchlib.BENCH_SAMPLE_SIZE,
                    help=f"Seeded random sample size (default {_benchlib.BENCH_SAMPLE_SIZE}; "
                    "0 = all IDs)")
    ap.add_argument("--sample-seed", type=int, default=_benchlib.BENCH_SEED)
    ap.add_argument("--first", "--max-complexes", dest="first", type=int, default=None,
                    help="Legacy 1.1.x selection: first N IDs in sorted order")
    ap.add_argument("--hetatm", choices=["none", "metals", "all"], default="metals",
                    help="HETATM records kept in the receptor (default metals)")
    ap.add_argument("--conformer-seed", type=int, default=_benchlib.BENCH_SEED,
                    help="ETKDG seed of the precomputed RDKit start conformer")
    args = ap.parse_args(argv)

    try:
        import rdkit  # noqa: F401
    except ImportError:
        ap.error("fetch_casf_rcsb.py needs RDKit (bond orders, protonation, conformers)")

    all_rows = _benchlib.read_pdbid_rows(args.ids)
    ids = _benchlib.select_complexes(all_rows, args.sample, args.sample_seed, args.first)
    modified = _benchlib.pipeline_modified_residues()
    os.makedirs(args.outdir, exist_ok=True)
    manifest = []
    for pdbid in ids:
        # The core-set CSV stores tokenised sequences and SMILES separated by </s>.
        smiles = all_rows[pdbid]["smiles"].split("</s>")[-1].replace(" ", "")
        try:
            info, reason = build(pdbid, smiles, args.outdir, args.hetatm,
                                 args.conformer_seed, modified_residues=modified)
        except Exception as e:  # network or conversion failure: record and go on
            info, reason = None, f"error: {e}"
        if info:
            manifest.append([pdbid, "ok", info["comp"], info["chain"], info["resnum"],
                             info["instances"], info["inchikey_pdb"],
                             info["inchikey_ligand"], info["inchikey_core_set"],
                             info["identity_check"], info["full_key_match"],
                             (f"rdkit conformer not written: {info['conformer_error']}"
                              if info.get("conformer_error") else "")])
            print(f"{pdbid}: {info['comp']} {info['chain']}{info['resnum']} "
                  f"({info['instances']} candidate(s), {info['identity_check']})",
                  flush=True)
        else:
            manifest.append([pdbid, "skipped", "", "", "", 0, "", "", "", "", "", reason])
            print(f"{pdbid}: skipped - {reason}", flush=True)
    with open(os.path.join(args.outdir, "manifest.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pdbid", "status", "comp_id", "chain", "resnum", "instances",
                    "inchikey_pdb", "inchikey_ligand", "inchikey_core_set",
                    "identity_check", "full_inchikey_match", "reason"])
        w.writerows(manifest)
    ok = sum(1 for m in manifest if m[1] == "ok")
    print(f"{ok}/{len(manifest)} complexes written to {args.outdir}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
