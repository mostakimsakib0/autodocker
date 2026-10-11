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
                            --sample-seed 42``, the same selection as
                            ``benchmark_redock.py``) or ``--first N``
  * Alternate locations   : per residue, the blank/'A' conformer is kept
                            (the highest-occupancy one when there is no 'A')
  * Receptor              : ``{pdbid}_protein.pdb`` = ATOM/TER records plus,
                            with ``--hetatm metals`` (default), metal-ion
                            HETATM records (Zn, Mg, Ca, Mn, Fe, ...);
                            ``--hetatm all`` keeps every non-water HETATM
                            except the chosen ligand instance
  * Ligand                : ``{pdbid}_ligand.mol2`` = one crystal instance of
                            the HET group whose heavy-atom formula matches the
                            core-set SMILES, with bond orders from the RCSB
                            ModelServer (Chemical Component Dictionary). The
                            identity is verified against the core-set SMILES:
                            the InChIKey connectivity block must match, or
                            (for SMILES whose charges/hydrogens differ, e.g.
                            charge-separated phosphonates) the heavy-atom
                            graph must be isomorphic. Stereo and protonation
                            are reported (full_inchikey_match), not required
  * Start conformer       : with RDKit available, ``{pdbid}_ligand_rdkit.mol2``
                            (+ ``.json``) = seeded ETKDGv3 + MMFF conformer for
                            ``benchmark_redock.py --start-conformer rdkit``
                            (lets the re-docking run where RDKit is missing)

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


def het_groups(records):
    """Heavy-atom element counts of every candidate ligand HET residue."""
    groups = defaultdict(Counter)
    order = []
    for line in records:
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        if resname in SKIP_RESNAMES:
            continue
        element = _element(line)
        if element in ("H", "D"):
            continue
        key = (resname, line[21], line[22:26].strip())
        if key not in groups:
            order.append(key)
        groups[key][element] += 1
    return [(k, groups[k]) for k in order]


def receptor_records(records, hetatm="metals", ligand_key=None):
    """Receptor lines: ATOM/TER plus the HETATM records selected by
    ``hetatm`` ('none', 'metals', or 'all' = every non-water HETATM except
    the residue ``ligand_key`` = (resname, chain, resnum))."""
    out = []
    for line in records:
        if line.startswith(("ATOM", "TER")):
            out.append(line)
        elif line.startswith("HETATM"):
            resname = line[17:20].strip()
            key = (resname, line[21], line[22:26].strip())
            if hetatm == "metals" and _element(line) in METAL_ELEMENTS \
                    and resname not in WATERS:
                out.append(line)
            elif hetatm == "all" and resname not in WATERS and key != ligand_key:
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


def inchikey(source, is_smiles=False):
    """Standard InChIKey via Open Babel (None when it cannot be computed)."""
    args = ["obabel", f"-:{source}"] if is_smiles else ["obabel", source]
    proc = subprocess.run(args + ["-oinchikey"], capture_output=True, text=True,
                          timeout=_benchlib.TOOL_TIMEOUT)
    key = proc.stdout.strip().split()
    return key[0] if key and len(key[0]) == 27 else None


def same_connectivity(key_a, key_b):
    """True when two InChIKeys share the connectivity (first) block."""
    return bool(key_a and key_b) and key_a.split("-")[0] == key_b.split("-")[0]


def sdf_heavy_graph(sdf_text):
    """Heavy-atom graph of the first V2000 molecule: ``(elements, edges)``.

    Bond orders, charges and hydrogens are ignored, so protonation states
    and charge-separated resonance forms describe the same skeleton.
    """
    lines = sdf_text.split("\n")
    n_atoms, n_bonds = int(lines[3][:3]), int(lines[3][3:6])
    elements = [lines[4 + i][31:34].strip().upper() for i in range(n_atoms)]
    heavy = {i: e for i, e in enumerate(elements) if e not in ("H", "D")}
    edges = set()
    for j in range(n_bonds):
        line = lines[4 + n_atoms + j]
        a, b = int(line[:3]) - 1, int(line[3:6]) - 1
        if a in heavy and b in heavy:
            edges.add((min(a, b), max(a, b)))
    index = {old: new for new, old in enumerate(sorted(heavy))}
    return ([heavy[i] for i in sorted(heavy)],
            sorted((index[a], index[b]) for a, b in edges))


def skeleton_hash(elements, edges, iterations=6):
    """Weisfeiler-Lehman hash of an element-labelled graph. Isomorphic
    graphs always hash equal; for drug-sized molecules a collision between
    different skeletons of the same formula is practically excluded."""
    neighbours = defaultdict(list)
    for a, b in edges:
        neighbours[a].append(b)
        neighbours[b].append(a)
    labels = list(elements)
    for _ in range(iterations):
        labels = [labels[i] + "(" + ",".join(sorted(labels[j] for j in neighbours[i])) + ")"
                  for i in range(len(labels))]
        # compress labels to keep strings short
        table = {lab: f"{k:x}" for k, lab in enumerate(sorted(set(labels)))}
        labels = [elements[i] + table[lab] for i, lab in enumerate(labels)]
    return "|".join(sorted(labels)) + f"#{len(edges)}"


def smiles_sdf(smiles):
    """Connection table of a SMILES string (Open Babel, no coordinates)."""
    return subprocess.run(["obabel", f"-:{smiles}", "-osdf"], capture_output=True,
                          text=True, timeout=_benchlib.TOOL_TIMEOUT).stdout


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


def build(pdbid, smiles, outdir, hetatm="metals", conformer_seed=42):
    pdb_text = fetch(PDB_URL.format(pdbid=pdbid.upper())).decode()
    records = filter_altlocs(first_model_records(pdb_text))
    target = smiles_formula(smiles)
    matches = [key for key, counts in het_groups(records) if counts == target]
    if not matches:
        return None, "no HET group matches the ligand formula (peptide or multi-residue ligand)"
    ref_key = inchikey(smiles, is_smiles=True)
    ref_skeleton = skeleton_hash(*sdf_heavy_graph(smiles_sdf(smiles)))
    cdir = os.path.join(outdir, pdbid)
    os.makedirs(cdir, exist_ok=True)
    sdf_path = os.path.join(cdir, f"{pdbid}_ligand.sdf")
    chosen, last_key, identity = None, None, None
    for comp, chain, resnum in matches:
        sdf = fetch(LIGAND_URL.format(pdbid=pdbid.lower(), comp=comp, chain=chain,
                                      resnum=resnum)).decode()
        sdf = normalize_sdf_elements(sdf.split("$$$$")[0] + "$$$$\n")
        with open(sdf_path, "w", encoding="utf-8") as f:
            f.write(sdf)
        last_key = inchikey(sdf_path)
        if same_connectivity(last_key, ref_key):
            chosen, identity = (comp, chain, resnum), "inchikey-connectivity"
            break
        if skeleton_hash(*sdf_heavy_graph(sdf)) == ref_skeleton:
            chosen, identity = (comp, chain, resnum), "heavy-atom-graph"
            break
    if chosen is None:
        os.remove(sdf_path)
        return None, (f"formula matches {len(matches)} HET group(s) but neither the "
                      f"InChIKey connectivity block nor the heavy-atom graph matches "
                      f"the core-set SMILES (core set {ref_key}, PDB {last_key})")
    mol2 = os.path.join(cdir, f"{pdbid}_ligand.mol2")
    subprocess.run(["obabel", sdf_path, "-omol2", "-O", mol2],
                   check=True, capture_output=True)
    with open(os.path.join(cdir, f"{pdbid}_protein.pdb"), "w", encoding="utf-8") as f:
        for line in receptor_records(records, hetatm, ligand_key=chosen):
            f.write(line + "\n")
        f.write("END\n")
    info = {"comp": chosen[0], "chain": chosen[1], "resnum": chosen[2],
            "instances": len(matches), "inchikey_pdb": last_key,
            "identity_check": identity,
            "inchikey_core_set": ref_key,
            "full_key_match": last_key == ref_key}
    try:
        rd_mol2 = os.path.join(cdir, f"{pdbid}_ligand_rdkit.mol2")
        meta = _benchlib.rdkit_start_conformer(sdf_path, rd_mol2, conformer_seed)
        with open(rd_mol2[:-len(".mol2")] + ".json", "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)
    except ImportError:
        pass
    except Exception as e:  # the crystal-start benchmark still works
        info["conformer_error"] = str(e)[:200]
    return info, None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ids", default=DEFAULT_IDS, help="CSV with pdbid,smiles columns")
    ap.add_argument("--sample", type=int, default=40,
                    help="Seeded random sample size (default 40; 0 = all IDs)")
    ap.add_argument("--sample-seed", type=int, default=42)
    ap.add_argument("--first", "--max-complexes", dest="first", type=int, default=None,
                    help="Legacy 1.1.x selection: first N IDs in sorted order")
    ap.add_argument("--hetatm", choices=["none", "metals", "all"], default="metals",
                    help="HETATM records kept in the receptor (default metals)")
    ap.add_argument("--conformer-seed", type=int, default=42,
                    help="ETKDG seed of the precomputed RDKit start conformer")
    args = ap.parse_args(argv)

    with open(args.ids, encoding="utf-8") as f:
        all_rows = {r["pdbid"].lower(): r for r in csv.DictReader(f)}
    if args.first:
        ids = sorted(all_rows)[:args.first]
    else:
        ids = _benchlib.select_sample(all_rows, args.sample, args.sample_seed)
    os.makedirs(args.outdir, exist_ok=True)
    manifest = []
    for pdbid in ids:
        # The core-set CSV stores tokenised sequences and SMILES separated by </s>.
        smiles = all_rows[pdbid]["smiles"].split("</s>")[-1].replace(" ", "")
        try:
            info, reason = build(pdbid, smiles, args.outdir, args.hetatm,
                                 args.conformer_seed)
        except Exception as e:  # network or conversion failure: record and go on
            info, reason = None, f"error: {e}"
        if info:
            manifest.append([pdbid, "ok", info["comp"], info["chain"], info["resnum"],
                             info["instances"], info["inchikey_pdb"],
                             info["identity_check"], info["full_key_match"], ""])
            print(f"{pdbid}: {info['comp']} {info['chain']}{info['resnum']} "
                  f"({info['instances']} candidate(s))", flush=True)
        else:
            manifest.append([pdbid, "skipped", "", "", "", 0, "", "", "", reason])
            print(f"{pdbid}: skipped - {reason}", flush=True)
    with open(os.path.join(args.outdir, "manifest.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["pdbid", "status", "comp_id", "chain", "resnum", "instances",
                    "inchikey", "identity_check", "full_inchikey_match", "reason"])
        w.writerows(manifest)
    ok = sum(1 for m in manifest if m[1] == "ok")
    print(f"{ok}/{len(manifest)} complexes written to {args.outdir}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
