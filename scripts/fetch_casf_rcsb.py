#!/usr/bin/env python3
"""Build the re-docking benchmark set from the public RCSB PDB.

The PDBbind/CASF-2016 files need a registered download. This script rebuilds
the same complexes from the RCSB PDB instead, so anyone can regenerate the
benchmark input without an account:

  * IDs and ligand SMILES : ``benchmarks/data/casf2016_core_set.csv``
                            (first ``--max-complexes`` IDs in sorted order)
  * Receptor              : ``{pdbid}_protein.pdb`` = the entry's ATOM records
                            (first model)
  * Ligand                : ``{pdbid}_ligand.mol2`` = the crystal instance of
                            the HET group whose heavy-atom formula matches the
                            core-set SMILES, with bond orders from the RCSB
                            ModelServer (Chemical Component Dictionary)

Ligands that are polymers (peptides) or span several residues have no single
matching HET group and are skipped. ``manifest.csv`` records the chosen
component, chain and residue, or the reason for skipping.

Usage:
  python3 scripts/fetch_casf_rcsb.py --outdir DATA_DIR [--max-complexes 40]
then
  python3 scripts/benchmark_redock.py --data-dir DATA_DIR ...
"""

import argparse
import csv
import os
import subprocess
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_IDS = os.path.join(HERE, "..", "benchmarks", "data", "casf2016_core_set.csv")
PDB_URL = "https://files.rcsb.org/download/{pdbid}.pdb"
LIGAND_URL = ("https://models.rcsb.org/v1/{pdbid}/ligand?label_comp_id={comp}"
              "&auth_asym_id={chain}&auth_seq_id={resnum}&encoding=sdf")
SKIP_RESNAMES = {"HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "ACT", "DMS",
                 "PEG", "CL", "NA", "K", "MG", "ZN", "CA", "MN", "FE", "CU",
                 "NI", "CO", "CD", "IOD", "BR"}


def fetch(url, retries=4):
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "autodocker-bench"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise
        except (urllib.error.URLError, ConnectionError, TimeoutError):
            pass
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"could not download {url}")


def smiles_formula(smiles):
    """Heavy-atom element counts of a SMILES string."""
    out = subprocess.run(["obabel", f"-:{smiles}", "-d", "-otxt", "--append", "formula"],
                         capture_output=True, text=True).stdout.strip()
    return formula_counts(out)


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


def het_groups(pdb_text):
    """Heavy-atom element counts of every HET residue in the first model."""
    groups = defaultdict(Counter)
    order = []
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        if resname in SKIP_RESNAMES:
            continue
        element = (line[76:78].strip() or line[12:14].strip()).upper()
        if element in ("H", "D"):
            continue
        key = (resname, line[21], line[22:26].strip())
        if key not in groups:
            order.append(key)
        groups[key][element] += 1
    return [(k, groups[k]) for k in order]


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


def build(pdbid, smiles, outdir):
    pdb_text = fetch(PDB_URL.format(pdbid=pdbid.upper())).decode()
    target = smiles_formula(smiles)
    matches = [key for key, counts in het_groups(pdb_text) if counts == target]
    if not matches:
        return None, "no HET group matches the ligand formula (peptide or multi-residue ligand)"
    comp, chain, resnum = matches[0]
    sdf = fetch(LIGAND_URL.format(pdbid=pdbid.lower(), comp=comp, chain=chain,
                                  resnum=resnum)).decode()
    sdf = normalize_sdf_elements(sdf.split("$$$$")[0] + "$$$$\n")
    cdir = os.path.join(outdir, pdbid)
    os.makedirs(cdir, exist_ok=True)
    sdf_path = os.path.join(cdir, f"{pdbid}_ligand.sdf")
    with open(sdf_path, "w") as f:
        f.write(sdf)
    mol2 = os.path.join(cdir, f"{pdbid}_ligand.mol2")
    subprocess.run(["obabel", sdf_path, "-omol2", "-O", mol2],
                   check=True, capture_output=True)
    with open(os.path.join(cdir, f"{pdbid}_protein.pdb"), "w") as f:
        for line in pdb_text.splitlines():
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "TER")):
                f.write(line + "\n")
        f.write("END\n")
    return (comp, chain, resnum, len(matches)), None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ids", default=DEFAULT_IDS, help="CSV with pdbid,smiles columns")
    ap.add_argument("--max-complexes", type=int, default=40,
                    help="first N IDs in sorted order (the benchmark's default subset)")
    args = ap.parse_args()

    with open(args.ids) as f:
        rows = sorted(csv.DictReader(f), key=lambda r: r["pdbid"])[:args.max_complexes]
    os.makedirs(args.outdir, exist_ok=True)
    manifest = []
    for row in rows:
        pdbid = row["pdbid"].lower()
        # The core-set CSV stores tokenised sequences and SMILES separated by </s>.
        smiles = row["smiles"].split("</s>")[-1].replace(" ", "")
        try:
            chosen, reason = build(pdbid, smiles, args.outdir)
        except Exception as e:  # network or conversion failure: record and go on
            chosen, reason = None, f"error: {e}"
        if chosen:
            comp, chain, resnum, n = chosen
            manifest.append([pdbid, "ok", comp, chain, resnum, n, ""])
            print(f"{pdbid}: {comp} {chain}{resnum} ({n} instance(s))", flush=True)
        else:
            manifest.append([pdbid, "skipped", "", "", "", 0, reason])
            print(f"{pdbid}: skipped - {reason}", flush=True)
    with open(os.path.join(args.outdir, "manifest.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pdbid", "status", "comp_id", "chain", "resnum", "instances", "reason"])
        w.writerows(manifest)
    ok = sum(1 for m in manifest if m[1] == "ok")
    print(f"{ok}/{len(manifest)} complexes written to {args.outdir}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
