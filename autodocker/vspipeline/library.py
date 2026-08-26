#!/usr/bin/env python3
"""LIBRARY GENERATION (PubChem + ADMET).

Split out of runner.py (behavior-preserving refactor). All cross-module
references go through the ``runner`` namespace at call time so the public
entry point and its monkeypatch surface stay unchanged.
"""
import os
import re
import csv
import json
import math
import time
import shlex
import uuid
import shutil
import fnmatch
import glob
import html
import statistics
import subprocess
import logging
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import List, Tuple, Optional, Dict
from urllib.parse import quote
from multiprocessing import Pool, cpu_count as mp_cpu_count

# The runner module owns the tool-path globals (OBABEL, VINA*, SMINA,
# COMMAND_TIMEOUT) and every public name; import it lazily-safe: during the
# initial import of runner this binds the partially-initialized module
# object, which is fine because attributes are only accessed at call time.
import runner  # noqa: E402

logger = logging.getLogger(__name__)

class LibraryManager:
    """Manages compound library sourcing and preparation."""

    PUBCHEM_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    MAX_COMPOUNDS = 20

    # Curated FDA-approved drug names (well-known actives) resolved through
    # PubChem PUG REST name->CID lookup.
    FDA_DRUG_NAMES = [
        "aspirin", "ibuprofen", "acetaminophen", "naproxen", "diclofenac",
        "celecoxib", "atorvastatin", "simvastatin", "rosuvastatin", "lovastatin",
        "metformin", "glibenclamide", "glipizide", "metoprolol", "propranolol",
        "amlodipine", "nifedipine", "verapamil", "diltiazem", "losartan",
        "valsartan", "candesartan", "furosemide", "hydrochlorothiazide", "spironolactone",
        "digoxin", "warfarin", "clopidogrel", "ticlopidine", "prasugrel",
        "omeprazole", "lansoprazole", "pantoprazole", "ranitidine", "cimetidine",
        "fluoxetine", "sertraline", "paroxetine", "citalopram", "escitalopram",
        "diazepam", "alprazolam", "lorazepam", "clonazepam", "zolpidem",
        "sildenafil", "tadalafil", "vardenafil", "finasteride", "dutasteride",
        "tamoxifen", "anastrozole", "letrozole", "methotrexate", "cyclophosphamide",
        "prednisone", "dexamethasone", "hydrocortisone", "levothyroxine", "amoxicillin",
        "ciprofloxacin", "azithromycin", "doxycycline", "fluconazole", "acyclovir",
        "valacyclovir", "allopurinol", "colchicine", "gabapentin", "pregabalin",
        "donepezil", "rivastigmine", "memantine", "levodopa", "carbidopa",
    ]

    FALLBACK_LIBRARY_NAMES = list(FDA_DRUG_NAMES)

    def __init__(self, outdir: str, ligands_input_dir: str,
                 minimize: bool = True, minimize_steps: int = 250):
        self.outdir = outdir
        self.ligands_input_dir = ligands_input_dir
        self.lib_dir = os.path.join(outdir, "ligands")
        self.minimized_dir = os.path.join(outdir, "ligands_minimized")
        self.metadata_file = os.path.join(outdir, "ligand_metadata.json")
        self.metadata = {}
        self.minimize = minimize
        self.minimize_steps = minimize_steps
        os.makedirs(self.lib_dir, exist_ok=True)
        os.makedirs(self.minimized_dir, exist_ok=True)
        # Allow overriding the PubChem REST base for tests/mirrors.
        self.pubchem_rest = os.environ.get('PUBCHEM_REST') or self.PUBCHEM_REST

    def create_fda_library(self, apply_admet: bool = True) -> List[str]:
        """Download FDA-approved drugs from PubChem.

        Uses PubChem PUG REST: ``compound/name/{name}/cids`` then
        ``compound/cid/{cid}/SDF`` (the old ZINC endpoint is dead).
        """
        logger.info(
            "[*] FDA library mode - downloading approved drugs from PubChem...")

        downloaded = []
        for name in self.FDA_DRUG_NAMES:
            if len(downloaded) >= self.MAX_COMPOUNDS:
                logger.info(
                    f"[*] Reached {self.MAX_COMPOUNDS} compounds, stopping download")
                break
            cids = self._pubchem_cid_for_name(name)
            if not cids:
                logger.debug(f"  No PubChem CID for {name}")
                continue
            pdbqt_file = self._prepare_pubchem_compound(
                str(cids[0]), runner._safe_ligand_id(name), apply_admet)
            if pdbqt_file:
                downloaded.append(pdbqt_file)
                logger.info(f"  [+] {name} (CID {cids[0]})")

        if not downloaded:
            logger.warning(
                "No FDA compounds could be downloaded; falling back to local SDF")
            return self._prepare_local_sdf(apply_admet=apply_admet)

        self._save_metadata()
        logger.info(
            f"[✔] Downloaded and prepared {len(downloaded)} FDA compounds from PubChem")
        return downloaded

    def create_custom_library(self, mw_min: int = 200, mw_max: int = 500,
                              logp_max: float = 5, apply_admet: bool = True) -> List[str]:
        """Download a custom library filtered by molecular properties.

        PubChem computed properties (MolecularWeight, XLogP) are used for
        the MW/LogP filters instead of the dead ZINC query API.
        """
        logger.info(
            f"[*] Custom library mode (MW: {mw_min}-{mw_max}, LogP: <={logp_max}) - PubChem")

        downloaded = []
        for name in self.FALLBACK_LIBRARY_NAMES:
            if len(downloaded) >= self.MAX_COMPOUNDS:
                logger.info(
                    f"[*] Reached {self.MAX_COMPOUNDS} compounds, stopping download")
                break
            cids = self._pubchem_cid_for_name(name)
            if not cids:
                continue
            cid = str(cids[0])
            props = self._pubchem_properties(cid)
            if props is None:
                continue
            try:
                mw = props.get("MolecularWeight")
                logp = props.get("XLogP")
                if mw is None:
                    continue
                if not (mw_min <= float(mw) <= mw_max):
                    logger.debug(f"  {name}: MW {mw} outside range; skipping")
                    continue
                if logp is not None and float(logp) > logp_max:
                    logger.debug(
                        f"  {name}: LogP {logp} > {logp_max}; skipping")
                    continue
            except (TypeError, ValueError):
                continue
            pdbqt_file = self._prepare_pubchem_compound(
                cid, runner._safe_ligand_id(name), apply_admet)
            if pdbqt_file:
                downloaded.append(pdbqt_file)
                logger.info(f"  [+] {name} (CID {cid})")

        if not downloaded:
            logger.warning(
                "No matching compounds found; falling back to local SDF")
            return self._prepare_local_sdf(apply_admet=apply_admet)

        self._save_metadata()
        logger.info(
            f"[✔] Downloaded and prepared {len(downloaded)} custom compounds from PubChem")
        return downloaded

    def _pubchem_cid_for_name(self, name: str) -> List[str]:
        """Resolve a drug name to PubChem Compound IDs (CIDs)."""
        url = (f"{self.pubchem_rest}/compound/name/"
               f"{quote(name)}/cids/JSON")
        try:
            data = json.loads(runner._http_get_bytes(url).decode("utf-8"))
        except Exception as e:
            logger.debug(f"PubChem name lookup failed for '{name}': {e}")
            return []
        cids = ((data.get("IdentifierList") or {}).get("CID")) or []
        return [str(c) for c in cids]

    def _pubchem_properties(self, cid: str) -> Optional[Dict]:
        """Fetch computed MolecularWeight and XLogP for a CID."""
        url = (f"{self.pubchem_rest}/compound/cid/{cid}/property/"
               f"MolecularWeight,XLogP/JSON")
        try:
            data = json.loads(runner._http_get_bytes(url).decode("utf-8"))
        except Exception as e:
            logger.debug(f"PubChem property lookup failed for CID {cid}: {e}")
            return None
        try:
            props = data["PropertyTable"]["Properties"][0]
        except (KeyError, IndexError):
            return None
        return props

    def _pubchem_download_sdf(self, cid: str, outpath: str) -> str:
        """Download the 3D SDF for a CID to ``outpath``."""
        url = (f"{self.pubchem_rest}/compound/cid/{cid}/SDF"
               f"?record_type=3d")
        content = runner._http_get_bytes(url)
        with open(outpath, "wb") as f:
            f.write(content)
        return outpath

    def _prepare_pubchem_compound(self, cid: str, name: str,
                                  apply_admet: bool) -> Optional[str]:
        """Download a PubChem SDF, apply ADMET and convert to PDBQT."""
        ligand_id = f"{name}_{cid}"
        sdf_file = os.path.join(self.lib_dir, f"{ligand_id}.sdf")
        pdbqt_file = os.path.join(self.lib_dir, f"{ligand_id}.pdbqt")

        try:
            self._pubchem_download_sdf(cid, sdf_file)
        except Exception as e:
            logger.debug(f"PubChem SDF download failed for {ligand_id}: {e}")
            return None

        if apply_admet:
            props = runner.ADMETFilter.parse_sdf_properties(sdf_file)
            passes, violations = runner.ADMETFilter.check_lipinski(props)
            if not passes:
                logger.debug(
                    f"  Skipped {ligand_id}: {', '.join(violations)}")
                try:
                    os.remove(sdf_file)
                except OSError:
                    pass
                return None

        try:
            pdbqt_file = self._prepare_sdf_to_pdbqt(
                sdf_file, pdbqt_file, ligand_id)
        except Exception as e:
            logger.debug(f"PDBQT preparation failed for {ligand_id}: {e}")
            return None

        self.metadata.setdefault(ligand_id, {})["source_url"] = (
            f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}")
        return pdbqt_file

    def _save_metadata(self):
        try:
            with open(self.metadata_file, "w") as f:
                json.dump(self.metadata, f, indent=2)
            logger.info(f"[✔] Ligand metadata saved: {self.metadata_file}")
        except IOError as e:
            logger.warning(f"Could not save ligand metadata: {e}")

    def _prepare_sdf_to_pdbqt(self, sdf_file: str, pdbqt_file: str, ligand_id: str) -> str:
        """Optionally minimize an SDF with MMFF94, then convert to PDBQT."""
        source_for_conversion = sdf_file
        minimized_file = os.path.join(self.minimized_dir, f"{ligand_id}.sdf")

        # Parse properties to get real data
        props = runner.ADMETFilter.parse_sdf_properties(sdf_file)
        if props is None:
            props = {}

        if self.minimize:
            try:
                runner.run([
                    runner.OBABEL, "-isdf", sdf_file, "-osdf", "-O", minimized_file,
                    "--gen3d", "--minimize", "--ff", "MMFF94",
                    "--steps", str(self.minimize_steps)
                ])
                source_for_conversion = minimized_file
            except Exception as e:
                logger.warning(
                    f"MMFF94 minimization failed for {ligand_id}: {e}. Converting original SDF.")

        # Convert to PDBQT
        try:
            runner.run([runner.OBABEL, "-isdf", source_for_conversion, "-opdbqt", "-O", pdbqt_file])
        except Exception as e:
            logger.error(f"Failed to convert {ligand_id} to PDBQT: {e}")
            # Create empty file to mark failure
            with open(pdbqt_file, 'w') as f:
                f.write(f"# Error converting {ligand_id}: {e}\n")
            raise

        # Verify output file was created and has content
        if not os.path.exists(pdbqt_file) or os.path.getsize(pdbqt_file) == 0:
            logger.error(
                f"PDBQT conversion failed for {ligand_id} - output file is empty")
            raise ValueError(
                f"PDBQT file is empty after conversion for {ligand_id}")

        # Verify charges were computed
        if not runner._ensure_pdbqt_has_charges(pdbqt_file):
            logger.debug(
                f"[*] Ligand {ligand_id} PDBQT has no charges - computing Gasteiger...")
            runner._fix_pdbqt_charges(pdbqt_file, source_for_conversion)

        # Store metadata with real values from properties
        self.metadata[ligand_id] = {
            "source_sdf": os.path.abspath(sdf_file),
            "prepared_sdf": os.path.abspath(source_for_conversion),
            "pdbqt": os.path.abspath(pdbqt_file),
            "mw": props.get("mw", 0),
            "logp": props.get("logp", 0),
            "tpsa": props.get("tpsa", 0),
            "rotors": props.get("rotors", 0),
            "hbd": props.get("hbd", 0),
            "hba": props.get("hba", 0),
            "source_url": "",
            "minimized": self.minimize and source_for_conversion == minimized_file,
            "properties_source": "obabel_descriptors" if props else "default",
        }

        return pdbqt_file

    def _prepare_local_sdf(self, apply_admet: bool = True) -> List[str]:
        """Convert local SDF/PDBQT files to ready docking inputs with smart auto-detection."""

        ligands_input = self.ligands_input_dir

        if not os.path.exists(ligands_input):
            raise FileNotFoundError(
                f"Ligands directory not found: {ligands_input}"
            )

        # -----------------------------
        # STEP 1: COLLECT FILES
        # -----------------------------
        sdf_files = []
        pdbqt_files = []
        mol2_files = []
        pdb_files = []

        if os.path.isfile(ligands_input):
            lname = ligands_input.lower()
            if lname.endswith(".sdf"):
                sdf_files = [ligands_input]
            elif lname.endswith(".pdbqt"):
                pdbqt_files = [ligands_input]
            elif lname.endswith(".mol2"):
                mol2_files = [ligands_input]
            elif lname.endswith(".pdb"):
                pdb_files = [ligands_input]
        else:
            for root, _, files in os.walk(ligands_input):
                for f in files:
                    path = os.path.join(root, f)

                    if f.lower().endswith(".sdf"):
                        sdf_files.append(path)
                    elif f.lower().endswith(".pdbqt"):
                        pdbqt_files.append(path)
                    elif f.lower().endswith(".mol2"):
                        mol2_files.append(path)
                    elif f.lower().endswith(".pdb"):
                        pdb_files.append(path)

        # -----------------------------
        # STEP 2: DECIDE MODE
        # -----------------------------
        # Diagnostic: log what we found to help debug empty-result cases
        try:
            logger.info(
                f"[DEBUG] Local ligands found - SDF: {len(sdf_files)}, PDBQT: {len(pdbqt_files)}, MOL2: {len(mol2_files)}, PDB: {len(pdb_files)}")
            sample = (sdf_files or pdbqt_files or mol2_files or pdb_files)[:5]
            if sample:
                logger.info(
                    f"[DEBUG] Sample files: {', '.join([os.path.basename(s) for s in sample])}")
        except Exception:
            pass
        if sdf_files:
            mode = "sdf"
            ligands = sdf_files
            logger.info(f"[*] SDF mode detected: {len(sdf_files)} ligands")

        elif pdbqt_files:
            mode = "pdbqt"
            ligands = pdbqt_files
            logger.info(f"[*] PDBQT mode detected: {len(pdbqt_files)} ligands")

        elif mol2_files:
            mode = "mol2"
            ligands = mol2_files
            logger.info(f"[*] MOL2 mode detected: {len(mol2_files)} ligands")

        elif pdb_files:
            mode = "pdb"
            ligands = pdb_files
            logger.info(f"[*] PDB mode detected: {len(pdb_files)} ligands")

        else:
            raise FileNotFoundError(
                "No ligands found (.sdf/.pdbqt/.mol2/.pdb)")

        # -----------------------------
        # STEP 3: PROCESS
        # -----------------------------
        admet = runner.ADMETFilter()
        out_files = []
        failed_ligands = []

        logger.info(f"[*] Preparing {len(ligands)} ligands...")

        for lig in ligands:

            try:
                # -------------------------
                # CASE A: SDF pipeline
                # -------------------------
                if mode == "sdf":
                    inp = lig
                    name = Path(lig).stem
                    out = os.path.join(self.lib_dir, f"{name}.pdbqt")

                    if apply_admet:
                        props = admet.parse_sdf_properties(inp)
                        if props is None:
                            failed_ligands.append(
                                (lig, "Invalid SDF properties"))
                            continue

                        ok, violations = admet.check_lipinski(props)
                        if not ok:
                            failed_ligands.append(
                                (lig, f"ADMET: {violations}"))
                            continue

                    self._prepare_sdf_to_pdbqt(inp, out, name)
                    out_files.append(out)

                    logger.info(f"  [✓] {name} (SDF) → PDBQT")

                # -------------------------
                # CASE B: PDBQT direct
                # -------------------------
                elif mode == "pdbqt":
                    out_files.append(lig)
                    logger.info(f"  [✓] {Path(lig).stem} (PDBQT direct)")

                # -------------------------
                # CASE C: MOL2 -> convert to PDBQT
                # -------------------------
                elif mode == "mol2":
                    inp = lig
                    name = Path(lig).stem
                    out = os.path.join(self.lib_dir, f"{name}.pdbqt")
                    try:
                        runner.run([runner.OBABEL, "-imol2", inp, "-opdbqt", "-O", out])
                        if not os.path.exists(out) or os.path.getsize(out) == 0:
                            raise RuntimeError(
                                "Conversion produced empty file")
                        if not runner._ensure_pdbqt_has_charges(out):
                            runner._fix_pdbqt_charges(out, inp)
                        out_files.append(out)
                        logger.info(f"  [✓] {name} (MOL2 -> PDBQT)")
                    except Exception as e:
                        failed_ligands.append((lig, f"Conversion failed: {e}"))
                        logger.error(f"  [✗] {lig}: {e}")

                # -------------------------
                # CASE D: PDB -> convert to PDBQT
                # -------------------------
                elif mode == "pdb":
                    inp = lig
                    name = Path(lig).stem
                    out = os.path.join(self.lib_dir, f"{name}.pdbqt")
                    try:
                        runner.run([runner.OBABEL, "-ipdb", inp, "-opdbqt", "-O", out])
                        if not os.path.exists(out) or os.path.getsize(out) == 0:
                            raise RuntimeError(
                                "Conversion produced empty file")
                        if not runner._ensure_pdbqt_has_charges(out):
                            runner._fix_pdbqt_charges(out, inp)
                        out_files.append(out)
                        logger.info(f"  [✓] {name} (PDB -> PDBQT)")
                    except Exception as e:
                        failed_ligands.append((lig, f"Conversion failed: {e}"))
                        logger.error(f"  [✗] {lig}: {e}")

            except Exception as e:
                failed_ligands.append((lig, str(e)))
                logger.error(f"  [✗] {lig}: {e}")

        # -----------------------------
        # STEP 4: FINAL REPORT
        # -----------------------------
        self._save_metadata()

        logger.info(
            f"[✔] Success: {len(out_files)}/{len(ligands)} ligands ready"
        )

        if failed_ligands:
            logger.warning(f"[!] Skipped {len(failed_ligands)} ligands:")
            for lig, reason in failed_ligands:
                logger.warning(f"   - {Path(lig).name}: {reason}")

        if not out_files:
            raise RuntimeError("No valid ligands prepared")

        return out_files
