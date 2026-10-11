#!/usr/bin/env python3
"""Ligand library preparation: local files and PubChem downloads.

Local mode prepares every SDF/PDBQT/MOL2/PDB file of ``-l``; the PubChem
modes download a fixed list of approved drugs (``--library fda``) or the same
list filtered by a property window (``--library custom``). Every prepared
ligand is recorded in ``ligand_metadata.json`` (source file, descriptors,
minimization, and for PubChem compounds the CID, download URL and retrieval
time).

Cross-module references go through the ``runner`` namespace at call time so
monkeypatching ``runner.<name>`` (e.g. ``runner.run`` in tests) keeps working.
"""
import json
import logging
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import quote

from vspipeline import defaults  # noqa: E402
from vspipeline import outputs  # noqa: E402
import runner  # noqa: E402  (tool paths and shared helpers, resolved at call time)
from vspipeline import conformers, net  # noqa: E402
from vspipeline.admet import (  # noqa: E402
    LIPINSKI_HBA_MAX,
    LIPINSKI_HBD_MAX,
    LIPINSKI_LOGP_MAX,
    LIPINSKI_MW_MAX,
    sdf_coordinates,
)
from vspipeline.charges import _read_coordinates  # noqa: E402
from vspipeline.hydrogens import run_obabel_deterministic  # noqa: E402
from vspipeline.utils import (  # noqa: E402
    ObabelNoOutput,
    _strip_molfile_timestamps,
    elements_without_gasteiger,
    molecule_elements,
)

logger = logging.getLogger(__name__)

# Local ligand formats in order of precedence: when several files share a
# name (e.g. X.mol2 and X.pdbqt in a DUD-E folder), the first format in this
# order is docked and the others are reported as superseded.
FORMAT_PRECEDENCE = ("sdf", "pdbqt", "mol2", "pdb")
_FORMAT_LABEL = {"sdf": "SDF", "pdbqt": "PDBQT", "mol2": "MOL2", "pdb": "PDB"}

# An SDF atom with |z| above this (SDF coordinates have 4 decimals) makes the
# molecule 3D; a molblock with every z at 0 is a 2D depiction.
Z_EPSILON = 1e-4
# A prepared ligand whose coordinates do not spread by more than this along
# some axis (PDBQT coordinates have 3 decimals) is planar-degenerate.
PLANAR_SPREAD_EPSILON = 1e-3
# Up to this many atoms a molecule is planar whatever its geometry (three
# points define a plane), so planarity is no sign of a 2D depiction.
MAX_ALWAYS_PLANAR_ATOMS = 3
# Directories of -l holding one of these files are AutoDocker output
# directories (an earlier run's -o) and are never read as ligands.
OUTPUT_DIR_MARKERS = (outputs.LIGAND_METADATA_JSON, outputs.RUN_MANIFEST_JSON)
# ... as are directories holding both ligand-preparation subdirectories (an
# earlier run interrupted before it wrote ligand_metadata.json).
OUTPUT_SUBDIR_MARKERS = (outputs.LIGANDS_DIR, outputs.LIGANDS_MINIMIZED_DIR)

DEFAULT_MINIMIZE_STEPS = defaults.MINIMIZE_STEPS
DEFAULT_MAX_COMPOUNDS = defaults.MAX_COMPOUNDS
# Seed of the 3D embedding of 2D inputs when the caller gives none (the
# runner passes --seed).
DEFAULT_EMBED_SEED = conformers.EMBED_SEED


class LigandRejected(Exception):
    """A ligand that cannot or must not be docked; the message is the reason."""


class FilterRejected(LigandRejected):
    """A ligand removed by the property filter (not a preparation failure)."""


class PubChemUnreachable(RuntimeError):
    """Several PubChem requests in a row failed at the network level."""


def _ligand_pdbqt_problem(path: str) -> Optional[str]:
    """Return why a PDBQT cannot be docked as a ligand, or None if it looks usable.

    Vina requires exactly one ligand torsion tree (one ROOT ... TORSDOF) per
    file. A receptor-style PDBQT has atoms and charges but no tree, several
    ROOT blocks come from several molecules or disconnected fragments (a
    salt) written to one file, and MODEL records from a multi-model file;
    Vina rejects all of them at parse time ("Unknown or inappropriate tag"),
    so they are caught here with a clear reason instead. Checked for every
    ligand PDBQT, given or prepared.
    """
    counts: Dict[str, int] = {}
    try:
        with open(path, errors="replace", encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    tag = line.split()[0]
                    counts[tag] = counts.get(tag, 0) + 1
    except OSError as e:
        return f"Unreadable PDBQT: {e}"
    if not counts.keys() & {"ATOM", "HETATM"}:
        return "PDBQT has no ATOM/HETATM records"
    missing = [t for t in ("ROOT", "TORSDOF") if t not in counts]
    if missing:
        return ("Not a ligand PDBQT (missing " + "/".join(missing) +
                "; looks like a receptor-style file)")
    if "MODEL" in counts:
        return ("PDBQT contains several models (MODEL records, e.g. a Vina "
                "output file); Vina docks one ligand per file")
    if counts["ROOT"] > 1 or counts["TORSDOF"] > 1:
        return (f"PDBQT holds {counts['ROOT']} ROOT blocks (several molecules or "
                "disconnected fragments, e.g. a salt); Vina docks one molecule per file")
    return None


def _count_models(path: str) -> int:
    """Number of MODEL records of a PDB file."""
    try:
        with open(path, errors="replace", encoding="utf-8") as fh:
            return sum(1 for line in fh if line.startswith("MODEL"))
    except OSError:
        return 0


def _count_molecules(path: str, fmt: str) -> int:
    """Number of molecule records in an SDF ($$$$) or MOL2 (@<TRIPOS>MOLECULE) file."""
    marker = "$$$$" if fmt == "sdf" else "@<TRIPOS>MOLECULE"
    try:
        with open(path, errors="replace", encoding="utf-8") as fh:
            n = sum(1 for line in fh if line.strip() == marker)
    except OSError:
        return 0
    # A single-molecule SDF may omit the trailing $$$$.
    return max(n, 1) if fmt == "sdf" else n


def _sdf_is_usable(path: str) -> bool:
    """True if ``path`` exists and holds at least one complete molblock."""
    try:
        with open(path, errors="replace", encoding="utf-8") as fh:
            return "M  END" in fh.read()
    except OSError:
        return False


def _sdf_has_3d_coords(sdf_file: str) -> bool:
    """True if the first molecule of an SDF (V2000 or V3000) has any atom
    with a non-zero Z coordinate, i.e. it is not a flat 2D depiction."""
    coords = sdf_coordinates(sdf_file)
    return bool(coords) and any(abs(z) > Z_EPSILON for _, _, z in coords)


def _is_flat(coords: Optional[List[Tuple[float, float, float]]]) -> bool:
    """True for coordinates with every z at 0 (a 2D depiction); False when
    there are none."""
    return bool(coords) and all(abs(z) <= Z_EPSILON for _, _, z in coords)


def _is_2d_input(path: str) -> bool:
    """A MOL2/PDB ligand drawn in 2D: more than MAX_ALWAYS_PLANAR_ATOMS
    atoms, every z at 0 (fewer atoms are planar whatever their geometry)."""
    coords = _read_coordinates(path)
    return len(coords) > MAX_ALWAYS_PLANAR_ATOMS and _is_flat(coords)


def _generated_3d_problem(sdf_file: str) -> Optional[str]:
    """Why freshly generated 3D coordinates are unusable, or None.

    Any one, two or three points lie in a plane, and RDKit's embedding puts
    such a molecule (water, CO2, HCN with its hydrogen) exactly in z = 0,
    so for up to MAX_ALWAYS_PLANAR_ATOMS atoms only coincident atoms are a
    problem; a larger molecule must leave the z = 0 plane.
    """
    coords = sdf_coordinates(sdf_file)
    if not coords:
        return "3D coordinate generation wrote no atoms"
    if len(coords) <= MAX_ALWAYS_PLANAR_ATOMS:
        if len(set(coords)) < len(coords):
            return "3D coordinate generation put atoms on top of each other"
        return None
    if _is_flat(coords):
        return "3D coordinate generation produced no 3D coordinates for a 2D input"
    return None


def _pdbqt_coordinate_spread(path: str) -> Optional[Tuple[float, float, float]]:
    """Extent (max - min) of the atom coordinates of a PDBQT along x, y, z."""
    xs, ys, zs = [], [], []
    try:
        with open(path, errors="replace", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
                    except ValueError:
                        continue
    except OSError:
        return None
    if not xs:
        return None
    return max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)


def _flat_axis(path: str) -> Optional[str]:
    """Axis along which a prepared ligand of more than
    MAX_ALWAYS_PLANAR_ATOMS atoms has no extent, else None."""
    spread = _pdbqt_coordinate_spread(path)
    if spread is None:
        return None
    with open(path, errors="replace", encoding="utf-8") as fh:
        n_atoms = sum(1 for line in fh if line.startswith(("ATOM", "HETATM")))
    if n_atoms <= MAX_ALWAYS_PLANAR_ATOMS:
        return None  # any three points are coplanar
    for axis, extent in zip("xyz", spread):
        if extent <= PLANAR_SPREAD_EPSILON:
            return axis
    return None


# PDBQT columns 71-76 (0-based slice 70:76): the partial charge.
PDBQT_CHARGE_COLUMNS = (70, 76)
# Coordinate tolerance (A) when matching the atoms of two conversions.
_ATOM_MATCH_TOLERANCE = 2e-3


def _pdbqt_charge(line: str) -> Optional[float]:
    try:
        return float(line[PDBQT_CHARGE_COLUMNS[0]:PDBQT_CHARGE_COLUMNS[1]])
    except ValueError:
        return None


def _zero_charge_polar_hydrogens(path: str) -> bool:
    """True if a PDBQT has a polar hydrogen (HD) with charge exactly 0:
    Gasteiger never gives a hydrogen on N/O/S a zero charge, so it was added
    after the charges were assigned."""
    with open(path, errors="replace", encoding="utf-8") as fh:
        return any(line.startswith(("ATOM", "HETATM")) and line[77:79].strip() == "HD"
                   and _pdbqt_charge(line) == 0.0 for line in fh)


def _same_pdbqt_atom(a: str, b: str) -> bool:
    """Same atom name and type, and (for heavy atoms) the same position, in
    two PDBQT records. Hydrogen positions may differ (they do not affect
    Gasteiger charges, which depend on the bonds only)."""
    if a[12:16] != b[12:16] or a[77:79].strip() != b[77:79].strip():
        return False
    if a[77:79].strip() in ("H", "HD", "HS"):
        return True
    try:
        return all(abs(float(a[k:k + 8]) - float(b[k:k + 8])) <= _ATOM_MATCH_TOLERANCE
                   for k in (30, 38, 46))
    except ValueError:
        return False


def _unsupported_elements_hint(fmt: str, path: str, seed: Optional[int] = None) -> str:
    """' (element(s) X have no ...)' naming the elements of ``path`` that
    Open Babel cannot assign Gasteiger charges to, or ''."""
    bad = elements_without_gasteiger(molecule_elements(path, fmt, seed=seed))
    if not bad:
        return ""
    return (f" (element(s) {', '.join(bad)} have no Gasteiger charge parameters in "
            "Open Babel, which then converts nothing; AutoDock Vina has no atom type "
            "for most of them either)")


def _is_output_dir(path: str) -> bool:
    """True if ``path`` looks like an AutoDocker output directory (see
    OUTPUT_DIR_MARKERS and OUTPUT_SUBDIR_MARKERS)."""
    if any(os.path.isfile(os.path.join(path, m)) for m in OUTPUT_DIR_MARKERS):
        return True
    return all(os.path.isdir(os.path.join(path, d)) for d in OUTPUT_SUBDIR_MARKERS)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class LibraryManager:
    """Sources ligands (local files or PubChem) and prepares docking PDBQTs."""

    PUBCHEM_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    # Human-readable compound page, recorded as each compound's source_url.
    PUBCHEM_COMPOUND_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/"
    MAX_COMPOUNDS = DEFAULT_MAX_COMPOUNDS
    MAX_NETWORK_FAILURES = 3

    # A fixed, hand-picked list of well-known approved drugs, resolved by
    # name through PubChem PUG REST (the first CID PubChem returns for the
    # name is used). It is NOT the FDA approved-drug list and not a
    # representative sample of it; screens take the first --max-compounds
    # entries that can be prepared, in this order. Ranitidine was removed in
    # 1.2.0 (withdrawn from the US market in April 2020, NDMA contamination);
    # the list is otherwise unchanged.
    FDA_DRUG_NAMES = [
        "aspirin", "ibuprofen", "acetaminophen", "naproxen", "diclofenac",
        "celecoxib", "atorvastatin", "simvastatin", "rosuvastatin", "lovastatin",
        "metformin", "glibenclamide", "glipizide", "metoprolol", "propranolol",
        "amlodipine", "nifedipine", "verapamil", "diltiazem", "losartan",
        "valsartan", "candesartan", "furosemide", "hydrochlorothiazide", "spironolactone",
        "digoxin", "warfarin", "clopidogrel", "ticlopidine", "prasugrel",
        "omeprazole", "lansoprazole", "pantoprazole", "cimetidine",
        "fluoxetine", "sertraline", "paroxetine", "citalopram", "escitalopram",
        "diazepam", "alprazolam", "lorazepam", "clonazepam", "zolpidem",
        "sildenafil", "tadalafil", "vardenafil", "finasteride", "dutasteride",
        "tamoxifen", "anastrozole", "letrozole", "methotrexate", "cyclophosphamide",
        "prednisone", "dexamethasone", "hydrocortisone", "levothyroxine", "amoxicillin",
        "ciprofloxacin", "azithromycin", "doxycycline", "fluconazole", "acyclovir",
        "valacyclovir", "allopurinol", "colchicine", "gabapentin", "pregabalin",
        "donepezil", "rivastigmine", "memantine", "levodopa", "carbidopa",
    ]

    # The custom library is the same list with a user-defined property window.
    FALLBACK_LIBRARY_NAMES = list(FDA_DRUG_NAMES)

    def __init__(self, outdir: str, ligands_input_dir: Optional[str],
                 minimize: bool = True, minimize_steps: int = DEFAULT_MINIMIZE_STEPS,
                 seed: int = DEFAULT_EMBED_SEED):
        self.outdir = outdir
        self.seed = seed
        self.ligands_input_dir = ligands_input_dir
        self.lib_dir = os.path.join(outdir, outputs.LIGANDS_DIR)
        self.minimized_dir = os.path.join(outdir, outputs.LIGANDS_MINIMIZED_DIR)
        self.metadata_file = os.path.join(outdir, outputs.LIGAND_METADATA_JSON)
        self.metadata: Dict[str, Dict] = {}
        self.minimize = minimize
        self.minimize_steps = minimize_steps
        # (name, reason) of ligands that could not be prepared, and of
        # same-name files that were not used (see FORMAT_PRECEDENCE).
        self.failed: List[Tuple[str, str]] = []
        self.superseded: List[Tuple[str, str]] = []
        os.makedirs(self.lib_dir, exist_ok=True)
        os.makedirs(self.minimized_dir, exist_ok=True)
        # Allow overriding the PubChem REST base for tests/mirrors.
        self.pubchem_rest = os.environ.get('PUBCHEM_REST') or self.PUBCHEM_REST
        self._network_failures = 0

    def _out_relpath(self, path: str) -> str:
        """Path relative to the output directory when it lies inside it (as
        ranking.csv's Docked_File), so metadata does not depend on where the
        output directory is; absolute otherwise (e.g. a PDBQT used as given)."""
        out = os.path.abspath(self.outdir)
        full = os.path.abspath(path)
        if os.path.commonpath([out, full]) == out:
            return os.path.relpath(full, out)
        return full

    # ------------------------------------------------------------------
    # PubChem libraries
    # ------------------------------------------------------------------

    def create_fda_library(self, apply_lipinski: bool = True, *,
                           apply_admet: Optional[bool] = None,
                           max_compounds: Optional[int] = None) -> List[str]:
        """Download and prepare compounds of ``FDA_DRUG_NAMES`` from PubChem.

        With ``apply_lipinski`` the Lipinski filter (zero violations) is
        applied to each downloaded SDF. ``apply_admet`` is a deprecated alias.
        Raises RuntimeError if no compound could be prepared; local ligands
        are never used as a substitute.
        """
        if apply_admet is not None:
            apply_lipinski = apply_admet
        limits = ({} if not apply_lipinski else
                  dict(mw_max=LIPINSKI_MW_MAX, logp_max=LIPINSKI_LOGP_MAX,
                       hbd_max=LIPINSKI_HBD_MAX, hba_max=LIPINSKI_HBA_MAX))
        logger.info("[*] FDA library mode - downloading approved drugs from PubChem...")
        return self._build_pubchem_library("fda", self.FDA_DRUG_NAMES, limits,
                                           max_compounds)

    def create_custom_library(self, mw_min: float = defaults.CUSTOM_MW_MIN,
                              mw_max: float = defaults.CUSTOM_MW_MAX,
                              logp_max: float = defaults.CUSTOM_LOGP_MAX,
                              apply_lipinski: bool = True, *,
                              apply_admet: Optional[bool] = None,
                              max_compounds: Optional[int] = None) -> List[str]:
        """Download ``FALLBACK_LIBRARY_NAMES`` from PubChem and keep the
        compounds inside a property window.

        The MW window and logP maximum REPLACE Lipinski's MW <= 500 and
        logP <= 5 cut-offs; with ``apply_lipinski`` Lipinski's HBD <= 5 and
        HBA <= 10 are applied as well. All four descriptors are computed by
        Open Babel from the downloaded structure (one logP model, recorded in
        ligand_metadata.json as ``logp_model``).
        """
        if apply_admet is not None:
            apply_lipinski = apply_admet
        limits = dict(mw_min=mw_min, mw_max=mw_max, logp_max=logp_max,
                      hbd_max=LIPINSKI_HBD_MAX if apply_lipinski else None,
                      hba_max=LIPINSKI_HBA_MAX if apply_lipinski else None)
        logger.info(
            f"[*] Custom library mode (MW {mw_min:g}-{mw_max:g}, logP <= {logp_max:g}"
            f"{', HBD <= 5, HBA <= 10' if apply_lipinski else ''}) - PubChem")
        return self._build_pubchem_library("custom", self.FALLBACK_LIBRARY_NAMES,
                                           limits, max_compounds)

    def _build_pubchem_library(self, mode: str, names: List[str], limits: Dict,
                               max_compounds: Optional[int]) -> List[str]:
        limit = self.MAX_COMPOUNDS if max_compounds is None else max_compounds
        downloaded: List[str] = []
        rejected: List[Tuple[str, str]] = []
        self.failed = []
        for name in names:
            if len(downloaded) >= limit:
                logger.info(f"[*] Reached {limit} compounds, stopping download")
                break
            cids = self._pubchem_cid_for_name(name)
            if not cids:
                self._fail(name, "no PubChem CID found for this name")
                continue
            try:
                pdbqt = self._prepare_pubchem_compound(
                    cids[0], runner._safe_ligand_id(name), limits, query=name)
            except PubChemUnreachable:
                raise
            except FilterRejected as e:
                rejected.append((name, str(e)))
                logger.info(f"  [-] {name} (CID {cids[0]}) rejected: {e}")
                continue
            except (LigandRejected, RuntimeError, OSError, ValueError) as e:
                self._fail(name, f"CID {cids[0]}: {e}")
                continue
            downloaded.append(pdbqt)
            logger.info(f"  [+] {name} (CID {cids[0]})")

        self._save_metadata()
        self._log_failures()
        if not downloaded:
            raise RuntimeError(
                f"No compound of the --library {mode} list could be prepared from "
                f"PubChem ({len(rejected)} rejected by the property filter, "
                f"{len(self.failed)} failed; see the warnings above). Ligands in "
                f"-l are not used in --library {mode} mode.")
        logger.info(f"[*] Downloaded and prepared {len(downloaded)} {mode} "
                    f"compounds from PubChem ({len(rejected)} rejected by the "
                    "property filter)")
        return downloaded

    def _pubchem_fetch(self, url: str) -> bytes:
        """GET a PubChem URL, counting consecutive network-level failures.

        After ``MAX_NETWORK_FAILURES`` in a row (connection errors, timeouts
        or persistent HTTP 429/5xx) PubChemUnreachable is raised instead of
        walking the whole list offline. Other errors (e.g. HTTP 404) are
        re-raised for the caller without changing the count.
        """
        try:
            content = runner._http_get_bytes(url)
        except Exception as e:
            if net.is_network_failure(e):
                self._network_failures += 1
                if self._network_failures >= self.MAX_NETWORK_FAILURES:
                    raise PubChemUnreachable(
                        f"PubChem is unreachable ({self._network_failures} failed "
                        f"requests in a row; last error: {e}). Check the network "
                        "connection or proxy, or use --library local with "
                        "ligand files in -l") from e
            raise
        self._network_failures = 0
        return content

    def _pubchem_cid_for_name(self, name: str) -> List[str]:
        """Resolve a drug name to PubChem Compound IDs (CIDs); [] if unknown."""
        url = f"{self.pubchem_rest}/compound/name/{quote(name)}/cids/JSON"
        try:
            data = json.loads(self._pubchem_fetch(url).decode("utf-8"))
        except PubChemUnreachable:
            raise
        except Exception as e:
            logger.warning(f"  [!] PubChem name lookup failed for '{name}': {e}")
            return []
        cids = ((data.get("IdentifierList") or {}).get("CID")) or []
        return [str(c) for c in cids]

    def _pubchem_download_sdf(self, cid: str, outpath: str) -> Dict[str, str]:
        """Download the SDF of a CID to ``outpath``; return its provenance.

        The 3D conformer record is preferred; PubChem has none for some
        (large or very flexible) compounds, whose 2D record is used instead
        and gets 3D coordinates during preparation.
        """
        url = f"{self.pubchem_rest}/compound/cid/{cid}/SDF?record_type=3d"
        record_type = "3d"
        try:
            content = self._pubchem_fetch(url)
        except Exception as e:
            if net._http_status(e) != 404:
                raise
            logger.info(f"  [*] CID {cid}: no 3D conformer in PubChem; using the 2D record")
            url = f"{self.pubchem_rest}/compound/cid/{cid}/SDF"
            record_type = "2d"
            content = self._pubchem_fetch(url)
        retrieved_at = _utc_now()
        with open(outpath, "wb") as f:
            f.write(content)
        return {"download_url": url, "pubchem_record_type": record_type,
                "retrieved_at": retrieved_at}

    def _prepare_pubchem_compound(self, cid: str, name: str, limits: Dict,
                                  query: Optional[str] = None) -> str:
        """Download a PubChem SDF, apply the property filter, convert to PDBQT.

        Raises FilterRejected when the filter rejects the compound, and
        LigandRejected/RuntimeError/OSError when it cannot be downloaded or
        prepared.
        """
        ligand_id = f"{name}_{cid}"
        sdf_file = os.path.join(self.lib_dir, f"{ligand_id}.sdf")
        pdbqt_file = os.path.join(self.lib_dir, f"{ligand_id}.pdbqt")
        try:
            provenance = self._pubchem_download_sdf(cid, sdf_file)
        except PubChemUnreachable:
            raise
        except Exception as e:
            raise RuntimeError(f"PubChem SDF download failed ({e})") from e

        source, removed = self._largest_fragment("sdf", sdf_file, ligand_id)
        props = None
        if limits:
            props = runner.ADMETFilter.parse_sdf_properties(source, seed=self.seed)
            passes, violations = runner.ADMETFilter.check_lipinski(props, **limits)
            if not passes:
                for path in {sdf_file, source}:
                    try:
                        os.remove(path)
                    except OSError:
                        pass
                raise FilterRejected("property filter: " + ", ".join(violations))

        self._prepare_sdf_to_pdbqt(source, pdbqt_file, ligand_id, props=props,
                                   original=sdf_file, fragments_removed=removed)
        self.metadata[ligand_id].update({
            "source": "pubchem",
            "pubchem_query": query or name,
            "pubchem_cid": cid,
            # The first CID PubChem returned for the name; records can
            # change over time, so the retrieval time is recorded too.
            "source_url": f"{self.PUBCHEM_COMPOUND_URL}{cid}",
            **provenance,
        })
        return pdbqt_file

    def _save_metadata(self):
        try:
            with open(self.metadata_file, "w", encoding="utf-8") as f:
                json.dump(self.metadata, f, indent=2, sort_keys=True)
            logger.info(f"[*] Ligand metadata saved: {self.metadata_file}")
        except IOError as e:
            logger.warning(f"Could not save ligand metadata: {e}")

    def _fail(self, name: str, reason: str) -> None:
        self.failed.append((name, reason))
        logger.warning(f"  [!] {name}: {reason}")

    def _log_failures(self) -> None:
        if self.failed:
            logger.warning(f"[!] {len(self.failed)} ligand(s) could not be prepared:")
            for name, reason in self.failed:
                logger.warning(f"   - {name}: {reason}")

    # ------------------------------------------------------------------
    # Conversion
    # ------------------------------------------------------------------

    def _convert_to_pdbqt(self, fmt: str, src: str, out: str,
                          extra_args: Tuple[str, ...] = ()) -> str:
        """Convert ``src`` (format ``fmt``) to a ligand PDBQT with Open Babel.

        The conversion runs through ``run_obabel_deterministic``: hydrogens
        Open Babel places at random are re-placed deterministically, so the
        same input always gives the same PDBQT. A failed conversion, a
        missing or empty output, or an output Vina cannot dock as one
        ligand (see ``_ligand_pdbqt_problem``) raises (and the file is
        removed); an element Open Babel has no Gasteiger parameters for is
        named. If every atom of the result carries a zero charge the
        charges are rebuilt from ``src``, which must therefore be a 3D
        structure (``--gen3d`` is never used on this path). Polar hydrogens
        added by ``-h`` after the charges were assigned get their Gasteiger
        charges (see ``_charge_added_hydrogens``).
        """
        if os.path.exists(out):
            os.remove(out)
        try:
            run_obabel_deterministic([f"-i{fmt}", src, "-opdbqt", "-O", out, *extra_args],
                                     out, kind="ligand", seed=self.seed)
        except RuntimeError as e:
            raise RuntimeError(f"{e}{_unsupported_elements_hint(fmt, src, self.seed)}") from e
        if not os.path.exists(out) or os.path.getsize(out) == 0:
            if os.path.exists(out):
                os.remove(out)
            raise RuntimeError("Open Babel produced an empty PDBQT"
                               + _unsupported_elements_hint(fmt, src, self.seed))
        if not runner._ensure_pdbqt_has_charges(out):
            try:
                runner._fix_pdbqt_charges(out, src, gen3d=False, seed=self.seed)
            except RuntimeError as e:
                raise RuntimeError(f"{e}{_unsupported_elements_hint(fmt, src, self.seed)}") from e
        elif _zero_charge_polar_hydrogens(out) and fmt in ("sdf", "mol2"):
            self._charge_added_hydrogens(fmt, src, out)
        problem = _ligand_pdbqt_problem(out)
        if problem:
            os.remove(out)
            raise LigandRejected(problem)
        return out

    def _charge_added_hydrogens(self, fmt: str, src: str, out: str) -> None:
        """Gasteiger charges for a PDBQT whose polar hydrogens were added
        (``-h``) after its charges: Open Babel assigns the charges before it
        adds hydrogens, so they kept 0.000. The molecule is protonated first
        and charged in a separate conversion, and only the charge column of
        ``out`` is replaced (coordinates, order, names, types and the
        torsion tree are unchanged). With MOL2 input the charges given in
        the file (computed without those hydrogens) are replaced by Gasteiger
        charges of the protonated molecule. If the two conversions do not
        match atom for atom, ``out`` is left as it is with a warning.
        """
        with tempfile.TemporaryDirectory(prefix="autodocker_") as tmp:
            protonated = os.path.join(tmp, f"protonated.{fmt}")
            charged = os.path.join(tmp, "charged.pdbqt")
            runner.run_obabel([f"-i{fmt}", src, f"-o{fmt}", "-O", protonated, "-h"],
                              outputs=[protonated], seed=self.seed)
            runner.run_obabel([f"-i{fmt}", protonated, "-opdbqt", "-O", charged,
                               "--partialcharge", "gasteiger"], outputs=[charged],
                              seed=self.seed)
            with open(charged, encoding="utf-8") as f:
                reference = [line for line in f if line.startswith(("ATOM", "HETATM"))]
        with open(out, encoding="utf-8") as f:
            lines = f.readlines()
        atoms = [i for i, line in enumerate(lines) if line.startswith(("ATOM", "HETATM"))]
        if len(atoms) != len(reference) or not all(
                _same_pdbqt_atom(lines[i], ref) for i, ref in zip(atoms, reference)):
            logger.warning(f"[!] {os.path.basename(out)}: could not match the atoms of "
                           "the separately charged molecule; added polar hydrogens "
                           "keep charge 0.000")
            return
        for i, ref in zip(atoms, reference):
            lines[i] = lines[i][:PDBQT_CHARGE_COLUMNS[0]] + \
                ref[PDBQT_CHARGE_COLUMNS[0]:PDBQT_CHARGE_COLUMNS[1]] + \
                lines[i][PDBQT_CHARGE_COLUMNS[1]:]
        with open(out, "w", encoding="utf-8") as f:
            f.writelines(lines)
        logger.info(f"[*] {os.path.basename(out)}: Gasteiger charges computed after "
                    "adding hydrogens (added polar hydrogens had 0.000)")

    def _generate_3d(self, sdf_file: str, ligand_id: str) -> Tuple[str, str]:
        """3D coordinates for the 2D molecule of ``sdf_file`` with
        ``conformers.generate_3d`` (RDKit ETKDGv3 seeded with the run seed +
        MMFF94 when RDKit is installed, else ``obabel --gen3d``, which is not
        reproducible). Returns (generated SDF, method); raises
        LigandRejected when no usable 3D structure results."""
        gen3d_file = os.path.join(self.minimized_dir, f"{ligand_id}_gen3d.sdf")
        logger.info(f"[*] {ligand_id}: input has no 3D coordinates; generating them")
        if os.path.exists(gen3d_file):
            os.remove(gen3d_file)
        try:
            method = conformers.generate_3d(
                sdf_file, gen3d_file, seed=self.seed, run_obabel=runner.run_obabel)
        except (ValueError, RuntimeError) as e:
            raise LigandRejected(f"3D coordinate generation failed: {e}") from e
        problem = _generated_3d_problem(gen3d_file)
        if problem:
            raise LigandRejected(problem)
        return gen3d_file, method

    def _reject_if_flat(self, pdbqt_file: str) -> None:
        """Reject a prepared ligand built from a 2D input that is still
        planar (see ``_flat_axis``)."""
        axis = _flat_axis(pdbqt_file)
        if axis:
            os.remove(pdbqt_file)
            raise LigandRejected(
                f"prepared ligand is planar (no extent along {axis}) although "
                "the 2D input should have received 3D coordinates; refusing "
                "to dock a flat structure")

    def _prepare_sdf_to_pdbqt(self, sdf_file: str, pdbqt_file: str, ligand_id: str,
                              props: Optional[Dict] = None, *,
                              original: Optional[str] = None,
                              fragments_removed: Optional[List[str]] = None) -> str:
        """Prepare an SDF ligand: 3D coordinates, optional MMFF94, PDBQT.

        * A 2D input (every z = 0) first gets 3D coordinates (see
          ``_generate_3d``); the method is recorded as
          ``coordinates_method``. The planar input is never converted or
          used as a charge source, and a result that is still planar is
          rejected.
        * MMFF94 minimization (``--minimize --ff MMFF94``) of the 3D structure
          runs when enabled; a failure is logged as a warning and recorded in
          the metadata.
        * The PDBQT gets hydrogens (``-h``) and Gasteiger charges.

        ``props`` are the descriptors if the caller already computed them.
        ``original`` is the input file when ``sdf_file`` is its largest
        fragment (see ``_largest_fragment``), whose removed fragments are
        ``fragments_removed``.
        """
        if props is None:
            props = runner.ADMETFilter.parse_sdf_properties(sdf_file, seed=self.seed) or {}
        input_is_2d = not _sdf_has_3d_coords(sdf_file)
        source = sdf_file
        coordinates = "input"
        minimized = False
        minimization_error = None

        coordinates_method = None
        if input_is_2d:
            source, coordinates_method = self._generate_3d(sdf_file, ligand_id)
            coordinates = "gen3d"

        if self.minimize:
            minimized_file = os.path.join(self.minimized_dir, f"{ligand_id}.sdf")
            cmd = ["-isdf", source, "-osdf", "-O", minimized_file,
                   "--minimize", "--ff", "MMFF94", "--steps", str(self.minimize_steps)]
            try:
                if os.path.exists(minimized_file):
                    os.remove(minimized_file)
                # Open Babel exits 0 even when MMFF94 cannot be set up (e.g.
                # boron compounds) and then writes nothing, so check the file.
                try:
                    runner.run_obabel(cmd, outputs=[minimized_file], seed=self.seed)
                except ObabelNoOutput:
                    raise RuntimeError("Open Babel wrote no minimized molecule "
                                       "(MMFF94 setup failed?)") from None
                if not _sdf_is_usable(minimized_file):
                    raise RuntimeError("Open Babel wrote no minimized molecule "
                                       "(MMFF94 setup failed?)")
                source = minimized_file
                minimized = True
            except (RuntimeError, OSError) as e:
                minimization_error = str(e)
                if os.path.exists(minimized_file):
                    os.remove(minimized_file)
                logger.warning(f"[!] {ligand_id}: MMFF94 minimization failed ({e}); "
                               "continuing with the unminimized structure")

        # Gasteiger charges are requested here: Open Babel's SDF reader marks
        # charges as known, so the PDBQT writer would otherwise write 0.000
        # everywhere (and the charge rebuild from the source would follow).
        self._convert_to_pdbqt("sdf", source, pdbqt_file,
                               ("-h", "--partialcharge", "gasteiger"))
        if input_is_2d:
            self._reject_if_flat(pdbqt_file)

        self.metadata[ligand_id] = {
            "source_sdf": self._out_relpath(original or sdf_file),
            "prepared_sdf": self._out_relpath(source),
            "pdbqt": self._out_relpath(pdbqt_file),
            "input_format": "sdf",
            "coordinates": coordinates,
            "coordinates_method": coordinates_method,
            "mw": props.get("mw"),
            "logp": props.get("logp"),
            "tpsa": props.get("tpsa"),
            "rotors": props.get("rotors"),
            "hbd": props.get("hbd"),
            "hba": props.get("hba"),
            "logp_model": props.get("logp_model"),
            "source_url": "",
            "minimized": minimized,
            "minimization_error": minimization_error,
            "properties_source": props.get("properties_source") or "unavailable",
        }
        if fragments_removed:
            self.metadata[ligand_id]["fragments_removed"] = list(fragments_removed)
        return pdbqt_file

    def _prepare_pdb_to_pdbqt(self, src: str, out: str, title: Optional[str] = None) -> str:
        """PDB ligand -> PDBQT with hydrogens added BEFORE Gasteiger charges.

        In one ``obabel -h --partialcharge gasteiger`` call Open Babel
        computes the charges first, so the added polar hydrogens kept 0.000.
        Writing the protonated molecule to a temporary PDB first gives every
        atom its charge; coordinates, atom types and the torsion tree are
        unchanged. ``--title`` sets the name record to ``title`` (default:
        the file name of ``src``, never its directory, so identical files in
        different folders give identical PDBQTs).
        """
        with tempfile.TemporaryDirectory(prefix="autodocker_") as tmp:
            protonated = os.path.join(tmp, "protonated.pdb")
            run_obabel_deterministic(["-ipdb", src, "-opdb", "-O", protonated, "-h"],
                                     protonated, kind="ligand", seed=self.seed)
            return self._convert_to_pdbqt(
                "pdb", protonated, out,
                ("--partialcharge", "gasteiger", "--title",
                 title if title is not None else os.path.basename(src)))

    def _prepare_flat_to_pdbqt(self, fmt: str, src: str, out: str, ligand_id: str,
                               title: Optional[str] = None) -> str:
        """A 2D (every z = 0) MOL2 or PDB ligand -> PDBQT through the same
        deterministic 3D generation as a 2D SDF (``_generate_3d``; never
        ``--gen3d``). Returns the coordinates method. The flat molecule is
        written as ``<ligand_id>_2d.sdf`` in the minimized directory first.
        Charges are Gasteiger charges of the generated structure (MOL2
        charges are not kept: they belong to the flat input). A PDB has no
        bond orders, and Open Babel perceives them from the flat geometry;
        when RDKit cannot read the result the ligand is rejected."""
        flat_sdf = os.path.join(self.minimized_dir, f"{ligand_id}_2d.sdf")
        title_args = ["--title", title] if title is not None else []
        runner.run_obabel([f"-i{fmt}", src, "-osdf", "-O", flat_sdf, *title_args],
                          outputs=[flat_sdf], seed=self.seed)
        try:
            gen3d_file, method = self._generate_3d(flat_sdf, ligand_id)
        except LigandRejected as e:
            hint = (" (bond orders perceived from a flat PDB are unreliable; supply "
                    "3D coordinates or an SDF/MOL2 with bond orders)" if fmt == "pdb" else "")
            raise LigandRejected(f"2D {fmt.upper()} input: {e}{hint}") from e
        self._convert_to_pdbqt("sdf", gen3d_file, out, ("-h", "--partialcharge", "gasteiger"))
        self._reject_if_flat(out)
        return method

    def _largest_fragment(self, fmt: str, path: str, ligand_id: str) -> Tuple[str, List[str]]:
        """Strip counter-ions and other disconnected fragments.

        Returns (file to prepare, canonical SMILES of the removed fragments).
        A single-fragment input is returned unchanged. Otherwise the
        fragment with the most heavy atoms is written (same format, input
        title kept) to ``<ligand_id>_fragment.<fmt>`` in the minimized
        directory and the others are removed, with a logged note; ties for
        the largest fragment (e.g. CCO.CCO, a 1:1 co-crystal) are rejected,
        as there is no way to tell which one is the ligand. Fragments are
        found by Open Babel (bonds as written, or perceived for PDB).
        """
        stdout = runner.run_obabel([f"-i{fmt}", path, "-ocan"], seed=self.seed)
        first = (stdout.splitlines() or [""])[0]
        smiles, _, title = first.partition("\t")
        if "." not in smiles.strip():
            return path, []
        stdout = runner.run_obabel([f"-i{fmt}", path, "--separate", "-d", "-ocan",
                                    "--append", "atoms"], seed=self.seed)
        fragments = []
        for line in stdout.splitlines():
            tokens = line.split()
            if len(tokens) >= 2:
                try:
                    fragments.append((tokens[0], int(tokens[-1])))
                except ValueError:
                    pass
        if len(fragments) < 2:
            raise LigandRejected(
                f"input has disconnected fragments ({smiles.strip()}) that Open Babel "
                "could not separate")
        largest = max(n for _, n in fragments)
        winners = [i for i, (_, n) in enumerate(fragments) if n == largest]
        if len(winners) > 1:
            raise LigandRejected(
                f"input holds {len(fragments)} disconnected fragments and "
                f"{len(winners)} of them share the largest size ({largest} heavy atoms: "
                f"{', '.join(sorted(fragments[i][0] for i in winners))}); cannot tell "
                "which one is the ligand. Keep one molecule per file")
        keep = winners[0]
        removed = sorted(smi for i, (smi, _) in enumerate(fragments) if i != keep)
        target = os.path.join(self.minimized_dir, f"{ligand_id}_fragment.{fmt}")
        with tempfile.TemporaryDirectory(prefix="autodocker_") as tmp:
            base = os.path.join(tmp, f"frag.{fmt}")
            runner.run_obabel([f"-i{fmt}", path, f"-o{fmt}", "-O", base, "--separate",
                               "-m", "--title", title.strip()], seed=self.seed)
            written = os.path.join(tmp, f"frag{keep + 1}.{fmt}")
            if not os.path.exists(written) or os.path.getsize(written) == 0:
                raise LigandRejected("Open Babel could not write the largest fragment")
            shutil.copyfile(written, target)
        if fmt == "sdf":
            _strip_molfile_timestamps(target)
        logger.warning(
            f"[!] {ligand_id}: kept the largest of {len(fragments)} disconnected "
            f"fragments ({fragments[keep][0]}, {largest} heavy atoms); removed "
            f"{', '.join(removed)} (counter-ions/solvent; recorded as fragments_removed)")
        return target, removed

    # ------------------------------------------------------------------
    # Local ligands
    # ------------------------------------------------------------------

    def _collect_local_files(self) -> List[Tuple[str, str]]:
        """(format, path) of every ligand file of ``-l``, grouped by format
        in FORMAT_PRECEDENCE order and sorted within each format."""
        ligands_input = self.ligands_input_dir
        if not ligands_input or not os.path.exists(ligands_input):
            raise FileNotFoundError(f"Ligands directory not found: {ligands_input}")

        by_format: Dict[str, List[str]] = {fmt: [] for fmt in FORMAT_PRECEDENCE}

        def add(path: str) -> None:
            fmt = os.path.splitext(path)[1].lower().lstrip(".")
            if fmt in by_format:
                by_format[fmt].append(path)

        if os.path.isfile(ligands_input):
            add(ligands_input)
        elif os.path.realpath(ligands_input) == os.path.realpath(self.outdir):
            raise ValueError(
                f"The output directory ({self.outdir}) is the ligand directory; "
                "choose a different -o (outputs would be read back as ligands)")
        elif _is_output_dir(ligands_input):
            raise ValueError(
                f"The ligand directory {ligands_input} is an AutoDocker output directory "
                f"(it holds {' or '.join(OUTPUT_DIR_MARKERS)} or "
                f"{'/ and '.join(OUTPUT_SUBDIR_MARKERS)}/); pass the original ligand "
                f"directory, or its {outputs.LIGANDS_DIR}/ subdirectory to re-dock the "
                "prepared PDBQTs")
        else:
            # Never prepare this run's own outputs, nor those of an earlier run
            # (an -o inside -l would otherwise feed ligands/*.pdbqt,
            # docked/*.pdbqt and receptor.pdbqt back in as ligands).
            out_real = os.path.realpath(self.outdir)
            skipped_outdir = False
            skipped_earlier = []
            for root, dirs, files in os.walk(ligands_input):
                kept = []
                for d in sorted(dirs):
                    full = os.path.join(root, d)
                    if os.path.realpath(full) == out_real:
                        skipped_outdir = True
                    elif _is_output_dir(full):
                        skipped_earlier.append(full)
                    else:
                        kept.append(d)
                dirs[:] = kept
                for f in sorted(files):
                    add(os.path.join(root, f))
            if skipped_outdir:
                logger.warning(
                    f"[!] The output directory {self.outdir} is inside the ligand "
                    f"directory {ligands_input}; its contents are not used as ligands")
            for d in skipped_earlier:
                logger.warning(
                    f"[!] {d} inside the ligand directory is the output directory of an "
                    "earlier AutoDocker run; its contents are not used as ligands")
        return [(fmt, path) for fmt in FORMAT_PRECEDENCE for path in by_format[fmt]]

    def _prepare_local_file(self, fmt: str, path: str, apply_lipinski: bool) -> str:
        """Prepare one local ligand file; return the PDBQT to dock.

        Raises LigandRejected (or RuntimeError/OSError from the conversion)
        with the reason when the file cannot be used. Files with several
        molecules (SDF/MOL2 records, PDB models) are rejected; counter-ions
        and other disconnected fragments are removed first (see
        ``_largest_fragment``).
        """
        name = Path(path).stem
        out = os.path.join(self.lib_dir, f"{name}.pdbqt")

        if fmt in ("sdf", "mol2"):
            n_mols = _count_molecules(path, fmt)
            if n_mols > 1:
                raise LigandRejected(
                    f"File holds {n_mols} molecules; Vina docks one per file (split "
                    f"it, e.g. obabel {Path(path).name} -O {name}_.{fmt} -m)")
        elif fmt == "pdb":
            n_models = _count_models(path)
            if n_models > 1:
                raise LigandRejected(
                    f"File holds {n_models} models (MODEL records); Vina docks one "
                    f"molecule per file (split it, e.g. obabel {Path(path).name} "
                    f"-O {name}_.pdb -m)")

        if fmt == "pdbqt":
            problem = _ligand_pdbqt_problem(path)
            if problem:
                raise LigandRejected(problem)
            source, removed = path, []
        else:
            source, removed = self._largest_fragment(fmt, path, name)

        if fmt == "sdf":
            props = None
            if apply_lipinski:
                props = runner.ADMETFilter.parse_sdf_properties(source, seed=self.seed)
                if props is None:
                    raise LigandRejected("Invalid SDF or incomplete descriptors")
                ok, violations = runner.ADMETFilter.check_lipinski(props)
                if not ok:
                    raise LigandRejected("Lipinski: " + ", ".join(violations))
            return self._prepare_sdf_to_pdbqt(source, out, name, props=props,
                                              original=path, fragments_removed=removed)

        coordinates, coordinates_method = "input", None
        if fmt == "pdbqt":
            out = path
        elif _is_2d_input(source):
            # A 2D MOL2/PDB gets 3D coordinates exactly like a 2D SDF.
            title = os.path.basename(path) if fmt == "pdb" else None
            coordinates_method = self._prepare_flat_to_pdbqt(fmt, source, out, name, title)
            coordinates = "gen3d"
        elif fmt == "mol2":
            self._convert_to_pdbqt("mol2", source, out, ("-h",))
        else:  # pdb
            self._prepare_pdb_to_pdbqt(source, out, title=os.path.basename(path))

        self.metadata[name] = {
            "source_file": os.path.abspath(path),
            "input_format": fmt,
            "pdbqt": self._out_relpath(out),
            "coordinates": coordinates,
            "coordinates_method": coordinates_method,
            "mw": None, "logp": None, "tpsa": None, "rotors": None,
            "hbd": None, "hba": None, "logp_model": None,
            "source_url": "",
            "minimized": False,
            "minimization_error": None,
            "properties_source": "unavailable",
        }
        if removed:
            self.metadata[name]["fragments_removed"] = removed
        return out

    def prepare_local(self, apply_lipinski: bool = True, *,
                      apply_admet: Optional[bool] = None) -> List[str]:
        """Prepare the local SDF/PDBQT/MOL2/PDB ligands of ``-l`` for docking.

        Files sharing a name (stem) are resolved by FORMAT_PRECEDENCE on
        first sight: the highest-precedence file is used (even if it then
        fails) and the others are reported as superseded, not as failures
        and not in the success count. The Lipinski filter and MMFF94
        minimization apply to SDF inputs only. ``apply_admet`` is a
        deprecated alias of ``apply_lipinski``.
        """
        if apply_admet is not None:
            apply_lipinski = apply_admet
        ligands = self._collect_local_files()
        if not ligands:
            raise FileNotFoundError("No ligands found (.sdf/.pdbqt/.mol2/.pdb)")

        counts = {fmt: sum(1 for f, _ in ligands if f == fmt) for fmt in FORMAT_PRECEDENCE}
        logger.info("[*] Ligand formats detected: " + ", ".join(
            f"{_FORMAT_LABEL[fmt]} {n}" for fmt, n in counts.items() if n))
        logger.debug("Sample ligand files: " + ", ".join(
            os.path.basename(p) for _, p in ligands[:5]))
        non_sdf = any(fmt != "sdf" for fmt, _ in ligands)
        if apply_lipinski and non_sdf:
            logger.info("[*] The Lipinski filter applies to SDF inputs only; "
                        "PDBQT, MOL2 and PDB ligands are docked unfiltered")
        if self.minimize and non_sdf:
            logger.info("[*] MMFF94 minimization applies to SDF inputs only; "
                        "PDBQT, MOL2 and PDB ligands keep their input coordinates "
                        "(a 2D MOL2/PDB gets seeded RDKit 3D coordinates)")

        self.failed = []
        self.superseded = []
        out_files: List[str] = []
        owners: Dict[str, Tuple[str, str]] = {}  # stem -> (format, path)
        for fmt, path in ligands:
            stem = Path(path).stem
            if stem in owners:
                owner_fmt, owner = owners[stem]
                if owner_fmt == fmt:
                    self.failed.append((path, f"Duplicate ligand name '{stem}' "
                                              f"(already taken by {owner})"))
                else:
                    self.superseded.append((path, owner))
                continue
            owners[stem] = (fmt, path)
            try:
                out_files.append(self._prepare_local_file(fmt, path, apply_lipinski))
            except (LigandRejected, RuntimeError, OSError, ValueError) as e:
                self.failed.append((path, str(e)))
                logger.error(f"  [x] {Path(path).name}: {e}")
                continue
            logger.info(f"  [ok] {stem} ({_FORMAT_LABEL[fmt]})")

        self._save_metadata()
        candidates = len(ligands) - len(self.superseded)
        logger.info(f"[*] Success: {len(out_files)}/{candidates} ligands ready")
        if self.superseded:
            logger.warning(
                f"[!] {len(self.superseded)} file(s) superseded by a same-name file "
                "of higher precedence (SDF > PDBQT > MOL2 > PDB) and not used:")
            for path, owner in self.superseded:
                logger.info(f"   - {Path(path).name}: superseded (same name as "
                            f"{Path(owner).name})")
        if self.failed:
            logger.warning(f"[!] Skipped {len(self.failed)} ligands:")
            for path, reason in self.failed:
                logger.warning(f"   - {Path(path).name}: {reason}")
        if not out_files:
            raise RuntimeError("No valid ligands prepared")
        return out_files
