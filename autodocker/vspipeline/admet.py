#!/usr/bin/env python3
"""SDF validation and the Lipinski rule-of-five property filter.

This is a drug-likeness filter, not an ADMET model: a ligand passes when it
has ZERO violations of Lipinski's four cut-offs (MW <= 500 Da, logP <= 5,
H-bond donors <= 5, H-bond acceptors <= 10). Lipinski et al. (1997) only
flagged compounds with more than one violation, so this filter is stricter.

Descriptors come from Open Babel (``obabel --append``):

* ``MW``   molecular weight;
* ``logP`` Open Babel's atom-contribution logP (Wildman & Crippen, 1999);
* ``HBD``  JoelLib donor count, SMARTS ``[!#6;!H0]`` (heavy non-carbon atoms
  carrying at least one H);
* ``HBA2`` JoelLib acceptor count, SMARTS
  ``[$([$([#8,#16]);!$(*=N~O);!$(*~N=O);X1,X2]),$([#7;v3;!$([nH]);!$(*(-a)-a)])]``
  (O and S with at most two connections, excluding nitro/nitroso oxygens, and
  trivalent N other than pyrrole-type [nH] and diaryl amines).

The donor/acceptor counts are typed SMARTS counts, not Lipinski's original
"number of N + O" and "number of NH + OH" sums; they are kept because they
are what every earlier AutoDocker release used, so a filter decision never
changes between releases. When Open Babel cannot compute them, PubChem-style
SDF data items (e.g. ``PUBCHEM_XLOGP3``) are used and the source is recorded.

Cross-module references go through the ``runner`` namespace at call time so
monkeypatching ``runner.<name>`` keeps working.
"""
import logging
import os
import re
from typing import Dict, List, Optional, Tuple

import runner  # noqa: E402  (tool paths and run_obabel, resolved at call time)

logger = logging.getLogger(__name__)

# Lipinski rule-of-five cut-offs (a value above the cut-off is a violation).
LIPINSKI_MW_MAX = 500.0
LIPINSKI_LOGP_MAX = 5.0
LIPINSKI_HBD_MAX = 5
LIPINSKI_HBA_MAX = 10

# Human-readable name of the logP model, recorded in ligand_metadata.json.
LOGP_MODEL_OBABEL = "Open Babel logP (Wildman-Crippen atom contributions)"
LOGP_MODEL_SDF = "SDF data item (as supplied, e.g. PubChem XLogP3)"


def sdf_coordinates(sdf_file: str) -> Optional[List[Tuple[float, float, float]]]:
    """Atom coordinates of the first molblock of an SDF (V2000 or V3000).

    Returns None when the file cannot be read or holds no parsable atom
    block. A V2000 counts line gives the atom count; a V3000 CTAB lists the
    atoms between ``M  V30 BEGIN ATOM`` and ``M  V30 END ATOM``.
    """
    try:
        with open(sdf_file, errors="replace", encoding="utf-8") as f:
            lines = f.read().split("$$$$")[0].splitlines()
    except OSError:
        return None
    if len(lines) < 4:
        return None
    counts = lines[3]
    coords: List[Tuple[float, float, float]] = []
    try:
        if "V3000" in counts:
            in_atoms = False
            for line in lines[4:]:
                body = line[6:].strip() if line.startswith("M  V30") else None
                if body is None:
                    continue
                if body.startswith("BEGIN ATOM"):
                    in_atoms = True
                elif body.startswith("END ATOM"):
                    break
                elif in_atoms:
                    # "index type x y z aamap [props]"
                    fields = body.split()
                    coords.append((float(fields[2]), float(fields[3]), float(fields[4])))
        else:
            n_atoms = int(counts[:3])
            for line in lines[4:4 + n_atoms]:
                coords.append((float(line[0:10]), float(line[10:20]), float(line[20:30])))
            if len(coords) != n_atoms:
                return None
    except (ValueError, IndexError):
        return None
    return coords or None


def _is_valid_sdf(sdf_file: str) -> Tuple[bool, str]:
    """Check that an SDF holds a complete molblock with atom coordinates.

    Returns:
        (is_valid, error_message)
    """
    if not os.path.exists(sdf_file):
        return False, f"File not found: {sdf_file}"

    # No size threshold: a complete molblock can be tiny (methane with its
    # hydrogens is ~90 bytes). Truncated downloads and error pages are caught
    # by the parse checks below (no "M  END", no atom block).
    if os.path.getsize(sdf_file) == 0:
        return False, "File is empty"
    try:
        with open(sdf_file, errors="replace", encoding="utf-8") as f:
            content = f.read()
    except OSError as e:
        return False, f"Error reading SDF: {e}"

    # A molblock always ends with "M  END"; an HTML/JSON error page saved as
    # .sdf (e.g. a failed PubChem download) does not.
    if "M  END" not in content:
        return False, "Invalid SDF format - missing M  END marker (not a molblock)"
    if sdf_coordinates(sdf_file) is None:
        return False, "Invalid SDF format - no atom coordinates found"
    return True, ""


def _parse_sdf_data_items(sdf_file: str, props: Dict) -> bool:
    """Fill descriptors from SDF data items ("> <TAG>" then the value line).

    Returns True if all four Lipinski descriptors were found.
    """
    with open(sdf_file, errors="replace", encoding="utf-8") as f:
        lines = f.read().split("$$$$")[0].splitlines()
    if lines:
        props['name'] = lines[0].strip() or props['name']
    tags = {}
    for i, line in enumerate(lines[:-1]):
        m = re.match(r"^>\s*.*<([^>]+)>", line)
        if m:
            tags[m.group(1).strip().upper()] = lines[i + 1].strip()
    for key, names in _SDF_PROPERTY_TAGS.items():
        for tag in names:
            if tag in tags:
                try:
                    props[key] = float(tags[tag])
                    break
                except ValueError:
                    pass
    return all(props[k] is not None for k in ("mw", "logp", "hbd", "hba"))


class LipinskiFilter:
    """Lipinski rule-of-five filter (zero violations allowed)."""

    @staticmethod
    def parse_sdf_properties(sdf_file: str, seed: Optional[int] = None) -> Optional[Dict]:
        """Descriptors of an SDF ligand, or None if it is invalid or the
        four Lipinski descriptors cannot all be determined. ``seed`` is
        Open Babel's OB_RANDOM_SEED (see ``utils.run_obabel``).

        ``properties_source`` records where they came from and
        ``logp_model`` which logP model produced the logP value.
        """
        is_valid, error_msg = _is_valid_sdf(sdf_file)
        if not is_valid:
            logger.warning(
                f"Invalid SDF file {os.path.basename(sdf_file)}: {error_msg}")
            return None

        props = {
            'mw': None,
            'logp': None,
            'hbd': None,
            'hba': None,
            'tpsa': None,
            'rotors': None,
            'name': os.path.basename(sdf_file),
            'properties_source': None,
            'logp_model': None,
        }

        obabel_props = _obabel_descriptors(sdf_file, seed=seed)
        if obabel_props:
            props.update(obabel_props)
            props['properties_source'] = "obabel_descriptors"
            props['logp_model'] = LOGP_MODEL_OBABEL
        else:
            try:
                if _parse_sdf_data_items(sdf_file, props):
                    props['properties_source'] = "sdf_data_items"
                    props['logp_model'] = LOGP_MODEL_SDF
            except OSError as e:
                logger.warning(
                    f"Could not parse SDF properties from {os.path.basename(sdf_file)}: {e}")

        if any(props[k] is None for k in ("mw", "logp", "hba", "hbd")):
            logger.warning(
                f"Incomplete properties for {os.path.basename(sdf_file)} - "
                "could not calculate all descriptors")
            return None
        return props

    @staticmethod
    def check_lipinski(props: Optional[Dict],
                       mw_max: Optional[float] = LIPINSKI_MW_MAX,
                       logp_max: Optional[float] = LIPINSKI_LOGP_MAX,
                       hbd_max: Optional[float] = LIPINSKI_HBD_MAX,
                       hba_max: Optional[float] = LIPINSKI_HBA_MAX,
                       mw_min: Optional[float] = None) -> Tuple[bool, List[str]]:
        """Return (passes, violations); a ligand passes with zero violations.

        The defaults are Lipinski's cut-offs; the custom PubChem library
        passes its own MW window and logP maximum instead. A cut-off of None
        is not checked.
        """
        if props is None:
            return False, ["Invalid SDF file - no properties available"]

        mw = props.get('mw')
        logp = props.get('logp')
        hbd = props.get('hbd')
        hba = props.get('hba')
        if mw is None or logp is None or hba is None or hbd is None:
            return False, ["Missing molecular properties"]

        violations = []
        if mw_min is not None and mw < mw_min:
            violations.append(f"MW={mw:.1f} (<{mw_min:g})")
        if mw_max is not None and mw > mw_max:
            violations.append(f"MW={mw:.1f} (>{mw_max:g})")
        if logp_max is not None and logp > logp_max:
            violations.append(f"LogP={logp:.2f} (>{logp_max:g})")
        if hbd_max is not None and hbd > hbd_max:
            violations.append(f"HBD={hbd:.0f} (>{hbd_max:g})")
        if hba_max is not None and hba > hba_max:
            violations.append(f"HBA={hba:.0f} (>{hba_max:g})")
        return len(violations) == 0, violations


# Deprecated alias: the filter was called "ADMET" before 2.0.0 although it
# only applies the four Lipinski cut-offs.
ADMETFilter = LipinskiFilter


# SDF data-item tags (upper case) that hold each descriptor, in order of
# preference; covers common exporters and PubChem downloads.
_SDF_PROPERTY_TAGS = {
    "mw": ["MW", "MOLWT", "MOLECULAR_WEIGHT", "MOLECULARWEIGHT", "PUBCHEM_MOLECULAR_WEIGHT"],
    "logp": ["LOGP", "XLOGP3", "XLOGP", "PUBCHEM_XLOGP3", "PUBCHEM_XLOGP3_AA"],
    "hbd": ["HBD", "HBOND_DONOR", "PUBCHEM_CACTVS_HBOND_DONOR"],
    "hba": ["HBA", "HBOND_ACCEPTOR", "PUBCHEM_CACTVS_HBOND_ACCEPTOR"],
}

# Open Babel descriptors appended to the SMILES output, in this order.
_OBABEL_DESCRIPTORS = ("MW", "logP", "TPSA", "rotors", "HBD", "HBA2")
_OBABEL_KEYS = ("mw", "logp", "tpsa", "rotors", "hbd", "hba")


def _obabel_descriptors(sdf_file: str, seed: Optional[int] = None) -> Optional[Dict]:
    """Descriptors of the first molecule computed by Open Babel, or None."""
    name = os.path.basename(sdf_file)
    try:
        stdout = runner.run_obabel(
            [sdf_file, "-osmi", "--append", " ".join(_OBABEL_DESCRIPTORS)], seed=seed)
    except (RuntimeError, OSError) as e:
        logger.debug(f"Open Babel descriptor calculation failed for {name}: {e}")
        return None

    first = stdout.strip().splitlines()[0] if stdout.strip() else ""
    tokens = first.split()
    # SMILES, (title,) then one token per descriptor.
    if len(tokens) < len(_OBABEL_DESCRIPTORS) + 1:
        logger.debug(f"Open Babel output incomplete for {name}: {first!r}")
        return None
    try:
        values = [float(t) for t in tokens[-len(_OBABEL_DESCRIPTORS):]]
    except ValueError as e:
        logger.debug(f"Could not parse Open Babel descriptors for {name}: {e}")
        return None
    return dict(zip(_OBABEL_KEYS, values))
