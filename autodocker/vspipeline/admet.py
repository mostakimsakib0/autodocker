#!/usr/bin/env python3
"""SDF VALIDATION + ADMET FILTERING (Lipinski's Rule).

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

def _is_valid_sdf(sdf_file: str) -> Tuple[bool, str]:
    """Validate if SDF file contains valid molecular data.

    Returns:
        (is_valid, error_message)
    """
    if not os.path.exists(sdf_file):
        return False, f"File not found: {sdf_file}"

    file_size = os.path.getsize(sdf_file)
    if file_size < 100:
        return False, f"File is too small ({file_size} bytes) - likely incomplete or error message"

    try:
        with open(sdf_file, 'r', errors='ignore') as f:
            content = f.read()

            # Check for error messages
            if '404' in content or 'NotFound' in content or 'Status:' in content or 'Message:' in content:
                if len(content) < 500:  # Likely an error response, not real molecule data
                    return False, "SDF file contains error message (likely failed download from PubChem)"

            # Check for molecule block end marker (essential for SDF format)
            if 'M  END' not in content:
                return False, "Invalid SDF format - missing M  END marker"

            # Check for actual atom/bond data (these are required in a valid SDF)
            lines = content.split('\n')

            # Look for the counts line (typically line 3 in SDF after header)
            has_atom_bond_data = False
            # Check first ~100 lines
            for i, line in enumerate(lines[2:min(100, len(lines))]):
                # Counts line format: " 92 96..." (atom count space bond count space...)
                # OR check for coordinate data lines (should have float coordinates)
                if len(line) >= 30:
                    try:
                        # Try to parse coordinates from typical PDB/SDF line
                        coords = line[:30]
                        float(coords[0:10])     # x coordinate
                        float(coords[10:20])    # y coordinate
                        float(coords[20:30])    # z coordinate
                        has_atom_bond_data = True
                        break
                    except (ValueError, IndexError):
                        pass

            if not has_atom_bond_data:
                return False, "Invalid SDF format - no coordinates or atom data found"

    except Exception as e:
        return False, f"Error reading SDF: {str(e)}"

    return True, ""


# =============================
# ADMET FILTERING (Lipinski's Rule)
# =============================
class ADMETFilter:
    """ADMET filtering based on Lipinski's rule of five."""

    @staticmethod
    def parse_sdf_properties(sdf_file: str) -> Optional[Dict]:
        """Extract molecular properties from SDF file.

        Returns:
            Dict with properties or None if SDF is invalid
        """
        # Validate SDF first
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
            'name': os.path.basename(sdf_file)
        }

        # Try to get descriptors from OpenBabel first (most reliable)
        obabel_props = _obabel_descriptors(sdf_file)
        if obabel_props:
            props.update(obabel_props)

        # If OpenBabel failed, try parsing from SDF file
        if not obabel_props:
            try:
                with open(sdf_file) as f:
                    lines = f.readlines()
                    if len(lines) > 0:
                        props['name'] = lines[0].strip() or props['name']

                        for i, line in enumerate(lines):
                            if 'MW' in line or 'MolWt' in line:
                                try:
                                    props['mw'] = float(line.split()[-1])
                                except ValueError:
                                    pass
                            if 'LogP' in line or 'LOGP' in line:
                                try:
                                    props['logp'] = float(line.split()[-1])
                                except ValueError:
                                    pass
                            if 'HBD' in line:
                                try:
                                    props['hbd'] = float(line.split()[-1])
                                except ValueError:
                                    pass
                            if 'HBA' in line or 'HBCount' in line:
                                try:
                                    if not props.get('hba'):
                                        props['hba'] = float(line.split()[-1])
                                except ValueError:
                                    pass
            except Exception as e:
                logger.warning(
                    f"Could not parse SDF properties from {os.path.basename(sdf_file)}: {e}")

        # Return None if critical properties are missing
        if any(v is None for v in [props['mw'], props['logp'], props['hba'], props['hbd']]):
            logger.warning(
                f"Incomplete properties for {os.path.basename(sdf_file)} - could not calculate all descriptors")
            return None

        return props

    @staticmethod
    def check_lipinski(props: Optional[Dict]) -> Tuple[bool, List[str]]:
        """Check Lipinski's rule of five compliance."""
        if props is None:
            return False, ["Invalid SDF file - no properties available"]

        violations = []

        mw = props.get('mw')
        logp = props.get('logp')
        hbd = props.get('hbd')
        hba = props.get('hba')

        if mw is None or logp is None or hba is None or hbd is None:
            return False, ["Missing molecular properties"]

        if mw > 500:
            violations.append(f"MW={mw:.1f} (>500)")
        if logp > 5:
            violations.append(f"LogP={logp:.2f} (>5)")
        if hbd > 5:
            violations.append(f"HBD={hbd:.0f} (>5)")
        if hba > 10:
            violations.append(f"HBA={hba:.0f} (>10)")

        return len(violations) == 0, violations


def _obabel_descriptors(sdf_file: str) -> Optional[Dict]:
    """Calculate useful ligand descriptors using Open Babel when available.

    Returns:
        Dict with descriptors or None if calculation failed
    """
    try:
        result = subprocess.run(
            [runner.OBABEL, sdf_file, "-osmi", "--append",
                "MW logP TPSA rotors HBD HBA2"],
            capture_output=True,
            text=True,
            timeout=runner.COMMAND_TIMEOUT,
        )

        if result.returncode != 0:
            logger.debug(
                f"Open Babel failed for {os.path.basename(sdf_file)}: {result.stderr}")
            return None

        if not result.stdout.strip():
            logger.debug(
                f"Open Babel returned empty output for {os.path.basename(sdf_file)}")
            return None

        tokens = result.stdout.strip().split()
        if len(tokens) < 7:
            logger.debug(
                f"Open Babel output incomplete for {os.path.basename(sdf_file)}: expected >=7 tokens, got {len(tokens)}")
            return None

        try:
            mw, logp, tpsa, rotors, hbd, hba = [float(t) for t in tokens[-6:]]
            return {
                "mw": mw,
                "logp": logp,
                "tpsa": tpsa,
                "rotors": rotors,
                "hbd": hbd,
                "hba": hba,
            }
        except ValueError as e:
            logger.debug(
                f"Could not parse Open Babel descriptors for {os.path.basename(sdf_file)}: {e}")
            return None

    except subprocess.TimeoutExpired:
        logger.debug(f"Open Babel timed out for {os.path.basename(sdf_file)}")
        return None
    except Exception as e:
        logger.debug(
            f"Open Babel descriptor calculation failed for {os.path.basename(sdf_file)}: {e}")
        return None
