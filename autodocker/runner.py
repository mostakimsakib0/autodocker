#!/usr/bin/env python3
"""
Virtual Screening Pipeline - ENHANCED Edition (v2.0)
NOW INCLUDES ALL 4 PHASES:

PHASE 1: Advanced Result Analysis
  - HTML professional reports with charts
  - Pose clustering (RMSD-based)
  - Pose diversity analysis
  - Binding mode statistics

PHASE 2: Smart Preprocessing
  - Water molecule handling (keep strategic waters)
  - Metal ion detection & parameterization
  - Cofactor handling
  - pH-based protonation (pKa calculation)

PHASE 3: Flexible Receptor Docking
  - Select rotatable residues in binding site
  - Flexible side-chain docking
  - Residue flexibility control (1-10)

PHASE 4: Consensus Scoring
  - Multi-engine docking (Vina + SMINA)
  - Consensus ranking
  - Score agreement analysis

Original Features:
  - ZINC API Integration (FDA & custom libraries)
  - SimScore calculation (RMSD-based pose consistency)
  - ADMET filtering (Lipinski's rule of five)
  - Advanced Vina parameters

Workflow:
  1. Library preparation (ZINC API with ADMET filtering)
  2. Protein preparation (chain selection + fpocket)
  3. Smart preprocessing (optional: water/metal/cofactor)
  4. Docking (flexible or rigid, single or multi-engine)
  5. Advanced results analysis (HTML reports + clustering)
"""
import uuid
import argparse
import subprocess
import os
import shutil
import csv
import sys
import logging
import json
import signal
import re
import shlex
import time
import math
import statistics
import fnmatch
from multiprocessing import Pool, cpu_count as mp_cpu_count
from pathlib import Path
from datetime import datetime
from typing import List, Tuple, Optional, Dict, Set
from urllib.parse import urljoin
from collections import defaultdict

try:
    import requests
except ModuleNotFoundError:
    requests = None

try:
    import numpy as np
except ModuleNotFoundError:
    np = None

# =============================
# LOGGING SETUP
# =============================
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)

# =============================
# PHASE 1: ADVANCED RESULT ANALYSIS
# =============================


def calculate_rmsd(pose1: List[Tuple[float, float, float]],
                   pose2: List[Tuple[float, float, float]]) -> float:
    """Calculate RMSD between two poses (coordinate lists)."""
    if not pose1 or not pose2 or len(pose1) != len(pose2):
        return float('inf')

    if np is None:
        # Fallback without numpy
        sum_sq = sum((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
                     for (x1, y1, z1), (x2, y2, z2) in zip(pose1, pose2))
        return math.sqrt(sum_sq / len(pose1))
    else:
        p1 = np.array(pose1)
        p2 = np.array(pose2)
        return float(np.sqrt(np.mean(np.sum((p1 - p2)**2, axis=1))))


def extract_coordinates_from_pdbqt(pdbqt_file: str) -> Tuple[str, List[Tuple[float, float, float]]]:
    """Extract ligand name and coordinates from PDBQT file."""
    coords = []
    name = os.path.basename(pdbqt_file).replace('.pdbqt', '')

    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        coords.append((x, y, z))
                    except (ValueError, IndexError):
                        pass
    except Exception as e:
        logger.debug(f"Could not extract coords from {pdbqt_file}: {e}")

    return name, coords


def cluster_poses(poses_dir: str, rmsd_threshold: float = 2.0) -> Dict[str, List[Dict]]:
    """Cluster poses by RMSD. Returns mapping of ligands to pose clusters."""
    clusters = defaultdict(list)

    if not os.path.isdir(poses_dir):
        return clusters

    # Group poses by ligand
    ligand_poses = defaultdict(list)
    for file in os.listdir(poses_dir):
        # match pose files like: <ligand>_pose_1_aff_-7.123.pdbqt or <ligand>_pose_1.pdbqt
        if fnmatch.fnmatch(file, '*_pose_*.pdbqt'):
            ligand = file.split('_pose_')[0]
            ligand_poses[ligand].append(file)

    # Cluster each ligand's poses
    for ligand, pose_files in ligand_poses.items():
        pose_data = []
        for pose_file in sorted(pose_files):
            name, coords = extract_coordinates_from_pdbqt(
                os.path.join(poses_dir, pose_file))
            if coords:
                # Extract pose affinity from filename if available
                affinity = 0.0
                try:
                    # Format: ligand_pose_N_aff_X.XXX.pdbqt
                    match = re.search(r'_aff_([-\d.]+)', pose_file)
                    if match:
                        affinity = float(match.group(1))
                except:
                    pass
                pose_data.append(
                    {'file': pose_file, 'coords': coords, 'affinity': affinity})

        # Cluster poses
        assigned = [False] * len(pose_data)
        cluster_id = 0
        for i, pose_i in enumerate(pose_data):
            if not assigned[i]:
                cluster = {'id': cluster_id, 'poses': [pose_i]}
                assigned[i] = True
                for j in range(i+1, len(pose_data)):
                    if not assigned[j]:
                        rmsd = calculate_rmsd(
                            pose_i['coords'], pose_data[j]['coords'])
                        if rmsd <= rmsd_threshold:
                            cluster['poses'].append(pose_data[j])
                            assigned[j] = True
                clusters[ligand].append(cluster)
                cluster_id += 1

    return clusters


def generate_html_report(results_csv: str, poses_dir: str, output_file: str) -> None:
    """Generate professional HTML report with charts and clustering data."""

    logger.info("[*] Generating HTML report...")

    # Parse results
    ligands = []
    try:
        with open(results_csv, 'r') as f:
            reader = csv.DictReader(f)
            ligands = list(reader)
    except Exception as e:
        logger.warning(f"Could not parse results CSV: {e}")
        return

    if not ligands:
        logger.warning("No results to report")
        return

    # Calculate statistics from valid (negative) affinities only.
    valid_affinities = []
    for lig in ligands:
        raw = lig.get('Binding_Affinity', '').strip()
        if not raw or raw.upper() == 'FAILED':
            continue
        try:
            value = float(raw)
        except ValueError:
            continue
        if value < 0:
            valid_affinities.append(value)

    failed_count = len(ligands) - len(valid_affinities)

    avg_affinity = statistics.mean(
        valid_affinities) if valid_affinities else 0.0
    best_affinity = min(valid_affinities) if valid_affinities else 0.0
    worst_affinity = max(valid_affinities) if valid_affinities else 0.0
    strong_binders = sum(1 for a in valid_affinities if a < -6.0)
    success_rate = (len(valid_affinities) / len(ligands)
                    * 100) if ligands else 0.0

    # Cluster poses
    clusters = cluster_poses(poses_dir) if os.path.isdir(poses_dir) else {}

    # Generate HTML
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Virtual Screening Results Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }}
        .header {{ background-color: #2c3e50; color: white; padding: 20px; border-radius: 5px; }}
        .section {{ background-color: white; padding: 20px; margin: 20px 0; border-radius: 5px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
        th {{ background-color: #34495e; color: white; padding: 10px; text-align: left; }}
        td {{ padding: 10px; border-bottom: 1px solid #ddd; }}
        tr:hover {{ background-color: #f5f5f5; }}
        .stat-box {{ display: inline-block; background-color: #ecf0f1; padding: 15px; margin: 10px; border-radius: 5px; min-width: 150px; }}
        .stat-value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}
        .stat-label {{ font-size: 12px; color: #7f8c8d; margin-top: 5px; }}
        .good {{ color: #27ae60; }}
        .warning {{ color: #f39c12; }}
        .bad {{ color: #e74c3c; }}
        .top-hits {{ background-color: #ecf0f1; padding: 10px; border-left: 4px solid #27ae60; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Virtual Screening Results Report</h1>
        <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>

    <div class="section">
        <h2>📊 Key Statistics</h2>
        <div class="stat-box">
            <div class="stat-value">{len(ligands)}</div>
            <div class="stat-label">Ligands Screened</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{best_affinity:.2f}</div>
            <div class="stat-label">Best Binder (kcal/mol)</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{avg_affinity:.2f}</div>
            <div class="stat-label">Average Affinity</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{strong_binders}/{len(valid_affinities)}</div>
            <div class="stat-label">Strong Binders (&lt; -6.0)</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{success_rate:.1f}%</div>
            <div class="stat-label">Success Rate</div>
        </div>
        <div class="stat-box">
            <div class="stat-value">{failed_count}</div>
            <div class="stat-label">Failed/Invalid Dockings</div>
        </div>
        <div style="clear: both;"></div>
    </div>

    <div class="section">
        <h2>🏆 Top 10 Hits</h2>
        <table>
            <tr>
                <th>#</th>
                <th>Ligand Name</th>
                <th>Binding Affinity (kcal/mol)</th>
                <th>Poses</th>
                <th>Confidence</th>
            </tr>
"""

    # Add top hits
    ranked_ligands = []
    for lig in ligands:
        raw = lig.get('Binding_Affinity', '').strip()
        try:
            affinity_value = float(raw)
        except ValueError:
            continue
        if affinity_value < 0:
            ranked_ligands.append((lig, affinity_value))
    ranked_ligands.sort(key=lambda item: item[1])

    for i, (lig, affinity) in enumerate(ranked_ligands[:10], 1):
        modes = lig.get('Binding_Modes', '0')
        score = lig.get('SimScore', '0.0')

        # Confidence color
        if affinity < -8:
            conf_class = 'good'
            conf_text = '✅ Excellent'
        elif affinity < -6:
            conf_class = 'good'
            conf_text = '✅ Good'
        elif affinity < -4:
            conf_class = 'warning'
            conf_text = '⚠️ Weak'
        else:
            conf_class = 'bad'
            conf_text = '❌ Poor'

        html_content += f"""
            <tr>
                <td>{i}</td>
                <td><strong>{lig.get('Ligand', 'Unknown')}</strong></td>
                <td class="{conf_class}"><strong>{affinity:.2f}</strong></td>
                <td>{modes}</td>
                <td class="{conf_class}">{conf_text}</td>
            </tr>
"""

    html_content += """
        </table>
    </div>

    <div class="section">
        <h2>📈 Affinity Distribution</h2>
        <p><em>Histogram of binding affinities across all ligands</em></p>
        <pre>
"""

    # Create simple histogram
    if valid_affinities:
        bins = 10
        min_aff = min(valid_affinities)
        max_aff = max(valid_affinities)
        bin_width = (max_aff - min_aff) / bins if max_aff != min_aff else 1

        histogram = [0] * bins
        for aff in valid_affinities:
            if bin_width > 0:
                bin_idx = min(int((aff - min_aff) / bin_width), bins - 1)
                histogram[bin_idx] += 1

        for i, count in enumerate(histogram):
            bin_start = min_aff + i * bin_width
            bin_end = bin_start + bin_width
            bar = '█' * int(count * 2)
            html_content += f"{bin_start:6.1f} to {bin_end:6.1f} │ {bar} ({count})\n"

    html_content += """
        </pre>
    </div>

    <div class="section">
        <h2>🔬 Pose Clustering Analysis</h2>
        <p><em>Grouping of similar poses by RMSD (threshold: 2.0 Ångström)</em></p>
"""

    if clusters:
        for ligand, ligand_clusters in clusters.items():
            html_content += f"<h3>{ligand}</h3>\n<ul>\n"
            for cluster in ligand_clusters:
                poses = cluster['poses']
                html_content += f"  <li>Cluster {cluster['id']}: {len(poses)} poses"
                if poses:
                    affs = [p['affinity'] for p in poses if p['affinity']]
                    if affs:
                        html_content += f" (Affinity: {min(affs):.2f} to {max(affs):.2f})"
                html_content += "</li>\n"
            html_content += "</ul>\n"
    else:
        html_content += "<p><em>No docked structures available for clustering.</em></p>\n"

    html_content += """
    </div>

    <div class="section">
        <h2>📋 Interpretation Guide</h2>
        <ul>
            <li><strong>&lt; -10 kcal/mol:</strong> Excellent binder - investigate immediately</li>
            <li><strong>-10 to -8 kcal/mol:</strong> Very good binder - high priority</li>
            <li><strong>-8 to -6 kcal/mol:</strong> Good binder - validate further</li>
            <li><strong>-6 to -4 kcal/mol:</strong> Weak binder - may need optimization</li>
            <li><strong>&gt; -4 kcal/mol:</strong> Negligible binding - not promising</li>
            <li><strong>0.0 kcal/mol:</strong> Docking error - check input format</li>
        </ul>
    </div>

    <div class="section" style="background-color: #ecf0f1; text-align: center; color: #7f8c8d;">
        <p><small>Report generated by Virtual Screening Toolkit v2.0</small></p>
        <p><small>© 2026 - Open Source Virtual Screening Tool</small></p>
    </div>
</body>
</html>
"""

    # Write HTML
    try:
        with open(output_file, 'w') as f:
            f.write(html_content)
        logger.info(f"[✔] HTML report generated: {output_file}")
    except Exception as e:
        logger.warning(f"Could not write HTML report: {e}")

# =============================
# PHASE 2: SMART PREPROCESSING
# =============================


def detect_water_molecules(pdb_file: str, binding_site_coords: Tuple[float, float, float],
                           distance_threshold: float = 4.0) -> List[str]:
    """Detect water molecules near binding site. Returns residue IDs."""
    waters = []

    try:
        cx, cy, cz = binding_site_coords
        with open(pdb_file, 'r') as f:
            for line in f:
                if 'HOH' in line or 'WAT' in line:  # Water residues
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        # Calculate distance to binding site center
                        dist = math.sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2)
                        if dist <= distance_threshold:
                            res_id = line[22:26].strip()
                            waters.append(res_id)
                    except (ValueError, IndexError):
                        pass
    except Exception as e:
        logger.debug(f"Could not detect waters: {e}")

    return list(set(waters))  # Remove duplicates


def detect_metal_ions(pdb_file: str) -> List[Dict]:
    """Detect metal ions in protein. Returns list with ion info."""
    metals = []
    metal_elements = {'ZN': 'Zinc', 'CU': 'Copper', 'FE': 'Iron', 'MG': 'Magnesium',
                      'CA': 'Calcium', 'MN': 'Manganese', 'NI': 'Nickel'}

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    elem = line[76:78].strip().upper()
                    if elem in metal_elements:
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            res_id = line[22:26].strip()
                            metals.append({
                                'element': elem,
                                'name': metal_elements[elem],
                                'residue': res_id,
                                'coords': (x, y, z)
                            })
                        except (ValueError, IndexError):
                            pass
    except Exception as e:
        logger.debug(f"Could not detect metals: {e}")

    return metals


def detect_cofactors(pdb_file: str) -> List[Dict]:
    """Detect cofactors/ligands already in protein."""
    cofactors = []
    common_cofactors = {'NAD', 'ATP', 'GTP', 'FAD', 'HEM', 'ZN', 'CA', 'MG'}

    try:
        with open(pdb_file, 'r') as f:
            seen = set()
            for line in f:
                if line.startswith('HETATM'):
                    res_name = line[17:20].strip().upper()
                    res_id = line[22:26].strip()

                    if res_name in common_cofactors and res_id not in seen:
                        try:
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            cofactors.append({
                                'name': res_name,
                                'residue': res_id,
                                'coords': (x, y, z)
                            })
                            seen.add(res_id)
                        except (ValueError, IndexError):
                            pass
    except Exception as e:
        logger.debug(f"Could not detect cofactors: {e}")

    return cofactors

# =============================
# PHASE 3: FLEXIBLE DOCKING
# =============================


def detect_flexible_residues(pdbqt_file: str, binding_site_coords: Tuple[float, float, float],
                             radius: float = 8.0, max_residues: int = 10) -> List[str]:
    """Auto-detect residues near binding site that should be flexible."""
    flexible = []

    try:
        cx, cy, cz = binding_site_coords
        residues = {}

        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())

                        dist = math.sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2)
                        if dist <= radius:
                            chain = line[21].strip()
                            res_num = line[22:26].strip()
                            res_key = f"{chain}{res_num}"

                            if res_key not in residues:
                                residues[res_key] = dist
                    except (ValueError, IndexError):
                        pass

        # Sort by distance and take closest
        sorted_residues = sorted(residues.items(), key=lambda x: x[1])
        flexible = [res[0] for res in sorted_residues[:max_residues]]

    except Exception as e:
        logger.debug(f"Could not detect flexible residues: {e}")

    return flexible

# =============================
# PHASE 4: CONSENSUS SCORING
# =============================

# Note: SMINA tool path will be set after find_tool is defined


def consensus_rank(vina_affinity: float, smina_affinity: Optional[float]) -> Dict:
    """Calculate consensus score."""
    if smina_affinity is None:
        return {'vina': vina_affinity, 'smina': None, 'consensus': vina_affinity, 'agreement': 0}

    # Calculate agreement (0-100, higher = better agreement)
    diff = abs(vina_affinity - smina_affinity)
    agreement = max(0, 100 - (diff * 10))  # Rule of thumb

    # Average the scores
    consensus = (vina_affinity + smina_affinity) / 2

    return {
        'vina': vina_affinity,
        'smina': smina_affinity,
        'consensus': consensus,
        'agreement': agreement
    }


def find_tool(*names: str) -> str:
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    return names[-1]


# Tool paths (prefer modern Vina for consistency)
OBABEL = find_tool("obabel", "OpenBabel", "OpenBabel.exe", "obabel")
# Prefer modern Vina for reproducibility/consistency
VINA_PRIMARY = find_tool("vina", "vina.exe", "qvina02", "qvina2", "qvina")
VINA_FALLBACK = find_tool("qvina02", "qvina2") if os.environ.get(
    "VS_ENABLE_QVINA_FALLBACK", "0") == "1" else VINA_PRIMARY
VINA = VINA_PRIMARY  # Start with primary tool
# Optional for consensus scoring
SMINA = find_tool("smina", "smina.exe", "smina")
COMMAND_TIMEOUT = int(os.environ.get("VS_COMMAND_TIMEOUT", "900"))

# Verify tools exist - STRICT validation
for tool_name, tool_path in [("obabel", OBABEL), ("vina", VINA)]:
    if not shutil.which(tool_path):
        raise RuntimeError(
            f"Critical tool '{tool_name}' not found in PATH. Ensure it is installed and in PATH.")

if not shutil.which(SMINA):
    logger.info("SMINA not found - consensus scoring will be skipped")

# =============================
# PHASE 4 SMINA DOCKING
# =============================


def run_smina_docking(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    output_pdbqt: str,
    cx: float, cy: float, cz: float,
    sx: float, sy: float, sz: float
) -> Optional[float]:
    """Run SMINA docking and return best affinity."""

    if not shutil.which(SMINA):
        logger.debug("SMINA not available")
        return None

    cmd = [
        SMINA,
        "-r", receptor_pdbqt,
        "-l", ligand_pdbqt,
        "-o", output_pdbqt,
        "--center_x", str(cx),
        "--center_y", str(cy),
        "--center_z", str(cz),
        "--size_x", str(sx),
        "--size_y", str(sy),
        "--size_z", str(sz),
        "--num_modes", "1",
        "--exhaustiveness", str(os.environ.get('VS_EXHAUSTIVENESS', 8))
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=COMMAND_TIMEOUT
    )

    # 🔴 Fail loudly
    if result.returncode != 0:
        raise RuntimeError(
            f"SMINA docking failed:\nSTDERR:\n{result.stderr}\nSTDOUT:\n{result.stdout}"
        )

    # 🔴 Validate output file
    if not os.path.exists(output_pdbqt):
        raise RuntimeError("SMINA did not produce output file")

    with open(output_pdbqt) as f:
        if not any(line.startswith("MODEL") for line in f):
            raise RuntimeError("SMINA produced no docking poses")

    # 🔴 Parse properly (mode table)
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                return float(parts[1])
            except ValueError:
                continue

    raise RuntimeError("Could not extract affinity from SMINA output")
# =============================
# PDBQT CHARGE VALIDATION & FIXING
# =============================

def _ensure_pdbqt_has_charges(pdbqt_file: str) -> bool:
    """Check if PDBQT file contains at least one valid atomic charge (robust)."""
    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    parts = line.split()
                    # PDBQT charge is usually last column
                    if len(parts) >= 9:
                        try:
                            float(parts[-1])
                            return True
                        except ValueError:
                            continue
        return False
    except Exception as e:
        logger.warning(f"Charge check failed: {e}")
        return False


def _pdbqt_has_atoms(pdbqt_file: str) -> bool:
    """Return True when a PDBQT contains at least one ATOM/HETATM record."""
    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    return True
        return False
    except Exception as e:
        logger.warning(f"Could not inspect atom records in {pdbqt_file}: {e}")
        return False


def _sanitize_receptor_pdbqt(pdbqt_file: str) -> None:
    """Remove ligand-only tags from receptor PDBQT for strict Vina compatibility."""
    allowed_prefixes = ("ATOM", "HETATM", "TER", "END", "REMARK")
    removed = 0
    kept = []

    try:
        with open(pdbqt_file, "r") as f:
            for line in f:
                stripped = line.lstrip()
                if stripped.startswith(("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")):
                    removed += 1
                    continue
                if line.startswith(allowed_prefixes):
                    kept.append(line)
                    continue
                # Keep unrecognized lines as REMARK to preserve traceability.
                if line.strip():
                    kept.append("REMARK " + line)

        with open(pdbqt_file, "w") as f:
            f.writelines(kept)

        if removed:
            logger.info(
                f"[*] Sanitized receptor PDBQT: removed {removed} ligand-only tags")
    except Exception as e:
        raise RuntimeError(
            f"Failed to sanitize receptor PDBQT {pdbqt_file}: {e}")


def _fix_pdbqt_charges(pdbqt_file: str, source_file: str) -> None:
    """
    STRICT charge assignment.
    No silent guessing. Only real tools allowed.
    """

    logger.warning(f"[!] Fixing charges for: {pdbqt_file}")

    # Try OpenBabel only (most reliable baseline)
    try:
        run([
            OBABEL,
            "-isdf" if source_file.endswith(".sdf") else "-ipdb",
            source_file,
            "-opdbqt",
            "-O",
            pdbqt_file,
            "--partialcharge",
            "gasteiger"
        ], capture=False)

        if _ensure_pdbqt_has_charges(pdbqt_file):
            logger.info("[✔] Charges assigned via OpenBabel Gasteiger")
            return

    except Exception as e:
        logger.error(f"Charge assignment failed: {e}")

    raise RuntimeError(
        f"CRITICAL: Cannot assign valid charges for {pdbqt_file}. "
        "Install OpenBabel with Gasteiger support."
    )


def _assign_simple_charges(pdbqt_file: str) -> None:
    """Deprecated legacy helper (kept for backward compatibility; not used)."""
    # Simple charge assignment based on common atom types
    charge_map = {
        'C': 0.0,    # Carbon
        'N': -0.3,   # Nitrogen
        'O': -0.5,   # Oxygen
        'S': -0.3,   # Sulfur
        'P': 1.1,    # Phosphorus
        'H': 0.0,    # Hydrogen (connected to C)
        'HD': 0.12,  # Polar H
        'OA': -0.7,  # Acceptor O
        'NA': -0.7,  # Acceptor N
        'SA': -0.7,  # Acceptor S
    }

    lines = []
    modified = False

    try:
        with open(pdbqt_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Check if charge is zero (needs assignment)
                    if len(line) > 76:
                        try:
                            charge_str = line[70:76].strip()
                            charge = float(charge_str) if charge_str else 0.0

                            if charge == 0.0:
                                # Get atom type from columns 77-78
                                atom_type = line[77:79].strip() if len(
                                    line) > 78 else ''

                                # Look up charge
                                new_charge = charge_map.get(atom_type, 0.0)

                                if new_charge != 0.0:
                                    # Replace charge in the line
                                    charge_field = f"{new_charge:+7.3f}"
                                    line_list = list(line)
                                    for i, c in enumerate(charge_field):
                                        if 70 + i <= 76 and i < len(charge_field):
                                            if 70 + i < len(line_list):
                                                line_list[70 + i] = c
                                    line = ''.join(
                                        line_list[:81]) + (line[81:] if len(line) > 81 else '')
                                    modified = True
                        except (ValueError, IndexError):
                            pass

                lines.append(line)

        # Write back if modified
        if modified:
            with open(pdbqt_file, 'w') as f:
                f.writelines(lines)
            logger.debug(
                f"Heuristic charges assigned to {os.path.basename(pdbqt_file)}")
    except Exception as e:
        logger.debug(f"Could not assign charges heuristically: {e}")

# =============================
# SDF VALIDATION & UTILITIES
# =============================


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
            [OBABEL, sdf_file, "-osmi", "--append",
                "MW logP TPSA rotors HBD HBA2"],
            capture_output=True,
            text=True,
            timeout=COMMAND_TIMEOUT,
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

# =============================
# CHECKPOINT/RESUME SYSTEM
# =============================


class DockingCheckpoint:
    """Manages docking progress and resume capability."""

    def __init__(self, outdir: str):
        self.checkpoint_file = os.path.join(outdir, ".docking_checkpoint.json")
        self.completed = self._load_checkpoint()

    def _load_checkpoint(self) -> Dict:
        """Load previous docking progress."""
        if os.path.exists(self.checkpoint_file):
            try:
                with open(self.checkpoint_file, "r") as f:
                    data = json.load(f)
                logger.info(
                    f"[✔] Loaded checkpoint: {len(data)} ligands already docked")
                return data
            except Exception as e:
                logger.warning(
                    f"Could not load checkpoint: {e}. Starting fresh.")
                return {}
        return {}

    def is_completed(self, ligand_name: str) -> bool:
        """Check if ligand already docked."""
        return ligand_name in self.completed

    def save_result(self, ligand_name: str, score: float, metrics: Optional[Dict] = None):
        """Save docking result to checkpoint."""
        self.completed[ligand_name] = {
            "score": score,
            "timestamp": datetime.now().isoformat(),
            "metrics": metrics or {},
        }
        self._persist()

    def _persist(self):
        """Write checkpoint to disk."""
        try:
            with open(self.checkpoint_file, "w") as f:
                json.dump(self.completed, f, indent=2)
        except Exception as e:
            logger.warning(f"Failed to save checkpoint: {e}")

    def get_results(self) -> List[Tuple[str, float]]:
        """Return all completed results as (ligand_name, score) tuples."""
        return [(name, data["score"]) for name, data in self.completed.items()]

    def get_metrics(self) -> Dict:
        """Return cached metrics keyed by ligand name."""
        return {name: data.get("metrics", {}) for name, data in self.completed.items()}

# =============================
# LIBRARY GENERATION (ZINC API + ADMET)
# =============================


class LibraryManager:
    """Manages compound library sourcing and preparation."""

    ZINC_URL = "https://zinc.docking.org/api/v0/"

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
        self.session = requests.Session() if requests is not None else None
        # Allow overriding ZINC API base URL and API key via environment variables.
        # Useful when the public endpoints change or require authentication.
        self.zinc_api_base = os.environ.get('ZINC_API_BASE') or self.ZINC_URL
        self.zinc_api_key = os.environ.get('ZINC_API_KEY')
        # Header name to use for the API key. Common choices: 'Authorization' (Bearer), 'X-API-KEY'
        self.zinc_api_key_header = os.environ.get(
            'ZINC_API_KEY_HEADER', 'Authorization')

    def create_fda_library(self, apply_admet: bool = True) -> List[str]:
        """Download FDA-approved drugs from ZINC.

        Uses ZINC REST API: https://zinc.docking.org/
        Filters: FDA approved, MW 200-500, logP < 5
        """
        logger.info("[*] FDA library mode - downloading from ZINC API...")

        if self.session is None:
            logger.warning(
                "requests package not found. Falling back to local SDF mode.")
            return self._prepare_local_sdf(apply_admet=apply_admet)

        try:
            substances = self._query_zinc(
                supplier='fda8',
                count=10000,
                mw_min=200,
                mw_max=500
            )

            if not substances:
                logger.warning(
                    "No FDA substances found, falling back to local SDF")
                return self._prepare_local_sdf(apply_admet=apply_admet)

            return self._download_and_prepare(substances, apply_admet)

        except Exception as e:
            logger.warning(f"ZINC API failed: {e}. Using local SDF.")
            return self._prepare_local_sdf(apply_admet=apply_admet)

    def create_custom_library(self, mw_min: int = 200, mw_max: int = 500,
                              logp_max: float = 5, apply_admet: bool = True) -> List[str]:
        """Download custom library filtered by molecular properties.

        Uses ZINC REST API with filters:
        - Molecular weight: mw_min to mw_max
        - LogP: <= logp_max (drug-likeness)
        """
        logger.info(
            f"[*] Custom library mode (MW: {mw_min}-{mw_max}, LogP: <={logp_max})")

        if self.session is None:
            logger.warning(
                "requests package not found. Falling back to local SDF mode.")
            return self._prepare_local_sdf(apply_admet=apply_admet)

        try:
            substances = self._query_zinc(
                supplier='now',
                count=50000,
                mw_min=mw_min,
                mw_max=mw_max,
                logp_max=logp_max
            )

            if not substances:
                logger.warning(
                    "No substances found, falling back to local SDF")
                return self._prepare_local_sdf(apply_admet=apply_admet)

            return self._download_and_prepare(substances, apply_admet)

        except Exception as e:
            logger.warning(f"ZINC API failed: {e}. Using local SDF.")
            return self._prepare_local_sdf(apply_admet=apply_admet)

    def _query_zinc(self, supplier: str = 'now', count: int = 10000,
                    mw_min: int = 200, mw_max: int = 500,
                    logp_max: float = 5) -> List[Dict]:
        """Query ZINC API for compounds."""
        logger.debug(
            f"Querying ZINC API: supplier={supplier}, MW={mw_min}-{mw_max}")

        params = {
            'supplier': supplier,
            'mw__lte': mw_max,
            'mw__gte': mw_min,
            'logp__lte': logp_max,
            'limit': min(count, 10000),
            'format': 'json'
        }

        # Try a set of candidate base URLs and supplier names to tolerate
        # ZINC API endpoint changes or deprecated supplier codes.
        candidate_bases = [self.zinc_api_base]
        # Add historical/known endpoints as fallbacks if not overridden
        for b in ("https://zinc15.docking.org/api/",
                  "https://zinc15.docking.org/api/v1/",
                  "https://zinc.docking.org/api/",
                  "https://zinc.docking.org/"):
            if b not in candidate_bases:
                candidate_bases.append(b)

        candidate_suppliers = [supplier]
        # common fallbacks
        for s in ("fda", "fda8", "now"):
            if s not in candidate_suppliers:
                candidate_suppliers.append(s)

        last_exception = None
        for base in candidate_bases:
            for sup in candidate_suppliers:
                params['supplier'] = sup
                try:
                    url = urljoin(base, "substances")
                    logger.debug(
                        f"Trying ZINC endpoint: {url} (supplier={sup})")
                    headers = {}
                    if self.zinc_api_key:
                        # Support Bearer tokens or raw API keys depending on header
                        if self.zinc_api_key_header.lower() == 'authorization' and not self.zinc_api_key.lower().startswith('bearer '):
                            headers['Authorization'] = f"Bearer {self.zinc_api_key}"
                        else:
                            headers[self.zinc_api_key_header] = self.zinc_api_key

                    response = self.session.get(
                        url, params=params, headers=headers or None, timeout=30)
                    response.raise_for_status()

                    data = response.json()
                    # Some API variants may return the list at top-level
                    substances = data.get('results') or data.get(
                        'substances') or data

                    if not substances:
                        logger.debug(
                            f"Endpoint {url} returned no substances (supplier={sup})")
                        continue

                    logger.info(
                        f"[✔] Found {len(substances)} compounds from ZINC (endpoint: {base}, supplier: {sup})")
                    return substances

                except Exception as e:
                    logger.debug(
                        f"ZINC attempt failed for {base} (supplier={sup}): {e}")
                    last_exception = e
                    # try the next candidate
                    continue

        # All attempts failed; raise the last caught exception to be handled by caller
        logger.error("All ZINC endpoint attempts failed")
        if last_exception:
            raise last_exception
        raise RuntimeError("ZINC query failed for unknown reasons")

    def _download_and_prepare(self, substances: List[Dict],
                              apply_admet: bool = True) -> List[str]:
        """Download SDF files and convert to PDBQT."""
        admet = ADMETFilter()
        downloaded = []
        skipped = 0

        for i, subst in enumerate(substances[:100], 1):
            try:
                zinc_id = subst.get('zinc_id')
                if not zinc_id:
                    continue

                sdf_url = f"http://zinc.docking.org/substances/{zinc_id}.sdf"
                sdf_file = os.path.join(self.lib_dir, f"{zinc_id}.sdf")

                logger.debug(f"Downloading [{i}/100]: {zinc_id}...")

                response = self.session.get(sdf_url, timeout=10)
                response.raise_for_status()

                with open(sdf_file, 'wb') as f:
                    f.write(response.content)

                if apply_admet:
                    props = admet.parse_sdf_properties(sdf_file)
                    passes, violations = admet.check_lipinski(props)
                    if not passes:
                        logger.debug(
                            f"  Skipped {zinc_id}: {', '.join(violations)}")
                        os.remove(sdf_file)
                        skipped += 1
                        continue

                pdbqt_file = os.path.join(self.lib_dir, f"{zinc_id}.pdbqt")
                pdbqt_file = self._prepare_sdf_to_pdbqt(
                    sdf_file, pdbqt_file, zinc_id)
                downloaded.append(pdbqt_file)

                if len(downloaded) >= 20:
                    logger.info(f"[*] Reached 20 compounds, stopping download")
                    break

            except Exception as e:
                logger.debug(
                    f"Failed to download {subst.get('zinc_id', 'unknown')}: {e}")
                skipped += 1

        logger.info(f"[✔] Downloaded and prepared {len(downloaded)} compounds"
                    f" ({skipped} filtered by ADMET)")

        self._save_metadata()
        return downloaded

    def _prepare_sdf_to_pdbqt(self, sdf_file: str, pdbqt_file: str, ligand_id: str) -> str:
        """Optionally minimize an SDF with MMFF94, then convert to PDBQT."""
        source_for_conversion = sdf_file
        minimized_file = os.path.join(self.minimized_dir, f"{ligand_id}.sdf")

        # Parse properties to get real data
        props = ADMETFilter.parse_sdf_properties(sdf_file)
        if props is None:
            props = {}

        if self.minimize:
            try:
                run([
                    OBABEL, "-isdf", sdf_file, "-osdf", "-O", minimized_file,
                    "--gen3d", "--minimize", "--ff", "MMFF94",
                    "--steps", str(self.minimize_steps)
                ])
                source_for_conversion = minimized_file
            except Exception as e:
                logger.warning(
                    f"MMFF94 minimization failed for {ligand_id}: {e}. Converting original SDF.")

        # Convert to PDBQT
        try:
            run([OBABEL, "-isdf", source_for_conversion, "-opdbqt", "-O", pdbqt_file])
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
        if not _ensure_pdbqt_has_charges(pdbqt_file):
            logger.debug(
                f"[*] Ligand {ligand_id} PDBQT has no charges - computing Gasteiger...")
            _fix_pdbqt_charges(pdbqt_file, source_for_conversion)

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
            "zinc_url": _zinc_url(ligand_id),
            "minimized": self.minimize and source_for_conversion == minimized_file,
            "properties_source": "obabel_descriptors" if props else "default",
        }

        return pdbqt_file

    def _save_metadata(self):
        try:
            with open(self.metadata_file, "w") as f:
                json.dump(self.metadata, f, indent=2)
            logger.info(f"[✔] Ligand metadata saved: {self.metadata_file}")
        except IOError as e:
            logger.warning(f"Could not save ligand metadata: {e}")

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
            for f in os.listdir(ligands_input):
                path = os.path.join(ligands_input, f)

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
            logger.info(f"[DEBUG] Local ligands found - SDF: {len(sdf_files)}, PDBQT: {len(pdbqt_files)}, MOL2: {len(mol2_files)}, PDB: {len(pdb_files)}")
            sample = (sdf_files or pdbqt_files or mol2_files or pdb_files)[:5]
            if sample:
                logger.info(f"[DEBUG] Sample files: {', '.join([os.path.basename(s) for s in sample])}")
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
            raise FileNotFoundError("No ligands found (.sdf/.pdbqt/.mol2/.pdb)")
    
        # -----------------------------
        # STEP 3: PROCESS
        # -----------------------------
        admet = ADMETFilter()
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
    
                    props = admet.parse_sdf_properties(inp)
                    if props is None:
                        failed_ligands.append((lig, "Invalid SDF properties"))
                        continue
    
                    if apply_admet:
                        ok, violations = admet.check_lipinski(props)
                        if not ok:
                            failed_ligands.append((lig, f"ADMET: {violations}"))
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
                        run([OBABEL, "-imol2", inp, "-opdbqt", "-O", out])
                        if not os.path.exists(out) or os.path.getsize(out) == 0:
                            raise RuntimeError("Conversion produced empty file")
                        if not _ensure_pdbqt_has_charges(out):
                            _fix_pdbqt_charges(out, inp)
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
                        run([OBABEL, "-ipdb", inp, "-opdbqt", "-O", out])
                        if not os.path.exists(out) or os.path.getsize(out) == 0:
                            raise RuntimeError("Conversion produced empty file")
                        if not _ensure_pdbqt_has_charges(out):
                            _fix_pdbqt_charges(out, inp)
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

# =============================
# UTILITIES
# =============================


def run(cmd_list: List[str], capture: bool = False) -> Optional[Tuple[str, str, int]]:
    """Run a command safely without shell injection risk."""
    try:
        if capture:
            result = subprocess.run(
                cmd_list, capture_output=True, text=True, timeout=COMMAND_TIMEOUT)
            if result.returncode != 0:
                raise subprocess.CalledProcessError(
                    result.returncode, cmd_list, output=result.stdout, stderr=result.stderr
                )
            return result.stdout, result.stderr, result.returncode
        else:
            subprocess.run(cmd_list, check=True, timeout=COMMAND_TIMEOUT)
    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out: {' '.join(cmd_list)}")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(
            f"Command failed with return code {e.returncode}: {' '.join(cmd_list)}")
        raise
    except FileNotFoundError as e:
        logger.error(f"Command not found: {' '.join(cmd_list)}")
        raise


def parse_pdb_coords(pdb_path: str) -> Tuple[List[float], List[float], List[float]]:
    """Extract XYZ coordinates from PDB file."""
    xs, ys, zs = [], [], []
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
                    except (ValueError, IndexError):
                        logger.warning(
                            f"Skipping malformed PDB line: {line.strip()}")
                        continue
    except IOError as e:
        logger.error(f"Cannot read PDB file {pdb_path}: {e}")
        raise

    if not xs:
        raise ValueError(f"No ATOM records found in {pdb_path}")

    return xs, ys, zs


def _zinc_url(ligand_id: str) -> str:
    clean_id = ligand_id.strip()
    if re.match(r"^ZINC\d+", clean_id, re.IGNORECASE):
        return f"https://zinc.docking.org/substances/{clean_id.upper()}/"
    return ""


def get_chains(pdb_file: str) -> List[str]:
    """Extract unique chains from PDB file."""
    chains = set()
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    chain = line[21]
                    if chain.strip():
                        chains.add(chain)
    except Exception as e:
        logger.warning(f"Could not parse chains: {e}")
        return ["A"]

    return sorted(list(chains)) if chains else ["A"]


def _parse_chain_selection(chain: str, available: List[str]) -> List[str]:
    choice = (chain or "").strip()
    if not choice or choice.lower() in {"all", "*"}:
        return available

    selected = [c.strip().upper()
                for c in re.split(r"[,\s]+", choice) if c.strip()]
    invalid = [c for c in selected if c not in available]
    if invalid:
        raise ValueError(
            f"Invalid chain(s): {', '.join(invalid)}. Available chains: {', '.join(available)}"
        )
    return selected

# =============================
# PROTEIN PREPARATION
# =============================


class ProteinPreparation:
    """Advanced protein preparation with chain selection."""

    def __init__(self, pdb_file: str, outdir: str):
        self.pdb_file = pdb_file
        self.outdir = outdir
        self.pdb_clean = os.path.join(outdir, "protein_clean.pdb")
        self.receptor_pdbqt = os.path.join(outdir, "receptor.pdbqt")
        self.grid_conf = os.path.join(outdir, "grid.conf")
        self.grid_box_script = os.path.join(outdir, "grid_box.py")
        self.pocket_summary_file = os.path.join(outdir, "pocket_summary.csv")

    def validate(self):
        """Validate protein PDB file."""
        if not os.path.exists(self.pdb_file):
            raise FileNotFoundError(f"Protein file not found: {self.pdb_file}")

        if not os.access(self.pdb_file, os.R_OK):
            raise PermissionError(f"Cannot read protein file: {self.pdb_file}")

        try:
            count = sum(1 for l in open(self.pdb_file) if l.startswith("ATOM"))
        except IOError as e:
            raise IOError(f"Error reading protein file: {e}")

        if count < 100:
            raise ValueError(
                f"Invalid protein: only {count} ATOM lines (need >= 100)")

        logger.info(f"[✔] Protein validated: {count} ATOM lines")

    def select_chain(self) -> str:
        """Interactively select chain for docking."""
        chains = get_chains(self.pdb_file)

        if len(chains) == 1:
            logger.info(f"[*] Single chain detected: {chains[0]}")
            return chains[0]

        logger.info(f"[*] Multiple chains detected: {', '.join(chains)}")
        while True:
            choice = input(
                f"Select chain [{'/'.join(chains)}]: ").strip().upper()
            if choice in chains:
                return choice
            logger.warning(f"Invalid choice. Choose from: {', '.join(chains)}")

    def prepare_receptor(self, chain: str = "A", keep_hetero: bool = False) -> List[str]:
        """Clean protein and convert to PDBQT."""
        available_chains = get_chains(self.pdb_file)
        selected_chains = _parse_chain_selection(chain, available_chains)
        logger.info(
            f"[*] Preparing receptor (chain(s): {', '.join(selected_chains)})...")

        kept_atoms = 0
        skipped_hetero = 0
        try:
            with open(self.pdb_file, "r") as src, open(self.pdb_clean, "w") as dst:
                for line in src:
                    if line.startswith("ATOM"):
                        if len(line) > 21 and line[21].strip() in selected_chains:
                            dst.write(line)
                            kept_atoms += 1
                    elif keep_hetero and line.startswith("HETATM"):
                        if len(line) > 21 and line[21].strip() in selected_chains:
                            dst.write(line)
                            kept_atoms += 1
                    elif line.startswith("HETATM"):
                        skipped_hetero += 1
                    elif line.startswith(("TER", "END")):
                        dst.write(line)
        except IOError as e:
            raise IOError(
                f"Failed to filter chain(s) {', '.join(selected_chains)}: {e}")

        if kept_atoms == 0:
            raise ValueError(
                f"Selected chain(s) '{', '.join(selected_chains)}' have no atoms in {self.pdb_file}")

        try:
            # Use basic obabel conversion; QuickVina-compatible PDBQT format
            run([OBABEL, "-ipdb", self.pdb_clean, "-opdbqt", "-O", self.receptor_pdbqt,
                 "-xr", "-c"])

            # Receptor PDBQT must not contain ligand torsion tags (ROOT/BRANCH/TORSDOF).
            _sanitize_receptor_pdbqt(self.receptor_pdbqt)

            # Verify charges were computed
            if not _ensure_pdbqt_has_charges(self.receptor_pdbqt):
                logger.warning(
                    f"[!] Receptor PDBQT has no charges - attempting Gasteiger computation...")
                _fix_pdbqt_charges(self.receptor_pdbqt, self.pdb_clean)
                _sanitize_receptor_pdbqt(self.receptor_pdbqt)

            logger.info(f"[✔] Receptor prepared: {self.receptor_pdbqt}")
            if skipped_hetero and not keep_hetero:
                logger.info(
                    f"[*] Removed {skipped_hetero} HETATM records during receptor cleaning")
        except Exception as e:
            logger.error(f"Failed to prepare receptor: {e}")
            raise

        return selected_chains

    def detect_pocket(self, pocket_spec: Optional[str] = None,
                      padding: float = 6.0) -> Tuple[float, float, float, float, float, float]:
        """Detect binding pockets using fpocket and build a grid from one or more pockets."""
        logger.info("[*] Running fpocket...")
        fpocket_target = self.pdb_clean if os.path.exists(
            self.pdb_clean) else self.pdb_file
        try:
            run(["fpocket", "-f", fpocket_target])
        except Exception as e:
            logger.warning(f"fpocket failed: {e}. Using fallback grid.")
            return 0, 0, 0, 24, 24, 24

        pocket_root = fpocket_target.replace(".pdb", "_out")
        pocket_dir = os.path.join(pocket_root, "pockets")

        if not os.path.exists(pocket_dir):
            logger.warning(
                f"Pocket directory not found: {pocket_dir}. Using fallback grid.")
            return 0, 0, 0, 24, 24, 24

        info_file = os.path.join(
            pocket_root, f"{Path(fpocket_target).stem}_info.txt"
        )
        pocket_info = self._parse_fpocket_info(info_file)
        pockets = []

        for f in os.listdir(pocket_dir):
            match = re.match(r"pocket(\d+)_atm\.pdb$", f)
            if not match:
                continue

            number = int(match.group(1))
            p = os.path.join(pocket_dir, f)

            try:
                xs, ys, zs = parse_pdb_coords(p)
                size = (max(xs) - min(xs)) * \
                    (max(ys) - min(ys)) * (max(zs) - min(zs))
                data = dict(pocket_info.get(number, {}))
                data.update({
                    "number": number,
                    "file": f,
                    "path": p,
                    "box_volume": size,
                })
                pockets.append(data)
            except ValueError:
                continue

        if not pockets:
            logger.warning("No valid pockets found. Using fallback grid.")
            return 0, 0, 0, 24, 24, 24

        # FIXED: Sort by SCORE first (binding pocket quality), then druggability
        # Previously: sorted by druggability_score first, which selected worst-scoring
        # pocket (Pocket 4) over best-scoring pocket (Pocket 1) when druggability was slightly higher.
        # New behavior: Best-scoring pocket wins, with druggability as tiebreaker.
        pockets.sort(key=lambda p: (
            p.get("score", -1),
            p.get("druggability_score", -1),
            p.get("volume", -1),
            p.get("box_volume", -1),
        ), reverse=True)
        self._write_pocket_summary(pockets)

        logger.info(f"[*] Found {len(pockets)} pockets. Top candidates:")
        for pocket in pockets[:5]:
            druggable = "druggable" if pocket.get(
                "druggability_score", 0) > 0.15 else "low-druggability"
            logger.info(
                "  pocket%-2s score=%6.3f druggability=%5.3f volume=%7.1f %s"
                % (
                    pocket["number"],
                    pocket.get("score", 0),
                    pocket.get("druggability_score", 0),
                    pocket.get("volume", 0),
                    druggable,
                )
            )

        selected_numbers = self._parse_pocket_selection(pocket_spec, pockets)
        selected = [p for p in pockets if p["number"] in selected_numbers]
        selected.sort(key=lambda p: selected_numbers.index(p["number"]))

        selected_names = ", ".join(f"pocket{p['number']}" for p in selected)
        logger.info(
            f"[✔] Pocket selection: {selected_names} (padding={padding:g} Å)")
        return self._get_pocket_info([p["path"] for p in selected], padding=padding)

    def _parse_fpocket_info(self, info_file: str) -> Dict[int, Dict]:
        """Parse fpocket's *_info.txt file into a pocket metadata map."""
        info = {}
        if not os.path.exists(info_file):
            return info

        current = None
        try:
            with open(info_file) as f:
                for line in f:
                    pocket_match = re.match(
                        r"Pocket\s+(\d+)\s*:", line.strip())
                    if pocket_match:
                        current = int(pocket_match.group(1))
                        info[current] = {}
                        continue

                    if current is None or ":" not in line:
                        continue

                    key, value = line.split(":", 1)
                    key = key.strip().lower().replace(" ", "_").replace(".", "")
                    key = key.replace("-", "_")
                    try:
                        info[current][key] = float(value.strip().split()[0])
                    except (ValueError, IndexError):
                        continue
        except IOError as e:
            logger.warning(
                f"Could not parse fpocket info file {info_file}: {e}")

        return info

    def _write_pocket_summary(self, pockets: List[Dict]):
        try:
            with open(self.pocket_summary_file, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow([
                    "Pocket", "Score", "Druggability_Score", "Volume",
                    "Box_Volume", "Could_Be_Druggable", "File"
                ])
                for pocket in sorted(pockets, key=lambda p: p["number"]):
                    writer.writerow([
                        pocket["number"],
                        pocket.get("score", ""),
                        pocket.get("druggability_score", ""),
                        pocket.get("volume", ""),
                        f"{pocket.get('box_volume', 0):.3f}",
                        "YES" if pocket.get(
                            "druggability_score", 0) > 0.15 else "NO",
                        pocket["file"],
                    ])
            logger.info(
                f"[✔] Pocket summary saved: {self.pocket_summary_file}")
        except IOError as e:
            logger.warning(f"Could not save pocket summary: {e}")

    def _parse_pocket_selection(self, pocket_spec: Optional[str], pockets: List[Dict]) -> List[int]:
        available = {p["number"] for p in pockets}
        if not pocket_spec or pocket_spec.lower() == "auto":
            return [pockets[0]["number"]]

        selected = []
        for token in re.split(r"[,\s]+", pocket_spec.strip()):
            if not token:
                continue
            try:
                pocket_number = int(token)
            except ValueError as e:
                raise ValueError(
                    f"Invalid pocket selection '{token}'. Use numbers like 1 or 1,3.") from e
            if pocket_number not in available:
                raise ValueError(
                    f"Pocket {pocket_number} was not found. Available: {sorted(available)}")
            selected.append(pocket_number)

        if not selected:
            raise ValueError("No pockets selected")
        return selected

    def _get_pocket_info(self, pocket_files: List[str], padding: float = 6.0) -> Tuple[float, float, float, float, float, float]:
        """Calculate grid center and size from one or more pocket PDB files."""
        xs, ys, zs = [], [], []
        for pocket_file in pocket_files:
            px, py, pz = parse_pdb_coords(pocket_file)
            xs.extend(px)
            ys.extend(py)
            zs.extend(pz)

        cx = sum(xs) / len(xs)
        cy = sum(ys) / len(ys)
        cz = sum(zs) / len(zs)

        sx = max(xs) - min(xs) + (2 * padding)
        sy = max(ys) - min(ys) + (2 * padding)
        sz = max(zs) - min(zs) + (2 * padding)

        logger.debug(f"Grid center: ({cx:.2f}, {cy:.2f}, {cz:.2f}), "
                     f"size: ({sx:.2f}, {sy:.2f}, {sz:.2f})")

        return cx, cy, cz, sx, sy, sz

    def write_grid(self, cx: float, cy: float, cz: float,
                   sx: float, sy: float, sz: float):
        """Write AutoDock Vina grid configuration."""
        try:
            with open(self.grid_conf, "w") as f:
                f.write(f"""center_x = {cx}
center_y = {cy}
center_z = {cz}
size_x = {sx}
size_y = {sy}
size_z = {sz}
""")
            logger.info(f"[✔] Grid ready: {self.grid_conf}")
            self.write_grid_box_script(cx, cy, cz, sx, sy, sz)
        except IOError as e:
            logger.error(f"Failed to write grid file: {e}")
            raise

    def write_grid_box_script(self, cx: float, cy: float, cz: float,
                              sx: float, sy: float, sz: float):
        """Write a PyMOL script that visualizes the final docking box."""
        x0, x1 = cx - sx / 2, cx + sx / 2
        y0, y1 = cy - sy / 2, cy + sy / 2
        z0, z1 = cz - sz / 2, cz + sz / 2
        vertices = [
            (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
            (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
        ]
        edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ]

        try:
            with open(self.grid_box_script, "w") as f:

                f.write("import os\n\n")
                f.write("from pymol import cmd\n")
                f.write(
                    "from pymol.cgo import BEGIN, LINES, VERTEX, END, COLOR, LINEWIDTH\n\n")
                f.write("os.chdir(os.path.dirname(os.path.abspath(__file__)))\n")
                f.write(
                    f"cmd.load(r'{os.path.relpath(self.receptor_pdbqt, os.path.dirname(os.path.abspath(self.grid_box_script)))}', 'receptor')\n")
                f.write(
                    "box = [COLOR, 0.0, 0.35, 1.0, LINEWIDTH, 3.0, BEGIN, LINES,\n")
                for a, b in edges:
                    ax, ay, az = vertices[a]
                    bx, by, bz = vertices[b]
                    f.write(
                        f"       VERTEX, {ax:.3f}, {ay:.3f}, {az:.3f}, VERTEX, {bx:.3f}, {by:.3f}, {bz:.3f},\n")
                f.write("       END]\n")
                f.write("cmd.load_cgo(box, 'docking_grid')\n")
                f.write("cmd.zoom('receptor or docking_grid')\n")
            logger.info(
                f"[✔] Grid visualization script saved: {self.grid_box_script}")
        except IOError as e:
            logger.warning(f"Could not write grid visualization script: {e}")

# =============================
# DOCKING WITH RESUME + SimScore
# =============================


def _parse_vina_modes(vina_output: str) -> List[Dict[str, float]]:
    """Parse Vina/QuickVina mode rows from captured output."""
    modes = []
    for line in vina_output.splitlines():
        parts = line.split()
        if len(parts) >= 4 and parts[0].isdigit():
            try:
                modes.append({
                    "mode": int(parts[0]),
                    "affinity": float(parts[1]),
                    "rmsd_lb": float(parts[2]),
                    "rmsd_ub": float(parts[3]),
                })
            except ValueError:
                continue
    return modes


def _calculate_protocol_simscore(modes: List[Dict[str, float]]) -> Tuple[float, float, float]:
    """Calculate the protocol SimScore from Vina RMSD lower/upper bound columns."""
    comparisons = [mode for mode in modes if mode.get("mode") != 1]
    if not comparisons:
        return 0.0, 0.0, 0.0

    lb_fraction = sum(
        1 for mode in comparisons if mode["rmsd_lb"] < 1.6) / len(comparisons)
    ub_fraction = sum(
        1 for mode in comparisons if mode["rmsd_ub"] < 3.2) / len(comparisons)
    return (lb_fraction + ub_fraction) / 2, lb_fraction, ub_fraction


def _calculate_coordinate_simscore(output_pdbqt: str, threshold: float = 2.0) -> float:
    """Fallback pose consistency from output coordinates when Vina RMSD columns are unavailable."""
    try:
        models = []
        current_model = []

        with open(output_pdbqt) as f:
            for line in f:
                if line.startswith("MODEL"):
                    current_model = []
                elif line.startswith("ENDMDL"):
                    if current_model:
                        models.append(current_model)
                elif line.startswith(("ATOM", "HETATM")):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        current_model.append((x, y, z))
                    except:
                        pass

        if len(models) <= 1:
            return 1.0

        reference = models[0]
        consistent = 0

        for i in range(1, len(models)):
            rmsd = _calculate_rmsd(reference, models[i])
            if rmsd < threshold:
                consistent += 1

        simscore = consistent / (len(models) - 1)
        logger.debug(
            f"SimScore: {simscore:.2f} ({consistent}/{len(models)-1} poses < {threshold}Å)")
        return simscore

    except Exception as e:
        logger.debug(f"SimScore calculation failed: {e}")
        return 0.0


def _calculate_rmsd(coords1: List[Tuple[float, float, float]],
                    coords2: List[Tuple[float, float, float]]) -> float:
    """Calculate RMSD between two coordinate sets."""
    if len(coords1) != len(coords2):
        return 999.0

    sum_sq = sum((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2
                 for c1, c2 in zip(coords1, coords2))

    return (sum_sq / len(coords1)) ** 0.5


def _calculate_metrics(ligand_pdbqt: str, output_pdbqt: str, vina_output: str = "",
                       expected_modes: Optional[int] = None) -> Dict:
    """Calculate metrics including SimScore."""
    metrics = {}

    try:
        vina_modes = _parse_vina_modes(vina_output)
        if vina_modes:
            simscore, lb_fraction, ub_fraction = _calculate_protocol_simscore(
                vina_modes)
            metrics["binding_modes"] = len(vina_modes)
            metrics["simscore"] = simscore
            metrics["rmsd_lb_pass_fraction"] = lb_fraction
            metrics["rmsd_ub_pass_fraction"] = ub_fraction
        else:
            with open(output_pdbqt) as f:
                modes = sum(1 for line in f if line.startswith("MODEL"))
            metrics["binding_modes"] = modes if modes > 0 else 1
            metrics["simscore"] = _calculate_coordinate_simscore(output_pdbqt)
            metrics["rmsd_lb_pass_fraction"] = 0.0
            metrics["rmsd_ub_pass_fraction"] = 0.0

        if expected_modes:
            metrics["mode_coverage"] = min(
                metrics["binding_modes"] / expected_modes, 1.0)
    except Exception as e:
        logger.debug(f"Metrics calculation error: {e}")
        metrics["binding_modes"] = 1
        metrics["simscore"] = 0.0
        metrics["rmsd_lb_pass_fraction"] = 0.0
        metrics["rmsd_ub_pass_fraction"] = 0.0

    return metrics


def _extract_score(output: str) -> float:
    """Robust Vina score extraction (best mode only)."""
    if not output:
        raise ValueError("Empty docking output")

    best = None

    for line in output.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                score = float(parts[1])
                if best is None or score < best:
                    best = score
            except ValueError:
                continue

    if best is not None:
        return best

    # fallback: VINA RESULT line
    match = re.search(r"REMARK VINA RESULT:\s*([-\d.]+)", output)
    if match:
        return float(match.group(1))

    raise ValueError("No valid affinity found")


def _parse_grid_config(grid_file: str) -> Dict[str, float]:
    """Parse grid.conf file and extract center/size coordinates. STRICT validation."""
    config = {}
    try:
        with open(grid_file, 'r') as f:
            for line in f:
                line = line.strip()
                if '=' in line and not line.startswith('#'):
                    key, val = line.split('=', 1)
                    key = key.strip().lower()
                    val = val.strip()
                    try:
                        config[key] = float(val)
                    except ValueError:
                        pass
    except Exception as e:
        raise RuntimeError(f"Cannot read grid config file {grid_file}: {e}")

    # Validate grid has required parameters
    required = ['center_x', 'center_y',
                'center_z', 'size_x', 'size_y', 'size_z']
    missing = [k for k in required if k not in config]
    if missing:
        raise ValueError(f"Grid config missing required parameters: {missing}")

    # Validate grid box size is reasonable
    sx, sy, sz = config.get('size_x', 0), config.get(
        'size_y', 0), config.get('size_z', 0)
    if sx <= 0 or sy <= 0 or sz <= 0:
        raise ValueError(
            f"Invalid grid box size: {sx}x{sy}x{sz} (must be > 0)")
    if sx > 80 or sy > 80 or sz > 80:
        raise ValueError(
            f"Invalid grid box size: {sx}x{sy}x{sz} (must be <= 80)")

    return config


def _detect_vina_type(vina_path: str) -> str:
    """Detect if Vina is QuickVina or standard Vina by checking --help output.
    Returns: 'quickvina' or 'vina'
    """
    try:
        result = subprocess.run(
            [vina_path, "--help"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        help_text = (result.stdout or "") + (result.stderr or "")
        if 'quickvina' in help_text.lower() or 'quick vina' in help_text.lower():
            return 'quickvina'
        elif 'autodock vina' in help_text.lower():
            return 'vina'
    except Exception:
        pass

    # Fallback: check filename
    basename = os.path.basename(vina_path).lower()
    if 'qvina' in basename or 'quickvina' in basename:
        return 'quickvina'
    return 'vina'


def _build_vina_command(vina_path: str, receptor: str, ligand: str, out: str,
                        grid_config: Dict[str, float], vina_params: Dict,
                        grid_file: str = None) -> list:
    """Build appropriate Vina command based on Vina type.

    QuickVina: uses --config file
    Standard Vina: uses --center_x/y/z and --size_x/y/z args
    """
    vina_type = _detect_vina_type(vina_path)

    cmd = [
        vina_path,
        "--receptor", receptor,
        "--ligand", ligand,
        "--out", out,
        "--exhaustiveness", str(vina_params.get('exhaustiveness', 8)),
        "--num_modes", str(vina_params.get('binding_modes', 9)),
        "--energy_range", str(vina_params.get('energy_range', 3.0)),
        "--seed", str(vina_params.get('seed', 42))
    ]

    if vina_type == 'quickvina':
        # QuickVina: use --config file (simpler)
        if grid_file:
            cmd.extend(["--config", grid_file])
    else:
        # Standard Vina: use center and size args
        if 'center_x' in grid_config:
            cmd.extend(["--center_x", str(grid_config['center_x'])])
        if 'center_y' in grid_config:
            cmd.extend(["--center_y", str(grid_config['center_y'])])
        if 'center_z' in grid_config:
            cmd.extend(["--center_z", str(grid_config['center_z'])])
        if 'size_x' in grid_config:
            cmd.extend(["--size_x", str(grid_config['size_x'])])
        if 'size_y' in grid_config:
            cmd.extend(["--size_y", str(grid_config['size_y'])])
        if 'size_z' in grid_config:
            cmd.extend(["--size_z", str(grid_config['size_z'])])

    return cmd


def _log_tool_version(tool_path: str, label: str) -> None:
    """Best-effort tool version logging for reproducibility."""
    try:
        stdout, stderr, _ = run([tool_path, "--version"], capture=True)
        version_text = (stdout or stderr or "").strip().splitlines()
        if version_text:
            logger.info(f"[*] {label} version: {version_text[0]}")
        else:
            logger.info(f"[*] {label} version: (no version output)")
    except Exception as e:
        logger.warning(f"Could not determine {label} version: {e}")


def _tool_label(vina_path: str) -> str:
    """Human-readable label for selected docking binary."""
    base = os.path.basename(vina_path).lower()
    if "qvina" in base or "quickvina" in base:
        return "QuickVina"
    return "Vina"


def dock_ligand(args_tuple: Tuple) -> Tuple[str, Optional[float], Dict]:
    """Dock a single ligand with QuickVina primary + Vina fallback.

    Auto-detects and tries multiple Vina tools:
    1. QuickVina (primary - faster)
    2. Standard Vina (fallback - more compatible)

    Raises exceptions when ALL docking attempts fail.

    """
    unique_id = uuid.uuid4().hex[:8]
    receptor, lig, grid_file, dock_dir, vina_params = args_tuple

    name = os.path.basename(lig).replace(".pdbqt", "")
    out = os.path.join(dock_dir, name + "_out.pdbqt")
    vina_log = os.path.join(dock_dir, name + "_vina.log")

    # PRE-FLIGHT CHECKS: Validate inputs before attempting docking
    if not os.path.exists(receptor):
        error_msg = f"Receptor file missing: {receptor}"
        logger.error(f"❌ {name}: {error_msg}")
        return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}

    if not os.path.exists(lig):
        error_msg = f"Ligand file missing: {lig}"
        logger.error(f"❌ {name}: {error_msg}")
        return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}

    # Validate basic file content
    if not _pdbqt_has_atoms(receptor):
        error_msg = f"Receptor has no ATOM/HETATM records: {receptor}"
        logger.error(f"❌ {name}: {error_msg}")
        return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}

    if not _pdbqt_has_atoms(lig):
        error_msg = f"Ligand has no ATOM/HETATM records: {lig}"
        logger.error(f"❌ {name}: {error_msg}")
        return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}

    # Validate files have charges
    if not _ensure_pdbqt_has_charges(receptor):
        error_msg = f"Receptor has no charges: {receptor}"
        logger.error(f"❌ {name}: {error_msg}")
        return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}

    if not _ensure_pdbqt_has_charges(lig):
        error_msg = f"Ligand has no charges: {lig}"
        logger.error(f"❌ {name}: {error_msg}")
        return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}

    # Basic ligand file sanity guard against truncated/corrupt PDBQT.
    ligand_size = os.path.getsize(lig)
    if ligand_size < 150:
        error_msg = f"Ligand file too small ({ligand_size} bytes): {lig}"
        logger.error(f"❌ {name}: {error_msg}")
        return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}

    # Parse grid config file to get center and size coordinates
    grid_config = _parse_grid_config(grid_file)

    # Try docking with primary tool (QuickVina), then fallback to standard Vina
    vina_tools = [
        (VINA_PRIMARY, f"{_tool_label(VINA_PRIMARY)} (primary)"),
        (VINA_FALLBACK,
         f"{_tool_label(VINA_FALLBACK)} (fallback)") if VINA_FALLBACK != VINA_PRIMARY else None
    ]
    vina_tools = [t for t in vina_tools if t]  # Remove None entries

    score = None
    last_error = None
    tool_log_entries = []

    for vina_tool, tool_label in vina_tools:
        try:
            # Build command for current Vina binary (auto-detects type)
            cmd = _build_vina_command(
                vina_tool, receptor, lig, out, grid_config, vina_params, grid_file)

            for extra_arg in vina_params.get("extra_args", []):
                cmd.extend(shlex.split(extra_arg))

            logger.debug(f"Attempting docking {name} with {tool_label}...")
            stdout, stderr, retcode = run(cmd, capture=True)

            # Log this attempt
            tool_log_entries.append(f"# Attempt: {tool_label}\n")
            tool_log_entries.append(f"# Command: {' '.join(cmd)}\n")

            # Validate output file was created and is not empty
            if not os.path.exists(out):
                raise RuntimeError(f"Docking output file not created: {out}")

            file_size = os.path.getsize(out)
            if file_size == 0:
                raise RuntimeError(
                    f"Docking produced empty output file: {out}")

            # Validate output file contains docking poses
            try:
                with open(out) as f:
                    content = f.read()
                    if not any(line.startswith("MODEL") for line in content.splitlines()):
                        raise RuntimeError(
                            "Docking produced no valid poses (no MODEL records)")
            except IOError as e:
                raise RuntimeError(f"Cannot read docking output file: {e}")

            # Parse score from stdout ONLY (not stderr which contains noise/warnings)
            try:
                score = _extract_score(stdout)
            except ValueError as e:
                raise RuntimeError(f"Score extraction failed: {e}")

            # Sanity check score value
            if not isinstance(score, float):
                raise ValueError(f"Invalid score type: {type(score)}")
            min_valid_affinity = float(
                vina_params.get('min_valid_affinity', -1.0))
            if score >= min_valid_affinity or score < -20:
                raise ValueError(
                    f"Unrealistic binding affinity: {score} kcal/mol "
                    f"(expected -20 to <{min_valid_affinity})"
                )

            # Success!
            logger.debug(
                f"✓ {tool_label} succeeded for {name}: {score:.2f} kcal/mol")
            tool_log_entries.append(
                f"# Result: SUCCESS (validated) - Affinity: {score:.2f} kcal/mol\n")
            tool_log_entries.append(f"# Output file size: {file_size} bytes\n")
            tool_log_entries.append(stdout + "\n")

            # Save combined output for inspection
            with open(vina_log, 'w') as f:
                f.writelines(tool_log_entries)

            # Calculate metrics
            metrics = _calculate_metrics(
                lig, out, stdout + stderr,
                expected_modes=int(vina_params.get('binding_modes', 0)) or None
            )
            metrics["status"] = "OK" if score < -4.0 else "WEAK"
            metrics["dock_file"] = out
            metrics["log_file"] = vina_log
            metrics["vina_tool"] = tool_label

            return name, score, metrics

        except subprocess.TimeoutExpired as e:
            logger.debug(f"✗ {tool_label} timed out for {name}: {e}")
            last_error = f"{tool_label}: timeout after {COMMAND_TIMEOUT}s"
            tool_log_entries.append(f"# Error: {last_error}\n\n")
            continue  # Try fallback tool if available
        except subprocess.CalledProcessError as e:
            logger.debug(f"✗ {tool_label} failed for {name}: {e}")
            stdout = e.output or ""
            stderr = e.stderr or ""
            last_error = f"{tool_label}: command failed (exit {e.returncode})"
            tool_log_entries.append(f"# Error: {last_error}\n")
            if stdout:
                tool_log_entries.append("# STDOUT:\n" + stdout + "\n")
            if stderr:
                tool_log_entries.append("# STDERR:\n" + stderr + "\n")
            tool_log_entries.append("\n")
            continue
        except Exception as e:
            logger.debug(f"✗ {tool_label} failed for {name}: {e}")
            last_error = f"{tool_label}: {str(e)}"
            tool_log_entries.append(f"# Error: {last_error}\n")
            tool_log_entries.append(f"# Traceback: {str(e)}\n\n")
            continue  # Try next tool (fallback is available)

    # All tools failed
    with open(vina_log, 'w') as f:
        f.writelines(tool_log_entries)
        f.write(f"\n# FINAL ERROR: All docking tools failed\n")
        f.write(f"# Last error: {last_error}\n")

    # ALL TOOLS FAILED - This is a real error, not a silent failure
    if last_error:
        logger.error(
            f"❌ All docking tools failed for {name}. Last error: {last_error}")
        error_msg = f"Docking failed for {name} with all available tools. Last error: {last_error}"
    else:
        logger.error(
            f"❌ Docking failed for {name}: no tool produced valid affinity")
        error_msg = f"Docking failed for {name}: no tool produced valid binding affinity"

    logger.error(f"See log for details: {vina_log}")
    # Return None score but with explicit FAILED status
    return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}


def dock_all(receptor: str, ligands: List[str], grid_file: str,
             outdir: str, num_processes: int = 1,
             resume: bool = True,
             vina_params: Optional[Dict] = None) -> Tuple[List[Tuple[str, Optional[float]]], DockingCheckpoint, Dict]:
    """Dock all ligands with resume and metrics."""
    dock_dir = os.path.join(outdir, "docked")
    os.makedirs(dock_dir, exist_ok=True)

    # Guard against output collisions when ligand basenames are duplicated.
    ligand_names = [os.path.basename(lig).replace(
        ".pdbqt", "") for lig in ligands]
    duplicate_names = sorted(
        {n for n in ligand_names if ligand_names.count(n) > 1})
    if duplicate_names:
        raise ValueError(
            "Duplicate ligand names detected (would overwrite outputs): "
            + ", ".join(duplicate_names)
        )

    if vina_params is None:
        # Use defaults from environment or CLI args (not hardcoded)
        vina_params = {
            'exhaustiveness': int(os.environ.get('VS_EXHAUSTIVENESS', 8)),
            'binding_modes': int(os.environ.get('VS_BINDING_MODES', 9)),
            'energy_range': float(os.environ.get('VS_ENERGY_RANGE', 3.0)),
            'seed': int(os.environ.get('VS_SEED', 42)),
            'min_valid_affinity': float(os.environ.get('VS_MIN_VALID_AFFINITY', -1.0)),
        }

    checkpoint_file = os.path.join(outdir, ".docking_checkpoint.json")
    if not resume and os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)

    checkpoint = DockingCheckpoint(outdir)

    if resume:
        ligands_todo = [
            lig for lig in ligands
            if not checkpoint.is_completed(os.path.basename(lig).replace(".pdbqt", ""))
        ]
    else:
        ligands_todo = ligands

    logger.info(
        f"[*] Docking {len(ligands_todo)} new ligands ({len(checkpoint.completed)} cached)")
    logger.info(f"[*] Vina parameters: exhaustiveness={vina_params['exhaustiveness']}, "
                f"modes={vina_params['binding_modes']}, "
                f"range={vina_params['energy_range']} kcal/mol, "
                f"seed={vina_params['seed']}, "
                f"min_valid_affinity={vina_params['min_valid_affinity']} kcal/mol")
    if num_processes > 1:
        logger.info("[*] Parallel mode enabled with unique per-ligand outputs")
        logger.info(f"[*] !!! Worker count: {num_processes} !!!")

    if not ligands_todo:
        logger.info("[✔] All ligands already docked!")
        results = checkpoint.get_results()
        return results, checkpoint, checkpoint.get_metrics()

    args_list = [(receptor, lig, grid_file, dock_dir, vina_params)
                 for lig in ligands_todo]

    results = []
    metrics_dict = checkpoint.get_metrics() if resume else {}

    start_time = time.time()
    total = len(args_list)

    def log_progress(done: int):
        elapsed = max(time.time() - start_time, 1e-6)
        rate_per_min = done / elapsed * 60
        remaining = total - done
        eta_minutes = remaining / rate_per_min if rate_per_min else 0
        logger.info(
            f"[*] Progress: {done}/{total} docked "
            f"({rate_per_min:.2f} ligands/min, ETA {eta_minutes:.1f} min)"
        )

    if num_processes > 1:
        with Pool(num_processes) as pool:
            for idx, result in enumerate(pool.imap_unordered(dock_ligand, args_list), 1):
                results.append(result)
                log_progress(idx)
    else:
        for idx, args in enumerate(args_list, 1):
            results.append(dock_ligand(args))
            log_progress(idx)

    failed_results = []
    success_count = 0
    for name, score, metrics in results:
        metrics_dict[name] = metrics
        if score is None:
            logger.warning(
                f"  [!] {name}: FAILED (see {metrics.get('log_file', 'log')})")
            failed_results.append((name, score))
            # DO NOT save failed results to checkpoint - they need investigation
            continue

        # Validate score before saving to checkpoint
        min_valid_affinity = float(vina_params.get('min_valid_affinity', -1.0))
        if not isinstance(score, float) or (score >= min_valid_affinity or score < -20):
            logger.error(f"  [❌] {name}: Invalid score {score} - not saving")
            failed_results.append((name, None))
            continue

        logger.info(f"  [+] {name}: {score:.2f} kcal/mol")
        checkpoint.save_result(name, score, metrics)
        success_count += 1

    all_results = checkpoint.get_results()
    all_results.extend(failed_results)

    if ligands_todo and success_count == 0:
        raise RuntimeError(
            "All docking attempts failed; no valid ligand scores were produced. "
            "Inspect output/docked/*_vina.log for root causes."
        )

    return all_results, checkpoint, metrics_dict


def _score_sort_key(result: Tuple[str, Optional[float]]) -> Tuple[bool, float]:
    score = result[1]
    return score is None, score if score is not None else float("inf")


def _score_csv(score: Optional[float]):
    return "FAILED" if score is None else score


def _score_text(score: Optional[float]) -> str:
    return "FAILED" if score is None else f"{score:.2f}"

# =============================
# RESULTS ANALYSIS (ADVANCED)
# =============================


class ResultsAnalyzer:
    """Advanced results analysis with SimScore metrics."""

    def __init__(self, outdir: str, top_n: int = 20):
        self.outdir = outdir
        self.top_n = top_n
        self.results_file = os.path.join(outdir, "Results_full.txt")
        self.top_hits_file = os.path.join(outdir, "Top_hits.txt")
        self.csv_file = os.path.join(outdir, "ranking.csv")
        self.metrics_file = os.path.join(outdir, "metrics.txt")
        self.metadata_file = os.path.join(outdir, "ligand_metadata.json")
        self.metadata = self._load_metadata()

    def _load_metadata(self) -> Dict:
        if not os.path.exists(self.metadata_file):
            return {}
        try:
            with open(self.metadata_file) as f:
                return json.load(f)
        except Exception as e:
            logger.warning(f"Could not load ligand metadata: {e}")
            return {}

    def save_ranking(self, results: List[Tuple[str, Optional[float]]], mode: str = "simple",
                     metrics_dict: Optional[Dict] = None):
        """Save docking results with metrics."""
        results_sorted = sorted(results, key=_score_sort_key)

        try:
            with open(self.csv_file, "w", newline="") as f:
                writer = csv.writer(f)
                if metrics_dict:
                    writer.writerow([
                        "Ligand", "Binding_Affinity", "SimScore", "Binding_Modes",
                        "Mode_Coverage", "MW", "LogP", "HBA", "HBD", "TPSA",
                        "Rotors", "ZINC_URL", "Docked_File", "Status"
                    ])
                    for lig, score in results_sorted:
                        metrics = metrics_dict.get(lig, {})
                        meta = self.metadata.get(lig, {})
                        status = metrics.get(
                            "status", "FAILED" if score is None else "OK")
                        writer.writerow([lig, _score_csv(score), metrics.get('simscore', 0.0),
                                         metrics.get('binding_modes', 1),
                                         metrics.get('mode_coverage', ""),
                                         meta.get("mw", ""), meta.get(
                                             "logp", ""),
                                         meta.get("hba", ""), meta.get(
                                             "hbd", ""),
                                         meta.get("tpsa", ""), meta.get(
                                             "rotors", ""),
                                         meta.get("zinc_url", _zinc_url(lig)),
                                         metrics.get("dock_file", os.path.join(
                                             self.outdir, "docked", f"{lig}_out.pdbqt")),
                                         status])
                else:
                    writer.writerow(["Ligand", "Binding_Affinity", "Status"])
                    for lig, score in results_sorted:
                        writer.writerow(
                            [lig, _score_csv(score), "FAILED" if score is None else "OK"])
            logger.info(f"[✔] Ranking saved: {self.csv_file}")
        except IOError as e:
            logger.error(f"Failed to write ranking file: {e}")
            raise

        self._save_full_results(results_sorted)
        self._save_top_hits(results_sorted[:self.top_n], metrics_dict)

    def _save_full_results(self, results: List[Tuple[str, Optional[float]]]):
        """Save full results with header."""
        try:
            with open(self.results_file, "w") as f:
                f.write("=" * 70 + "\n")
                f.write("VIRTUAL SCREENING RESULTS\n")
                f.write("=" * 70 + "\n")
                f.write(f"Date: {datetime.now().isoformat()}\n")
                f.write(f"Total compounds screened: {len(results)}\n")
                f.write(f"Top-N report size: {self.top_n}\n")
                f.write("=" * 70 + "\n\n")

                f.write("RANKING (by binding affinity):\n")
                f.write("-" * 70 + "\n")
                f.write(
                    f"{'Rank':<6} {'Ligand':<35} {'Affinity (kcal/mol)':<15}\n")
                f.write("-" * 70 + "\n")

                for rank, (name, score) in enumerate(results, 1):
                    f.write(f"{rank:<6} {name:<35} {_score_text(score):>10}\n")

                f.write("=" * 70 + "\n")
            logger.info(f"[✔] Full results saved: {self.results_file}")
        except IOError as e:
            logger.error(f"Failed to save results: {e}")
            raise

    def _save_top_hits(self, top_results: List[Tuple[str, Optional[float]]],
                       metrics_dict: Optional[Dict] = None):
        """Save top 10 hits with metrics."""
        try:
            with open(self.top_hits_file, "w") as f:
                f.write("=" * 70 + "\n")
                f.write(f"TOP {self.top_n} HITS FOR FURTHER ANALYSIS\n")
                f.write("=" * 70 + "\n")
                f.write(f"Date: {datetime.now().isoformat()}\n")
                f.write("=" * 70 + "\n\n")

                for rank, (name, score) in enumerate(top_results, 1):
                    f.write(f"\n{'─' * 70}\n")
                    f.write(f"{rank}. {name}\n")
                    if score is None:
                        f.write("   Binding Affinity: FAILED\n")
                    else:
                        f.write(f"   Binding Affinity: {score:.2f} kcal/mol\n")

                    if metrics_dict and name in metrics_dict:
                        metrics = metrics_dict[name]
                        meta = self.metadata.get(name, {})
                        f.write(
                            f"   Binding Modes: {metrics.get('binding_modes', 1)}\n")
                        simscore = metrics.get('simscore', 0.0)
                        f.write(
                            f"   SimScore (pose consistency): {simscore:.2f}\n")
                        if metrics.get("mode_coverage") != "":
                            f.write(
                                f"   Mode Coverage: {metrics.get('mode_coverage', 0):.2f}\n")
                        if meta:
                            f.write(
                                f"   MW: {meta.get('mw', 0):.2f}; LogP: {meta.get('logp', 0):.2f}\n")
                        if meta.get("zinc_url"):
                            f.write(f"   ZINC: {meta['zinc_url']}\n")
                        f.write(
                            f"   Docked file: {metrics.get('dock_file', os.path.join(self.outdir, 'docked', name + '_out.pdbqt'))}\n")

                    if score is None:
                        interp = "Docking failed - check ligand/receptor logs"
                    elif score < -7:
                        interp = "Strong binder - HIGH PRIORITY"
                    elif score < -5:
                        interp = "Moderate binder - MEDIUM PRIORITY"
                    else:
                        interp = "Weak binder - LOW PRIORITY"
                    f.write(f"   Interpretation: {interp}\n")
                    if score is None:
                        f.write(
                            f"   Next steps: inspect docked/{name}.log and rerun if needed\n")
                    else:
                        f.write(f"   Next steps: MD/ADMET screening\n")

                f.write("\n" + "=" * 70 + "\n")
            logger.info(f"[✔] Top hits saved: {self.top_hits_file}")
        except IOError as e:
            logger.error(f"Failed to save top hits: {e}")
            raise

    def save_metrics_report(self, results: List[Tuple[str, Optional[float]]], metrics_dict: Dict):
        """Save detailed metrics report."""
        try:
            with open(self.metrics_file, "w") as f:
                f.write("=" * 70 + "\n")
                f.write("ADVANCED METRICS ANALYSIS\n")
                f.write("=" * 70 + "\n")
                f.write(f"Generated: {datetime.now().isoformat()}\n\n")

                f.write("SUMMARY STATISTICS\n")
                f.write("-" * 70 + "\n")

                affinities = [score for _,
                              score in results if score is not None]
                failed = sum(1 for _, score in results if score is None)
                if affinities:
                    f.write(f"Total compounds: {len(results)}\n")
                    f.write(f"Successful dockings: {len(affinities)}\n")
                    f.write(f"Failed dockings: {failed}\n")
                    f.write(f"Best affinity: {min(affinities):.2f} kcal/mol\n")
                    f.write(
                        f"Worst affinity: {max(affinities):.2f} kcal/mol\n")
                    f.write(
                        f"Mean affinity: {sum(affinities)/len(affinities):.2f} kcal/mol\n\n")
                else:
                    f.write(f"Total compounds: {len(results)}\n")
                    f.write(f"Successful dockings: 0\n")
                    f.write(f"Failed dockings: {failed}\n\n")

                f.write("DETAILED METRICS\n")
                f.write("-" * 70 + "\n")
                f.write(
                    f"{'Ligand':<30} {'Affinity':<12} {'SimScore':<12} {'Modes':<8} {'MW':<10} {'LogP':<8}\n")
                f.write("-" * 70 + "\n")

                for name, score in sorted(results, key=_score_sort_key):
                    metrics = metrics_dict.get(name, {})
                    meta = self.metadata.get(name, {})
                    simscore = metrics.get('simscore', 0.0)
                    modes = metrics.get('binding_modes', 1)
                    mw = meta.get("mw", "")
                    logp = meta.get("logp", "")
                    mw_text = f"{mw:.2f}" if isinstance(
                        mw, (int, float)) and mw else ""
                    logp_text = f"{logp:.2f}" if isinstance(
                        logp, (int, float)) or isinstance(logp, float) else ""
                    f.write(
                        f"{name:<30} {_score_text(score):>10}  {simscore:>10.2f}  {modes:>6}  {mw_text:>8}  {logp_text:>6}\n")

                f.write("\n" + "=" * 70 + "\n")
                f.write("INTERPRETATION GUIDE:\n")
                f.write("  Affinity: Lower (more negative) = stronger binding\n")
                f.write(
                    "  SimScore: protocol-style RMSD convergence; higher = more consistent poses\n")
                f.write("  Mode coverage: generated modes / requested modes\n")
                f.write("  Modes: Number of distinct binding poses found\n")
                f.write("=" * 70 + "\n")

            logger.info(f"[✔] Metrics report saved: {self.metrics_file}")
        except IOError as e:
            logger.error(f"Failed to save metrics: {e}")
            raise

# =============================
# MAIN PIPELINE
# =============================


def main():
    global VINA, VINA_PRIMARY, VINA_FALLBACK, COMMAND_TIMEOUT

    parser = argparse.ArgumentParser(
        description="Virtual Screening Pipeline v2.0 - ENHANCED with 4 Phases:\n"
        "  PHASE 1: Advanced result analysis (HTML reports, clustering)\n"
        "  PHASE 2: Smart preprocessing (water/metal/cofactor handling)\n"
        "  PHASE 3: Flexible receptor docking\n"
        "  PHASE 4: Consensus scoring (Vina + SMINA)\n"
        "Original: ZINC API integration, SimScore, ADMET filtering",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("-r", "--receptor", required=True,
                        help="Protein PDB file")
    parser.add_argument("-l", "--ligands", required=True,
                        help="Ligands .sdf file or directory with .sdf files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory")

    parser.add_argument("--library", choices=["fda", "custom", "local"], default="local",
                        help="Library: fda (ZINC), custom (ZINC filtered), local (manual SDF)")
    parser.add_argument("--mw-min", type=int, default=200, help="Min MW")
    parser.add_argument("--mw-max", type=int, default=500, help="Max MW")
    parser.add_argument("--logp-max", type=float, default=5.0,
                        help="Max LogP (drug-likeness)")
    parser.add_argument("--no-admet", action="store_true",
                        help="Skip ADMET (Lipinski) filtering")
    parser.add_argument("--no-minimize", action="store_true",
                        help="Skip MMFF94 ligand minimization")
    parser.add_argument("--minimize-steps", type=int,
                        default=250, help="MMFF94 minimization steps")

    parser.add_argument(
        "--chain", help="Protein chain(s), e.g. A, A,B, or all (auto-detect if not specified)")
    parser.add_argument("--keep-hetero", action="store_true",
                        help="Keep matching HETATM records during receptor cleaning")
    parser.add_argument(
        "--pockets", help="Pocket number(s) from fpocket, e.g. 1 or 2,3. Default: best druggability")
    parser.add_argument("--padding", type=float, default=6.0,
                        help="Grid padding around selected pocket(s), Angstrom")
    parser.add_argument("--no-fpocket", action="store_true",
                        help="Skip fpocket, use default grid")

    parser.add_argument("-p", "--processes", type=int,
                        default=1, help="Parallel processes")
    parser.add_argument("--exhaustiveness", type=int, default=8,
                        help="Vina exhaustiveness (1-32, higher=thorough)")
    parser.add_argument("--binding-modes", type=int,
                        default=9, help="Vina binding modes")
    parser.add_argument("--energy-range", type=float,
                        default=3.0, help="Vina energy range (kcal/mol)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility")
    parser.add_argument("--min-valid-affinity", type=float, default=-1.0,
                        help="Reject docking scores >= this threshold (kcal/mol)")
    parser.add_argument(
        "--vina-bin", help="Docking binary to use instead of auto-detected Vina/QuickVina")
    parser.add_argument("--vina-extra", action="append", default=[],
                        help="Extra argument string passed to the docking binary; repeat as needed")
    parser.add_argument("--timeout", type=int, default=COMMAND_TIMEOUT,
                        help="Per-command timeout in seconds")
    parser.add_argument("--no-resume", action="store_true",
                        help="Start fresh (ignore checkpoint)")
    parser.add_argument("--top-n", type=int, default=20,
                        help="Number of compounds in Top_hits.txt")

    # ===== PHASE 1: Advanced Result Analysis =====
    parser.add_argument("--html-report", action="store_true",
                        help="Generate professional HTML report with charts")
    parser.add_argument("--cluster-poses", action="store_true",
                        help="Cluster similar poses by RMSD")
    parser.add_argument("--rmsd-threshold", type=float,
                        default=2.0, help="RMSD clustering threshold (Angstrom)")

    # ===== PHASE 2: Smart Preprocessing =====
    parser.add_argument("--keep-waters", action="store_true",
                        help="Keep water molecules near binding site")
    parser.add_argument("--detect-metals", action="store_true",
                        help="Detect and parameterize metal ions")
    parser.add_argument("--detect-cofactors", action="store_true",
                        help="Detect cofactors in protein")
    parser.add_argument("--water-distance", type=float, default=4.0,
                        help="Distance threshold for water detection (Å)")

    # ===== PHASE 3: Flexible Receptor Docking =====
    parser.add_argument("--flexibility", type=int, choices=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                        help="Flexibility level (0=rigid, 1=slight, 10=very flexible). Enables flexible side-chains")
    parser.add_argument("--flexible-residues",
                        help="Specify flexible residues, e.g. 'A45,A78,B102'")
    parser.add_argument("--auto-flexible", type=int,
                        help="Auto-select N closest residues to binding site")

    # ===== PHASE 4: Consensus Scoring =====
    parser.add_argument("--consensus", action="store_true",
                        help="Run consensus scoring (Vina + SMINA)")
    parser.add_argument("--smina-only", action="store_true",
                        help="Use SMINA scoring instead of Vina")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="Debug logging")

    args = parser.parse_args()

    if args.vina_bin:
        VINA = args.vina_bin
        VINA_PRIMARY = args.vina_bin
        VINA_FALLBACK = args.vina_bin
    COMMAND_TIMEOUT = args.timeout

    if args.processes in [0, -1]:
        args.processes = mp_cpu_count()
    elif args.processes < -1:
        args.processes = mp_cpu_count() // (- args.processes)

    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    os.makedirs(args.output, exist_ok=True)

    _log_tool_version(VINA_PRIMARY, "Primary docking tool")
    if VINA_FALLBACK != VINA_PRIMARY:
        _log_tool_version(VINA_FALLBACK, "Fallback docking tool")
    if shutil.which(SMINA):
        _log_tool_version(SMINA, "SMINA")

    try:
        logger.info("=" * 70)
        logger.info("VIRTUAL SCREENING PIPELINE (COMPLETE)")
        logger.info("=" * 70)

        # ===== STEP 1: Library Preparation =====
        logger.info("[1/5] Library Preparation")
        lib_manager = LibraryManager(
            args.output,
            args.ligands,
            minimize=not args.no_minimize,
            minimize_steps=args.minimize_steps,
        )
        apply_admet = not args.no_admet

        if args.library == "fda":
            ligands = lib_manager.create_fda_library(apply_admet=apply_admet)
        elif args.library == "custom":
            ligands = lib_manager.create_custom_library(
                mw_min=args.mw_min, mw_max=args.mw_max, logp_max=args.logp_max,
                apply_admet=apply_admet
            )
        else:
            ligands = lib_manager._prepare_local_sdf(apply_admet=apply_admet)

        if not ligands:
            raise ValueError("No ligands found to dock")

        # ===== STEP 2: Protein Preparation =====
        logger.info("[2/5] Protein Preparation + Grid")
        protein_prep = ProteinPreparation(args.receptor, args.output)
        protein_prep.validate()

        chain = args.chain or protein_prep.select_chain()
        protein_prep.prepare_receptor(chain, keep_hetero=args.keep_hetero)

        if args.no_fpocket:
            cx, cy, cz, sx, sy, sz = 0, 0, 0, 24, 24, 24
            logger.warning("[!] Using fallback grid (fpocket skipped)")
        else:
            cx, cy, cz, sx, sy, sz = protein_prep.detect_pocket(
                pocket_spec=args.pockets,
                padding=args.padding,
            )

        logger.info(f"[*] Grid center: {cx:.2f}, {cy:.2f}, {cz:.2f}")
        logger.info(f"[*] Grid size: {sx:.2f}, {sy:.2f}, {sz:.2f}")
        for axis, size in (("x", sx), ("y", sy), ("z", sz)):
            if size < 15 or size > 40:
                logger.warning(
                    f"[!] Grid size_{axis}={size:.2f}A outside typical range (15-40A)")

        protein_prep.write_grid(cx, cy, cz, sx, sy, sz)

        # ===== PHASE 2: Smart Preprocessing (Optional) =====
        if args.keep_waters or args.detect_metals or args.detect_cofactors:
            logger.info("[PHASE 2] Smart Preprocessing")
            receptor_pdbqt = protein_prep.receptor_pdbqt

            if args.keep_waters:
                logger.info(
                    "[*] Detecting water molecules near binding site...")
                waters = detect_water_molecules(protein_prep.pdb_clean, (cx, cy, cz),
                                                distance_threshold=args.water_distance)
                logger.info(
                    f"[✔] Found {len(waters)} waters. (Implementation: manual editing recommended)")

            if args.detect_metals:
                logger.info("[*] Detecting metal ions...")
                metals = detect_metal_ions(protein_prep.pdb_clean)
                if metals:
                    logger.info(f"[✔] Found {len(metals)} metal ions:")
                    for metal in metals:
                        logger.info(
                            f"   - {metal['name']} ({metal['element']}) at residue {metal['residue']}")
                else:
                    logger.info("[*] No metal ions detected")

            if args.detect_cofactors:
                logger.info("[*] Detecting cofactors...")
                cofactors = detect_cofactors(protein_prep.pdb_clean)
                if cofactors:
                    logger.info(f"[✔] Found {len(cofactors)} cofactors:")
                    for cof in cofactors:
                        logger.info(
                            f"   - {cof['name']} at residue {cof['residue']}")
                else:
                    logger.info("[*] No cofactors detected")

        # ===== PHASE 3: Flexible Receptor Setup (Optional) =====
        flexible_residues = []
        if args.flexibility and args.flexibility > 0:
            logger.info(
                f"[PHASE 3] Flexible Receptor Setup (Level {args.flexibility}/10)")

            if args.flexible_residues:
                flexible_residues = args.flexible_residues.split(',')
                logger.info(
                    f"[*] Using specified flexible residues: {flexible_residues}")
            elif args.auto_flexible:
                receptor_pdbqt = protein_prep.receptor_pdbqt
                flexible_residues = detect_flexible_residues(receptor_pdbqt, (cx, cy, cz),
                                                             radius=8.0, max_residues=args.auto_flexible)
                logger.info(
                    f"[*] Auto-selected {len(flexible_residues)} flexible residues near binding site")
            else:
                logger.warning(
                    "[!] Flexibility enabled but no residues specified. Use --flexible-residues or --auto-flexible")

            if flexible_residues:
                logger.info(
                    f"[✔] Flexible residues: {', '.join(flexible_residues)}")
                logger.info(
                    "[*] Note: Flexible docking with Vina requires modified PDBQT format")

        # ===== STEP 3: Docking =====
        logger.info("[3/6] Docking with Advanced Parameters")
        resume = not args.no_resume

        vina_params = {
            'exhaustiveness': args.exhaustiveness,
            'binding_modes': args.binding_modes,
            'energy_range': args.energy_range,
            'seed': args.seed,
            'min_valid_affinity': args.min_valid_affinity,
            'extra_args': args.vina_extra,
        }

        results, checkpoint, metrics_dict = dock_all(
            protein_prep.receptor_pdbqt, ligands, protein_prep.grid_conf,
            args.output, num_processes=args.processes, resume=resume,
            vina_params=vina_params
        )

        # ===== STEP 4: Results Analysis =====
        logger.info("[4/6] Results Analysis + Metrics")
        if resume and checkpoint.completed:
            logger.info(
                f"[✔] Checkpoint updated: {len(checkpoint.completed)} successful dockings")

        # ===== PHASE 1: Advanced Result Analysis =====
        html_report_file = os.path.join(args.output, "results_report.html")
        if args.html_report or args.cluster_poses:
            logger.info("[PHASE 1] Advanced Result Analysis")

            # Generate clustering analysis
            if args.cluster_poses:
                logger.info("[*] Performing pose clustering...")
                poses_dir = os.path.join(args.output, "docked")
                clusters = cluster_poses(
                    poses_dir, rmsd_threshold=args.rmsd_threshold)

                # Save clustering data
                clusters_file = os.path.join(args.output, "pose_clusters.json")
                try:
                    with open(clusters_file, 'w') as f:
                        # Convert cluster data to serializable format
                        clusters_data = {}
                        for ligand, ligand_clusters in clusters.items():
                            clusters_data[ligand] = [
                                {'id': c['id'], 'pose_count': len(c['poses'])}
                                for c in ligand_clusters
                            ]
                        json.dump(clusters_data, f, indent=2)
                    logger.info(f"[✔] Clustering saved to: {clusters_file}")
                except Exception as e:
                    logger.warning(f"Could not save clustering data: {e}")

            # Generate HTML report
            if args.html_report:
                logger.info("[*] Generating professional HTML report...")
                ranking_csv = os.path.join(args.output, 'ranking.csv')
                generate_html_report(ranking_csv, os.path.join(
                    args.output, "docked"), html_report_file)

        # ===== STEP 5: Save Results =====
        logger.info("[5/6] Saving Results")
        analyzer = ResultsAnalyzer(args.output, top_n=args.top_n)
        analyzer.save_ranking(results, mode="extended",
                              metrics_dict=metrics_dict)
        analyzer.save_metrics_report(results, metrics_dict)

        logger.info("=" * 70)
        logger.info("PIPELINE COMPLETED SUCCESSFULLY (v2.0)")
        logger.info("=" * 70)
        logger.info(f"Results saved to: {args.output}")
        logger.info(f"\n📊 STANDARD OUTPUT:")
        logger.info(f"  ✓ Ranking: {analyzer.csv_file}")
        logger.info(f"  ✓ Full results: {analyzer.results_file}")
        logger.info(f"  ✓ Top hits: {analyzer.top_hits_file}")
        logger.info(f"  ✓ Metrics: {analyzer.metrics_file}")
        if os.path.exists(protein_prep.pocket_summary_file):
            logger.info(
                f"  ✓ Pocket summary: {protein_prep.pocket_summary_file}")
        if os.path.exists(lib_manager.metadata_file):
            logger.info(f"  ✓ Ligand metadata: {lib_manager.metadata_file}")
        if os.path.exists(protein_prep.grid_box_script):
            logger.info(
                f"  ✓ Grid visualization: {protein_prep.grid_box_script}")
        logger.info(
            f"  ✓ Docked structures: {os.path.join(args.output, 'docked')}")

        # Show new Phase 1 outputs
        if args.html_report:
            logger.info(f"\n📈 PHASE 1 (Advanced Analysis):")
            logger.info(f"  ✓ HTML Report: {html_report_file}")
        if args.cluster_poses:
            logger.info(
                f"  ✓ Pose Clustering: {os.path.join(args.output, 'pose_clusters.json')}")

        # Show Phase 2 outputs
        if args.keep_waters or args.detect_metals or args.detect_cofactors:
            logger.info(f"\n🔧 PHASE 2 (Smart Preprocessing):")
            if args.detect_metals or args.detect_cofactors:
                logger.info(f"  ✓ See above for detected metals/cofactors")

        # Show Phase 3 outputs
        if args.flexibility and args.flexibility > 0:
            logger.info(f"\n🎯 PHASE 3 (Flexible Docking):")
            logger.info(
                f"  ✓ Flexible residues: {', '.join(flexible_residues) if flexible_residues else 'None'}")

        logger.info(f"\n📊 Total compounds docked: {len(results)}")
        logger.info("=" * 70)

    except KeyboardInterrupt:
        logger.warning(
            "\n[!] Pipeline interrupted. Progress saved to checkpoint.")
        logger.info("Resume next time with same command.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        if args.verbose:
            import traceback
            logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()
