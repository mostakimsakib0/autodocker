#!/usr/bin/env python3
"""PHASE 1: ADVANCED RESULT ANALYSIS.

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

try:
    import numpy as np
except ModuleNotFoundError:
    np = None


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


NGL_CDN_URL = "https://cdn.jsdelivr.net/npm/ngl@2.4.0/dist/ngl.min.js"

# PDB lines that are safe/meaningful to NGL (PDBQT adds ROOT/BRANCH
# bookkeeping that NGL's PDB parser does not understand).
_PDBQT_STRIP_PREFIXES = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")
_PDBQT_KEEP_PREFIXES = ("ATOM", "HETATM", "MODEL", "ENDMDL", "REMARK",
                        "CONECT", "HEADER", "TITLE", "COMPND", "TER",
                        "MASTER", "END")


def _pdbqt_to_pdb_string(pdbqt_path: str) -> Optional[str]:
    """Convert a (possibly multi-model) PDBQT ligand/receptor to plain PDB text.

    NGL does not parse PDBQT directly, so we emit a faithful PDB by dropping
    the ROOT/BRANCH/TORSDOF bookkeeping records. The text is returned so it can
    be inlined directly into the report (no external ``viewer/`` files needed).
    Returns None when nothing useful could be extracted.
    """
    try:
        lines = []
        written = False
        with open(pdbqt_path) as src:
            for line in src:
                s = line.strip()
                if s.startswith(_PDBQT_STRIP_PREFIXES):
                    continue
                if s.startswith(_PDBQT_KEEP_PREFIXES):
                    lines.append(line)
                    written = True
        if not written:
            return None
        return "".join(lines)
    except Exception as e:
        logger.warning(f"[!] Could not convert {pdbqt_path} for viewer: {e}")
        return None


def _fetch_ngljs() -> Optional[bytes]:
    """Return the NGL viewer library bytes for inlining into the report.

    Resolution order (offline-first):
      1. ``/ngl`` — the directory the Docker image populates from the pnpm
         ``tools/ngl`` package during build (contains NGL's dist bundle,
         e.g. ``ngl.umd.js``). This is the preferred, always-available copy.
      2. The repo checkout's pnpm install at
         ``tools/ngl/node_modules/ngl/dist`` — so local (non-Docker) runs
         work offline too, with no re-download.
      3. A download from the CDN, only if no local copy exists.

    Returns None only when none are available; the caller then omits the
    viewer and shows a clear notice instead.
    """
    # Three levels up: vspipeline/report.py -> autodocker/ -> repo root
    repo_root = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    candidates = [
        # Prefer the standalone bundle (ngl.js) which inlines three.js etc.
        # and defines a global NGL. ngl.umd.js needs external peer globals
        # and must NOT be inlined on its own.
        "/ngl/ngl.js",
        "/ngl/ngl.min.js",
        "/ngl/ngl.umd.js",
        os.path.join(repo_root, "tools", "ngl", "node_modules", "ngl",
                     "dist", "ngl.js"),
        os.path.join(repo_root, "tools", "ngl", "node_modules", "ngl",
                     "dist", "ngl.min.js"),
        os.path.join(repo_root, "tools", "ngl", "node_modules", "ngl",
                     "dist", "ngl.umd.js"),
    ]
    for path in candidates:
        if os.path.exists(path) and os.path.getsize(path) > 0:
            try:
                with open(path, "rb") as f:
                    data = f.read()
                if data:
                    logger.info(f"[✔] Using local NGL viewer ({path}, "
                                f"{len(data) // 1024} KB) — no network needed")
                    return data
            except OSError as e:
                logger.debug(f"Could not read NGL from {path}: {e}")

    for attempt in (1, 2):
        try:
            return runner._http_get_bytes(NGL_CDN_URL, timeout=60)
        except Exception as e:
            logger.debug(f"NGL download attempt {attempt} failed: {e}")
            time.sleep(1)
    logger.warning(
        "[!] Could not bundle NGL viewer (no local copy and no network). "
        "The report will show a notice instead of the interactive viewer.")
    return None


def _build_viewer_section(ranked_ligands: List[Tuple], poses_dir: str,
                          receptor_pdbqt: Optional[str], output_file=None) -> str:
    """Create the interactive NGL 3D viewer HTML section.

    Everything needed for the viewer (the NGL library and the receptor/pose
    PDB texts) is inlined directly into the report HTML. The result is a
    single self-contained file that opens from ``file://`` with a double-click
    — no ``python3 -m http.server`` and no ``viewer/`` directory required.

    When the NGL library cannot be bundled (no network at report time) the
    viewer is omitted entirely and replaced with a clear notice. A solid
    message is preferred over a half-working report that silently needs
    internet. An empty string is returned when there are no viewable hits.
    """
    top = ranked_ligands[:10]
    if not top:
        return ""

    # 1) Bundled NGL library (inline <script>). If we cannot bundle it,
    #    drop the viewer with a clear notice — no CDN fallback ambiguity.
    # Resolve via runner so test monkeypatches are honored.
    ngl_js = runner._fetch_ngljs()
    if ngl_js is None:
        return """
    <div class="section">
        <h2>🖥️ 3D Structure Viewer</h2>
        <p style="color:#e67e22;">
            ⚠️ 3D viewer unavailable — the NGL viewer library could not be
            downloaded (no internet during report generation). The rest of
            this report is unaffected. Regenerate the report online to get
            the interactive viewer.
        </p>
    </div>
"""

    # 2) Receptor PDB text (inlined, not written to disk).
    receptor_text = None
    if receptor_pdbqt and os.path.exists(receptor_pdbqt):
        receptor_text = _pdbqt_to_pdb_string(receptor_pdbqt)

    # 3) Best pose PDB text per top hit (inlined).
    hits = []
    for lig, affinity in top:
        name = str(lig.get('Ligand', lig) if isinstance(lig, dict) else lig)
        if not os.path.isdir(poses_dir):
            break
        # Docked output is <name>_out.pdbqt (multi-model). Prefer the
        # first model as the representative (best) pose.
        candidates = sorted(glob.glob(os.path.join(
            poses_dir, f"{name}_out.pdbqt")))
        if not candidates:
            candidates = sorted(glob.glob(os.path.join(
                poses_dir, f"{name}_pose_*.pdbqt")))
        if not candidates:
            continue
        pose_text = _pdbqt_to_pdb_string(candidates[0])
        if not pose_text:
            continue
        hits.append({
            "name": name,
            "affinity": affinity,
            "pose": pose_text,
        })

    if not hits:
        return """
    <div class="section">
        <h2>🖥️ 3D Structure Viewer</h2>
        <p style="color:#e67e22;">
            No docked structures were available to display in the 3D viewer.
        </p>
    </div>
"""

    # Inline all data as one JSON blob; escape closing-script tags so the
    # PDB text cannot break out of the <script> context.
    data = {
        "receptor": receptor_text or "",
        "hits": hits,
    }
    data_json = json.dumps(data).replace("</", "<\\/")
    ngl_script = ngl_js.decode("utf-8", errors="replace")

    section = f"""
    <div class="section">
        <h2>🖥️ 3D Structure Viewer</h2>
        <p><em>Interactive NGL viewer — receptor (cartoon) with the best docked pose of each top hit. This report is fully self-contained and opens by double-clicking (no web server required).</em></p>
        <div style="margin-bottom:10px;">
            <label for="viewer-select" style="font-weight:bold;">View hit:</label>
            <select id="viewer-select" onchange="showHit()">
            </select>
            <button type="button" onclick="toggleSpin()" style="margin-left:8px;">Spin</button>
            <button type="button" onclick="resetView()" style="margin-left:8px;">Reset view</button>
        </div>
        <div id="viewer-message" style="margin-bottom:8px;color:#e67e22;"></div>
        <div id="viewer-container" style="width:100%;height:520px;background:#0b1622;border-radius:5px;"></div>
    </div>

    <div class="section" id="viewer-log" style="display:none;">
        <h2>🔍 Viewer Diagnostics</h2>
        <pre id="viewer-error"></pre>
    </div>

    <script>{ngl_script}</script>
    <script>
    (function () {{
        var data = {data_json};
        var receptorText = data.receptor || "";
        var hits = data.hits || [];
        var stage = null;

        function showHit() {{
            var sel = document.getElementById("viewer-select");
            var name = sel.value;
            var msg = document.getElementById("viewer-message");
            var selected = null;
            for (var i = 0; i < hits.length; i++) {{
                if (hits[i].name === name) selected = hits[i];
            }}
            if (!selected) return;
            if (typeof NGL === "undefined") {{
                msg.textContent = "NGL viewer failed to load.";
                return;
            }}
            if (stage) {{ stage.dispose(); }}
            var holder = document.getElementById("viewer-container");
            stage = new NGL.Stage("viewer-container", {{ backgroundColor: "#0b1622" }});
            msg.textContent = "Loading " + name + "…";
            var jobs = [];
            if (receptorText) {{
                jobs.push(stage.loadFile(
                    new Blob([receptorText], {{ type: "text/plain" }}),
                    {{ ext: "pdb" }}).then(function (o) {{
                        o.addRepresentation("cartoon", {{ colorScheme: "chainid" }});
                }}));
            }}
            jobs.push(stage.loadFile(
                new Blob([selected.pose], {{ type: "text/plain" }}),
                {{ ext: "pdb" }}).then(function (o) {{
                o.addRepresentation("ball+stick", {{
                    colorScheme: "element",
                    aspectRatio: 2.2,
                    multipleBond: "symmetric",
                }});
            }}));
            Promise.all(jobs).then(function () {{
                stage.autoView();
                msg.textContent = "";
                // Firefox does not present the WebGL canvas until a reflow;
                // nudge a resize on the next frame so it draws without a click.
                requestAnimationFrame(function () {{
                    if (stage && stage.viewer) stage.viewer.requestResize();
                }});
            }}).catch(function (err) {{
                msg.textContent = "Viewer error: " + err.message;
                document.getElementById("viewer-log").style.display = "";
                document.getElementById("viewer-error").textContent =
                    err.stack || String(err);
            }});
        }}

        function toggleSpin() {{
            if (stage) stage.toggleSpin();
        }}
        function resetView() {{
            if (stage) stage.autoView();
        }}

        // Keep the canvas sized to its container (and help Firefox repaint).
        window.addEventListener("resize", function () {{
            if (stage && stage.viewer) stage.viewer.requestResize();
        }});
        document.getElementById("viewer-container").addEventListener("click", function () {{
            if (stage && stage.viewer) stage.viewer.requestResize();
        }});

        window.addEventListener("load", function () {{
            var sel = document.getElementById("viewer-select");
            hits.forEach(function (h) {{
                var opt = document.createElement("option");
                opt.value = h.name;
                opt.textContent = h.name + " (" + (h.affinity != null ? h.affinity.toFixed(2) : "?") + " kcal/mol)";
                sel.appendChild(opt);
            }});
            if (hits.length) {{ sel.value = hits[0].name; showHit(); }}
        }});
    }})();
    </script>
"""
    return section


def generate_html_report(results_csv: str, poses_dir: str, output_file: str,
                         receptor_pdbqt: Optional[str] = None) -> None:
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
                <td><strong>{html.escape(str(lig.get('Ligand', 'Unknown')))}</strong></td>
                <td class="{conf_class}"><strong>{affinity:.2f}</strong></td>
                <td>{modes}</td>
                <td class="{conf_class}">{conf_text}</td>
            </tr>
"""

    html_content += """
        </table>
    </div>
"""

    # Interactive 3D viewer (NGL) for the top hits + receptor.
    html_content += _build_viewer_section(
        ranked_ligands, poses_dir, receptor_pdbqt, output_file)

    html_content += """
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
            html_content += f"<h3>{html.escape(str(ligand))}</h3>\n<ul>\n"
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
