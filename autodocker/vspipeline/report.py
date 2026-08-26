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


_cluster_cache: Dict[Tuple[str, float], Dict] = {}


def cluster_poses(poses_dir: str, rmsd_threshold: float = 2.0) -> Dict[str, List[Dict]]:
    """Cluster poses by RMSD. Returns mapping of ligands to pose clusters.

    Results are cached per (poses_dir, threshold) so the same directory isn't
    re-scanned when the report and results phases both request clustering.
    """
    cache_key = (poses_dir, rmsd_threshold)
    if cache_key in _cluster_cache:
        return _cluster_cache[cache_key]
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

    _cluster_cache[cache_key] = dict(clusters)
    return dict(clusters)


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


def _confidence(affinity: float):
    """Map an affinity to a (css_class, label) for visual confidence."""
    if affinity < -8:
        return 'good', '✅ Excellent'
    if affinity < -6:
        return 'good', '✅ Good'
    if affinity < -4:
        return 'warning', '⚠️ Weak'
    return 'bad', '❌ Poor'


def _stat_card(value, label, cls=""):
    return (f'<div class="stat-box {cls}"><div class="stat-value">{value}'
            f'</div><div class="stat-label">{label}</div></div>')


def _svg_histogram(values, bins=12, width=720, height=240):
    """Interactive SVG histogram of affinity values (kcal/mol)."""
    if not values:
        return '<p class="muted">No valid affinities to chart.</p>'
    mn, mx = min(values), max(values)
    if mn == mx:
        mn -= 0.5
        mx += 0.5
    bw = (mx - mn) / bins
    counts = [0] * bins
    for v in values:
        i = min(int((v - mn) / bw), bins - 1)
        counts[i] += 1
    cmax = max(counts) or 1
    pad_l, pad_b = 42, 28
    plot_w = width - pad_l - 12
    plot_h = height - pad_b - 12
    bw_px = plot_w / bins
    parts = [('<svg class="chart-svg" viewBox="0 0 %d %d" role="img" '
             'aria-label="Affinity histogram" '
             'xmlns="http://www.w3.org/2000/svg">' % (width, height))]
    parts.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="currentColor"/>'
                 % (pad_l, height - pad_b, width - 12, height - pad_b))
    parts.append('<line x1="%d" y1="12" x2="%d" y2="%d" stroke="currentColor"/>'
                 % (pad_l, pad_l, height - pad_b))
    for i, c in enumerate(counts):
        x = pad_l + i * bw_px
        h = (c / cmax) * plot_h
        y = height - pad_b - h
        lo = mn + i * bw
        hi = lo + bw
        parts.append(
            '<rect class="chart-bar" x="%.1f" y="%.1f" width="%.1f" '
            'height="%.1f" data-tip="Affinity %.1f to %.1f kcal/mol: %d ligand(s)" '
            'fill="var(--accent)" rx="1"></rect>'
            % (x, y, bw_px - 1, h, lo, hi, c))
    for frac in (0.0, 0.5, 1.0):
        ax = pad_l + frac * plot_w
        av = mn + frac * (mx - mn)
        parts.append('<text class="chart-label" x="%.1f" y="%d" text-anchor="middle">%.1f</text>'
                     % (ax, height - pad_b + 16, av))
    parts.append('</svg>')
    return "".join(parts)


def _scatter_chart(data, width=720, height=260):
    """data: list of {name, mw, aff}. Returns a clickable SVG scatter plot."""
    if not data:
        return '<p class="muted">No ligand properties to chart.</p>'
    mws = [d["mw"] for d in data]
    affs = [d["aff"] for d in data]
    mw_min, mw_max = min(mws), max(mws)
    aff_min, aff_max = min(affs), max(affs)
    if mw_min == mw_max:
        mw_min -= 1.0
        mw_max += 1.0
    if aff_min == aff_max:
        aff_min -= 0.5
        aff_max += 0.5
    pad_l, pad_b = 52, 34
    plot_w = width - pad_l - 16
    plot_h = height - pad_b - 16
    parts = [('<svg class="chart-svg" viewBox="0 0 %d %d" role="img" '
             'aria-label="Molecular weight versus affinity" '
             'xmlns="http://www.w3.org/2000/svg">' % (width, height))]
    parts.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="currentColor"/>'
                 % (pad_l, height - pad_b, width - 16, height - pad_b))
    parts.append('<line x1="%d" y1="16" x2="%d" y2="%d" stroke="currentColor"/>'
                 % (pad_l, pad_l, height - pad_b))
    parts.append('<text class="chart-label" x="%d" y="%d" text-anchor="middle">MW (Da)</text>'
                 % (width // 2, height - 4))
    parts.append('<text class="chart-label" x="14" y="%d" transform="rotate(-90 14 %d)" text-anchor="middle">Affinity</text>'
                 % (height // 2, height // 2))
    for d in data:
        cx = pad_l + (d["mw"] - mw_min) / (mw_max - mw_min) * plot_w
        cy = height - pad_b - (d["aff"] - aff_min) / (aff_max - aff_min) * plot_h
        tip = "%s — MW %.1f Da, %.2f kcal/mol" % (d["name"], d["mw"], d["aff"])
        parts.append(
            '<circle class="scatter-pt" cx="%.1f" cy="%.1f" r="5" '
            'data-name="%s" data-tip="%s" fill="var(--accent)" '
            'style="cursor:pointer"></circle>'
            % (cx, cy, html.escape(str(d["name"])), html.escape(tip)))
    parts.append('</svg>')
    return "".join(parts)


def _REPORT_CSS(has_viewer: bool = False) -> str:
    base = """
:root { --bg:#f5f7fa; --card:#fff; --text:#1f2d3d; --muted:#7f8c8d;
  --head-bg:#2c3e50; --head-fg:#fff; --th-bg:#34495e; --th-fg:#fff;
  --row-alt:#f4f6f8; --border:#dfe6e9; --accent:#27ae60; --code-bg:#ecf0f1; }
html[data-theme="dark"] { --bg:#10151b; --card:#1a222c; --text:#e6edf3;
  --muted:#8b98a5; --head-bg:#16202b; --head-fg:#e6edf3; --th-bg:#223041;
  --th-fg:#e6edf3; --row-alt:#141b23; --border:#2c3947; --accent:#2ecc71;
  --code-bg:#0d141c; }
* { box-sizing:border-box; }
body { font-family:-apple-system,Segoe UI,Arial,sans-serif; margin:0;
  background:var(--bg); color:var(--text); line-height:1.55; }
.layout { display:flex; align-items:flex-start; gap:0; max-width:1280px; margin:0 auto; }
.toc { position:sticky; top:0; align-self:flex-start; width:210px; flex:0 0 210px;
  padding:20px 14px; height:100vh; overflow-y:auto; border-right:1px solid var(--border); }
.toc h3 { margin:0 0 8px; font-size:13px; text-transform:uppercase; letter-spacing:.06em; color:var(--muted); }
.toc ul { list-style:none; margin:0; padding:0; }
.toc li { margin:2px 0; }
.toc a { display:block; padding:6px 10px; border-radius:6px; color:var(--text); text-decoration:none; font-size:14px; }
.toc a:hover { background:var(--row-alt); }
main { flex:1 1 auto; min-width:0; padding:18px 22px 60px; }
.header { background:linear-gradient(135deg,var(--head-bg),#1a2530); color:var(--head-fg); padding:24px 22px; border-radius:8px; margin-bottom:16px; }
.header h1 { margin:0 0 4px; font-size:22px; }
.header p { margin:0; opacity:.85; }
.section { background:var(--card); padding:18px 20px; margin:16px 0; border-radius:8px; box-shadow:0 2px 6px rgba(0,0,0,.07); scroll-margin-top:14px; }
.section h2 { margin-top:0; border-bottom:2px solid var(--accent); padding-bottom:6px; font-size:18px; }
table { width:100%; border-collapse:collapse; margin:10px 0; }
th { background:var(--th-bg); color:var(--th-fg); padding:9px; text-align:left; }
td { padding:8px 9px; border-bottom:1px solid var(--border); }
tr:hover { background:var(--row-alt); }
tr.top-row td:first-child::before { content:"\\2605 "; color:#f39c12; }
.stat-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:12px; }
.stat-box { background:var(--row-alt); border:1px solid var(--border); padding:14px; border-radius:8px; text-align:center; }
.stat-value { font-size:22px; font-weight:700; color:var(--accent); }
.stat-label { font-size:12px; color:var(--muted); margin-top:5px; }
.good { color:#27ae60; } .warning { color:#f39c12; } .bad { color:#e74c3c; }
.muted { color:var(--muted); }
.meta-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(210px,1fr)); gap:10px; }
.meta-item { background:var(--row-alt); border:1px solid var(--border); border-radius:6px; padding:9px 12px; }
.meta-label { display:block; font-size:11px; text-transform:uppercase; letter-spacing:.05em; color:var(--muted); }
.meta-value { font-weight:700; }
.chart-svg { width:100%; max-width:760px; height:auto; color:var(--text); }
.chart-bar:hover { opacity:.75; cursor:pointer; }
.chart-label { font-size:12px; fill:currentColor; }
.chart-tip { position:fixed; display:none; pointer-events:none; background:#222; color:#fff; padding:5px 9px; border-radius:5px; font-size:12px; z-index:50; box-shadow:0 2px 8px rgba(0,0,0,.3); max-width:320px; }
.float-btn { position:fixed; top:12px; z-index:20; padding:6px 12px; border:1px solid var(--border); border-radius:6px; background:var(--card); color:var(--text); cursor:pointer; box-shadow:0 1px 4px rgba(0,0,0,.18); }
#theme-toggle { right:12px; } #print-btn { right:110px; }
.viewer-help { font-size:13px; color:var(--muted); }
.viewer-controls { margin:0 0 10px; display:flex; flex-wrap:wrap; gap:8px; align-items:center; }
.viewer-controls select, .viewer-controls button { padding:6px 10px; border:1px solid var(--border); border-radius:5px; background:var(--card); color:var(--text); cursor:pointer; }
#viewer-message { margin:0 0 8px; color:#e67e22; min-height:18px; }
.viewer-canvas { touch-action:none; width:100%; height:520px; background:#0b1622; border-radius:6px; }
@media (max-width:880px) { .layout { flex-direction:column; }
  .toc { position:static; width:100%; height:auto; flex:none; border-right:none; border-bottom:1px solid var(--border); }
  .toc ul { display:flex; flex-wrap:wrap; gap:4px; } main { padding:14px; }
  #print-btn { right:12px; } #theme-toggle { right:100px; } }
@media print { body { background:#fff; color:#000; margin:0; }
  .layout { display:block; max-width:none; }
  .toc, .float-btn, .no-print { display:none !important; }
  .section { box-shadow:none; break-inside:avoid-page; }
  a[href^="data:"]::after { content:" (ranking.csv embedded)"; } }
"""
    return base


_THEME_JS = """
(function(){ var t=document.getElementById('theme-toggle');
  var cur=localStorage.getItem('vs-theme');
  if(cur){ document.documentElement.setAttribute('data-theme',cur); }
  if(t){ t.onclick=function(){ var n=document.documentElement.getAttribute('data-theme')==='dark'?'light':'dark';
    document.documentElement.setAttribute('data-theme',n); localStorage.setItem('vs-theme',n); }; }
})();
"""
_TABLE_JS = """
(function(){ var s=document.getElementById('results-search');
  if(s){ s.oninput=function(){ var q=s.value.toLowerCase();
    document.querySelectorAll('#all-results tbody tr').forEach(function(tr){
      tr.style.display = tr.textContent.toLowerCase().indexOf(q)>=0 ? '' : 'none'; }); }; }
  document.querySelectorAll('#all-results th.sortable').forEach(function(th){
    th.onclick=function(){ var tb=th.closest('table').querySelector('tbody');
      var i=Array.prototype.indexOf.call(th.parentNode.children,th);
      var asc=!th.classList.contains('sort-asc');
      th.parentNode.querySelectorAll('th').forEach(function(c){c.classList.remove('sort-asc','sort-desc');});
      th.classList.add(asc?'sort-asc':'sort-desc');
      Array.prototype.slice.call(tb.children).sort(function(a,b){
        var x=a.children[i].textContent.trim(), y=b.children[i].textContent.trim();
        var nx=parseFloat(x), ny=parseFloat(y);
        if(!isNaN(nx)&&!isNaN(ny)) return asc?nx-ny:ny-nx;
        return asc? x.localeCompare(y): y.localeCompare(x); }).forEach(function(r){tb.appendChild(r);});
    }; });
})();
"""
_CHART_JS = """
(function(){ var tip=document.getElementById('chart-tip');
  if(!tip){ tip=document.createElement('div'); tip.id='chart-tip'; tip.className='chart-tip';
    document.body.appendChild(tip); }
  function bind(){ document.querySelectorAll('[data-tip]').forEach(function(el){
    el.addEventListener('mousemove',function(e){ tip.style.display='block';
      tip.textContent=el.getAttribute('data-tip'); tip.style.left=(e.clientX+12)+'px';
      tip.style.top=(e.clientY+12)+'px'; });
    el.addEventListener('mouseleave',function(){ tip.style.display='none'; }); }); }
  bind();
  window.addEventListener('load',function(){ if(window.selectHit){
    document.querySelectorAll('.scatter-pt').forEach(function(pt){
      pt.addEventListener('click',function(){ window.selectHit(pt.getAttribute('data-name')); }); }); } });
})();
"""


def _build_viewer_section(ranked_ligands: List[Tuple], poses_dir: str,
                          receptor_pdbqt: Optional[str] = None,
                          output_file=None, grid_box=None) -> str:
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
        "gridBox": grid_box,
    }
    data_json = json.dumps(data).replace("</", "<\\/")
    ngl_script = ngl_js.decode("utf-8", errors="replace")

    section = f"""
    <div class="section" id="viewer">
        <h2>🖥️ 3D Structure Viewer</h2>
        <p class="viewer-help">Receptor (cartoon) with the best docked pose of each top hit, plus the docking grid box when available. Fully self-contained — opens by double-clicking, no server needed. Click a point in the scatter plot to jump to that ligand.</p>
        <div class="viewer-controls no-print">
            <label for="viewer-select"><strong>View hit:</strong></label>
            <select id="viewer-select" onchange="vsShowHit()"></select>
            <button type="button" onclick="vsSpin()">Spin</button>
            <button type="button" onclick="vsReset()">Reset</button>
            <button type="button" onclick="vsOverlayAll()">Overlay all top hits</button>
            <button type="button" onclick="vsSnapshot()">📷 Snapshot</button>
            <label><input type="checkbox" id="vs-receptor" checked onchange="vsToggle('receptor')"> Receptor</label>
            <label><input type="checkbox" id="vs-ligand" checked onchange="vsToggle('ligand')"> Ligand</label>
            <label><input type="checkbox" id="vs-box" {'checked' if grid_box else ''} onchange="vsToggle('box')"> Grid box</label>
            <label><input type="checkbox" id="vs-labels" onchange="vsToggle('labels')"> Labels</label>
            <select id="vs-style" onchange="vsStyle()">
                <option value="ball+stick">Ball &amp; Stick</option>
                <option value="licorice">Licorice</option>
                <option value="hyperball">Hyperball</option>
                <option value="surface">Surface</option>
            </select>
        </div>
        <div id="viewer-message"></div>
        <div id="viewer-container" class="viewer-canvas"></div>
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
        var gridBox = data.gridBox || null;
        var stage = null;
        var components = {{ receptor: null, ligand: null, box: null }};
        var currentStyle = "ball+stick";

        function vsShowHit() {{
            var sel = document.getElementById("viewer-select");
            var name = sel.value;
            var msg = document.getElementById("viewer-message");
            var selected = null;
            for (var i = 0; i < hits.length; i++) {{ if (hits[i].name === name) selected = hits[i]; }}
            if (!selected) return;
            if (typeof NGL === "undefined") {{ msg.textContent = "NGL viewer failed to load."; return; }}
            if (stage) {{ stage.dispose(); stage = null; components = {{ receptor:null, ligand:null, box:null }}; }}
            var holder = document.getElementById("viewer-container");
            stage = new NGL.Stage("viewer-container", {{ backgroundColor: "#0b1622" }});
            msg.textContent = "Loading " + name + "…";
            var jobs = [];
            if (receptorText) {{
                jobs.push(stage.loadFile(new Blob([receptorText], {{ type: "text/plain" }}), {{ ext: "pdb" }}).then(function (o) {{
                    components.receptor = o;
                    o.addRepresentation("cartoon", {{ colorScheme: "chainid", name: "receptor" }});
                }}));
            }}
            jobs.push(stage.loadFile(new Blob([selected.pose], {{ type: "text/plain" }}), {{ ext: "pdb" }}).then(function (o) {{
                components.ligand = o;
                o.addRepresentation(currentStyle, {{ colorScheme: "element", aspectRatio: 2.2, multipleBond: "symmetric", name: "ligand" }});
            }}));
            if (gridBox && gridBox.center && gridBox.size && document.getElementById("vs-box").checked) {{
                try {{
                    var c = gridBox.center, s = gridBox.size, hx=s[0]/2, hy=s[1]/2, hz=s[2]/2;
                    var shape = new NGL.Shape("Grid Box");
                    shape.addBox([c[0]-hx,c[1]-hy,c[2]-hz], [c[0]+hx,c[1]+hy,c[2]+hz], [0.18,0.8,0.44]);
                    components.box = stage.addComponentFromObject(shape);
                }} catch (e) {{ /* box is optional */ }}
            }}
            Promise.all(jobs).then(function () {{
                stage.autoView();
                msg.textContent = "";
                requestAnimationFrame(function () {{ if (stage && stage.viewer) stage.viewer.requestResize(); }});
            }}).catch(function (err) {{
                msg.textContent = "Viewer error: " + err.message;
                document.getElementById("viewer-log").style.display = "";
                document.getElementById("viewer-error").textContent = err.stack || String(err);
            }});
        }}

        window.vsShowHit = vsShowHit;
        window.vsSpin = function () {{ if (stage) stage.toggleSpin(); }};
        window.vsReset = function () {{ if (stage) stage.autoView(); }};
        window.vsSnapshot = function () {{
            if (!stage) return;
            stage.viewer.render(function () {{
                var a = document.createElement("a");
                a.href = stage.viewer.canvas.toDataURL("image/png");
                a.download = "docking_view.png"; a.click();
            }});
        }};
        window.vsToggle = function (which) {{
            if (!stage) return;
            var show = document.getElementById("vs-" + which).checked;
            if (which === "labels") {{
                if (components.ligand) {{ show ? components.ligand.addRepresentation("label", {{ labelType:"atomname", name:"lbl" }}) : components.ligand.removeRepresentation("lbl"); }}
                return;
            }}
            if (which === "receptor" && components.receptor) components.receptor.setVisibility(show);
            if (which === "ligand" && components.ligand) components.ligand.setVisibility(show);
            if (which === "box" && components.box) components.box.setVisibility(show);
        }};
        window.vsStyle = function () {{
            currentStyle = document.getElementById("vs-style").value;
            if (components.ligand) {{
                components.ligand.removeRepresentation("ligand");
                components.ligand.addRepresentation(currentStyle, {{ colorScheme:"element", aspectRatio:2.2, multipleBond:"symmetric", name:"ligand" }});
            }}
        }};
        window.vsOverlayAll = function () {{
            if (!stage || !hits.length) return;
            for (var i = 0; i < hits.length; i++) {{
                (function (h) {{
                    stage.loadFile(new Blob([h.pose], {{ type: "text/plain" }}), {{ ext: "pdb" }}).then(function (o) {{
                        o.addRepresentation(currentStyle, {{ colorScheme: "chainname", aspectRatio: 2.2, name: "ov" + i }});
                    }});
                }})(hits[i]);
            }}
            stage.autoView();
        }};
        window.selectHit = function (name) {{
            var sel = document.getElementById("viewer-select");
            if (!sel) return;
            for (var i = 0; i < sel.options.length; i++) {{ if (sel.options[i].value === name) {{ sel.selectedIndex = i; break; }} }}
            vsShowHit();
            document.getElementById("viewer-container").scrollIntoView({{ behavior: "smooth", block: "center" }});
        }};

        window.addEventListener("resize", function () {{ if (stage && stage.viewer) stage.viewer.requestResize(); }});
        document.getElementById("viewer-container").addEventListener("click", function () {{ if (stage && stage.viewer) stage.viewer.requestResize(); }});
        window.addEventListener("load", function () {{
            var sel = document.getElementById("viewer-select");
            hits.forEach(function (h) {{
                var opt = document.createElement("option");
                opt.value = h.name;
                opt.textContent = h.name + " (" + (h.affinity != null ? h.affinity.toFixed(2) : "?") + " kcal/mol)";
                sel.appendChild(opt);
            }});
            if (hits.length) {{ sel.value = hits[0].name; vsShowHit(); }}
        }});
    }})();
    </script>
"""
    return section


def generate_html_report(results_csv: str, poses_dir: str, output_file: str,
                          receptor_pdbqt: Optional[str] = None,
                          meta: Optional[Dict] = None,
                          grid_file: Optional[str] = None,
                          scorer_agreement: Optional[Dict] = None) -> None:
    """Generate a self-contained, professional HTML screening report.

    Includes a table of contents, run metadata, summary stat cards, an
    interactive NGL 3D viewer (with docking-box overlay when ``grid_file`` is
    supplied), an affinity histogram, a molecular-weight vs affinity scatter
    plot, a sortable/filterable results table, pose clustering and an
    interpretation guide. The render is fully offline (NGL + PDB data inlined).
    """
    logger.info("[*] Generating HTML report...")

    ligands = []
    try:
        with open(results_csv, 'r') as f:
            ligands = list(csv.DictReader(f))
    except Exception as e:
        logger.warning(f"Could not parse results CSV: {e}")
        return
    if not ligands:
        logger.warning("No results to report")
        return

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
    avg_affinity = statistics.mean(valid_affinities) if valid_affinities else 0.0
    best_affinity = min(valid_affinities) if valid_affinities else 0.0
    worst_affinity = max(valid_affinities) if valid_affinities else 0.0
    strong_binders = sum(1 for a in valid_affinities if a < -6.0)
    success_rate = (len(valid_affinities) / len(ligands) * 100) if ligands else 0.0

    # Ranked ligands (negative affinity = valid dock)
    ranked_ligands = []
    for lig in ligands:
        try:
            av = float(lig.get('Binding_Affinity', '').strip())
        except ValueError:
            continue
        if av < 0:
            ranked_ligands.append((lig, av))
    ranked_ligands.sort(key=lambda item: item[1])

    # Docking grid box for the viewer
    grid_box = None
    if grid_file and os.path.exists(grid_file):
        try:
            gc = runner._parse_grid_config(grid_file)
            if gc and gc.get('center') and gc.get('size'):
                grid_box = {"center": list(gc['center']), "size": list(gc['size'])}
        except Exception as e:
            logger.debug(f"Could not parse grid config for viewer: {e}")

    # Scatter data (MW vs affinity)
    scatter_data = []
    for lig in ligands:
        try:
            mw = float(lig.get('MW', '').strip())
            av = float(lig.get('Binding_Affinity', '').strip())
        except ValueError:
            continue
        if av < 0 and mw > 0:
            scatter_data.append({"name": lig.get('Ligand', '?'),
                                 "mw": mw, "aff": av})

    clusters = cluster_poses(poses_dir) if os.path.isdir(poses_dir) else {}

    css = _REPORT_CSS(bool(ranked_ligands))
    generated = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    version = getattr(runner, "__version__", "2.0")

    meta_items = ""
    if meta:
        for k, v in meta.items():
            meta_items += (f'<div class="meta-item"><span class="meta-label">'
                           f'{html.escape(str(k))}</span>'
                           f'<span class="meta-value">{html.escape(str(v))}</span></div>')

    # Top hits rows
    top_rows = ""
    top_names = set()
    for i, (lig, affinity) in enumerate(ranked_ligands[:10], 1):
        top_names.add(str(lig.get('Ligand', 'Unknown')))
        conf_class, conf_text = _confidence(affinity)
        top_rows += f"""
            <tr class="top-row">
                <td>{i}</td>
                <td><strong>{html.escape(str(lig.get('Ligand', 'Unknown')))}</strong></td>
                <td class="{conf_class}"><strong>{affinity:.2f}</strong></td>
                <td>{html.escape(str(lig.get('Binding_Modes', '0')))}</td>
                <td class="{conf_class}">{conf_text}</td>
            </tr>"""

    # All results rows
    all_rows = ""
    for lig in ligands:
        raw = lig.get('Binding_Affinity', '').strip()
        try:
            av = float(raw)
        except ValueError:
            av = float('nan')
        if av == av and av < 0:
            conf_class, conf_text = _confidence(av)
        else:
            conf_class, conf_text = "bad", "❌ N/A"
        cls = "top-row" if str(lig.get('Ligand', '')) in top_names else ""
        all_rows += f"""
            <tr class="{cls}">
                <td>{html.escape(str(lig.get('Ligand', 'Unknown')))}</td>
                <td class="{conf_class}">{html.escape(raw)}</td>
                <td>{html.escape(str(lig.get('Binding_Modes', '-')))}</td>
                <td>{html.escape(str(lig.get('SimScore', '-')))}</td>
                <td>{html.escape(str(lig.get('MW', '-')))}</td>
                <td>{html.escape(str(lig.get('LogP', '-')))}</td>
                <td>{html.escape(str(lig.get('Status', '-')))}</td>
            </tr>"""

    # Global, citable scorer-concordance metric (Spearman rho) -- shown as the
    # defensible agreement statistic. The per-ligand "Agreement (heuristic)"
    # column in consensus_ranking.csv is only a quick UI indicator.
    agreement_section = ""
    if scorer_agreement and scorer_agreement.get("rho") is not None:
        _rho = float(scorer_agreement["rho"])
        _n = int(scorer_agreement.get("n", 0))
        _abs = abs(_rho)
        if _abs >= 0.7:
            _strength = "strong"
        elif _abs >= 0.4:
            _strength = "moderate"
        elif _abs >= 0.2:
            _strength = "weak"
        else:
            _strength = "negligible"
        _direction = ("concordant" if _rho > 0 else
                      "discordant" if _rho < 0 else "uncorrelated")
        agreement_section = f"""
            <div class="section" id="agreement">
                <h2>🤝 Scorer Agreement (Vina vs SMINA)</h2>
                <p class="muted"><em>Global, citable measure of scorer concordance
                across the screen (Spearman rank correlation).</em></p>
                <div class="stat-grid">
                    {_stat_card(f"{_rho:.3f}", "Spearman &rho;")}
                    {_stat_card(_n, "Paired ligands (n)")}
                </div>
                <p>The two scoring functions are <strong>{_strength}
                {_direction}</strong> (|&rho;| = {_abs:.3f}). This &rho; is the
                statistically grounded agreement metric and may be reported as
                such.</p>
                <p class="muted"><small>Note: the per-ligand
                <code>Agreement (heuristic)</code> column in
                <code>consensus_ranking.csv</code> is only a quick UI indicator
                (max(0, 100 &minus; |&Delta;|&times;10)) and is
                <strong>not</strong> a validated metric &mdash; do not cite it
                quantitatively.</small></p>
            </div>
"""

    html_content = f"""<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Virtual Screening Results Report</title>
    <style>{css}</style>
</head>
<body>
    <button id="theme-toggle" class="float-btn" type="button">🌓 Theme</button>
    <button id="print-btn" class="float-btn no-print" type="button" onclick="window.print()">🖨️ PDF</button>
    <div class="layout">
        <nav class="toc">
            <h3>Contents</h3>
            <ul>
                <li><a href="#summary">Summary</a></li>
                <li><a href="#top">Top Hits</a></li>
                <li><a href="#viewer">3D Viewer</a></li>
                <li><a href="#dist">Affinity Distribution</a></li>
                <li><a href="#scatter">MW vs Affinity</a></li>
                <li><a href="#all">All Results</a></li>
                <li><a href="#clusters">Pose Clustering</a></li>
                <li><a href="#guide">Interpretation</a></li>
            </ul>
        </nav>
        <main>
            <div class="header">
                <h1>🧬 Virtual Screening Results Report</h1>
                <p>Generated: {generated} · AutoDock Toolkit v{version}</p>
            </div>

            <div class="section" id="summary">
                <h2>📊 Run Metadata</h2>
                <div class="meta-grid">{meta_items or '<p class="muted">No run metadata supplied.</p>'}</div>
                <h2 style="margin-top:18px;">Key Statistics</h2>
                <div class="stat-grid">
                    {_stat_card(len(ligands), "Ligands Screened")}
                    {_stat_card(f"{best_affinity:.2f}", "Best Binder (kcal/mol)")}
                    {_stat_card(f"{avg_affinity:.2f}", "Average Affinity")}
                    {_stat_card(f"{strong_binders}/{len(valid_affinities)}", "Strong Binders (&lt; -6.0)")}
                    {_stat_card(f"{success_rate:.1f}%", "Success Rate")}
                    {_stat_card(failed_count, "Failed / Invalid")}
                </div>
            </div>

            {agreement_section}

            <div class="section" id="top">
                <h2>🏆 Top 10 Hits</h2>
                <table>
                    <tr><th>#</th><th>Ligand Name</th><th>Affinity (kcal/mol)</th><th>Poses</th><th>Confidence</th></tr>
                    {top_rows}
                </table>
            </div>

            {_build_viewer_section(ranked_ligands, poses_dir, receptor_pdbqt, output_file, grid_box=grid_box)}

            <div class="section" id="dist">
                <h2>📈 Affinity Distribution</h2>
                <p class="muted"><em>Histogram of binding affinities across all docked ligands (kcal/mol).</em></p>
                {_svg_histogram(valid_affinities)}
            </div>

            <div class="section" id="scatter">
                <h2>🔬 Molecular Weight vs Affinity</h2>
                <p class="muted"><em>Each point is a ligand; click a point to load it in the 3D viewer.</em></p>
                {_scatter_chart(scatter_data)}
            </div>

            <div class="section" id="all">
                <h2>📋 All Results</h2>
                <input id="results-search" type="text" placeholder="Filter ligands…" class="no-print">
                <table id="all-results">
                    <tr>
                        <th class="sortable">Ligand</th>
                        <th class="sortable">Affinity</th>
                        <th class="sortable">Poses</th>
                        <th class="sortable">SimScore</th>
                        <th class="sortable">MW</th>
                        <th class="sortable">LogP</th>
                        <th class="sortable">Status</th>
                    </tr>
                    {all_rows}
                </table>
            </div>

            <div class="section" id="clusters">
                <h2>🔬 Pose Clustering Analysis</h2>
                <p class="muted"><em>Grouping of similar poses by RMSD (threshold: 2.0 Å).</em></p>
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
        html_content += "<p class='muted'><em>No docked structures available for clustering.</em></p>\n"

    html_content += f"""
            </div>

            <div class="section" id="guide">
                <h2>📖 Interpretation Guide</h2>
                <ul>
                    <li><strong>&lt; -10 kcal/mol:</strong> Excellent binder — investigate immediately</li>
                    <li><strong>-10 to -8 kcal/mol:</strong> Very good binder — high priority</li>
                    <li><strong>-8 to -6 kcal/mol:</strong> Good binder — validate further</li>
                    <li><strong>-6 to -4 kcal/mol:</strong> Weak binder — may need optimization</li>
                    <li><strong>&gt; -4 kcal/mol:</strong> Negligible binding — not promising</li>
                    <li><strong>0.0 kcal/mol:</strong> Docking error — check input format</li>
                </ul>
            </div>

            <div class="section" style="text-align:center; color:var(--muted);">
                <p><small>Report generated by AutoDock Virtual Screening Toolkit v{version}</small></p>
            </div>
        </main>
    </div>
    <script>{_THEME_JS}</script>
    <script>{_TABLE_JS}</script>
    <script>{_CHART_JS}</script>
</body>
</html>
"""

    try:
        with open(output_file, 'w') as f:
            f.write(html_content)
        logger.info(f"[✔] HTML report generated: {output_file}")
    except Exception as e:
        logger.warning(f"Could not write HTML report: {e}")
