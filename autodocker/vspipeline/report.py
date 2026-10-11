#!/usr/bin/env python3
"""HTML report and pose clustering.

All cross-module references go through the ``runner`` namespace at call
time so the public entry point and its monkeypatch surface stay unchanged.
The report's CSS, JavaScript and viewer markup live in ``templates/``.

Note: the PDBQT model/REMARK parsing below (``_models_from_vina_output``)
should switch to the pipeline's shared parser in ``vspipeline.docking``
(``parse_pdbqt_results``) once that exists.
"""
import os
import re
import csv
import json
import math
import time
import hashlib
import html
import statistics
import logging
from collections import defaultdict
from datetime import datetime
from typing import List, Tuple, Optional, Dict

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
    """RMSD between two coordinate lists of the same molecule, pairing atoms
    by their order in the file. No superposition and no symmetry correction
    (equivalent atoms swapped by symmetry count as displaced), so it is a
    quick clustering distance, not a pose-accuracy metric."""
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


_cluster_cache: Dict[Tuple, Dict] = {}


def cluster_poses(poses_dir: str, rmsd_threshold: float = 2.0,
                  ligands: Optional[List[str]] = None) -> Dict[str, List[Dict]]:
    """Greedy clustering of each ligand's docked modes.

    The distance is the heavy-atom RMSD by atom order (``calculate_rmsd``:
    no superposition, not symmetry-corrected). Returns a mapping of ligand
    to pose clusters.

    ``ligands`` limits clustering to the ligands of the current screen, so
    output files left in ``docked/`` by an earlier run are ignored.
    Results are cached per (poses_dir, threshold, ligands) so the same
    directory isn't re-scanned when the report and results phases both
    request clustering.
    """
    keep = None if ligands is None else set(ligands)
    cache_key = (poses_dir, rmsd_threshold, None if keep is None else frozenset(keep))
    if cache_key in _cluster_cache:
        return _cluster_cache[cache_key]
    clusters = defaultdict(list)

    if not os.path.isdir(poses_dir):
        return clusters

    # Vina (<ligand>_out.pdbqt) and SMINA (<ligand>_smina.pdbqt) write every
    # mode of a ligand as a MODEL in one file; Vina wins if both exist.
    vina_outputs = {}
    for file in sorted(os.listdir(poses_dir)):
        for suffix in ('_out.pdbqt', '_smina.pdbqt'):
            if file.endswith(suffix):
                lig = file[:-len(suffix)]
                if suffix == '_out.pdbqt' or lig not in vina_outputs:
                    vina_outputs[lig] = file

    for ligand, out_file in vina_outputs.items():
        if keep is not None and ligand not in keep:
            continue
        pose_data = _models_from_vina_output(os.path.join(poses_dir, out_file))
        _cluster_into(clusters, ligand, pose_data, rmsd_threshold)

    _cluster_cache[cache_key] = dict(clusters)
    return dict(clusters)


_HYDROGEN_TYPES = ("H", "HD", "HS")


def _models_from_vina_output(path: str) -> List[Dict]:
    """One pose per MODEL of a Vina output file: heavy-atom coordinates
    (flexible side chains excluded) and the mode's affinity."""
    poses = []
    current, affinity, model, in_flex = [], 0.0, 0, False
    try:
        with open(path) as f:
            for line in f:
                # Flexible side chains (BEGIN_RES..END_RES) are receptor atoms.
                if line.startswith('BEGIN_RES'):
                    in_flex = True
                elif line.startswith('END_RES'):
                    in_flex = False
                elif in_flex:
                    continue
                elif line.startswith('MODEL'):
                    current, affinity = [], 0.0
                    model += 1
                elif line.startswith('REMARK VINA RESULT:'):
                    try:
                        affinity = float(line.split()[3])
                    except (IndexError, ValueError):
                        pass
                elif line.startswith(('ATOM', 'HETATM')):
                    # Heavy atoms only: the AutoDock type (last field) of
                    # hydrogens is H/HD/HS.
                    fields = line.split()
                    if fields and fields[-1] in _HYDROGEN_TYPES:
                        continue
                    try:
                        current.append((float(line[30:38]), float(line[38:46]),
                                        float(line[46:54])))
                    except ValueError:
                        pass
                elif line.startswith('ENDMDL') and current:
                    poses.append({'file': f"{os.path.basename(path)}#{model}",
                                  'coords': current, 'affinity': affinity})
    except OSError as e:
        logger.debug(f"Could not read poses from {path}: {e}")
    return poses


def _cluster_into(clusters: Dict, ligand: str, pose_data: List[Dict],
                  rmsd_threshold: float) -> None:
    """Greedy RMSD clustering of one ligand's poses (best pose first)."""
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


# NGL is inlined into the report. A downloaded copy must be exactly this
# file (the standalone bundle that defines a global ``NGL``); anything else
# is rejected and the viewer is left out of the report.
NGL_VERSION = "2.4.0"
NGL_CDN_URL = f"https://cdn.jsdelivr.net/npm/ngl@{NGL_VERSION}/dist/ngl.js"
NGL_SHA256 = "be8d4b16f1d73458efd70d2a3319ef2888fb6fa821732c319a99ff43a8e4f7e0"
TEMPLATE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")

# PDB lines that are safe/meaningful to NGL (PDBQT adds ROOT/BRANCH
# bookkeeping that NGL's PDB parser does not understand).
_PDBQT_STRIP_PREFIXES = ("ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF")
_PDBQT_KEEP_PREFIXES = ("ATOM", "HETATM", "MODEL", "ENDMDL", "REMARK",
                        "CONECT", "HEADER", "TITLE", "COMPND", "TER",
                        "MASTER", "END")


def _pdbqt_to_pdb_string(pdbqt_path: str, first_model_only: bool = False) -> Optional[str]:
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
                if first_model_only and s.startswith("ENDMDL"):
                    lines.append(line)
                    break
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


def _ngl_local_candidates() -> List[str]:
    """Local NGL bundles, in order of preference."""
    # Three levels up: vspipeline/report.py -> autodocker/ -> repo root
    repo_root = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__))))
    return ["/ngl/ngl.js",
            os.path.join(repo_root, "tools", "ngl", "node_modules", "ngl",
                         "dist", "ngl.js")]


def _fetch_ngljs() -> Optional[bytes]:
    """Return the NGL viewer library bytes for inlining into the report.

    Resolution order:
      1. ``/ngl/ngl.js``: the copy the Docker image installs at build time
         from ``tools/ngl`` (no network needed at report time).
      2. The repo checkout's ``tools/ngl/node_modules/ngl/dist/ngl.js``, for
         local (non-Docker) runs.
      3. A download of ``NGL_CDN_URL`` (needs network). The download is used
         only if its SHA-256 equals ``NGL_SHA256``; otherwise it is
         discarded.

    Only the standalone ``ngl.js`` bundle is used: ``ngl.umd.js`` needs
    external peer globals and does not work when inlined on its own. Local
    copies are logged with their SHA-256 but not rejected when they differ
    from the pinned file (they come from the image or checkout).

    Returns None when no usable copy exists; the caller then omits the
    viewer and shows a notice instead.
    """
    for path in _ngl_local_candidates():
        if os.path.exists(path) and os.path.getsize(path) > 0:
            try:
                with open(path, "rb") as f:
                    data = f.read()
            except OSError as e:
                logger.debug(f"Could not read NGL from {path}: {e}")
                continue
            digest = hashlib.sha256(data).hexdigest()
            note = "" if digest == NGL_SHA256 else f" (not the pinned {NGL_VERSION} build)"
            logger.info(f"[*] Using local NGL viewer {path} "
                        f"(sha256 {digest[:12]}{note})")
            return data

    for attempt in (1, 2):
        try:
            data = runner._http_get_bytes(NGL_CDN_URL, timeout=60)
        except Exception as e:
            logger.debug(f"NGL download attempt {attempt} failed: {e}")
            time.sleep(1)
            continue
        digest = hashlib.sha256(data).hexdigest()
        if digest != NGL_SHA256:
            logger.warning(
                f"[!] Downloaded NGL {NGL_VERSION} has SHA-256 {digest}, expected "
                f"{NGL_SHA256}; not embedding it. The report is written without "
                "the 3D viewer.")
            return None
        return data
    logger.warning(
        "[!] Could not bundle the NGL viewer (no local copy and the download "
        "failed). The report is written without the 3D viewer.")
    return None


# Heuristic score bands (kcal/mol) for colouring and the interpretation
# guide. A band applies when the score is below its upper bound. Vina scores
# are not calibrated binding free energies; these bands only sort a screen.
AFFINITY_BANDS = (
    (-8.0, "good", "Excellent"),
    (-6.0, "good", "Good"),
    (-4.0, "warning", "Weak"),
    (None, "bad", "Poor"),
)
STRONG_BINDER_CUTOFF = -6.0  # "Good" or better

# Status values of ranking.csv rows that carry no usable score.
_UNSCORED_STATUS = ("FAILED", "REJECTED", "ERROR", "OUT_OF_WINDOW")


def _confidence(affinity: float):
    """Map an affinity to a (css_class, label) from ``AFFINITY_BANDS``."""
    for upper, css, label in AFFINITY_BANDS:
        if upper is None or affinity < upper:
            return css, label
    return AFFINITY_BANDS[-1][1], AFFINITY_BANDS[-1][2]


def _interpretation_items() -> str:
    """Guide entries generated from ``AFFINITY_BANDS``."""
    items, lower = [], None
    for upper, _css, label in AFFINITY_BANDS:
        if lower is None:
            rng = f"&lt; {upper:g}"
        elif upper is None:
            rng = f"&ge; {lower:g}"
        else:
            rng = f"{lower:g} to {upper:g}"
        items.append(f"<li><strong>{rng} kcal/mol:</strong> {label}</li>")
        lower = upper
    return "\n                    ".join(items)


def _valid_score(row: Dict) -> Optional[float]:
    """The row's affinity when it is a usable docking score, else None."""
    status = str(row.get('Status', '') or '').strip().upper()
    if status.startswith(_UNSCORED_STATUS):
        return None
    try:
        value = float(str(row.get('Binding_Affinity', '')).strip())
    except ValueError:
        return None
    if value != value or value >= 0:
        return None
    return value


def _template(name: str) -> str:
    with open(os.path.join(TEMPLATE_DIR, name), encoding="utf-8") as f:
        return f.read()


def _fill(template: str, values: Dict[str, str]) -> str:
    """Replace ``__VS_*__`` placeholders in one pass, so inserted text is
    never scanned for further placeholders."""
    return re.sub(r"__VS_[A-Z_]+__", lambda m: values.get(m.group(0), m.group(0)),
                  template)


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
    """Report stylesheet (``templates/report.css``)."""
    return _template("report.css")


def _viewer_notice(message: str) -> str:
    return f"""
    <div class="section" id="viewer">
        <h2>3D Structure Viewer</h2>
        <p class="viewer-note">{message}</p>
    </div>
"""


def _pose_file(lig, name: str, poses_dir: str) -> Optional[str]:
    """Docked pose file of a ranking row: the recorded Docked_File (absolute,
    or relative to the output directory), then ``docked/<basename>``, then
    ``<name>_out.pdbqt``."""
    recorded = lig.get('Docked_File', '') if isinstance(lig, dict) else ''
    candidates = []
    if recorded:
        if os.path.isabs(recorded):
            candidates.append(recorded)
        else:
            candidates.append(os.path.join(os.path.dirname(os.path.abspath(poses_dir)),
                                           recorded))
        candidates.append(os.path.join(poses_dir, os.path.basename(recorded)))
    candidates.append(os.path.join(poses_dir, f"{name}_out.pdbqt"))
    return next((p for p in candidates if os.path.isfile(p)), None)


def _build_viewer_section(ranked_ligands: List[Tuple], poses_dir: str,
                          receptor_pdbqt: Optional[str] = None,
                          output_file=None, grid_box=None) -> str:
    """Create the interactive NGL 3D viewer HTML section.

    The NGL library and the receptor/pose PDB texts are inlined into the
    report, so the written file opens from ``file://`` without a server or
    network access. Producing it may need the network once (see
    ``_fetch_ngljs``); when no verified NGL copy is available the viewer is
    omitted and replaced by a notice. Returns "" when there are no hits.
    """
    top = ranked_ligands[:10]
    if not top:
        return ""

    # Resolve via runner so test monkeypatches are honored.
    ngl_js = runner._fetch_ngljs()
    if ngl_js is None:
        return _viewer_notice(
            "3D viewer unavailable: no verified copy of the NGL viewer library "
            f"(version {NGL_VERSION}) was available when this report was "
            "generated (no local copy, and the download failed or did not "
            "match the pinned checksum). The rest of this report is unaffected.")

    receptor_text = None
    if receptor_pdbqt and os.path.exists(receptor_pdbqt):
        receptor_text = _pdbqt_to_pdb_string(receptor_pdbqt)

    # Best pose (first MODEL) of each top hit.
    hits = []
    if os.path.isdir(poses_dir):
        for lig, affinity in top:
            name = str(lig.get('Ligand', lig) if isinstance(lig, dict) else lig)
            path = _pose_file(lig, name, poses_dir)
            pose_text = _pdbqt_to_pdb_string(path, first_model_only=True) if path else None
            if pose_text:
                hits.append({"name": name, "affinity": affinity, "pose": pose_text})

    if not hits:
        return _viewer_notice(
            "No docked structures were available to display in the 3D viewer.")

    # Inline all data as one JSON blob; escape every "<" so the PDB text
    # cannot close or confuse the <script> element (</script>, <!--).
    data_json = json.dumps({"receptor": receptor_text or "", "hits": hits,
                            "gridBox": grid_box}).replace("<", "\\u003c")
    return _fill(_template("viewer.html"), {
        "__VS_BOX_CHECKED__": "checked" if grid_box else "",
        "__VS_DATA_JSON__": data_json,
        "__VS_NGL_SCRIPT__": ngl_js.decode("utf-8", errors="replace"),
    })


def generate_html_report(results_csv: str, poses_dir: str, output_file: str,
                          receptor_pdbqt: Optional[str] = None,
                          meta: Optional[Dict] = None,
                          grid_file: Optional[str] = None,
                          scorer_agreement: Optional[Dict] = None,
                          rmsd_threshold: float = 2.0) -> bool:
    """Write a single-file HTML screening report.

    Includes a table of contents, run metadata, summary statistics, an
    interactive NGL 3D viewer (with the docking box when ``grid_file`` is
    supplied), an affinity histogram, a molecular-weight vs affinity scatter
    plot, a sortable/filterable results table, pose clustering and an
    interpretation guide. The written file needs no network or server (NGL
    and the structures are inlined); see ``_fetch_ngljs`` for where NGL
    comes from while the report is generated.

    A row counts as scored when its Binding_Affinity is a negative number
    and its Status is not FAILED/REJECTED/ERROR/OUT_OF_WINDOW.
    """
    logger.info("[*] Generating HTML report...")

    ligands = []
    try:
        with open(results_csv, 'r') as f:
            ligands = list(csv.DictReader(f))
    except Exception as e:
        logger.warning(f"Could not parse results CSV: {e}")
        return False
    if not ligands:
        logger.warning("No results to report")
        return False

    scores = [_valid_score(lig) for lig in ligands]
    valid_affinities = [v for v in scores if v is not None]

    failed_count = len(ligands) - len(valid_affinities)
    avg_affinity = statistics.mean(valid_affinities) if valid_affinities else 0.0
    best_affinity = min(valid_affinities) if valid_affinities else 0.0
    worst_affinity = max(valid_affinities) if valid_affinities else 0.0
    strong_binders = sum(1 for a in valid_affinities if a < STRONG_BINDER_CUTOFF)
    success_rate = (len(valid_affinities) / len(ligands) * 100) if ligands else 0.0

    # Scored ligands, best first (name breaks ties).
    ranked_ligands = sorted(
        ((lig, v) for lig, v in zip(ligands, scores) if v is not None),
        key=lambda item: (item[1], str(item[0].get('Ligand', ''))))

    # Docking grid box for the viewer
    grid_box = None
    if grid_file and os.path.exists(grid_file):
        try:
            gc = runner._parse_grid_config(grid_file)
            if gc:
                grid_box = {"center": [gc['center_x'], gc['center_y'], gc['center_z']],
                            "size": [gc['size_x'], gc['size_y'], gc['size_z']]}
        except Exception as e:
            logger.debug(f"Could not parse grid config for viewer: {e}")

    # Scatter data (MW vs affinity)
    scatter_data = []
    for lig, av in zip(ligands, scores):
        try:
            mw = float(str(lig.get('MW', '')).strip())
        except ValueError:
            continue
        if av is not None and mw > 0:
            scatter_data.append({"name": lig.get('Ligand', '?'),
                                 "mw": mw, "aff": av})

    clusters = (cluster_poses(poses_dir, rmsd_threshold=rmsd_threshold,
                              ligands=[l.get('Ligand') for l in ligands])
                if os.path.isdir(poses_dir) else {})

    css = _REPORT_CSS(bool(ranked_ligands))
    generated = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    version = getattr(runner, "__version__", None) or "unknown"

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
    for lig, av in zip(ligands, scores):
        raw = str(lig.get('Binding_Affinity', '')).strip()
        if av is not None:
            conf_class, conf_text = _confidence(av)
        else:
            conf_class, conf_text = "bad", "N/A"
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
                <h2>Scorer Agreement (Vina vs SMINA)</h2>
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
    <button id="theme-toggle" class="float-btn" type="button">Theme</button>
    <button id="print-btn" class="float-btn no-print" type="button" onclick="window.print()">PDF</button>
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
                <h1>Virtual Screening Results Report</h1>
                <p>Generated: {generated} · AutoDocker {version}</p>
            </div>

            <div class="section" id="summary">
                <h2>Run Metadata</h2>
                <div class="meta-grid">{meta_items or '<p class="muted">No run metadata supplied.</p>'}</div>
                <h2 style="margin-top:18px;">Key Statistics</h2>
                <div class="stat-grid">
                    {_stat_card(len(ligands), "Ligands Screened")}
                    {_stat_card(f"{best_affinity:.2f}", "Best Binder (kcal/mol)")}
                    {_stat_card(f"{avg_affinity:.2f}", "Average Affinity")}
                    {_stat_card(f"{strong_binders}/{len(valid_affinities)}", f"Score &lt; {STRONG_BINDER_CUTOFF:g} kcal/mol")}
                    {_stat_card(f"{success_rate:.1f}%", "Scored ligands")}
                    {_stat_card(failed_count, "Failed / not scored")}
                </div>
            </div>

            {agreement_section}

            <div class="section" id="top">
                <h2>Top 10 Hits</h2>
                <table>
                    <thead><tr><th>#</th><th>Ligand Name</th><th>Affinity (kcal/mol)</th><th>Poses</th><th>Score band</th></tr></thead>
                    <tbody>{top_rows}
                    </tbody>
                </table>
            </div>

            {_build_viewer_section(ranked_ligands, poses_dir, receptor_pdbqt, output_file, grid_box=grid_box)}

            <div class="section" id="dist">
                <h2>Affinity Distribution</h2>
                <p class="muted"><em>Histogram of binding affinities across all docked ligands (kcal/mol).</em></p>
                {_svg_histogram(valid_affinities)}
            </div>

            <div class="section" id="scatter">
                <h2>Molecular Weight vs Affinity</h2>
                <p class="muted"><em>Each point is a ligand; click a point to load it in the 3D viewer.</em></p>
                {_scatter_chart(scatter_data)}
            </div>

            <div class="section" id="all">
                <h2>All Results</h2>
                <input id="results-search" type="text" placeholder="Filter ligands…" class="no-print">
                <table id="all-results">
                    <thead>
                    <tr>
                        <th class="sortable">Ligand</th>
                        <th class="sortable">Affinity</th>
                        <th class="sortable">Poses</th>
                        <th class="sortable">SimScore</th>
                        <th class="sortable">MW</th>
                        <th class="sortable">LogP</th>
                        <th class="sortable">Status</th>
                    </tr>
                    </thead>
                    <tbody>{all_rows}
                    </tbody>
                </table>
            </div>

            <div class="section" id="clusters">
                <h2>Pose Clustering</h2>
                <p class="muted"><em>Greedy grouping of each ligand's docked modes by heavy-atom
                RMSD (threshold {rmsd_threshold:g} Å). The RMSD pairs atoms by file order without
                superposition and is not symmetry-corrected, so symmetric groups can split clusters.</em></p>
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
                <h2>Interpretation Guide</h2>
                <p class="muted">Score bands used for the colours in this report. They are
                heuristics for sorting a screen: docking scores are not calibrated binding free
                energies, and hits need experimental or orthogonal validation.</p>
                <ul>
                    {_interpretation_items()}
                </ul>
            </div>

            <div class="section" style="text-align:center; color:var(--muted);">
                <p><small>Report generated by AutoDocker {version}</small></p>
            </div>
        </main>
    </div>
    <script>{_template("report.js")}</script>
</body>
</html>
"""

    try:
        with open(output_file, 'w') as f:
            f.write(html_content)
        logger.info(f"[✔] HTML report generated: {output_file}")
        return True
    except Exception as e:
        logger.warning(f"Could not write HTML report: {e}")
        if os.path.exists(output_file):
            os.remove(output_file)  # never leave a stale report from an earlier run
        return False
