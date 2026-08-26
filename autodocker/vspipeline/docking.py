#!/usr/bin/env python3
"""DOCKING WITH RESUME + SimScore + METRICS.

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


def _parse_grid_triplet(value: str, label: str) -> Tuple[float, float, float]:
    """Parse an explicit 'X,Y,Z' grid triplet. STRICT validation: any malformed
    value raises, never silently falls back to auto-detection."""
    parts = [p.strip() for p in value.split(',')]
    if len(parts) != 3:
        raise ValueError(
            f"--{label} must be three comma-separated values 'X,Y,Z', got '{value}'")
    try:
        result = tuple(float(p) for p in parts)
    except ValueError:
        raise ValueError(
            f"--{label} values must be floats, got '{value}'")
    for v in result:
        if v != v:  # NaN
            raise ValueError(
                f"--{label} values must be finite floats, got '{value}'")
    return result


def _parse_grid_config(grid_file: str) -> Dict[str, float]:
    """Parse grid.conf file and extract center/size coordinates. STRICT validation.
    Results are cached by (path, mtime) so the same file isn't re-parsed per ligand."""
    try:
        mtime = os.path.getmtime(grid_file)
    except OSError:
        mtime = -1.0
    cache_key = (grid_file, mtime)
    if cache_key in _GRID_CACHE:
        return _GRID_CACHE[cache_key]
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

    _GRID_CACHE[cache_key] = config
    return config


# Cache repeated, expensive, per-ligand inspections so they run once per
# unique input instead of N times across the ligand set.
_VINA_THREADS_CACHE: Dict[str, bool] = {}
_VINA_TYPE_CACHE: Dict[str, str] = {}
_GRID_CACHE: Dict[Tuple[str, float], Dict] = {}
_RECEPTOR_CACHE: Dict[str, Tuple[bool, bool]] = {}


def _cached_receptor_check(receptor: str) -> Tuple[bool, bool]:
    """Return (has_atoms, has_charges) for the receptor, cached by path."""
    if receptor in _RECEPTOR_CACHE:
        return _RECEPTOR_CACHE[receptor]
    res = (runner._pdbqt_has_atoms(receptor),
           runner._ensure_pdbqt_has_charges(receptor))
    _RECEPTOR_CACHE[receptor] = res
    return res


def _vina_supports_threads(vina_path: str) -> bool:
    """Probe whether the Vina binary accepts --threads. Cached per binary."""
    if vina_path in _VINA_THREADS_CACHE:
        return _VINA_THREADS_CACHE[vina_path]
    supported = False
    try:
        result = subprocess.run(
            [vina_path, "--help_advanced"],
            capture_output=True, text=True, timeout=5,
        )
        help_text = (result.stdout or "") + (result.stderr or "")
        supported = "--threads" in help_text
    except Exception:
        supported = False
    _VINA_THREADS_CACHE[vina_path] = supported
    return supported


def _detect_vina_type(vina_path: str) -> str:
    """Detect if Vina is QuickVina or standard Vina by checking --help output.
    Returns: 'quickvina' or 'vina'. Cached per binary path.
    """
    if vina_path in _VINA_TYPE_CACHE:
        return _VINA_TYPE_CACHE[vina_path]
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
    result = 'quickvina' if ('qvina' in basename or 'quickvina' in basename) else 'vina'
    _VINA_TYPE_CACHE[vina_path] = result
    return result


def _build_vina_command(vina_path: str, receptor: str, ligand: str, out: str,
                        grid_config: Dict[str, float], vina_params: Dict,
                        grid_file: str = None, flex_file: str = None) -> list:
    """Build appropriate Vina command based on Vina type.

    QuickVina: uses --config file
    Standard Vina: uses --center_x/y/z and --size_x/y/z args
    Flexible docking (--flex) is only supported by standard Vina.
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

    if vina_params.get('threads'):
        # Resolve via runner so test monkeypatches are honored.
        if runner._vina_supports_threads(vina_path):
            cmd.extend(["--threads", str(vina_params['threads'])])
        else:
            logger.warning(
                f"[!] Vina binary {vina_path} does not support --threads; "
                "ignoring requested thread limit")

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

    if flex_file and vina_type != 'quickvina':
        cmd.extend(["--flex", flex_file])

    return cmd


def _log_tool_version(tool_path: str, label: str) -> None:
    """Best-effort tool version logging for reproducibility."""
    try:
        stdout, stderr, _ = runner.run([tool_path, "--version"], capture=True)
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
    receptor, lig, grid_file, dock_dir, vina_params, flex_file = args_tuple

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

    if flex_file and not os.path.exists(flex_file):
        error_msg = f"Flexible receptor file missing: {flex_file}"
        logger.error(f"❌ {name}: {error_msg}")
        return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}

    # Validate basic file content (ligand is per-ligand; receptor is cached)
    if not runner._pdbqt_has_atoms(lig):
        error_msg = f"Ligand has no ATOM/HETATM records: {lig}"
        logger.error(f"❌ {name}: {error_msg}")
        return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}

    # Validate files have charges (receptor check cached per path)
    _rec_has_atoms, _rec_has_charges = _cached_receptor_check(receptor)
    if not _rec_has_atoms:
        error_msg = f"Receptor has no ATOM/HETATM records: {receptor}"
        logger.error(f"❌ {name}: {error_msg}")
        return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}

    if not _rec_has_charges:
        error_msg = f"Receptor has no charges: {receptor}"
        logger.error(f"❌ {name}: {error_msg}")
        return name, None, {"status": "FAILED", "dock_file": out, "log_file": vina_log, "error": error_msg}

    if not runner._ensure_pdbqt_has_charges(lig):
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
        (runner.VINA_PRIMARY, f"{_tool_label(runner.VINA_PRIMARY)} (primary)"),
        (runner.VINA_FALLBACK,
         f"{_tool_label(runner.VINA_FALLBACK)} (fallback)") if runner.VINA_FALLBACK != runner.VINA_PRIMARY else None
    ]
    vina_tools = [t for t in vina_tools if t]  # Remove None entries

    score = None
    last_error = None
    tool_log_entries = []

    for vina_tool, tool_label in vina_tools:
        try:
            # Build command for current Vina binary (auto-detects type)
            cmd = _build_vina_command(
                vina_tool, receptor, lig, out, grid_config, vina_params,
                grid_file, flex_file)

            for extra_arg in vina_params.get("extra_args", []):
                cmd.extend(shlex.split(extra_arg))

            logger.debug(f"Attempting docking {name} with {tool_label}...")
            stdout, stderr, retcode = runner.run(cmd, capture=True)

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
            last_error = f"{tool_label}: timeout after {runner.COMMAND_TIMEOUT}s"
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
             vina_params: Optional[Dict] = None,
             flex_file: Optional[str] = None) -> Tuple[List[Tuple[str, Optional[float]]], runner.DockingCheckpoint, Dict]:
    """Dock all ligands with resume and metrics.

    ``flex_file`` optionally points to a Vina flexible-receptor PDBQT
    (Phase 3 flexible docking).
    """
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
            'threads': (int(os.environ['VS_THREADS'])
                        if os.environ.get('VS_THREADS') else None),
        }
    else:
        # Merge caller-provided params over defaults so a partial dict never
        # causes a KeyError downstream.
        defaults = {
            'exhaustiveness': int(os.environ.get('VS_EXHAUSTIVENESS', 8)),
            'binding_modes': int(os.environ.get('VS_BINDING_MODES', 9)),
            'energy_range': float(os.environ.get('VS_ENERGY_RANGE', 3.0)),
            'seed': int(os.environ.get('VS_SEED', 42)),
            'min_valid_affinity': float(os.environ.get('VS_MIN_VALID_AFFINITY', -1.0)),
            'threads': (int(os.environ['VS_THREADS'])
                        if os.environ.get('VS_THREADS') else None),
        }
        defaults.update(vina_params)
        vina_params = defaults

    checkpoint_file = os.path.join(outdir, ".docking_checkpoint.json")
    journal_file = os.path.join(outdir, ".docking_checkpoint.jsonl")
    if not resume:
        for _f in (checkpoint_file, journal_file):
            if os.path.exists(_f):
                os.remove(_f)

    checkpoint = runner.DockingCheckpoint(outdir)

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

    args_list = [(receptor, lig, grid_file, dock_dir, vina_params, flex_file)
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
    checkpoint.flush()  # consolidate O(1) journal writes into the JSON file

    # Fail loud only when the *whole* screen produced no valid score.
    # On a resumed run, cached successes may already exist while the only
    # remaining ligands are permanent outliers (unrealistic affinities,
    # malformed inputs) that will never succeed; that is a partial success,
    # not a pipeline failure.
    if not any(score is not None for _, score in all_results):
        raise RuntimeError(
            "All docking attempts failed; no valid ligand scores were produced. "
            "Inspect output/docked/*_vina.log for root causes."
        )

    return all_results, checkpoint, metrics_dict
