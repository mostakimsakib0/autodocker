#!/usr/bin/env python3
"""Docking: per-ligand Vina runs, run fingerprint, SimScore and metrics.

Split out of runner.py. ``runner`` is imported lazily (see _LazyRunner) so
this module can be imported first, e.g. by a spawned worker. The engine
paths and the per-command timeout reach the worker processes explicitly,
through the :class:`EngineConfig` stored in each :class:`DockJob`, so they
survive the ``spawn`` start method (a spawned worker re-imports ``runner``
and would otherwise see its import-time defaults). Helpers that are only
called, not configured (``ligand_name``, ``_pdbqt_has_atoms``,
``_vina_cpu_flag`` ...), are resolved through the ``runner`` namespace at
call time so the tests' monkeypatches keep working.
"""
from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import re
import shlex
import subprocess
import time
import traceback
from collections import Counter
from dataclasses import dataclass, field
from enum import Enum
from multiprocessing import Pool
from typing import Dict, List, Optional, Sequence, Tuple, Union
from vspipeline import defaults  # noqa: E402



class _LazyRunner:
    """Attribute access to the ``runner`` module, imported on first use.

    runner.py owns the tool-path globals and the public names, and imports
    this module while it initializes. Importing runner lazily keeps this
    module importable on its own: a ``spawn``-started worker unpickles
    ``dock_ligand`` by importing vspipeline.docking first, and an eager
    ``import runner`` there would re-enter runner's ``from vspipeline.docking
    import ...`` on a half-initialized module. Every access reads the live
    module, so ``monkeypatch.setattr(runner, ...)`` in tests is honored.
    """

    def __getattr__(self, name):
        import runner as _runner
        return getattr(_runner, name)


runner = _LazyRunner()

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Named constants
# ---------------------------------------------------------------------------

# Accepted affinity window (kcal/mol). A Vina score is *in the window* iff
#
#     AFFINITY_FLOOR <= score < min_valid_affinity
#
# i.e. the lower bound is inclusive and the upper bound (``--min-valid-
# affinity``, default DEFAULT_MIN_VALID_AFFINITY) is exclusive: a score equal
# to -1.0 is out of the window, a score equal to -20.0 is in it. Scores below
# the floor are physically implausible for drug-like ligands and usually point
# to a broken input (clashes, wrong atom types). Out-of-window scores are still
# valid engine results: they are reported with status OUT_OF_WINDOW, never
# re-docked and never counted as failures.
AFFINITY_FLOOR = -20.0
DEFAULT_MIN_VALID_AFFINITY = defaults.MIN_VALID_AFFINITY

# In-window scores at or above this value get status WEAK instead of OK.
WEAK_AFFINITY_THRESHOLD = -4.0

# A ligand PDBQT smaller than this is treated as truncated/corrupt.
MIN_LIGAND_PDBQT_BYTES = 150

# Largest accepted grid-box edge (Angstrom).
MAX_GRID_EDGE = 80.0

# SimScore cut-offs (Angstrom) on Vina's RMSD columns; see
# _calculate_protocol_simscore.
SIMSCORE_RMSD_LB_CUTOFF = 1.6
SIMSCORE_RMSD_UB_CUTOFF = 3.2

# Timeout (s) for the short --help/--version probes of engine binaries.
_PROBE_TIMEOUT = 10


def affinity_in_window(score: float, min_valid_affinity: float) -> bool:
    """True iff ``AFFINITY_FLOOR <= score < min_valid_affinity``."""
    return AFFINITY_FLOOR <= score < min_valid_affinity


def window_text(min_valid_affinity: float) -> str:
    """Human-readable form of the accepted window, for logs and reports."""
    return f"{AFFINITY_FLOOR:g} <= score < {min_valid_affinity:g} kcal/mol"


class Status(str, Enum):
    """Per-ligand result status, as written to ranking.csv.

    OK / WEAK: in-window score (WEAK when >= WEAK_AFFINITY_THRESHOLD).
    OUT_OF_WINDOW: the engine produced a valid result whose score is outside
    the accepted window; the score is reported but not ranked with the hits.
    FAILED: no score (bad input, engine crash, timeout, unreadable output).
    """
    OK = "OK"
    WEAK = "WEAK"
    OUT_OF_WINDOW = "OUT_OF_WINDOW"
    FAILED = "FAILED"

    def __str__(self) -> str:
        return self.value

    @classmethod
    def of(cls, value, score: Optional[float] = None) -> "Status":
        """Normalize a stored status string (including the pre-2.0
        ``REJECTED (x)`` form) to a Status; missing -> OK/FAILED by score."""
        if isinstance(value, cls):
            return value
        text = str(value or "").strip()
        if text.startswith("REJECTED"):
            return cls.OUT_OF_WINDOW
        try:
            return cls(text)
        except ValueError:
            return cls.FAILED if score is None else cls.OK


def status_for_score(score: float) -> Status:
    """Status of an in-window score."""
    return Status.OK if score < WEAK_AFFINITY_THRESHOLD else Status.WEAK


# ---------------------------------------------------------------------------
# Engine configuration passed to workers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class EngineConfig:
    """Docking engines and limits, passed explicitly to every worker."""
    vina_primary: Optional[str]
    vina_fallback: Optional[str] = None
    smina: Optional[str] = None
    timeout: int = 900

    @classmethod
    def from_runner(cls) -> "EngineConfig":
        """Snapshot of the current runner globals (set by runner.main)."""
        return cls(vina_primary=runner.VINA_PRIMARY,
                   vina_fallback=runner.VINA_FALLBACK,
                   smina=runner.SMINA,
                   timeout=int(runner.COMMAND_TIMEOUT))

    def docking_engines(self) -> List[Tuple[str, str]]:
        """[(path, label)] in the order they are tried. The fallback is only
        present when it differs from the primary (runner enables it with
        VS_ENABLE_QVINA_FALLBACK=1)."""
        engines = []
        if self.vina_primary:
            engines.append((self.vina_primary,
                            f"{_tool_label(self.vina_primary)} (primary)"))
        if self.vina_fallback and self.vina_fallback != self.vina_primary:
            engines.append((self.vina_fallback,
                            f"{_tool_label(self.vina_fallback)} (fallback)"))
        return engines


@dataclass(frozen=True)
class DockJob:
    """Everything one worker needs to dock one ligand."""
    receptor: str
    ligand: str
    grid_file: str
    dock_dir: str
    vina_params: Dict
    flex_file: Optional[str] = None
    engine: EngineConfig = field(default_factory=EngineConfig.from_runner)


def _as_job(job: Union[DockJob, Sequence]) -> DockJob:
    """Accept a DockJob or the legacy 6-tuple
    ``(receptor, ligand, grid_file, dock_dir, vina_params, flex_file)``."""
    if isinstance(job, DockJob):
        return job
    receptor, lig, grid_file, dock_dir, vina_params, flex_file = job
    return DockJob(receptor, lig, grid_file, dock_dir, vina_params or {},
                   flex_file, EngineConfig.from_runner())


# ---------------------------------------------------------------------------
# Parsing engine output
# ---------------------------------------------------------------------------

_VINA_RESULT_RE = re.compile(r"^REMARK VINA RESULT:\s+(\S+)\s+(\S+)\s+(\S+)")
_MINIMIZED_AFFINITY_RE = re.compile(r"^REMARK minimizedAffinity\s+(\S+)")


def _parse_pdbqt_result_tokens(path: str) -> List[Tuple[str, Optional[str], Optional[str]]]:
    """Raw (affinity, rmsd_lb, rmsd_ub) tokens, one per MODEL (see
    parse_pdbqt_results)."""
    records: List[Tuple[str, Optional[str], Optional[str]]] = []
    vina = smina = None

    def close():
        if smina is not None:
            records.append((smina, None, None))
        elif vina is not None:
            records.append(vina)

    saw_model = False
    with open(path) as f:
        for line in f:
            if line.startswith("MODEL"):
                if saw_model:
                    close()
                saw_model = True
                vina = smina = None
            elif line.startswith("ENDMDL"):
                close()
                saw_model = False
                vina = smina = None
            elif vina is None and line.startswith("REMARK VINA RESULT:"):
                m = _VINA_RESULT_RE.match(line)
                if m:
                    vina = (m.group(1), m.group(2), m.group(3))
            elif smina is None and line.startswith("REMARK minimizedAffinity"):
                m = _MINIMIZED_AFFINITY_RE.match(line)
                if m:
                    smina = m.group(1)
    close()  # a file without MODEL records holds a single pose
    return records


def parse_pdbqt_results(path: str) -> List[Tuple[float, Optional[float], Optional[float]]]:
    """Scores written into a docking output PDBQT, one tuple per MODEL in
    file order: ``(affinity, rmsd_lb, rmsd_ub)``.

    This is the authoritative score source. The stdout table is rounded
    (Vina 1.2 prints 4 significant digits, -11.521 -> -11.52; QuickVina and
    SMINA print one decimal), while the PDBQT keeps Vina's 3 decimals.

    * Vina / QuickVina: ``REMARK VINA RESULT: affinity rmsd_lb rmsd_ub``.
    * SMINA: ``REMARK minimizedAffinity x``; its RMSD fields are None. When a
      model carries both (``smina --score_only`` copies the input's Vina
      remark), minimizedAffinity wins because it is the score SMINA computed.

    Models without either remark are skipped. Raises OSError if the file
    cannot be read and ValueError on a non-numeric score token.
    """
    out = []
    for aff, lb, ub in _parse_pdbqt_result_tokens(path):
        out.append((float(aff),
                    None if lb is None else float(lb),
                    None if ub is None else float(ub)))
    return out


def _decimals(token: str) -> int:
    """Number of decimals printed in a numeric token ('-11.52' -> 2)."""
    token = token.strip().lower()
    if "e" in token:
        return 6
    return len(token.split(".", 1)[1]) if "." in token else 0


def _scores_agree(a_token: str, b_token: str) -> bool:
    """True when two printed roundings of one value are consistent."""
    tol = 0.5 * 10 ** -_decimals(a_token) + 0.5 * 10 ** -_decimals(b_token) + 1e-9
    return abs(float(a_token) - float(b_token)) <= tol


def _mode_table_lines(output: str) -> List[str]:
    """Lines of the Vina/QuickVina/SMINA results table.

    The table follows a ``-----+------`` separator and ends at the first line
    that is not a mode row; parsing only those lines keeps other digit-led
    output (warnings, verbose messages) out of the scores. Output without a
    separator is returned whole, for older or abbreviated formats.
    """
    lines = output.splitlines()
    for i, line in enumerate(lines):
        if line.strip().startswith("-----+"):
            rows = []
            for row in lines[i + 1:]:
                parts = row.split()
                if not parts or not parts[0].isdigit():
                    break
                rows.append(row)
            return rows
    return lines


def _stdout_best_affinity_token(output: str) -> Optional[str]:
    """Affinity token of the first mode row of a stdout table, or None."""
    for line in _mode_table_lines(output or ""):
        parts = line.replace("|", " ").split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                float(parts[1])
            except ValueError:
                continue
            return parts[1]
    return None


def _parse_vina_modes(vina_output: str) -> List[Dict[str, float]]:
    """Parse Vina/QuickVina mode rows from captured stdout (rounded values;
    the pipeline itself uses parse_pdbqt_results)."""
    modes = []
    for line in _mode_table_lines(vina_output):
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


def _extract_score(output: str) -> float:
    """Best (lowest) affinity in a stdout mode table, or in a
    ``REMARK VINA RESULT`` line of the given text. Rounded to print precision;
    only used to cross-check the PDBQT score."""
    if not output:
        raise ValueError("Empty docking output")

    best = None
    for line in _mode_table_lines(output):
        parts = line.split()
        if len(parts) >= 2 and parts[0].isdigit():
            try:
                score = float(parts[1])
            except ValueError:
                continue
            if best is None or score < best:
                best = score
    if best is not None:
        return best

    match = re.search(r"REMARK VINA RESULT:\s*([-\d.]+)", output)
    if match:
        return float(match.group(1))

    raise ValueError("No valid affinity found")


# ---------------------------------------------------------------------------
# SimScore and metrics
# ---------------------------------------------------------------------------

def _calculate_protocol_simscore(
        modes: List[Dict[str, float]]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """Pose-consistency SimScore from Vina's RMSD columns.

    Definition: let M be the modes other than mode 1 (the best pose), and
    rmsd_lb / rmsd_ub the lower / upper-bound RMSD of each mode to mode 1 as
    reported by Vina (``REMARK VINA RESULT``). Then

        f_lb = |{m in M : rmsd_lb(m) < SIMSCORE_RMSD_LB_CUTOFF (1.6 A)}| / |M|
        f_ub = |{m in M : rmsd_ub(m) < SIMSCORE_RMSD_UB_CUTOFF (3.2 A)}| / |M|
        SimScore = (f_lb + f_ub) / 2,  in [0, 1]

    Higher means the alternative poses stay close to the best one. Returns
    ``(simscore, f_lb, f_ub)``; all three are None (undefined) when there is
    no alternative mode (|M| = 0, e.g. ``--binding-modes 1``) or a mode has no
    RMSD data (SMINA output).
    """
    comparisons = [mode for mode in modes if mode.get("mode") != 1]
    if not comparisons or any(mode.get("rmsd_lb") is None or mode.get("rmsd_ub") is None
                              for mode in comparisons):
        return None, None, None
    lb_fraction = sum(
        1 for mode in comparisons if mode["rmsd_lb"] < SIMSCORE_RMSD_LB_CUTOFF) / len(comparisons)
    ub_fraction = sum(
        1 for mode in comparisons if mode["rmsd_ub"] < SIMSCORE_RMSD_UB_CUTOFF) / len(comparisons)
    return (lb_fraction + ub_fraction) / 2, lb_fraction, ub_fraction


def _calculate_metrics(output_pdbqt: str, expected_modes: Optional[int] = None,
                       results: Optional[List[Tuple[float, Optional[float], Optional[float]]]] = None
                       ) -> Dict:
    """Mode count, SimScore and mode coverage of one docking output.

    ``results`` are the parsed ``parse_pdbqt_results(output_pdbqt)`` (read
    from the file when not given). ``simscore`` is None when undefined.
    """
    if results is None:
        results = parse_pdbqt_results(output_pdbqt)
    modes = [{"mode": i, "affinity": aff, "rmsd_lb": lb, "rmsd_ub": ub}
             for i, (aff, lb, ub) in enumerate(results, 1)]
    simscore, lb_fraction, ub_fraction = _calculate_protocol_simscore(modes)
    metrics = {
        "binding_modes": len(modes),
        "simscore": simscore,
        "rmsd_lb_pass_fraction": lb_fraction,
        "rmsd_ub_pass_fraction": ub_fraction,
    }
    if expected_modes:
        metrics["mode_coverage"] = min(len(modes) / expected_modes, 1.0)
    return metrics


# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------

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
        if not math.isfinite(v):
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

    required = ['center_x', 'center_y',
                'center_z', 'size_x', 'size_y', 'size_z']
    missing = [k for k in required if k not in config]
    if missing:
        raise ValueError(f"Grid config missing required parameters: {missing}")

    sx, sy, sz = config['size_x'], config['size_y'], config['size_z']
    if sx <= 0 or sy <= 0 or sz <= 0:
        raise ValueError(
            f"Invalid grid box size: {sx}x{sy}x{sz} (must be > 0)")
    if sx > MAX_GRID_EDGE or sy > MAX_GRID_EDGE or sz > MAX_GRID_EDGE:
        raise ValueError(
            f"Invalid grid box size: {sx}x{sy}x{sz} (must be <= {MAX_GRID_EDGE:g})")

    _GRID_CACHE[cache_key] = config
    return config


# ---------------------------------------------------------------------------
# Caches of repeated per-ligand inspections
# ---------------------------------------------------------------------------

_VINA_CPU_FLAG_CACHE: Dict[str, Optional[str]] = {}
_VINA_TYPE_CACHE: Dict[str, str] = {}
_FLEX_SUPPORT_CACHE: Dict[str, bool] = {}
_GRID_CACHE: Dict[Tuple[str, float], Dict] = {}
_RECEPTOR_CACHE: Dict[Tuple[str, int, int], Tuple[bool, bool]] = {}
_ENGINE_IDENTITY_CACHE: Dict[Tuple[str, int, int], Dict[str, str]] = {}


def _file_key(path: str) -> Tuple[str, int, int]:
    """Cache key that changes when the file is rewritten."""
    st = os.stat(path)
    return (os.path.abspath(path), st.st_mtime_ns, st.st_size)


def _cached_receptor_check(receptor: str) -> Tuple[bool, bool]:
    """Return (has_atoms, has_nonzero_charges) for the receptor, cached by
    (path, mtime_ns, size)."""
    from vspipeline.charges import _pdbqt_has_nonzero_charges
    key = _file_key(receptor)
    if key not in _RECEPTOR_CACHE:
        _RECEPTOR_CACHE[key] = (runner._pdbqt_has_atoms(receptor),
                                _pdbqt_has_nonzero_charges(receptor))
    return _RECEPTOR_CACHE[key]


def _help_text(path: str) -> Optional[str]:
    try:
        result = subprocess.run([path, "--help"], capture_output=True,
                                text=True, timeout=_PROBE_TIMEOUT)
    except Exception:
        return None
    return (result.stdout or "") + (result.stderr or "")


def _vina_cpu_flag(vina_path: str) -> Optional[str]:
    """Return the option that limits the binary's CPU threads, or None.

    Vina 1.2.x and QuickVina 2 call it --cpu; some forks use --threads.
    Cached per binary.
    """
    if vina_path in _VINA_CPU_FLAG_CACHE:
        return _VINA_CPU_FLAG_CACHE[vina_path]
    flag = None
    for help_opt in ("--help_advanced", "--help"):
        try:
            result = subprocess.run(
                [vina_path, help_opt],
                capture_output=True, text=True, timeout=_PROBE_TIMEOUT,
            )
        except Exception:
            continue
        help_text = (result.stdout or "") + (result.stderr or "")
        if re.search(r"--cpu\b", help_text):
            flag = "--cpu"
            break
        if re.search(r"--threads\b", help_text):
            flag = "--threads"
            break
    _VINA_CPU_FLAG_CACHE[vina_path] = flag
    return flag


def _detect_vina_type(vina_path: str) -> str:
    """Detect if Vina is QuickVina or standard Vina by checking --help output.
    Returns: 'quickvina' or 'vina'. Cached per binary path.
    """
    if vina_path in _VINA_TYPE_CACHE:
        return _VINA_TYPE_CACHE[vina_path]
    help_text = (_help_text(vina_path) or "").lower()
    if 'quickvina' in help_text or 'quick vina' in help_text:
        result = 'quickvina'
    elif 'autodock vina' in help_text:
        result = 'vina'
    else:
        basename = os.path.basename(vina_path).lower()
        result = 'quickvina' if ('qvina' in basename or 'quickvina' in basename) else 'vina'
    _VINA_TYPE_CACHE[vina_path] = result
    return result


def _engine_supports_flex(vina_path: str) -> bool:
    """Whether the engine accepts ``--flex``: Vina and QuickVina 2 do. When
    the binary's --help can be read it must list --flex; otherwise the
    detected engine type decides. Cached per binary."""
    if vina_path in _FLEX_SUPPORT_CACHE:
        return _FLEX_SUPPORT_CACHE[vina_path]
    supported = runner._detect_vina_type(vina_path) in ("vina", "quickvina")
    help_text = _help_text(vina_path)
    if supported and help_text and "--" in help_text:
        supported = re.search(r"--flex\b", help_text) is not None
    _FLEX_SUPPORT_CACHE[vina_path] = supported
    return supported


def _engine_version(path: str) -> str:
    """First line of ``<engine> --version`` (empty when unavailable)."""
    try:
        result = subprocess.run([path, "--version"], capture_output=True,
                                text=True, timeout=_PROBE_TIMEOUT)
    except Exception:
        return ""
    for line in ((result.stdout or "") + "\n" + (result.stderr or "")).splitlines():
        if line.strip():
            return line.strip()
    return ""


def _engine_identity(path: Optional[str]) -> Dict[str, str]:
    """{path, sha256, version} of an engine binary, cached per
    (path, mtime_ns, size). Paths that do not exist hash to ''."""
    if not path:
        return {"path": "", "sha256": "", "version": ""}
    resolved = os.path.realpath(path) if os.path.exists(path) else path
    try:
        key = _file_key(resolved)
    except OSError:
        return {"path": path, "sha256": "", "version": ""}
    if key not in _ENGINE_IDENTITY_CACHE:
        _ENGINE_IDENTITY_CACHE[key] = {
            "path": path,
            "sha256": _file_sha256(resolved),
            "version": _engine_version(path),
        }
    return _ENGINE_IDENTITY_CACHE[key]


# ---------------------------------------------------------------------------
# Command construction
# ---------------------------------------------------------------------------

def _build_vina_command(vina_path: str, receptor: str, ligand: str, out: str,
                        grid_config: Dict[str, float], vina_params: Dict,
                        grid_file: str = None, flex_file: str = None) -> list:
    """Build the docking command for a Vina-family binary.

    QuickVina: box from ``--config grid_file``. Standard Vina: explicit
    ``--center_*``/``--size_*``. ``--flex`` is passed whenever ``flex_file``
    is set, for every engine (dock_ligand refuses engines that cannot take it).
    """
    vina_type = runner._detect_vina_type(vina_path)

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
        cpu_flag = runner._vina_cpu_flag(vina_path)
        if cpu_flag:
            cmd.extend([cpu_flag, str(vina_params['threads'])])
        else:
            logger.warning(
                f"[!] Vina binary {vina_path} has no --cpu/--threads option; "
                "ignoring requested thread limit")

    if vina_type == 'quickvina':
        if grid_file:
            cmd.extend(["--config", grid_file])
    else:
        for key in ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z"):
            if key in grid_config:
                cmd.extend([f"--{key}", str(grid_config[key])])

    if flex_file:
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


# ---------------------------------------------------------------------------
# Docking one ligand
# ---------------------------------------------------------------------------

class _DockingError(Exception):
    """A per-ligand problem that makes the result FAILED."""


def _failed(name: str, out: str, log: str, error: str, **extra) -> Tuple[str, None, Dict]:
    logger.error(f"[x] {name}: {error}")
    return name, None, {"status": Status.FAILED.value, "dock_file": out,
                        "log_file": log, "error": error, **extra}


def _preflight(job: DockJob) -> Dict[str, float]:
    """Validate the inputs of one job; return the parsed grid config.

    Raises _DockingError for missing/empty/truncated inputs, an unreadable
    grid, or an engine that cannot take ``--flex`` in flexible mode. All-zero
    partial charges only produce a warning: Vina, QuickVina and SMINA scoring
    functions do not use partial charges.
    """
    if not os.path.exists(job.receptor):
        raise _DockingError(f"Receptor file missing: {job.receptor}")
    if not os.path.exists(job.ligand):
        raise _DockingError(f"Ligand file missing: {job.ligand}")
    if job.flex_file and not os.path.exists(job.flex_file):
        raise _DockingError(f"Flexible receptor file missing: {job.flex_file}")
    if not runner._pdbqt_has_atoms(job.ligand):
        raise _DockingError(f"Ligand has no ATOM/HETATM records: {job.ligand}")
    rec_has_atoms, _ = _cached_receptor_check(job.receptor)
    if not rec_has_atoms:
        raise _DockingError(f"Receptor has no ATOM/HETATM records: {job.receptor}")
    ligand_size = os.path.getsize(job.ligand)
    if ligand_size < MIN_LIGAND_PDBQT_BYTES:
        raise _DockingError(f"Ligand file too small ({ligand_size} bytes): {job.ligand}")
    from vspipeline.charges import _pdbqt_has_nonzero_charges
    if not _pdbqt_has_nonzero_charges(job.ligand):
        logger.warning(
            f"[!] {runner.ligand_name(job.ligand)}: all ligand partial charges are zero "
            "(not used by the Vina scoring function; docking continues)")
    if not job.engine.docking_engines():
        raise _DockingError("No docking engine configured")
    if job.flex_file:
        unsupported = [p for p, _ in job.engine.docking_engines()
                       if not _engine_supports_flex(p)]
        if unsupported:
            raise _DockingError(
                "Flexible docking requested but engine(s) do not accept --flex: "
                + ", ".join(unsupported)
                + " (refusing to dock against a receptor without the flexible side chains)")
    try:
        return _parse_grid_config(job.grid_file)
    except (ValueError, RuntimeError) as e:
        raise _DockingError(str(e))


def _run_engine(cmd: List[str], timeout: int, log) -> Tuple[str, str]:
    """Run one engine command, recording it in the per-ligand log first so
    the log names the command even when it crashes or times out. Returns
    (stdout, stderr); raises CalledProcessError / TimeoutExpired."""
    log.write(f"# Command: {shlex.join(cmd)}\n")
    log.flush()
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode, cmd, output=result.stdout, stderr=result.stderr)
    return result.stdout or "", result.stderr or ""


def _validate_output(name: str, out: str, stdout: str
                     ) -> Tuple[float, List[Tuple[float, Optional[float], Optional[float]]]]:
    """Check the output PDBQT and return (best score, parsed results).

    The score is the first ``REMARK VINA RESULT`` of the PDBQT; the stdout
    table only cross-checks it (a warning when they disagree beyond print
    precision).
    """
    if not os.path.exists(out):
        raise _DockingError(f"Docking output file not created: {out}")
    if os.path.getsize(out) == 0:
        raise _DockingError(f"Docking produced empty output file: {out}")
    try:
        with open(out) as f:
            has_model = any(line.startswith("MODEL") for line in f)
        tokens = _parse_pdbqt_result_tokens(out)
    except OSError as e:
        raise _DockingError(f"Cannot read docking output file: {e}")
    if not has_model:
        raise _DockingError("Docking produced no valid poses (no MODEL records)")
    if not tokens:
        raise _DockingError(f"No REMARK VINA RESULT line in docking output: {out}")
    try:
        results = [(float(a), None if lb is None else float(lb),
                    None if ub is None else float(ub)) for a, lb, ub in tokens]
    except ValueError as e:
        raise _DockingError(f"Unreadable score in docking output {out}: {e}")
    score = results[0][0]
    if not math.isfinite(score):
        raise _DockingError(f"Non-finite score {score} in docking output {out}")
    stdout_token = _stdout_best_affinity_token(stdout)
    if stdout_token is not None and not _scores_agree(tokens[0][0], stdout_token):
        logger.warning(
            f"[!] {name}: output PDBQT score {tokens[0][0]} disagrees with the stdout "
            f"table ({stdout_token}); using the PDBQT value")
    return score, results


def _build_metrics(job: DockJob, out: str, log_path: str, score: float,
                   results, engine_path: str, label: str) -> Dict:
    expected = int(job.vina_params.get('binding_modes', 0)) or None
    metrics = _calculate_metrics(out, expected_modes=expected, results=results)
    min_valid = float(job.vina_params.get('min_valid_affinity', DEFAULT_MIN_VALID_AFFINITY))
    metrics.update({
        "dock_file": out,
        "log_file": log_path,
        "engine": _tool_label(engine_path),
        "engine_path": engine_path,
        "vina_tool": label,
        "fallback_used": label.endswith("(fallback)"),
        "raw_score": score,
    })
    if affinity_in_window(score, min_valid):
        metrics["status"] = status_for_score(score).value
    else:
        metrics["status"] = Status.OUT_OF_WINDOW.value
        metrics["error"] = (f"score {score:.3f} kcal/mol outside the accepted window "
                            f"({window_text(min_valid)})")
    return metrics


def dock_ligand(job: Union[DockJob, Sequence]) -> Tuple[str, Optional[float], Dict]:
    """Dock one ligand; never raises for per-ligand problems.

    ``job`` is a :class:`DockJob` (the legacy 6-tuple ``(receptor, ligand,
    grid_file, dock_dir, vina_params, flex_file)`` is still accepted and uses
    the current runner globals).

    The primary engine (Vina by default) is tried first. The fallback engine
    (QuickVina, only when runner was started with VS_ENABLE_QVINA_FALLBACK=1)
    is tried only if the primary crashes, times out or writes unusable
    output; that is logged at WARNING and recorded in ``metrics["engine"]``
    so a ranking that mixes engines is visible. A valid result is never
    re-docked, whatever its score.

    Returns ``(name, score, metrics)``:

    * in-window result: score is the PDBQT affinity, status OK or WEAK;
    * out-of-window result: score None, ``metrics["raw_score"]`` holds the
      affinity and status is OUT_OF_WINDOW;
    * failure: score None, status FAILED and ``metrics["error"]`` the reason.

    Every attempt is recorded in ``docked/<name>_vina.log``.
    """
    job = _as_job(job)
    name = runner.ligand_name(job.ligand)
    out = os.path.join(job.dock_dir, name + "_out.pdbqt")
    log_path = os.path.join(job.dock_dir, name + "_vina.log")

    try:
        grid_config = _preflight(job)
    except _DockingError as e:
        return _failed(name, out, log_path, str(e))

    engines = job.engine.docking_engines()
    last_error = None
    with open(log_path, "w") as log:
        for idx, (engine_path, label) in enumerate(engines):
            # Never let a previous attempt's (or run's) poses pass for this one.
            if os.path.exists(out):
                os.remove(out)
            log.write(f"# Attempt: {label}\n")
            try:
                cmd = _build_vina_command(
                    engine_path, job.receptor, job.ligand, out, grid_config,
                    job.vina_params, job.grid_file, job.flex_file)
                for extra_arg in job.vina_params.get("extra_args", []) or []:
                    cmd.extend(shlex.split(extra_arg))
                stdout, stderr = _run_engine(cmd, job.engine.timeout, log)
                score, results = _validate_output(name, out, stdout)
                metrics = _build_metrics(job, out, log_path, score, results,
                                         engine_path, label)
                log.write(f"# Result: {metrics['status']} - Affinity: {score:.3f} kcal/mol\n")
                if metrics["status"] == Status.OUT_OF_WINDOW.value:
                    log.write(f"# Note: {metrics['error']}\n")
                    logger.warning(f"  [!] {name}: {metrics['error']}")
                log.write(f"# Output file size: {os.path.getsize(out)} bytes\n")
                log.write(stdout + "\n")
                if stderr.strip():
                    log.write("# STDERR:\n" + stderr + "\n")
                if idx > 0:
                    logger.warning(
                        f"[!] {name}: docked with fallback engine {label} "
                        f"({engine_path}); the ranking mixes engines")
                score_out = score if metrics["status"] != Status.OUT_OF_WINDOW.value else None
                return name, score_out, metrics
            except subprocess.TimeoutExpired:
                last_error = f"{label}: timeout after {job.engine.timeout}s"
                log.write(f"# Error: {last_error}\n\n")
            except subprocess.CalledProcessError as e:
                last_error = f"{label}: command failed (exit {e.returncode})"
                log.write(f"# Error: {last_error}\n")
                if e.output:
                    log.write("# STDOUT:\n" + e.output + "\n")
                if e.stderr:
                    log.write("# STDERR:\n" + e.stderr + "\n")
                log.write("\n")
            except _DockingError as e:
                last_error = f"{label}: {e}"
                log.write(f"# Error: {last_error}\n\n")
            except Exception as e:
                last_error = f"{label}: {e}"
                log.write(f"# Error: {last_error}\n")
                log.write("# Traceback:\n" + traceback.format_exc() + "\n")
            if idx + 1 < len(engines):
                logger.warning(
                    f"[!] {name}: {last_error}; retrying with {engines[idx + 1][1]}")
        log.write("\n# FINAL ERROR: All docking tools failed\n")
        log.write(f"# Last error: {last_error}\n")

    logger.error(f"See log for details: {log_path}")
    return _failed(name, out, log_path,
                   f"Docking failed for {name} with all available tools. "
                   f"Last error: {last_error}")


# ---------------------------------------------------------------------------
# Run fingerprint
# ---------------------------------------------------------------------------

def _file_sha256(path: Optional[str]) -> str:
    if not path or not os.path.exists(path):
        return ""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _receptor_sha256(path: Optional[str]) -> str:
    """SHA-256 of the complete receptor (or flex) PDBQT, hydrogens included.

    Hydrogen coordinates matter: Vina perceives donor hydrogens by bond
    distance (model.cpp, assign_bonds / bonded_to_HD), so moving a polar H
    can change which heavy atom is a donor, and different Open Babel H
    placements were observed to give different scores.
    """
    return _file_sha256(path)


def _run_fingerprint(receptor: str, grid_file: str, flex_file: Optional[str],
                     vina_params: Dict, engine: Optional[EngineConfig] = None) -> str:
    """Hash of everything that determines a ligand's docking result apart from
    the ligand itself: receptor and flex files (byte-exact), grid file,
    docking parameters and the engine binaries (path, SHA-256, version). The
    thread count is left out: seeded Vina output does not depend on it."""
    engine = engine or EngineConfig.from_runner()
    params = {k: v for k, v in vina_params.items() if k != "threads"}
    payload = json.dumps({
        "receptor": _file_sha256(receptor),
        "grid": _file_sha256(grid_file),
        "flex": _file_sha256(flex_file),
        "params": params,
        "engines": [_engine_identity(engine.vina_primary),
                    _engine_identity(engine.vina_fallback)],
    }, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()


# ---------------------------------------------------------------------------
# Docking a library
# ---------------------------------------------------------------------------

def _default_vina_params() -> Dict:
    """Defaults for keys a caller leaves out ($VS_* environment, else built-in)."""
    return {
        'exhaustiveness': int(os.environ.get('VS_EXHAUSTIVENESS', defaults.EXHAUSTIVENESS)),
        'binding_modes': int(os.environ.get('VS_BINDING_MODES', defaults.BINDING_MODES)),
        'energy_range': float(os.environ.get('VS_ENERGY_RANGE', defaults.ENERGY_RANGE)),
        'seed': int(os.environ.get('VS_SEED', defaults.SEED)),
        'min_valid_affinity': float(os.environ.get('VS_MIN_VALID_AFFINITY',
                                                   DEFAULT_MIN_VALID_AFFINITY)),
        'threads': (int(os.environ['VS_THREADS'])
                    if os.environ.get('VS_THREADS') else None),
    }


def dock_all(receptor: str, ligands: List[str], grid_file: str,
             outdir: str, num_processes: int = 1,
             resume: bool = True,
             vina_params: Optional[Dict] = None,
             flex_file: Optional[str] = None,
             engine: Optional[EngineConfig] = None
             ) -> Tuple[List[Tuple[str, Optional[float]]], "runner.DockingCheckpoint", Dict]:
    """Dock all ligands with resume and metrics.

    ``vina_params`` keys given by the caller always win; only missing keys
    are filled from ``$VS_*`` / built-in defaults. ``flex_file`` optionally
    points to a Vina flexible-receptor PDBQT. ``engine`` defaults to a
    snapshot of the runner globals.

    Returns ``(results, checkpoint, metrics)``. ``results`` holds
    ``(name, score)`` with score None for FAILED and OUT_OF_WINDOW ligands
    (the latter keep their affinity in ``metrics[name]["raw_score"]``).
    Raises RuntimeError only when no ligand produced any score at all.
    """
    engine = engine or EngineConfig.from_runner()
    dock_dir = os.path.join(outdir, "docked")
    os.makedirs(dock_dir, exist_ok=True)

    # Guard against output collisions when ligand basenames are duplicated.
    ligand_names = [runner.ligand_name(lig) for lig in ligands]
    duplicate_names = sorted(
        n for n, k in Counter(ligand_names).items() if k > 1)
    if duplicate_names:
        raise ValueError(
            "Duplicate ligand names detected (would overwrite outputs): "
            + ", ".join(duplicate_names)
        )

    vina_params = {**_default_vina_params(), **(vina_params or {})}

    # Each Vina process uses every core by default, so N workers would
    # oversubscribe the CPU N-fold. Split the cores between workers unless
    # the user chose a limit. Seeded Vina output does not depend on --cpu.
    if not vina_params.get('threads') and num_processes > 1:
        vina_params['threads'] = max(1, runner.available_cpus() // num_processes)

    ligand_hashes = {runner.ligand_name(lig): _file_sha256(lig) for lig in ligands}
    checkpoint = runner.DockingCheckpoint(
        outdir, fingerprint=_run_fingerprint(receptor, grid_file, flex_file,
                                             vina_params, engine),
        resume=resume)

    if resume:
        ligands_todo = [
            lig for lig in ligands
            if not checkpoint.is_completed(runner.ligand_name(lig),
                                           ligand_hashes[runner.ligand_name(lig)])
        ]
    else:
        ligands_todo = ligands
    # A ligand docked again (its file changed) must not keep its old result.
    for lig in ligands_todo:
        checkpoint.completed.pop(runner.ligand_name(lig), None)

    min_valid_affinity = float(vina_params['min_valid_affinity'])
    logger.info(
        f"[*] Docking {len(ligands_todo)} new ligands "
        f"({len(ligands) - len(ligands_todo)} cached)")
    logger.info(f"[*] Vina parameters: exhaustiveness={vina_params['exhaustiveness']}, "
                f"modes={vina_params['binding_modes']}, "
                f"range={vina_params['energy_range']} kcal/mol, "
                f"seed={vina_params['seed']}, "
                f"cpu_per_process={vina_params.get('threads') or 'all'}, "
                f"accepted window: {window_text(min_valid_affinity)}")
    logger.info("[*] Engines: " + ", ".join(
        f"{label}={path}" for path, label in engine.docking_engines())
        + f"; timeout={engine.timeout}s")
    if num_processes > 1:
        logger.info(f"[*] Parallel mode: {num_processes} workers, unique per-ligand outputs")

    if os.path.exists(receptor):
        has_atoms, has_charges = _cached_receptor_check(receptor)
        if has_atoms and not has_charges:
            logger.warning(
                f"[!] All receptor partial charges are zero in {receptor} "
                "(not used by the Vina scoring function; docking continues)")

    if not ligands_todo:
        logger.info("[*] All ligands already docked (checkpoint)")
        results = checkpoint.get_results(ligand_names)
        return results, checkpoint, checkpoint.get_metrics(ligand_names)

    jobs = [DockJob(receptor, lig, grid_file, dock_dir, vina_params, flex_file, engine)
            for lig in ligands_todo]

    metrics_dict = checkpoint.get_metrics(ligand_names) if resume else {}
    failed_results = []

    def record(result):
        """Store one ligand's result as soon as it arrives, so an interrupted
        screen keeps everything docked so far."""
        name, score, metrics = result
        metrics_dict[name] = metrics
        status = Status.of(metrics.get("status"), score)
        if status is Status.OUT_OF_WINDOW:
            # A seeded run that scored outside the window gives the same
            # score again, so it is cached like an in-window result.
            checkpoint.save_result(name, None, metrics,
                                   ligand_hash=ligand_hashes.get(name))
            return
        if score is None:
            logger.warning(
                f"  [!] {name}: FAILED (see {metrics.get('log_file', 'log')})")
            # Failures are not checkpointed; they are retried on resume.
            failed_results.append((name, None))
            return
        logger.info(f"  [+] {name}: {score:.3f} kcal/mol")
        checkpoint.save_result(name, score, metrics, ligand_hash=ligand_hashes.get(name))

    start_time = time.time()
    total = len(jobs)

    def log_progress(done: int):
        elapsed = max(time.time() - start_time, 1e-6)
        rate_per_min = done / elapsed * 60
        remaining = total - done
        eta_minutes = remaining / rate_per_min if rate_per_min else 0
        logger.info(
            f"[*] Progress: {done}/{total} docked "
            f"({rate_per_min:.2f} ligands/min, ETA {eta_minutes:.1f} min)"
        )

    try:
        if num_processes > 1:
            with Pool(num_processes) as pool:
                for idx, result in enumerate(pool.imap_unordered(dock_ligand, jobs), 1):
                    record(result)
                    log_progress(idx)
        else:
            for idx, job in enumerate(jobs, 1):
                record(dock_ligand(job))
                log_progress(idx)
    finally:
        checkpoint.flush()  # consolidate O(1) journal writes into the JSON file

    all_results = checkpoint.get_results(ligand_names)
    all_results.extend(failed_results)

    # Fail loud only when no ligand produced any score (in or out of the
    # window). A resumed run whose remaining ligands keep failing is a partial
    # success, and an all-out-of-window screen is still a reportable result.
    any_score = any(score is not None for _, score in all_results) or any(
        Status.of(metrics_dict.get(name, {}).get("status")) is Status.OUT_OF_WINDOW
        for name, _ in all_results)
    if not any_score:
        raise RuntimeError(
            "All docking attempts failed; no ligand produced a score. "
            "Inspect output/docked/*_vina.log for root causes."
        )

    return all_results, checkpoint, metrics_dict
