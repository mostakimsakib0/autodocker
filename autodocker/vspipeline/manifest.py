#!/usr/bin/env python3
"""run_manifest.json: versions, settings and checksums of a run.

runner.py writes the manifest at the end of every run that got past option
validation, whatever its outcome (see runner._main_locked); the entries that
describe the run's progress come from pipeline.RunState.
"""
import hashlib
import json
import logging
import math
import os
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple

from vspipeline import hydrogens, tools, utils
from vspipeline.version import __version__

logger = logging.getLogger(__name__)

# Timeout (seconds) for the "--version" probes of the external tools.
VERSION_PROBE_TIMEOUT = 30

_ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")
# (version flag, regular expression picking the version line or None).
_VERSION_PROBES = {
    "obabel": ("-V", None),
    "fpocket": ("-h", re.compile(r"fpocket\s+\d[\w.]*")),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _sha256(path: str) -> Optional[str]:
    try:
        digest = hashlib.sha256()
        with open(path, "rb") as f:
            for block in iter(lambda: f.read(1 << 20), b""):
                digest.update(block)
        return digest.hexdigest()
    except OSError:
        return None


def probe_tool(kind: str, path: Optional[str]) -> Dict:
    """Path, resolved path, version line and sha256 of an external tool."""
    info = {"path": path, "resolved_path": None, "version": None, "sha256": None}
    if not path:
        return info
    resolved = shutil.which(path)
    resolved = os.path.realpath(resolved) if resolved else None
    info["resolved_path"] = resolved
    info["sha256"] = _sha256(resolved) if resolved else None
    flag, pattern = _VERSION_PROBES.get(kind, ("--version", None))
    try:
        proc = subprocess.run([path, flag], capture_output=True, text=True,
                              timeout=VERSION_PROBE_TIMEOUT)
        text = _ANSI_ESCAPE.sub("", (proc.stdout or "") + "\n" + (proc.stderr or ""))
    except (OSError, subprocess.SubprocessError) as e:
        info["version"] = f"unknown ({type(e).__name__}: {e})"
        return info
    if pattern:
        match = pattern.search(text)
        info["version"] = match.group(0) if match else None
    else:
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        info["version"] = lines[0] if lines else None
    return info


def _version_text(version: Optional[Tuple[int, ...]]) -> Optional[str]:
    return ".".join(str(v) for v in version) if version else None


def check_openbabel(seed: int) -> Dict:
    """Open Babel version policy; returns the run_manifest.json entries.

    Open Babel >= 3.2.0 is recommended and 3.2.1 is what the container image
    builds and the 2.x results were produced with. An older Open Babel (e.g.
    3.1.1 from a Linux distribution) is accepted with a loud warning: it
    ignores OB_RANDOM_SEED and over-protonates pre-protonated Lys/Arg, but
    the input-hydrogen removal and the hydrogen re-placement still make its
    output deterministic. Results then differ from the image's (MMFF94
    minimization), and the manifest records the version so they stay
    traceable."""
    obabel = tools.settings().obabel
    version = utils.obabel_version(obabel)
    seeded = version is not None and version >= utils.OBABEL_SEEDABLE_VERSION
    validated = _version_text(utils.OBABEL_VALIDATED_VERSION)
    if version is None:
        logger.warning(f"[!] Could not read the Open Babel version from '{obabel} -V'; "
                       f"AutoDocker {__version__} is validated with Open Babel {validated}")
    elif not seeded:
        logger.warning("!" * 70)
        logger.warning(
            f"[!] Open Babel {_version_text(version)} found; AutoDocker {__version__} "
            f"needs Open Babel >= {_version_text(utils.OBABEL_SEEDABLE_VERSION)} for its "
            f"published behaviour (the container image ships {validated}).")
        logger.warning(
            "[!] This version ignores OB_RANDOM_SEED, and its MMFF94 minimization and "
            "PDBQT output differ: results are reproducible run to run, but NOT "
            f"comparable with results from Open Babel {validated}. See docs/openbabel.md.")
        logger.warning("!" * 70)
    elif version != utils.OBABEL_VALIDATED_VERSION:
        logger.info(f"[*] Open Babel {_version_text(version)} (AutoDocker {__version__} "
                    f"is validated with {validated})")
    return {
        "openbabel_version": _version_text(version),
        "openbabel_random_seed": {
            "variable": utils.OB_RANDOM_SEED_ENV,
            "value": seed,
            "hydrogen_probe_values": [hydrogens.probe_seed(seed, 1),
                                      hydrogens.probe_seed(seed, 2)],
            "honoured": seeded,
            "policy": (
                "Every Open Babel call runs with OB_RANDOM_SEED set to the run's --seed; "
                "the hydrogen-placement probe conversions use the next two seeds, and "
                "hydrogens that move between seeds are re-placed from heavy-atom "
                "geometry. Open Babel < 3.2.0 ignores the variable (clock seed)."),
        },
    }


def _environment_snapshot() -> Dict[str, str]:
    return {k: v for k, v in sorted(os.environ.items())
            if k.startswith("VS_") or k in ("OMP_NUM_THREADS", "AUTODOCKER_VERSION",
                                            "PUBCHEM_REST", "BABEL_LIBDIR",
                                            "BABEL_DATADIR")}


def _python_packages() -> Dict[str, Optional[str]]:
    """Versions of the optional Python packages that change results (None =
    not installed). Without RDKit, 2D ligands get 3D coordinates from obabel
    --gen3d, which is not reproducible."""
    packages: Dict[str, Optional[str]] = {}
    for name in ("rdkit", "numpy"):
        try:
            module = __import__(name)
            packages[name] = str(getattr(module, "__version__", "unknown"))
        except ImportError:
            packages[name] = None
    return packages


def _json_safe(value):
    """``value`` with every non-finite float replaced by its string ("nan",
    "inf"), so run_manifest.json is strict JSON (no NaN/Infinity tokens)."""
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


def write_run_manifest(path: str, manifest: Dict) -> None:
    """Write run_manifest.json with sorted keys (deterministic layout) as
    strict JSON, atomically: a reader (or a run killed mid-write) sees the
    previous file or the new one, never a truncated one."""
    tmp = f"{path}.tmp.{os.getpid()}"
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(_json_safe(manifest), f, indent=2, sort_keys=True, default=str,
                      allow_nan=False)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def build_manifest(args, argv: List[str], tools_used: Dict, started_at: str,
                   status: str, exit_code: int, error: Optional[str],
                   extra: Dict) -> Dict:
    """The manifest of a run: ``args`` (the resolved options), the tools
    probed (``tools_used``), the outcome and ``extra`` (e.g.
    pipeline.RunState.manifest_entries())."""
    resolved = {k: v for k, v in sorted(vars(args).items())}
    return {
        "autodocker_version": __version__,
        "argv": [sys.argv[0]] + list(argv),
        "args": resolved,
        "seed": getattr(args, "seed", None),
        "status": status,
        "exit_code": exit_code,
        "error": error,
        "python": {"version": platform.python_version(),
                   "implementation": platform.python_implementation(),
                   "executable": sys.executable},
        "platform": platform.platform(),
        "cpu_count": {"os": os.cpu_count(), "available": utils.available_cpus()},
        "tools": tools_used,
        "python_packages": _python_packages(),
        "environment": _environment_snapshot(),
        "command_timeout_s": tools.settings().command_timeout,
        "started_at": started_at,
        "finished_at": utc_now(),
        **extra,
    }
