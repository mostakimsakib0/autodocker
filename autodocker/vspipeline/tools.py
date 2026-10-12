#!/usr/bin/env python3
"""External tools: discovery on PATH and the run's tool settings.

The tools are looked up once, when this module is imported. The pipeline
modules read the paths and the per-command timeout through ``settings()``
at call time; runner.main() applies --vina-bin, --qvina-fallback and
--timeout with ``configure()`` once the options are validated. Docking and
scoring workers do not read these settings: they receive the engine paths
and the timeout explicitly (docking.EngineConfig, built from ``settings()``
when docking starts), so they survive the ``spawn`` start method.
"""
import os
import shutil
from dataclasses import dataclass, replace
from typing import List, Optional

from vspipeline import defaults

DEFAULT_COMMAND_TIMEOUT = defaults.COMMAND_TIMEOUT  # seconds

# Problems with environment variables read at import time; reported by
# runner.main() after argument parsing (so --version still works).
ENV_ERRORS: List[str] = []


def find_tool(*names: str) -> Optional[str]:
    """Return the path of the first *name* found on PATH, or None.

    Returns None when nothing is found instead of a placeholder name so
    callers can reliably distinguish 'found' from 'not found'.
    """
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    return None


def env_command_timeout() -> int:
    """$VS_COMMAND_TIMEOUT (seconds), else the default. An invalid value is
    recorded in ENV_ERRORS and the default is used."""
    raw = os.environ.get("VS_COMMAND_TIMEOUT", "").strip()
    if not raw:
        return DEFAULT_COMMAND_TIMEOUT
    try:
        value = int(raw)
        if value <= 0:
            raise ValueError
        return value
    except ValueError:
        ENV_ERRORS.append(
            f"VS_COMMAND_TIMEOUT must be a positive integer (seconds), got {raw!r}")
        return DEFAULT_COMMAND_TIMEOUT


@dataclass(frozen=True)
class ToolSettings:
    """Paths of the external tools (None = not found) and the per-command
    timeout of the run."""
    obabel: Optional[str]
    vina_primary: Optional[str]
    # Equals vina_primary (no fallback) unless --qvina-fallback is given;
    # runner.main() then sets it to QuickVina 2.
    vina_fallback: Optional[str]
    qvina: Optional[str]
    smina: Optional[str]
    fpocket: Optional[str]
    # reduce and AutoDock4/AutoGrid4 are not used by the pipeline (receptor
    # hydrogens come from Open Babel); the paths are kept for the real-tool
    # tests.
    reduce: Optional[str]
    autodock: Optional[str]
    autogrid: Optional[str]
    command_timeout: int


def discover() -> ToolSettings:
    """The tools found on PATH (prefer modern Vina for consistency) and the
    timeout of $VS_COMMAND_TIMEOUT."""
    vina = find_tool("vina", "vina.exe", "qvina02", "qvina2", "qvina")
    return ToolSettings(
        obabel=find_tool("obabel", "OpenBabel", "OpenBabel.exe"),
        vina_primary=vina,
        vina_fallback=vina,
        qvina=find_tool("qvina02", "qvina2"),
        smina=find_tool("smina", "smina.exe"),
        fpocket=find_tool("fpocket"),
        reduce=find_tool("reduce"),
        autodock=find_tool("autodock4", "autodock"),
        autogrid=find_tool("autogrid4", "autogrid"),
        command_timeout=env_command_timeout(),
    )


_settings = discover()


def settings() -> ToolSettings:
    """The current tool settings."""
    return _settings


def configure(**changes) -> ToolSettings:
    """Replace fields of the current settings (``ToolSettings`` field names);
    returns the new settings."""
    global _settings
    _settings = replace(_settings, **changes)
    return _settings
