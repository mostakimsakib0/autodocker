#!/usr/bin/env python3
"""AutoDocker's version string (single source of truth: autodocker/VERSION)."""
import os
from typing import Optional

FALLBACK_VERSION = "2.0.0"


def read_version(base_dir: Optional[str] = None) -> str:
    """Version string: autodocker/VERSION (``base_dir``/VERSION) if present,
    else $AUTODOCKER_VERSION, else the version this file was released with."""
    base_dir = base_dir or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(base_dir, "VERSION")
    try:
        with open(path, encoding="utf-8") as f:
            version = f.read().strip()
        if version:
            return version
    except OSError:
        pass
    return os.environ.get("AUTODOCKER_VERSION", "").strip() or FALLBACK_VERSION


__version__ = read_version()
