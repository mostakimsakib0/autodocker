#!/usr/bin/env python3
"""HTTP FETCHING.

Split out of runner.py (behavior-preserving refactor).
"""
import logging

try:
    import requests
except ModuleNotFoundError:
    requests = None

try:
    import urllib.request
except ImportError:  # pragma: no cover
    urllib = None

logger = logging.getLogger(__name__)


def _http_get_bytes(url: str, timeout: int = 30) -> bytes:
    """Fetch a URL, using requests when available and urllib otherwise."""
    if requests is not None:
        resp = requests.get(url, timeout=timeout)
        resp.raise_for_status()
        return resp.content
    if urllib is None:
        raise RuntimeError("No HTTP library available (requests/urllib)")
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read()
