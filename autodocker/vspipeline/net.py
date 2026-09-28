#!/usr/bin/env python3
"""HTTP FETCHING (with retry/backoff).

Split out of runner.py (behavior-preserving refactor).
"""
import logging
import time
from typing import Optional

try:
    import requests
except ModuleNotFoundError:
    requests = None

try:
    import urllib.request
    import urllib.error
except ImportError:  # pragma: no cover
    urllib = None

logger = logging.getLogger(__name__)

# Retry policy: transient network failures (connection resets, timeouts,
# HTTP 429/5xx) are retried with exponential backoff. Non-transient HTTP
# errors (4xx other than 429) fail immediately.
DEFAULT_RETRIES = 3
DEFAULT_BACKOFF = 1.0  # seconds; doubled per attempt
RETRYABLE_STATUS = {429, 500, 502, 503, 504}


class _RetryableHTTPError(Exception):
    """Internal marker for a transient failure worth retrying."""


def _http_get_once(url: str, timeout: int = 30) -> bytes:
    """Single fetch attempt, using requests when available and urllib otherwise."""
    if requests is not None:
        resp = requests.get(url, timeout=timeout)
        if resp.status_code in RETRYABLE_STATUS:
            raise _RetryableHTTPError(f"HTTP {resp.status_code} from {url}")
        resp.raise_for_status()
        return resp.content
    if urllib is None:
        raise RuntimeError("No HTTP library available (requests/urllib)")
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return resp.read()
    except urllib.error.HTTPError as e:
        if e.code in RETRYABLE_STATUS:
            raise _RetryableHTTPError(f"HTTP {e.code} from {url}") from e
        raise


def _http_get_bytes(url: str, timeout: int = 30,
                    retries: Optional[int] = None,
                    backoff: float = DEFAULT_BACKOFF) -> bytes:
    """Fetch a URL with retries and exponential backoff.

    Transient failures (connection errors, timeouts, HTTP 429/5xx) are
    retried up to ``retries`` times (default: ``DEFAULT_RETRIES``) with a
    backoff that doubles between attempts. Permanent client errors (e.g.
    404) raise immediately.
    """
    attempts = DEFAULT_RETRIES if retries is None else max(1, retries)
    last_error: Optional[BaseException] = None
    for attempt in range(1, attempts + 1):
        try:
            return _http_get_once(url, timeout=timeout)
        except _RetryableHTTPError as e:
            last_error = e
        except (OSError, RuntimeError) as e:
            # requests raises requests.exceptions.* subclasses of OSError;
            # urllib raises URLError (subclass of OSError).
            last_error = e
        if attempt < attempts:
            sleep_for = backoff * (2 ** (attempt - 1))
            logger.debug(
                f"HTTP attempt {attempt}/{attempts} failed for {url}: {last_error}; "
                f"retrying in {sleep_for:.1f}s")
            time.sleep(sleep_for)
    assert last_error is not None
    raise last_error