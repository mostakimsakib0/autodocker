#!/usr/bin/env python3
"""HTTP fetching with retry and exponential backoff (PubChem downloads)."""
import http.client
import logging
import time
import urllib.error
import urllib.request
from typing import Optional

try:
    import requests
except ModuleNotFoundError:  # the stdlib client below is used instead
    requests = None

logger = logging.getLogger(__name__)

# Per-request timeout in seconds (connect + read).
HTTP_TIMEOUT = 30

# Retry policy: transient network failures (connection resets, timeouts,
# HTTP 429/5xx) are retried with exponential backoff. Non-transient HTTP
# errors (4xx other than 429) fail immediately.
DEFAULT_RETRIES = 3
DEFAULT_BACKOFF = 1.0  # seconds; doubled per attempt
RETRYABLE_STATUS = {429, 500, 502, 503, 504}


class _RetryableHTTPError(Exception):
    """A transient HTTP failure (429/5xx); raised once retries are exhausted."""


def _http_status(error: BaseException) -> Optional[int]:
    """HTTP status of a permanent HTTP error (requests or urllib), else None."""
    status = getattr(error, "code", None)  # urllib.error.HTTPError
    response = getattr(error, "response", None)  # requests.HTTPError
    if status is None and response is not None:
        status = getattr(response, "status_code", None)
    if isinstance(status, int) and status not in RETRYABLE_STATUS:
        return status
    return None


def is_network_failure(error: BaseException) -> bool:
    """True if ``error`` means the server could not be reached or kept
    failing transiently (connection error, timeout, persistent 429/5xx),
    as opposed to an answer such as 404."""
    if isinstance(error, _RetryableHTTPError):
        return True
    return (isinstance(error, (OSError, http.client.HTTPException))
            and _http_status(error) is None)


def _http_get_once(url: str, timeout: float = HTTP_TIMEOUT) -> bytes:
    """Single fetch attempt, using requests when installed and urllib otherwise."""
    if requests is not None:
        resp = requests.get(url, timeout=timeout)
        if resp.status_code in RETRYABLE_STATUS:
            raise _RetryableHTTPError(f"HTTP {resp.status_code} from {url}")
        resp.raise_for_status()
        return resp.content
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return resp.read()
    except urllib.error.HTTPError as e:
        if e.code in RETRYABLE_STATUS:
            raise _RetryableHTTPError(f"HTTP {e.code} from {url}") from e
        raise


def _http_get_bytes(url: str, timeout: float = HTTP_TIMEOUT,
                    retries: Optional[int] = None,
                    backoff: float = DEFAULT_BACKOFF) -> bytes:
    """Fetch a URL with retries and exponential backoff.

    Transient failures (connection errors, timeouts, HTTP 429/5xx) are tried
    up to ``retries`` times in total (default: ``DEFAULT_RETRIES``) with a
    backoff that doubles between attempts. Permanent client errors (e.g.
    404) and programming errors raise immediately.
    """
    attempts = DEFAULT_RETRIES if retries is None else max(1, retries)
    last_error: Optional[BaseException] = None
    for attempt in range(1, attempts + 1):
        try:
            return _http_get_once(url, timeout=timeout)
        except (_RetryableHTTPError, OSError, http.client.HTTPException) as e:
            # requests' exceptions subclass OSError; urllib raises URLError
            # (an OSError). An HTTPError carrying a 4xx status other than 429
            # is permanent and is raised at once instead of being retried.
            if not is_network_failure(e):
                raise
            last_error = e
        if attempt < attempts:
            sleep_for = backoff * (2 ** (attempt - 1))
            logger.debug(
                f"HTTP attempt {attempt}/{attempts} failed for {url}: {last_error}; "
                f"retrying in {sleep_for:.1f}s")
            time.sleep(sleep_for)
    assert last_error is not None
    raise last_error
