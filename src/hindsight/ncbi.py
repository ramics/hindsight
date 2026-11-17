"""HTTP access to NCBI: the FTP archive over HTTPS, and E-utilities.

Two things here are load-bearing for correctness rather than convenience:

* E-utilities answers a rate-limit rejection with HTTP 200 and a JSON body, so the
  status code alone is not enough to tell a good response from a throttled one.
* Large archive files are fetched with HTTP Range resume, because a 1.3 GB catalog
  download that dies at 90% should not start over.
"""

from __future__ import annotations

import os
import random
import threading
import time
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import requests

from . import USER_AGENT

FTP_BASE = "https://ftp.ncbi.nlm.nih.gov"
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# E-utilities allows 3 requests/second anonymously and 10 with an API key. We stay
# under both: the published ceiling is per-IP, and on shared egress it is easy to be
# the one who trips it.
RATE_WITHOUT_KEY = 2.0
RATE_WITH_KEY = 8.0

RETRY_STATUS = frozenset({408, 425, 429, 500, 502, 503, 504})


class NcbiError(RuntimeError):
    pass


class RateLimitedError(NcbiError):
    """E-utilities rejected the request for exceeding the per-IP rate limit."""


class _RateLimiter:
    """Minimum-interval limiter, safe to share across threads."""

    def __init__(self, per_second: float) -> None:
        self._interval = 1.0 / per_second if per_second > 0 else 0.0
        self._lock = threading.Lock()
        self._next = 0.0

    def wait(self) -> None:
        if not self._interval:
            return
        with self._lock:
            now = time.monotonic()
            due = max(now, self._next)
            self._next = due + self._interval
        delay = due - now
        if delay > 0:
            time.sleep(delay)


@dataclass
class DownloadResult:
    path: Path
    bytes_written: int
    resumed_from: int
    from_cache: bool


class NcbiClient:
    """Thin client for the two NCBI surfaces hindsight needs."""

    def __init__(
        self,
        api_key: str | None = None,
        timeout: float = 120.0,
        max_attempts: int = 6,
        rate_per_second: float | None = None,
        session: requests.Session | None = None,
    ) -> None:
        self.api_key = api_key if api_key is not None else os.environ.get("NCBI_API_KEY") or None
        self.timeout = timeout
        self.max_attempts = max_attempts
        rate = rate_per_second
        if rate is None:
            rate = RATE_WITH_KEY if self.api_key else RATE_WITHOUT_KEY
        self._limiter = _RateLimiter(rate)
        self.session = session or requests.Session()
        self.session.headers["User-Agent"] = USER_AGENT

    # -- generic ---------------------------------------------------------------

    def _sleep_backoff(self, attempt: int) -> None:
        time.sleep(min(2.0**attempt, 60.0) * (0.5 + random.random()))

    def _request(self, method: str, url: str, **kwargs) -> requests.Response:
        last: Exception | None = None
        for attempt in range(self.max_attempts):
            try:
                response = self.session.request(method, url, timeout=self.timeout, **kwargs)
            except requests.RequestException as exc:  # connection reset, DNS, TLS
                last = exc
                self._sleep_backoff(attempt)
                continue
            if response.status_code in RETRY_STATUS:
                last = NcbiError(f"{response.status_code} from {url}")
                response.close()
                self._sleep_backoff(attempt)
                continue
            response.raise_for_status()
            return response
        raise NcbiError(f"giving up on {url} after {self.max_attempts} attempts") from last

    def get_text(self, url: str, max_bytes: int | None = None) -> str:
        response = self._request("GET", url)
        if max_bytes is None:
            return response.text
        return response.content[:max_bytes].decode("utf-8", "replace")

    def head(self, url: str) -> requests.structures.CaseInsensitiveDict:
        return self._request("HEAD", url).headers

    def download(
        self,
        url: str,
        dest: Path,
        *,
        resume: bool = True,
        expected_size: int | None = None,
        progress: bool = False,
        chunk_size: int = 1 << 20,
    ) -> DownloadResult:
        """Download ``url`` to ``dest``, resuming a partial transfer when possible.

        A completed download is left alone if its size matches the server's, so this
        doubles as the cache check.
        """
        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        part = dest.with_suffix(dest.suffix + ".part")

        if expected_size is None:
            try:
                expected_size = int(self.head(url).get("Content-Length", "")) or None
            except (NcbiError, ValueError, TypeError):
                expected_size = None

        if dest.exists() and (expected_size is None or dest.stat().st_size == expected_size):
            return DownloadResult(dest, dest.stat().st_size, 0, from_cache=True)

        offset = part.stat().st_size if (resume and part.exists()) else 0
        if offset and expected_size and offset >= expected_size:
            offset = 0
            part.unlink()

        headers = {"Range": f"bytes={offset}-"} if offset else {}
        response = self._request("GET", url, headers=headers, stream=True)
        if offset and response.status_code != 206:
            # Server ignored the range; start clean rather than concatenate garbage.
            offset = 0
            part.unlink(missing_ok=True)

        bar = None
        if progress:
            from tqdm import tqdm

            bar = tqdm(
                total=expected_size,
                initial=offset,
                unit="B",
                unit_scale=True,
                desc=dest.name,
            )
        written = 0
        mode = "ab" if offset else "wb"
        try:
            with open(part, mode) as handle:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if not chunk:
                        continue
                    handle.write(chunk)
                    written += len(chunk)
                    if bar is not None:
                        bar.update(len(chunk))
        finally:
            if bar is not None:
                bar.close()
            response.close()

        total = part.stat().st_size
        if expected_size is not None and total != expected_size:
            raise NcbiError(
                f"short download for {url}: got {total} bytes, expected {expected_size}"
            )
        part.replace(dest)
        return DownloadResult(dest, written, offset, from_cache=False)

    def stream_lines(self, url: str, chunk_size: int = 1 << 20) -> Iterator[bytes]:
        """Stream a URL as raw chunks without buffering it on disk."""
        response = self._request("GET", url, stream=True)
        try:
            yield from response.iter_content(chunk_size=chunk_size)
        finally:
            response.close()

    # -- E-utilities -----------------------------------------------------------

    def efetch(
        self,
        db: str,
        ids: list[str],
        rettype: str = "fasta",
        retmode: str = "text",
    ) -> str:
        """Fetch records by ``accession.version``.

        NCBI honours the version suffix, returning the record as it was at that
        revision. ``tests/test_fetch.py::test_version_pinning_is_honoured`` asserts
        this against the live service, because the whole leakage argument rests on it.
        """
        if not ids:
            return ""
        params = {
            "db": db,
            "id": ",".join(ids),
            "rettype": rettype,
            "retmode": retmode,
            "tool": "hindsight",
        }
        if self.api_key:
            params["api_key"] = self.api_key

        last: Exception | None = None
        for attempt in range(self.max_attempts):
            self._limiter.wait()
            try:
                response = self._request("POST", f"{EUTILS_BASE}/efetch.fcgi", data=params)
            except NcbiError as exc:
                last = exc
                continue
            text = response.text
            if _looks_like_rate_limit(text):
                last = RateLimitedError(text.strip()[:200])
                self._sleep_backoff(attempt)
                continue
            if not text.strip():
                last = NcbiError("empty response from efetch")
                self._sleep_backoff(attempt)
                continue
            return text
        raise NcbiError(f"efetch failed for {len(ids)} ids in db={db}") from last


def _looks_like_rate_limit(text: str) -> bool:
    head = text[:400].lstrip()
    return head.startswith("{") and "rate limit" in head.lower()


def ftp_url(path: str) -> str:
    return f"{FTP_BASE}/{path.lstrip('/')}"
