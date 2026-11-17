import pytest
import requests

from hindsight import ncbi


class FakeResponse:
    def __init__(self, status_code=200, text="", content=b"", headers=None):
        self.status_code = status_code
        self.text = text
        self._content = content or text.encode()
        self.headers = requests.structures.CaseInsensitiveDict(headers or {})

    @property
    def content(self):
        return self._content

    def iter_content(self, chunk_size=1):
        for i in range(0, len(self._content), chunk_size):
            yield self._content[i : i + chunk_size]

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code}")

    def close(self):
        pass


class FakeSession:
    def __init__(self, responses):
        self.responses = list(responses)
        self.requests_made: list[tuple[str, str, dict]] = []
        self.headers: dict[str, str] = {}

    def request(self, method, url, **kwargs):
        self.requests_made.append((method, url, kwargs))
        value = self.responses.pop(0)
        if isinstance(value, Exception):
            raise value
        return value


@pytest.fixture(autouse=True)
def no_sleeping(monkeypatch):
    monkeypatch.setattr(ncbi.time, "sleep", lambda _s: None)


def client(responses, **kwargs):
    return ncbi.NcbiClient(session=FakeSession(responses), api_key=None, **kwargs)


def test_rate_limit_body_is_detected_despite_http_200():
    """E-utilities answers a throttled request with 200 and a JSON body.

    Trusting the status code here would write the string
    '{"error":"API rate limit exceeded"}' into the database as a sequence.
    """
    body = '{"error":"API rate limit exceeded","api-key":"1.2.3.4","count":"4","limit":"3"}'
    good = ">NC_004718.3 ok\nACGT\n"
    c = client([FakeResponse(text=body), FakeResponse(text=good)])
    assert c.efetch("nuccore", ["NC_004718.3"]) == good


def test_efetch_gives_up_and_raises():
    body = '{"error":"API rate limit exceeded"}'
    c = client([FakeResponse(text=body)] * 3, max_attempts=3)
    with pytest.raises(ncbi.NcbiError):
        c.efetch("nuccore", ["NC_004718.3"])


def test_efetch_with_no_ids_makes_no_request():
    session = FakeSession([])
    c = ncbi.NcbiClient(session=session)
    assert c.efetch("nuccore", []) == ""
    assert session.requests_made == []


def test_efetch_retries_a_server_error():
    good = ">x ok\nACGT\n"
    c = client([FakeResponse(status_code=503), FakeResponse(text=good)])
    assert c.efetch("nuccore", ["x"]) == good


def test_efetch_retries_an_empty_body():
    good = ">x ok\nACGT\n"
    c = client([FakeResponse(text="   "), FakeResponse(text=good)])
    assert c.efetch("nuccore", ["x"]) == good


def test_efetch_sends_the_api_key_when_given():
    session = FakeSession([FakeResponse(text=">x\nA\n")])
    c = ncbi.NcbiClient(session=session, api_key="secret")
    c.efetch("nuccore", ["x"])
    _, _, kwargs = session.requests_made[0]
    assert kwargs["data"]["api_key"] == "secret"
    assert kwargs["data"]["id"] == "x"


def test_request_retries_connection_errors_then_raises():
    c = client([requests.ConnectionError("reset")] * 3, max_attempts=3)
    with pytest.raises(ncbi.NcbiError, match="giving up"):
        c.get_text("https://example.invalid/x")


def test_download_resumes_from_a_partial_file(tmp_path):
    dest = tmp_path / "catalog.gz"
    part = tmp_path / "catalog.gz.part"
    part.write_bytes(b"AAAA")
    session = FakeSession(
        [
            FakeResponse(headers={"Content-Length": "8"}),  # HEAD
            FakeResponse(status_code=206, content=b"BBBB"),
        ]
    )
    c = ncbi.NcbiClient(session=session)
    result = c.download("https://example.test/catalog.gz", dest)

    assert dest.read_bytes() == b"AAAABBBB"
    assert result.resumed_from == 4
    assert result.bytes_written == 4
    _, _, kwargs = session.requests_made[1]
    assert kwargs["headers"]["Range"] == "bytes=4-"


def test_download_starts_over_when_the_server_ignores_the_range(tmp_path):
    dest = tmp_path / "f.gz"
    (tmp_path / "f.gz.part").write_bytes(b"AAAA")
    session = FakeSession(
        [
            FakeResponse(headers={"Content-Length": "4"}),
            FakeResponse(status_code=200, content=b"CCCC"),
        ]
    )
    c = ncbi.NcbiClient(session=session)
    c.download("https://example.test/f.gz", dest)
    assert dest.read_bytes() == b"CCCC", "a 200 must replace, not append"


def test_download_uses_a_complete_cached_file(tmp_path):
    dest = tmp_path / "f.gz"
    dest.write_bytes(b"DDDD")
    session = FakeSession([FakeResponse(headers={"Content-Length": "4"})])
    c = ncbi.NcbiClient(session=session)
    result = c.download("https://example.test/f.gz", dest)
    assert result.from_cache
    assert len(session.requests_made) == 1, "only the HEAD; no body transferred"


def test_download_rejects_a_short_transfer(tmp_path):
    session = FakeSession(
        [
            FakeResponse(headers={"Content-Length": "10"}),
            FakeResponse(status_code=200, content=b"short"),
        ]
    )
    c = ncbi.NcbiClient(session=session)
    with pytest.raises(ncbi.NcbiError, match="short download"):
        c.download("https://example.test/f.gz", tmp_path / "f.gz")


def test_rate_limiter_spaces_requests(monkeypatch):
    slept: list[float] = []
    monkeypatch.setattr(ncbi.time, "sleep", slept.append)
    clock = iter([0.0, 0.0, 0.0, 0.0])
    monkeypatch.setattr(ncbi.time, "monotonic", lambda: next(clock))
    limiter = ncbi._RateLimiter(2.0)
    limiter.wait()
    limiter.wait()
    assert slept and slept[-1] == pytest.approx(0.5)


def test_rate_defaults_depend_on_the_api_key():
    assert ncbi.NcbiClient(api_key=None)._limiter._interval == pytest.approx(
        1 / ncbi.RATE_WITHOUT_KEY
    )
    assert ncbi.NcbiClient(api_key="k")._limiter._interval == pytest.approx(1 / ncbi.RATE_WITH_KEY)


def test_looks_like_rate_limit():
    assert ncbi._looks_like_rate_limit('{"error":"API rate limit exceeded"}')
    assert not ncbi._looks_like_rate_limit(">NC_004718.3 rate limit in a defline\nACGT\n")
