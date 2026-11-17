"""The SRA mirror path: locating a run, and refusing to convert what will not fit."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from hindsight.samples import sra


class FakeResponse:
    def __init__(self, payload, status=200):
        self._payload = payload
        self.status_code = status

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        return json.loads(json.dumps(self._payload))


class FakeSession:
    def __init__(self, payload, status=200):
        self.payload = payload
        self.status = status
        self.calls = []

    def get(self, url, params=None, timeout=None):
        self.calls.append((url, params))
        return FakeResponse(self.payload, self.status)


def sdl_payload(accession="SRR1", size=1000, md5="abc", service="s3"):
    return {
        "version": "2",
        "result": [
            {
                "bundle": accession,
                "status": 200,
                "msg": "ok",
                "files": [
                    {
                        "object": f"srapub|{accession}",
                        "accession": accession,
                        "type": "sra",
                        "name": accession,
                        "size": size,
                        "md5": md5,
                        "locations": [
                            {
                                "service": service,
                                "region": "us-east-1",
                                "link": f"https://example/{accession}",
                            }
                        ],
                    }
                ],
            }
        ],
    }


def test_locate_returns_the_virtual_host_style_odp_url():
    session = FakeSession(sdl_payload("SRR10971381", size=2743429212, md5="5496"))
    obj = sra.locate("SRR10971381", session=session)
    assert obj.size == 2743429212
    assert obj.md5 == "5496"
    assert obj.service == "sra-pub-run-odp"
    # Path style invites ambient AWS credentials into a request to a public bucket and
    # is then refused with AccessDenied; virtual-host style stays anonymous.
    assert obj.url == "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR10971381/SRR10971381"
    assert "s3.us-east-1.amazonaws.com" not in obj.url


def test_locate_reports_an_empty_result():
    with pytest.raises(RuntimeError, match="no result"):
        sra.locate("SRRnope", session=FakeSession({"version": "2", "result": []}))


def test_locate_reports_a_refusal_status():
    payload = sdl_payload()
    payload["result"][0]["status"] = 404
    payload["result"][0]["msg"] = "not found"
    with pytest.raises(RuntimeError, match="refused"):
        sra.locate("SRR1", session=FakeSession(payload))


def test_locate_requires_an_s3_hosted_object():
    with pytest.raises(RuntimeError, match="no S3-hosted"):
        sra.locate("SRR1", session=FakeSession(sdl_payload(service="gs")))


def test_conversion_space_rule():
    gib = 2**30
    # The index sample: a 2.7 GiB container needs far more than the 23 GiB that was
    # free on the machine this was built on, which is why it came from ENA instead.
    assert not sra.conversion_fits(int(2.7 * gib), int(23 * gib))
    assert sra.conversion_fits(int(2.7 * gib), int(40 * gib))
    # A typical negative: 0.15 GiB container, trivially convertible.
    assert sra.conversion_fits(int(0.15 * gib), int(23 * gib))


def test_conversion_space_factor_is_explicit():
    assert sra.CONVERSION_SPACE_FACTOR == 14
    assert sra.conversion_fits(10, 140)
    assert not sra.conversion_fits(10, 139)


def test_md5_of_a_known_file(tmp_path):
    path = tmp_path / "x"
    path.write_bytes(b"hindsight")
    assert sra.md5_file(path) == "2b85162ea1dc28002cc8fbebf051acd7"


def test_conversion_refuses_without_the_toolkit(tmp_path, monkeypatch):
    monkeypatch.setattr(sra.shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match="not on PATH"):
        sra.fasterq_dump_paired(tmp_path / "x.sra", tmp_path / "out", accession="SRR1", threads=1)


def test_conversion_is_skipped_when_the_pair_is_already_there(tmp_path, monkeypatch):
    out = tmp_path / "out"
    out.mkdir()
    (out / "SRR1_1.fastq.gz").write_bytes(b"")
    (out / "SRR1_2.fastq.gz").write_bytes(b"")

    def explode(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("fasterq-dump should not run when the pair exists")

    monkeypatch.setattr(sra.subprocess, "run", explode)
    r1, r2 = sra.fasterq_dump_paired(tmp_path / "x.sra", out, accession="SRR1", threads=1)
    assert r1.name == "SRR1_1.fastq.gz" and r2.name == "SRR1_2.fastq.gz"


def test_conversion_creates_its_own_log_directory(tmp_path, monkeypatch):
    """A sample whose first stage is the conversion has no logs/ directory yet.

    It cost an hour to learn that this surfaced as "SRA mirror unavailable" and a
    fallback to ENA at 1 MB/s, rather than as the missing directory it was.
    """
    out = tmp_path / "raw"
    monkeypatch.setattr(sra.shutil, "which", lambda name: f"/usr/bin/{name}")

    def fake_run(argv, **kwargs):
        if argv[0] == "fasterq-dump":
            for mate in (1, 2):
                (out / f"SRR1_{mate}.fastq").write_text("@a\nACGT\n+\nIIII\n")
        else:  # pigz
            target = Path(argv[-1])
            target.with_suffix(".fastq.gz").write_bytes(b"")
            target.unlink()
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(sra.subprocess, "run", fake_run)
    r1, r2 = sra.fasterq_dump_paired(
        tmp_path / "x.sra",
        out,
        accession="SRR1",
        threads=1,
        log_path=tmp_path / "logs" / "fasterq-dump.log",
    )
    assert r1.exists() and r2.exists()
    assert (tmp_path / "logs" / "fasterq-dump.log").exists()


def test_a_half_finished_pair_is_cleared_before_converting(tmp_path, monkeypatch):
    """One mate present means an attempt stopped partway, and pigz will not overwrite.

    Left in place it failed the conversion at its last step, which the fetch logic
    reported as the mirror being unavailable and answered with a 700 MB download from
    ENA at 1 MB/s.
    """
    out = tmp_path / "raw"
    out.mkdir()
    (out / "SRR1_1.fastq.gz").write_bytes(b"half a download")
    monkeypatch.setattr(sra.shutil, "which", lambda name: f"/usr/bin/{name}")

    def fake_run(argv, **kwargs):
        if argv[0] == "fasterq-dump":
            for mate in (1, 2):
                (out / f"SRR1_{mate}.fastq").write_text("@a\nACGT\n+\nIIII\n")
        else:
            target = Path(argv[-1])
            target.with_suffix(".fastq.gz").write_bytes(b"freshly converted")
            target.unlink()
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(sra.subprocess, "run", fake_run)
    r1, _ = sra.fasterq_dump_paired(tmp_path / "x.sra", out, accession="SRR1", threads=1)
    assert r1.read_bytes() == b"freshly converted"
