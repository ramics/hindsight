"""Chunked conversion: the resume bookkeeping and the stop condition.

``fastq-dump`` itself is not exercised here — these tests stand in a fake dumper that
writes a known number of spots per chunk, so the behaviour under test is the part that
can silently truncate a 28-million-read sample: when it stops, what it appends, and
what it does when a previous attempt died halfway.
"""

from __future__ import annotations

import gzip
import json

import pytest

from hindsight.samples import reads, sra


def test_chunked_conversion_covers_every_spot_once(fake_dumper, tmp_path):
    sra_path = tmp_path / "SRRX.sra"
    sra_path.write_bytes(b"")
    out = tmp_path / "raw"
    r1, r2 = sra.fastq_dump_chunked(sra_path, out, accession="SRRX", chunk_spots=3, threads=1)
    assert reads.count_reads(r1) == 7
    assert [reads.normalise_id(h) for h, *_ in reads.iter_records(r1)] == [
        f"spot{i}" for i in range(1, 8)
    ]
    assert reads.count_reads(r2) == 7
    # 3 + 3 + 1 spots, then one empty chunk that ends the loop.
    assert [c[2] for c in fake_dumper["calls"]] == [3, 3, 1, 0]


def test_concatenated_chunks_are_a_valid_single_gzip_stream(fake_dumper, tmp_path):
    sra_path = tmp_path / "SRRX.sra"
    sra_path.write_bytes(b"")
    r1, _ = sra.fastq_dump_chunked(
        sra_path, tmp_path / "raw", accession="SRRX", chunk_spots=2, threads=1
    )
    with gzip.open(r1, "rb") as handle:
        assert handle.read().count(b"@spot") == 7


def test_the_marker_records_the_spot_count_and_completion(fake_dumper, tmp_path):
    sra_path = tmp_path / "SRRX.sra"
    sra_path.write_bytes(b"")
    out = tmp_path / "raw"
    sra.fastq_dump_chunked(sra_path, out, accession="SRRX", chunk_spots=4, threads=1)
    marker = json.loads((out / ".SRRX.chunks.json").read_text())
    assert marker["spots"] == 7
    assert marker["complete"] is True


def test_a_finished_pair_is_not_converted_again(fake_dumper, tmp_path):
    out = tmp_path / "raw"
    out.mkdir()
    (out / "SRRX_1.fastq.gz").write_bytes(b"")
    (out / "SRRX_2.fastq.gz").write_bytes(b"")
    sra.fastq_dump_chunked(tmp_path / "SRRX.sra", out, accession="SRRX", threads=1)
    assert fake_dumper["calls"] == []


def test_bookkeeping_without_its_parts_starts_over(fake_dumper, tmp_path):
    """A marker that survived while its .part files did not must not be trusted.

    Resuming from spot 5 by appending to a file that no longer exists would produce a
    FASTQ missing its first four spots, with nothing to indicate it.
    """
    sra_path = tmp_path / "SRRX.sra"
    sra_path.write_bytes(b"")
    out = tmp_path / "raw"
    out.mkdir()
    (out / ".SRRX.chunks.json").write_text(json.dumps({"next_spot": 5, "spots": 4}))
    r1, _ = sra.fastq_dump_chunked(sra_path, out, accession="SRRX", chunk_spots=10, threads=1)
    assert reads.count_reads(r1) == 7
    assert fake_dumper["calls"][0][0] == 1


def test_a_container_with_no_spots_is_an_error(fake_dumper, tmp_path):
    fake_dumper["total"] = 0
    sra_path = tmp_path / "SRRX.sra"
    sra_path.write_bytes(b"")
    with pytest.raises(RuntimeError, match="no spots converted"):
        sra.fastq_dump_chunked(sra_path, tmp_path / "raw", accession="SRRX", threads=1)
