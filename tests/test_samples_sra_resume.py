"""Resuming a chunked conversion that was interrupted mid-append.

This is the failure that made the test: a conversion killed while appending a chunk
leaves a partial gzip member on the end of the part file. Without the byte lengths in
the marker, the resume appends a complete copy of the same chunk behind the truncated
one, and the result is a 2 GB FASTQ with duplicated reads that still opens.
"""

from __future__ import annotations

import gzip
import json

from hindsight.samples import reads, sra


def test_a_half_written_chunk_is_discarded_on_resume(fake_dumper, tmp_path):
    sra_path = tmp_path / "SRRX.sra"
    sra_path.write_bytes(b"")
    out = tmp_path / "raw"

    # Two chunks of three spots have been appended and recorded...
    fake_dumper["total"] = 6
    sra.fastq_dump_chunked(sra_path, out, accession="SRRX", chunk_spots=3, threads=1)
    r1 = out / "SRRX_1.fastq.gz"
    good1 = r1.read_bytes()
    good2 = (out / "SRRX_2.fastq.gz").read_bytes()
    marker = json.loads((out / ".SRRX.chunks.json").read_text())
    # ...and the seventh spot is still to come.
    fake_dumper["total"] = 7

    # ...and then the process died partway through appending the next one.
    part1 = out / "SRRX_1.fastq.gz.part"
    part2 = out / "SRRX_2.fastq.gz.part"
    truncated_at = {
        "next_spot": 7,
        "spots": 6,
        "part1_bytes": marker["part1_bytes"],
        "part2_bytes": marker["part2_bytes"],
    }
    part1.write_bytes(good1[: marker["part1_bytes"]] + gzip.compress(b"@spot6/1\nAC")[:20])
    part2.write_bytes(good2[: marker["part2_bytes"]] + gzip.compress(b"@spot6/2\nAC")[:20])
    r1.unlink()
    (out / "SRRX_2.fastq.gz").unlink()
    (out / ".SRRX.chunks.json").write_text(json.dumps(truncated_at))

    messages = []
    sra.fastq_dump_chunked(
        sra_path, out, accession="SRRX", chunk_spots=3, threads=1, log=messages.append
    )

    assert any("half-written chunk" in m for m in messages)
    ids = [reads.normalise_id(h) for h, *_ in reads.iter_records(r1)]
    assert ids == [f"spot{i}" for i in range(1, 8)]
    assert len(ids) == len(set(ids)), "no duplicated reads from the re-appended chunk"


def test_a_marker_without_part_lengths_is_not_trusted(fake_dumper, tmp_path):
    """An older marker carries no byte lengths, so its parts cannot be verified."""
    sra_path = tmp_path / "SRRX.sra"
    sra_path.write_bytes(b"")
    out = tmp_path / "raw"
    out.mkdir()
    (out / "SRRX_1.fastq.gz.part").write_bytes(gzip.compress(b"@spot1/1\nAC\n+\nII\n"))
    (out / "SRRX_2.fastq.gz.part").write_bytes(gzip.compress(b"@spot1/2\nAC\n+\nII\n"))
    (out / ".SRRX.chunks.json").write_text(json.dumps({"next_spot": 2, "spots": 1}))

    messages = []
    sra.fastq_dump_chunked(
        sra_path, out, accession="SRRX", chunk_spots=10, threads=1, log=messages.append
    )
    assert any("from the start" in m for m in messages)
    assert reads.count_reads(out / "SRRX_1.fastq.gz") == 7
