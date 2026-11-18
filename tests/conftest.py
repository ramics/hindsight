"""Fixtures shared across the test modules.

The fake dumper lives here rather than in one test module because two modules need
it, and importing a fixture from another test module needs ``tests`` to be a package
— which it is not, and making it one to share one fixture is the wrong trade.
"""

from __future__ import annotations

import gzip
import subprocess

import pytest

from hindsight.arms import diamond
from hindsight.samples import sra


def _blast_row(contig, subject, taxid, bitscore, *, pident=90.0, length=100, qlen=30000):
    values = {
        "qseqid": contig,
        "sseqid": subject,
        "staxids": str(taxid),
        "pident": f"{pident}",
        "length": str(length),
        "mismatch": "0",
        "gapopen": "0",
        "qstart": "1",
        "qend": "300",
        "sstart": "1",
        "send": "100",
        "evalue": "1e-50",
        "bitscore": f"{bitscore}",
        "qlen": str(qlen),
        "slen": "1000",
    }
    return "\t".join(values[field] for field in diamond.BLAST_FIELDS)


def _write_blast_hits(path, rows):
    with gzip.open(path, "wt") as handle:
        for row in rows:
            handle.write(row + "\n")
    return path


@pytest.fixture
def blast_row():
    """One row of DIAMOND tabular output, in the column order the arm reads.

    A fixture rather than a plain function because ``tests`` is deliberately not a
    package, so one test module cannot import a helper from another. Built from
    ``diamond.BLAST_FIELDS`` rather than by writing the columns out in order, so a change
    to that constant moves the tests with it instead of silently making them assert
    about the wrong column.
    """
    return _blast_row


@pytest.fixture
def write_blast_hits():
    """Write hit rows to a gzipped tabular file, the way the arm keeps them."""
    return _write_blast_hits


@pytest.fixture
def fake_dumper(monkeypatch, tmp_path):
    """Replace fastq-dump and pigz with a dumper that serves a fixed pool of spots."""
    state = {"total": 7, "calls": []}
    real_run = subprocess.run

    def fake_run(argv, **kwargs):
        if argv[0] == "fastq-dump":
            start = int(argv[argv.index("-N") + 1])
            stop = int(argv[argv.index("-X") + 1])
            out = argv[argv.index("-O") + 1]
            sra_path = argv[-1]
            stem = sra_path.rsplit("/", 1)[-1].rsplit(".", 1)[0]
            spots = [s for s in range(start, stop + 1) if s <= state["total"]]
            state["calls"].append((start, stop, len(spots)))
            for mate in (1, 2):
                target = f"{out}/{stem}_{mate}.fastq"
                with open(target, "w") as handle:
                    for spot in spots:
                        handle.write(f"@spot{spot}/{mate}\nACGT\n+\nIIII\n")
            return subprocess.CompletedProcess(argv, 0)
        if argv[0] == "pigz":
            source = argv[-1]
            with open(source, "rb") as handle:
                data = handle.read()
            kwargs["stdout"].write(gzip.compress(data))
            return subprocess.CompletedProcess(argv, 0)
        return real_run(argv, **kwargs)

    monkeypatch.setattr(sra.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(sra.subprocess, "run", fake_run)
    return state
