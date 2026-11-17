"""Minimal FASTA reading and writing.

Deliberately dependency-free and deliberately strict about headers: verifying that
the record we got back carries the accession *and version* we asked for is the only
check standing between a time-correct database and a quietly modernised one.
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

LINE_WIDTH = 70

_HEADER_ACCESSION_RE = re.compile(r"^\s*([A-Za-z0-9_.]+)")


class FastaError(ValueError):
    pass


@dataclass(frozen=True)
class Record:
    header: str
    sequence: str

    @property
    def accession(self) -> str:
        return accession_from_header(self.header)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.sequence.encode("ascii", "replace")).hexdigest()


def accession_from_header(header: str) -> str:
    """Pull the accession out of a defline, with or without a leading '>'.

    NCBI deflines look like ``>NC_004718.3 SARS coronavirus Tor2, complete genome``.
    """
    text = header[1:] if header.startswith(">") else header
    match = _HEADER_ACCESSION_RE.match(text)
    if not match:
        raise FastaError(f"cannot read an accession from defline: {header!r}")
    return match.group(1)


def parse(text: str) -> Iterator[Record]:
    header: str | None = None
    chunks: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield Record(header, "".join(chunks))
            header = line[1:]
            chunks = []
        else:
            if header is None:
                raise FastaError("sequence data before the first defline")
            chunks.append(line)
    if header is not None:
        yield Record(header, "".join(chunks))


def write(handle, records: Iterable[Record], width: int = LINE_WIDTH) -> int:
    count = 0
    for record in records:
        handle.write(f">{record.header}\n")
        sequence = record.sequence
        for start in range(0, len(sequence), width):
            handle.write(sequence[start : start + width])
            handle.write("\n")
        count += 1
    return count


def write_file(path: Path, records: Iterable[Record], width: int = LINE_WIDTH) -> int:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        return write(handle, records, width=width)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(Path(path), "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()
