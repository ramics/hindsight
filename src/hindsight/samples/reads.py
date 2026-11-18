"""FASTQ streaming: counting, labelling, and pulling out an exact set of reads.

Everything here streams. The index sample is 28 million read pairs and the machine
this has to run on cheaply is a spot instance with a few gigabytes of RAM free, so no
function in this module holds more than one read in memory at a time. The one thing
kept in memory is a set of read *identifiers* — at the top dilution level that is
50,000 strings, which is nothing.
"""

from __future__ import annotations

import gzip
import hashlib
import io
from collections.abc import Iterator
from pathlib import Path


def open_maybe_gzip(path: Path, mode: str = "rt"):
    path = Path(path)
    # Returns the handle for the caller to manage; every caller in this package uses
    # it in a `with`.
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)  # noqa: SIM115


def normalise_id(raw: str) -> str:
    """The fragment identifier shared by both mates.

    FASTQ headers carry the mate in several dialects: ``@SRR1/1`` in the old Illumina
    style, ``@SRR1 1:N:0:TAG`` in Casava 1.8, and sometimes both. Alignment output
    (SAM column 1) strips the mate marker. Labelling reads from a BAM and then
    subsetting a FASTQ only works if both sides agree on what the identifier is, so
    everything in this package passes through here.
    """
    name = raw[1:] if raw.startswith("@") else raw
    name = name.split()[0] if name.split() else name
    if name.endswith(("/1", "/2")):
        name = name[:-2]
    return name


def iter_records(path: Path) -> Iterator[tuple[str, str, str, str]]:
    """Yield ``(header, sequence, plus, quality)`` for each record, verbatim.

    Lines are returned with their newlines intact so a subset can be written back out
    byte-identically to the input, which keeps the dilution files comparable to the
    source file rather than merely equivalent to it.
    """
    with open_maybe_gzip(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"{path}: truncated FASTQ record at {header.strip()!r}")
            if not header.startswith("@"):
                raise ValueError(f"{path}: expected a header, got {header[:40]!r}")
            yield header, seq, plus, qual


def count_reads(path: Path) -> int:
    """Count records in a FASTQ, by reading it. No division of file sizes by four."""
    total = 0
    with open_maybe_gzip(path, "rb") as handle:
        buffered = io.BufferedReader(handle, buffer_size=1 << 20)
        for _ in buffered:
            total += 1
    if total % 4:
        raise ValueError(f"{path}: {total} lines is not a whole number of FASTQ records")
    return total // 4


def read_ids(path: Path) -> set[str]:
    """Every fragment identifier in a FASTQ."""
    return {normalise_id(header) for header, _, _, _ in iter_records(path)}


def ids_from_sam(path: Path) -> set[str]:
    """Fragment identifiers from a SAM/TSV listing produced by ``samtools view``.

    The caller is expected to have filtered the flags already (mapped reads only, and
    whatever mapping-quality threshold it wants). This reads column one.
    """
    ids: set[str] = set()
    with open_maybe_gzip(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("@"):
                continue
            first = line.split("\t", 1)[0]
            if first:
                ids.add(normalise_id(first))
    return ids


def write_subset(
    source: Path,
    dest: Path,
    keep: set[str],
    *,
    compresslevel: int = 6,
) -> int:
    """Copy the records whose identifier is in ``keep``, in source order.

    Source order rather than selection order, so the output is a function of the input
    file and the identifier set alone — no dependence on how the set was built or
    iterated. Returns the number of records written.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with gzip.open(dest, "wt", compresslevel=compresslevel) as out:
        for header, seq, plus, qual in iter_records(source):
            if normalise_id(header) in keep:
                out.write(header)
                out.write(seq)
                out.write(plus)
                out.write(qual)
                written += 1
    return written


def write_excluding(
    source: Path,
    dest: Path,
    drop: set[str],
    *,
    compresslevel: int = 6,
) -> int:
    """Copy every record whose identifier is *not* in ``drop``, in source order."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with gzip.open(dest, "wt", compresslevel=compresslevel) as out:
        for header, seq, plus, qual in iter_records(source):
            if normalise_id(header) not in drop:
                out.write(header)
                out.write(seq)
                out.write(plus)
                out.write(qual)
                written += 1
    return written


def keeps_fragment(read_id: str, *, seed: int, fraction: float) -> bool:
    """Whether a seeded hash thins this fragment in or out.

    Used to cap the negatives at a fixed sequencing depth. A two-pass seeded sample
    would need every identifier in memory — 24 million strings for the largest
    negative — and taking the first N pairs instead would take them from the head of
    the file, which on an Illumina run is a corner of the flowcell rather than a
    random draw.

    So the decision is made per fragment from its identifier alone: keep it when the
    leading bytes of ``sha256(f"{seed}:{id}")`` fall below ``fraction`` of the range.
    One pass, no memory, no order dependence, and reproducible from the seed. The
    count that results is binomial around the target rather than exact, so the
    *actual* number kept is what gets recorded.
    """
    if fraction >= 1.0:
        return True
    if fraction <= 0.0:
        return False
    digest = hashlib.sha256(f"{seed}:{read_id}".encode()).digest()
    # 8 bytes is ample precision for a depth cap and avoids a big-int conversion per
    # read; at 24 million reads the quantisation error is far below the binomial
    # spread.
    value = int.from_bytes(digest[:8], "big")
    return value < fraction * (1 << 64)


def write_hash_fraction(
    source: Path,
    dest: Path,
    *,
    seed: int,
    fraction: float,
    compresslevel: int = 6,
) -> int:
    """Thin a FASTQ to approximately ``fraction`` of its fragments, deterministically.

    Both mates thinned with the same seed keep the same fragments, so the pair stays
    in step without either file having to know what the other did.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with gzip.open(dest, "wt", compresslevel=compresslevel) as out:
        for header, seq, plus, qual in iter_records(source):
            if keeps_fragment(normalise_id(header), seed=seed, fraction=fraction):
                out.write(header)
                out.write(seq)
                out.write(plus)
                out.write(qual)
                written += 1
    return written


def concatenate(parts: list[Path], dest: Path) -> None:
    """Concatenate gzip members into one file.

    Valid because the gzip format is concatenable: ``cat a.gz b.gz`` decompresses to
    the concatenation of the two streams. This keeps building a dilution level to one
    cheap byte copy of the background plus a small viral part, rather than
    recompressing five gigabytes per level.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as out:
        for part in parts:
            with open(part, "rb") as handle:
                while block := handle.read(1 << 22):
                    out.write(block)


def sha256_file(path: Path, *, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def count_redacted(path: Path) -> int:
    """Count fragments whose sequence is entirely N.

    These are not a quirk of one sample. NCBI redacts reads it identifies as human in
    submissions of human-derived material, and a redacted read is distributed as a
    full-length run of N at quality zero. In the index sample that is 49% of the run:
    the archive has already host-depleted the data that the benchmark is about.

    It matters because it decides what "reads per million" means. A detection limit
    quoted against the raw spot count is quoted against a denominator half of which is
    placeholder, and is therefore twice as flattering as it looks.
    """
    total = 0
    for _, seq, _, _ in iter_records(path):
        stripped = seq.strip()
        if stripped and stripped.count("N") == len(stripped):
            total += 1
    return total
