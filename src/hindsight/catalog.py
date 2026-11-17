"""Reading an archived RefSeq release catalog and turning it into a viral manifest.

The catalog is a 6-column TSV, gzipped, around 1.3 GB compressed for a 2019-era
release:

    1. taxonomy ID
    2. species name
    3. accession.version
    4. release directories the accession appears in, '|' delimited
    5. RefSeq status
    6. sequence length

Column 4 is how we pick the viral set, and column 3 is why the whole exercise works:
it names the exact revision of every record in that release.
"""

from __future__ import annotations

import csv
import gzip
import io
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from .accessions import AccessionError, Molecule, classify, split_version
from .ncbi import NcbiClient
from .releases import Release

CATALOG_COLUMNS = 6
VIRAL_DIRECTORY = "viral"

MANIFEST_HEADER = [
    "accession",
    "base_accession",
    "version",
    "molecule",
    "taxid",
    "organism",
    "refseq_status",
    "length",
    "directories",
]


class CatalogError(RuntimeError):
    pass


@dataclass(frozen=True)
class CatalogRow:
    taxid: int
    organism: str
    accession: str
    directories: tuple[str, ...]
    refseq_status: str
    length: int

    @property
    def base_accession(self) -> str:
        return split_version(self.accession)[0]

    @property
    def version(self) -> int:
        return split_version(self.accession)[1]

    @property
    def molecule(self) -> Molecule:
        return classify(self.accession)


@dataclass
class ManifestStats:
    rows_scanned: int = 0
    rows_viral: int = 0
    rows_malformed: int = 0
    rows_unclassified: int = 0
    nucleotide: int = 0
    protein: int = 0
    nucleotide_bases: int = 0
    protein_residues: int = 0
    taxids: int = 0

    def as_dict(self) -> dict[str, int]:
        return dict(vars(self))


def parse_row(line: str) -> CatalogRow:
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) != CATALOG_COLUMNS:
        raise CatalogError(f"expected {CATALOG_COLUMNS} columns, got {len(fields)}: {line!r}")
    taxid, organism, accession, directories, status, length = fields
    return CatalogRow(
        taxid=int(taxid),
        organism=organism,
        accession=accession,
        directories=tuple(d for d in directories.split("|") if d),
        refseq_status=status,
        length=int(length),
    )


def is_viral(row: CatalogRow) -> bool:
    return VIRAL_DIRECTORY in row.directories


def iter_viral_rows(
    lines: Iterable[str], stats: ManifestStats | None = None
) -> Iterator[CatalogRow]:
    """Yield the viral rows of a catalog, counting what was skipped and why.

    Malformed rows are counted rather than fatal: a single corrupt line in a 200
    million row catalog should be reported in the provenance, not abort a two-hour
    build. The caller decides what counts as too many.
    """
    stats = stats if stats is not None else ManifestStats()
    for line in lines:
        if not line.strip():
            continue
        stats.rows_scanned += 1
        try:
            row = parse_row(line)
        except (CatalogError, ValueError):
            stats.rows_malformed += 1
            continue
        if not is_viral(row):
            continue
        stats.rows_viral += 1
        try:
            _ = row.molecule
        except AccessionError:
            stats.rows_unclassified += 1
            continue
        yield row


def catalog_path(cache_dir: Path, release: Release) -> Path:
    return Path(cache_dir) / "catalogs" / release.catalog_filename


def resolve_catalog_url(client: NcbiClient, release: Release) -> str:
    """The URL a release's catalog actually lives at.

    Every superseded release is under ``archive/``; the current one sits in the
    release-catalog directory itself and only moves when the next release lands. We
    check rather than guess, so a build run the week of a release does not 404.
    """
    archived = release.catalog_url(archived=True)
    try:
        client.head(archived)
        return archived
    except Exception:  # noqa: BLE001 - any failure here means "try the live directory"
        return release.catalog_url(archived=False)


def download_catalog(
    client: NcbiClient,
    release: Release,
    cache_dir: Path,
    *,
    progress: bool = False,
    url: str | None = None,
) -> tuple[Path, str]:
    """Download the release catalog into the cache, resuming a partial transfer."""
    dest = catalog_path(cache_dir, release)
    resolved = url or resolve_catalog_url(client, release)
    client.download(resolved, dest, progress=progress)
    return dest, resolved


def open_catalog_lines(path: Path) -> Iterator[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        yield from handle


def stream_catalog_lines(
    client: NcbiClient, release: Release, url: str | None = None
) -> Iterator[str]:
    """Decompress the catalog straight off the wire, never touching disk.

    Cheaper on disk, but a dropped connection means starting over, so this is opt-in.
    """
    raw = _ChunkReader(client.stream_lines(url or resolve_catalog_url(client, release)))
    with gzip.open(raw, "rt", encoding="utf-8", errors="replace") as handle:
        yield from handle


class _ChunkReader(io.RawIOBase):
    """Adapt an iterator of byte chunks to a readable file object."""

    def __init__(self, chunks: Iterator[bytes]) -> None:
        self._chunks = chunks
        self._buffer = b""

    def readable(self) -> bool:
        return True

    def readinto(self, target) -> int:  # type: ignore[override]
        while not self._buffer:
            try:
                self._buffer = next(self._chunks)
            except StopIteration:
                return 0
        take = min(len(target), len(self._buffer))
        target[:take] = self._buffer[:take]
        self._buffer = self._buffer[take:]
        return take


def manifest_sort_key(row: CatalogRow) -> tuple[str, str, int]:
    return (row.molecule.value, row.base_accession, row.version)


def write_manifest(
    rows: Iterable[CatalogRow],
    out_path: Path,
    stats: ManifestStats | None = None,
) -> ManifestStats:
    """Write the viral manifest, sorted so the output is byte-identical between runs.

    The whole viral set for a 2019 release is a few hundred thousand rows, so sorting
    in memory is fine and keeps the output deterministic without a temp file dance.
    """
    stats = stats if stats is not None else ManifestStats()
    collected = sorted(rows, key=manifest_sort_key)
    seen_taxids: set[int] = set()

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(MANIFEST_HEADER)
        for row in collected:
            molecule = row.molecule
            if molecule is Molecule.NUCLEOTIDE:
                stats.nucleotide += 1
                stats.nucleotide_bases += row.length
            else:
                stats.protein += 1
                stats.protein_residues += row.length
            seen_taxids.add(row.taxid)
            writer.writerow(
                [
                    row.accession,
                    row.base_accession,
                    row.version,
                    molecule.value,
                    row.taxid,
                    row.organism,
                    row.refseq_status,
                    row.length,
                    "|".join(row.directories),
                ]
            )
    stats.taxids = len(seen_taxids)
    return stats


def read_manifest(path: Path) -> list[dict[str, str]]:
    with open(Path(path), newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(MANIFEST_HEADER) - set(reader.fieldnames or [])
        if missing:
            raise CatalogError(f"manifest {path} is missing columns: {sorted(missing)}")
        return list(reader)


def manifest_accessions(path: Path, molecule: Molecule) -> list[str]:
    return [row["accession"] for row in read_manifest(path) if row["molecule"] == molecule.value]
