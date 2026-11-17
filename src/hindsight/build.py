"""The end-to-end build: cutoff date in, time-correct database out."""

from __future__ import annotations

import datetime as dt
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from . import catalog as catalog_mod
from . import fetch as fetch_mod
from . import provenance as provenance_mod
from . import releases as releases_mod
from . import taxonomy as taxonomy_mod
from .accessions import Molecule
from .ncbi import NcbiClient

MANIFEST_NAME = "manifest.tsv"
NUCLEOTIDE_FASTA = "viral.fna"
PROTEIN_FASTA = "viral.faa"
SEQID_MAP = "seqid2taxid.map"
WITHDRAWN_NAME = "withdrawn.tsv"
TAXONOMY_SUBDIR = "taxonomy"


@dataclass
class BuildResult:
    out_dir: Path
    release: releases_mod.Release
    manifest_path: Path
    provenance_path: Path | None
    outputs: dict[str, dict]


def _echo(message: str, quiet: bool = False) -> None:
    if not quiet:
        print(message, file=sys.stderr, flush=True)


def build_manifest(
    client: NcbiClient,
    release: releases_mod.Release,
    out_dir: Path,
    cache_dir: Path,
    *,
    stream: bool = False,
    progress: bool = False,
    quiet: bool = False,
) -> tuple[Path, catalog_mod.ManifestStats, int | None, str]:
    """Resolve the viral accession set for a release and write the manifest."""
    out_dir = Path(out_dir)
    stats = catalog_mod.ManifestStats()
    catalog_bytes: int | None = None
    catalog_url = catalog_mod.resolve_catalog_url(client, release)

    if stream:
        _echo(f"streaming catalog for {release} (nothing cached)", quiet)
        lines = catalog_mod.stream_catalog_lines(client, release, catalog_url)
    else:
        _echo(f"fetching catalog for {release}", quiet)
        path, catalog_url = catalog_mod.download_catalog(
            client, release, cache_dir, progress=progress, url=catalog_url
        )
        catalog_bytes = path.stat().st_size
        lines = catalog_mod.open_catalog_lines(path)

    manifest_path = out_dir / MANIFEST_NAME
    rows = catalog_mod.iter_viral_rows(lines, stats)
    catalog_mod.write_manifest(rows, manifest_path, stats)
    _echo(
        f"manifest: {stats.nucleotide} nucleotide + {stats.protein} protein records "
        f"across {stats.taxids} taxa",
        quiet,
    )
    return manifest_path, stats, catalog_bytes, catalog_url


def _write_withdrawn(out_dir: Path, reports: Sequence[fetch_mod.FetchReport]) -> int:
    rows = [
        (report.molecule.value, accession)
        for report in reports
        for accession in sorted(report.withdrawn)
    ]
    path = Path(out_dir) / WITHDRAWN_NAME
    with open(path, "w") as handle:
        handle.write("molecule\taccession\n")
        for molecule, accession in rows:
            handle.write(f"{molecule}\t{accession}\n")
    return len(rows)


def build(
    client: NcbiClient,
    cutoff: dt.date,
    out_dir: Path,
    cache_dir: Path,
    *,
    molecules: Sequence[Molecule] = (Molecule.NUCLEOTIDE, Molecule.PROTEIN),
    releases: list[releases_mod.Release] | None = None,
    with_taxonomy: bool = True,
    taxdump_flavour: str = "taxdmp",
    batch_size: int | None = None,
    limit: int | None = None,
    stream_catalog: bool = False,
    manifest_only: bool = False,
    progress: bool = False,
    quiet: bool = False,
    command: list[str] | None = None,
) -> BuildResult:
    out_dir = Path(out_dir)
    cache_dir = Path(cache_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    table = releases if releases is not None else releases_mod.load_releases()
    release = releases_mod.resolve(cutoff, table)
    following = releases_mod.next_release_after(release, table)
    _echo(
        f"cutoff {cutoff.isoformat()} -> {release}; "
        f"first excluded release: {following or 'none (this is the latest)'}",
        quiet,
    )

    manifest_path, stats, catalog_bytes, catalog_url = build_manifest(
        client,
        release,
        out_dir,
        cache_dir,
        stream=stream_catalog,
        progress=progress,
        quiet=quiet,
    )

    if manifest_only:
        return BuildResult(out_dir, release, manifest_path, None, {})

    manifest_rows = catalog_mod.read_manifest(manifest_path)
    fetch_reports: list[fetch_mod.FetchReport] = []
    fasta_names: list[str] = []

    for molecule in molecules:
        accessions = [r["accession"] for r in manifest_rows if r["molecule"] == molecule.value]
        if limit is not None:
            accessions = accessions[:limit]
        if not accessions:
            continue
        shard_dir = cache_dir / "records" / f"release{release.number}" / molecule.value
        results, report = fetch_mod.fetch_molecule(
            client,
            molecule,
            accessions,
            shard_dir,
            batch_size=batch_size,
            progress=progress,
        )
        name = NUCLEOTIDE_FASTA if molecule is Molecule.NUCLEOTIDE else PROTEIN_FASTA
        written = fetch_mod.assemble(results, out_dir / name)
        fasta_names.append(name)
        fetch_reports.append(report)
        _echo(
            f"{molecule.value}: wrote {written} records to {name} "
            f"({len(report.withdrawn)} withdrawn, "
            f"{len(report.version_mismatches)} version mismatches)",
            quiet,
        )
        if report.version_mismatches:
            _echo(
                "  WARNING: NCBI returned records at versions we did not ask for; "
                "they were rejected. See provenance.json.",
                quiet,
            )

    _write_withdrawn(out_dir, fetch_reports)

    taxonomy_record: dict | None = None
    if with_taxonomy:
        snapshots = taxonomy_mod.list_snapshots(client, taxdump_flavour)
        snapshot = taxonomy_mod.resolve_snapshot(release.date, snapshots)
        _echo(f"taxonomy: {snapshot}", quiet)
        archive = taxonomy_mod.download_snapshot(client, snapshot, cache_dir, progress=progress)
        tax_dir = out_dir / TAXONOMY_SUBDIR
        taxonomy_mod.extract_snapshot(archive, tax_dir)
        mapped = taxonomy_mod.write_seqid_map(manifest_rows, out_dir / SEQID_MAP)
        check = taxonomy_mod.check_taxids(
            manifest_rows,
            tax_dir / "nodes.dmp",
            tax_dir / "merged.dmp",
            tax_dir / "delnodes.dmp",
        )
        taxonomy_record = {
            "snapshot_date": snapshot.date.isoformat(),
            "flavour": snapshot.flavour,
            "url": snapshot.url,
            "seqid2taxid_entries": mapped,
            "lag_days_behind_release": (release.date - snapshot.date).days,
            "check": check.as_dict(),
        }
        if check.missing:
            unexplained = (
                set(check.missing)
                - set(check.resolvable_via_merged)
                - set(check.deleted_in_snapshot)
            )
            _echo(
                f"  note: {len(check.missing)} catalog taxids are absent from the "
                f"{snapshot.date.isoformat()} tree "
                f"({len(check.deleted_in_snapshot)} already deleted, "
                f"{len(check.resolvable_via_merged)} merged, "
                f"{len(unexplained)} unexplained). See provenance.json.",
                quiet,
            )

    output_names = [MANIFEST_NAME, WITHDRAWN_NAME, *fasta_names]
    if with_taxonomy:
        output_names.append(SEQID_MAP)
        output_names.extend(f"{TAXONOMY_SUBDIR}/{m}" for m in taxonomy_mod.REQUIRED_DUMPS)
    outputs = provenance_mod.checksum_outputs(out_dir, output_names)

    record = provenance_mod.build_record(
        cutoff=cutoff,
        release=release,
        next_release=following,
        catalog_url=catalog_url,
        catalog_bytes=catalog_bytes,
        manifest_stats=stats.as_dict() | {"limit_applied": limit},
        fetch_reports=[r.as_dict() for r in fetch_reports],
        taxonomy=taxonomy_record,
        outputs=outputs,
        command=command,
    )
    provenance_path = provenance_mod.write(out_dir, record)
    _echo(f"wrote {provenance_path}", quiet)

    return BuildResult(out_dir, release, manifest_path, provenance_path, outputs)
