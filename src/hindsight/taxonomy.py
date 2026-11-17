"""Time-correct taxonomy, from NCBI's dated taxonomy dump archive.

A database built from a 2019 catalog but labelled with today's taxonomy is only
half time-correct: taxon names, ranks and the shape of the tree have all moved since,
and *Severe acute respiratory syndrome-related coronavirus* in particular gained a
child in 2020. NCBI archives monthly snapshots, so we take the last snapshot at or
before the release date and use that tree.

Coverage of the archive is finite: ``taxdmp`` snapshots begin 2014-08-01 and
``new_taxdump`` snapshots begin 2018-12-01. For a cutoff earlier than that there is
no time-correct taxonomy to be had, and hindsight says so rather than silently
reaching for a later tree.
"""

from __future__ import annotations

import csv
import datetime as dt
import re
import zipfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from .ncbi import NcbiClient, ftp_url

TAXDUMP_DIR = "pub/taxonomy/taxdump_archive"
REQUIRED_DUMPS = ("nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp")

FLAVOURS = ("taxdmp", "new_taxdump")


class TaxonomyError(RuntimeError):
    pass


@dataclass(frozen=True, order=True)
class TaxdumpSnapshot:
    date: dt.date
    flavour: str

    @property
    def filename(self) -> str:
        return f"{self.flavour}_{self.date.isoformat()}.zip"

    @property
    def url(self) -> str:
        return ftp_url(f"{TAXDUMP_DIR}/{self.filename}")

    def __str__(self) -> str:  # pragma: no cover - display only
        return f"{self.flavour} snapshot {self.date.isoformat()}"


def _href_re(flavour: str) -> re.Pattern[str]:
    return re.compile(rf'href="{re.escape(flavour)}_(\d{{4}}-\d{{2}}-\d{{2}})\.zip"')


def list_snapshots(client: NcbiClient, flavour: str = "taxdmp") -> list[TaxdumpSnapshot]:
    if flavour not in FLAVOURS:
        raise TaxonomyError(f"unknown taxdump flavour {flavour!r}; expected one of {FLAVOURS}")
    index = client.get_text(ftp_url(TAXDUMP_DIR + "/"))
    dates = sorted({dt.date.fromisoformat(d) for d in _href_re(flavour).findall(index)})
    if not dates:
        raise TaxonomyError(f"no {flavour} snapshots found in {TAXDUMP_DIR}")
    return [TaxdumpSnapshot(date=d, flavour=flavour) for d in dates]


def resolve_snapshot(
    release_date: dt.date, snapshots: Iterable[TaxdumpSnapshot]
) -> TaxdumpSnapshot:
    """The last snapshot at or on the release date.

    On-or-before is right here, unlike release selection: the snapshot describes the
    tree as it stood on its own date, and a snapshot taken the same morning as the
    release is contemporaneous with it rather than ahead of it.
    """
    eligible = [s for s in snapshots if s.date <= release_date]
    if not eligible:
        earliest = min(snapshots)
        raise TaxonomyError(
            f"no archived taxonomy snapshot at or before {release_date.isoformat()}; "
            f"the earliest {earliest.flavour} snapshot is {earliest.date.isoformat()}. "
            "Re-run with --no-taxonomy, or accept that taxonomy cannot be time-correct "
            "for this cutoff."
        )
    return max(eligible)


def download_snapshot(
    client: NcbiClient,
    snapshot: TaxdumpSnapshot,
    cache_dir: Path,
    *,
    progress: bool = False,
) -> Path:
    dest = Path(cache_dir) / "taxdumps" / snapshot.filename
    client.download(snapshot.url, dest, progress=progress)
    return dest


def extract_snapshot(
    archive: Path, out_dir: Path, members: Iterable[str] = REQUIRED_DUMPS
) -> list[Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    with zipfile.ZipFile(archive) as zf:
        available = set(zf.namelist())
        wanted = list(members)
        missing = [m for m in wanted if m not in available]
        if missing:
            raise TaxonomyError(f"{archive.name} is missing expected members: {missing}")
        for member in wanted:
            target = out_dir / member
            with zf.open(member) as src, open(target, "wb") as dst:
                while chunk := src.read(1 << 20):
                    dst.write(chunk)
            written.append(target)
    return written


def read_taxids(nodes_dmp: Path) -> set[int]:
    """Every taxid present in a dated ``nodes.dmp``."""
    taxids: set[int] = set()
    with open(Path(nodes_dmp), encoding="utf-8", errors="replace") as handle:
        for line in handle:
            field, _, _ = line.partition("\t")
            if field.strip():
                taxids.add(int(field.strip()))
    return taxids


def read_deleted(delnodes_dmp: Path) -> set[int]:
    deleted: set[int] = set()
    with open(Path(delnodes_dmp), encoding="utf-8", errors="replace") as handle:
        for line in handle:
            field = line.split("|", 1)[0].strip()
            if field:
                deleted.add(int(field))
    return deleted


def read_merged(merged_dmp: Path) -> dict[int, int]:
    mapping: dict[int, int] = {}
    with open(Path(merged_dmp), encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 2 and parts[0] and parts[1]:
                mapping[int(parts[0])] = int(parts[1])
    return mapping


def write_seqid_map(manifest_rows: Iterable[dict[str, str]], out_path: Path) -> int:
    """Write a Kraken2-style ``accession.version<TAB>taxid`` map."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(out_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in manifest_rows:
            writer.writerow([row["accession"], row["taxid"]])
            count += 1
    return count


@dataclass
class TaxonomyCheck:
    manifest_taxids: int
    present: int
    missing: list[int]
    resolvable_via_merged: dict[int, int]
    deleted_in_snapshot: list[int]

    def as_dict(self) -> dict:
        return {
            "manifest_taxids": self.manifest_taxids,
            "present_in_dated_nodes": self.present,
            "missing_count": len(self.missing),
            "missing": sorted(self.missing)[:200],
            "resolvable_via_merged": {
                str(k): v for k, v in sorted(self.resolvable_via_merged.items())
            },
            "deleted_in_snapshot": sorted(self.deleted_in_snapshot)[:200],
            "unexplained": sorted(
                set(self.missing) - set(self.resolvable_via_merged) - set(self.deleted_in_snapshot)
            )[:200],
        }


def check_taxids(
    manifest_rows: Iterable[dict[str, str]],
    nodes_dmp: Path,
    merged_dmp: Path | None = None,
    delnodes_dmp: Path | None = None,
) -> TaxonomyCheck:
    """Confirm the catalog's taxids exist in the contemporaneous tree.

    Mostly they do, because the catalog and the snapshot are from the same era, but
    not always, and the exceptions are informative rather than alarming. Release 97
    cites taxid 10658 (Enterobacteria phage PRD1), which the 2019-11-01 taxonomy
    already lists as deleted: NCBI's own two products disagree by a few days. We split
    the misses into merged, deleted and unexplained so the provenance says which kind
    it was, instead of leaving a bare count for someone to chase later.
    """
    wanted = {int(row["taxid"]) for row in manifest_rows}
    present = read_taxids(nodes_dmp)
    missing = sorted(wanted - present)
    merged = read_merged(merged_dmp) if merged_dmp and Path(merged_dmp).exists() else {}
    deleted = read_deleted(delnodes_dmp) if delnodes_dmp and Path(delnodes_dmp).exists() else set()
    return TaxonomyCheck(
        manifest_taxids=len(wanted),
        present=len(wanted) - len(missing),
        missing=missing,
        resolvable_via_merged={t: merged[t] for t in missing if t in merged},
        deleted_in_snapshot=[t for t in missing if t in deleted],
    )
