"""The RefSeq release table: which releases exist, and when each one was issued.

Picking the release is the whole time control, so it gets its own module and its own
vendored, reviewable table. Two facts about the archive shaped this:

* Release notes exist for every archived release, and state the release date in their
  header. That is the authoritative date, not a file modification time — archived
  catalogs are re-stamped when they are moved into ``archive/``, so release 97's
  catalog carries an mtime of January 2020 despite being the November 2019 release.
* RefSeq numbering jumps from 99 (March 2020) straight to 200 (July 2020). Releases
  100-199 do not exist. Any code that assumes contiguous numbering, or that infers a
  date by interpolating between release numbers, is wrong.
"""

from __future__ import annotations

import datetime as dt
import re
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

from .ncbi import NcbiClient, ftp_url

CATALOG_DIR = "refseq/release/release-catalog"
CATALOG_ARCHIVE_DIR = f"{CATALOG_DIR}/archive"
NOTES_DIR = "refseq/release/release-notes"
NOTES_ARCHIVE_DIR = f"{NOTES_DIR}/archive"

TABLE_RESOURCE = "release_dates.tsv"

_MONTHS = {
    m: i
    for i, m in enumerate(
        [
            "january",
            "february",
            "march",
            "april",
            "may",
            "june",
            "july",
            "august",
            "september",
            "october",
            "november",
            "december",
        ],
        start=1,
    )
}

# Tolerant on purpose: the notes contain "May 2 , 2004" and "January 23,2010".
_HEADER_RE = re.compile(
    r"Release\s+(?P<release>\d+)\s*\n\s*"
    r"(?P<month>[A-Z][a-z]+)\s+(?P<day>\d{1,2})\s*,?\s*(?P<year>\d{4})",
)
_CATALOG_HREF_RE = re.compile(r'href="RefSeq-release(\d+)\.catalog\.gz"')
_NOTES_HREF_RE = re.compile(r'href="RefSeq-release(\d+)\.txt"')


class ReleaseError(RuntimeError):
    pass


@dataclass(frozen=True, order=True)
class Release:
    date: dt.date
    number: int

    @property
    def catalog_filename(self) -> str:
        return f"RefSeq-release{self.number}.catalog.gz"

    @property
    def notes_filename(self) -> str:
        return f"RefSeq-release{self.number}.txt"

    def catalog_url(self, *, archived: bool = True) -> str:
        directory = CATALOG_ARCHIVE_DIR if archived else CATALOG_DIR
        return ftp_url(f"{directory}/{self.catalog_filename}")

    def notes_url(self, *, archived: bool = True) -> str:
        directory = NOTES_ARCHIVE_DIR if archived else NOTES_DIR
        return ftp_url(f"{directory}/{self.notes_filename}")

    def __str__(self) -> str:  # pragma: no cover - display only
        return f"release {self.number} ({self.date.isoformat()})"


# -- the vendored table --------------------------------------------------------


def parse_table(text: str) -> list[Release]:
    releases: list[Release] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            raise ReleaseError(f"malformed release table row: {raw!r}")
        releases.append(Release(number=int(fields[0]), date=dt.date.fromisoformat(fields[1])))
    if not releases:
        raise ReleaseError("release table is empty")
    return sorted(releases)


def format_table(releases: list[Release]) -> str:
    lines = [
        "# RefSeq release -> release date, parsed from the release notes header.",
        "# Regenerate with: hindsight refresh-releases",
        "# release\tdate",
    ]
    lines.extend(f"{r.number}\t{r.date.isoformat()}" for r in sorted(releases))
    return "\n".join(lines) + "\n"


def load_releases(path: Path | None = None) -> list[Release]:
    """Load the release table, from ``path`` if given or the vendored copy otherwise."""
    if path is not None:
        return parse_table(Path(path).read_text())
    text = resources.files("hindsight.data").joinpath(TABLE_RESOURCE).read_text()
    return parse_table(text)


def check_monotonic(releases: list[Release]) -> None:
    """Assert that release numbers and dates agree on their ordering.

    If NCBI ever back-dates a release, or our parser picks up the wrong date, the
    cheapest place to notice is here rather than in a database that silently contains
    the future.
    """
    by_number = sorted(releases, key=lambda r: r.number)
    for earlier, later in zip(by_number, by_number[1:], strict=False):
        if later.date <= earlier.date:
            raise ReleaseError(
                f"release {later.number} ({later.date}) is not after "
                f"release {earlier.number} ({earlier.date})"
            )


# -- resolution ----------------------------------------------------------------


def resolve(cutoff: dt.date, releases: list[Release] | None = None) -> Release:
    """Return the last release issued *strictly before* ``cutoff``.

    Strictly before, not on-or-before: a cutoff date is the first day that must not
    be visible to the database, and a release issued that morning is already inside
    the excluded window.
    """
    table = releases if releases is not None else load_releases()
    eligible = [r for r in table if r.date < cutoff]
    if not eligible:
        earliest = min(table, key=lambda r: r.date)
        raise ReleaseError(
            f"no RefSeq release precedes {cutoff.isoformat()}; "
            f"the earliest archived release is {earliest}"
        )
    return max(eligible)


def next_release_after(release: Release, releases: list[Release] | None = None) -> Release | None:
    """The release immediately following ``release``, for reporting the blind window."""
    table = releases if releases is not None else load_releases()
    later = [r for r in table if r.date > release.date]
    return min(later) if later else None


# -- refreshing from the archive ----------------------------------------------


def parse_month(name: str) -> int:
    """Resolve a month name, tolerating the typos in NCBI's own release notes.

    Release 232's notes say "Septmeber". The names are hand-written prose in a header,
    so we accept an exact match or an unambiguous three-letter prefix, and nothing
    looser: a fuzzy match would happily read "Smarch" as March and mis-date a release
    by a season, which is the kind of quiet error this tool exists to avoid.
    """
    text = name.strip().lower()
    if text in _MONTHS:
        return _MONTHS[text]
    matches = {v for k, v in _MONTHS.items() if len(text) >= 3 and k.startswith(text[:3])}
    if len(matches) == 1:
        return next(iter(matches))
    raise ReleaseError(f"unrecognised month {name!r}")


def parse_notes_date(text: str, expected_release: int | None = None) -> dt.date:
    match = _HEADER_RE.search(text[:6000])
    if not match:
        raise ReleaseError("could not find a 'Release N / <date>' header in release notes")
    number = int(match.group("release"))
    if expected_release is not None and number != expected_release:
        raise ReleaseError(
            f"release notes header says release {number}, expected {expected_release}"
        )
    month = parse_month(match.group("month"))
    return dt.date(int(match.group("year")), month, int(match.group("day")))


def list_archived_releases(client: NcbiClient) -> list[int]:
    """Release numbers that have both an archived catalog and archived release notes."""
    catalogs = {
        int(n)
        for n in _CATALOG_HREF_RE.findall(client.get_text(ftp_url(CATALOG_ARCHIVE_DIR + "/")))
    }
    notes = {
        int(n) for n in _NOTES_HREF_RE.findall(client.get_text(ftp_url(NOTES_ARCHIVE_DIR + "/")))
    }
    return sorted(catalogs & notes)


def list_current_release(client: NcbiClient) -> int | None:
    """The live (not yet archived) release number, if the top-level directory has one."""
    numbers = {
        int(n) for n in _CATALOG_HREF_RE.findall(client.get_text(ftp_url(CATALOG_DIR + "/")))
    }
    return max(numbers) if numbers else None


def refresh_releases(client: NcbiClient, progress: bool = False) -> list[Release]:
    """Rebuild the release table by reading every release's notes header."""
    archived = list_archived_releases(client)
    current = list_current_release(client)
    wanted = [(n, True) for n in archived]
    if current is not None and current not in archived:
        wanted.append((current, False))

    iterator = wanted
    if progress:
        from tqdm import tqdm

        iterator = tqdm(wanted, desc="release notes", unit="release")

    releases: list[Release] = []
    for number, is_archived in iterator:
        placeholder = Release(number=number, date=dt.date(1970, 1, 1))
        text = client.get_text(placeholder.notes_url(archived=is_archived), max_bytes=20_000)
        releases.append(Release(number=number, date=parse_notes_date(text, number)))

    check_monotonic(releases)
    return sorted(releases)
