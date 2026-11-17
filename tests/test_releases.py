import datetime as dt

import pytest

from hindsight import releases


def test_vendored_table_loads_and_is_consistent():
    table = releases.load_releases()
    assert len(table) > 100
    releases.check_monotonic(table)


def test_numbering_gap_is_real():
    """RefSeq jumps from 99 to 200. Code must not assume contiguous numbering."""
    numbers = {r.number for r in releases.load_releases()}
    assert 99 in numbers
    assert 200 in numbers
    assert not numbers & set(range(100, 200))


@pytest.mark.parametrize(
    ("number", "date"),
    [
        (93, "2019-03-13"),
        (97, "2019-11-04"),
        (98, "2020-01-06"),
        (99, "2020-03-02"),
    ],
)
def test_known_release_dates(number, date):
    """Dates quoted in the benchmark's design, pinned so a bad refresh is visible."""
    table = {r.number: r for r in releases.load_releases()}
    assert table[number].date == dt.date.fromisoformat(date)


def test_resolve_picks_the_benchmark_release():
    """A December 2019 cutoff must land on release 97, November 2019."""
    chosen = releases.resolve(dt.date(2019, 12, 1))
    assert chosen.number == 97
    assert chosen.date == dt.date(2019, 11, 4)


def test_resolve_is_strictly_before_the_cutoff():
    """A cutoff falling on a release date excludes that release.

    The cutoff is the first day that must not be visible, so a release issued that
    morning is already inside the excluded window.
    """
    table = releases.load_releases()
    on_release_day = releases.resolve(dt.date(2019, 11, 4), table)
    assert on_release_day.number == 96
    day_after = releases.resolve(dt.date(2019, 11, 5), table)
    assert day_after.number == 97


def test_resolve_before_the_archive_begins():
    with pytest.raises(releases.ReleaseError, match="no RefSeq release precedes"):
        releases.resolve(dt.date(1999, 1, 1))


def test_next_release_after_reports_the_excluded_one():
    table = releases.load_releases()
    chosen = releases.resolve(dt.date(2019, 12, 1), table)
    following = releases.next_release_after(chosen, table)
    assert following is not None
    assert following.number == 98


def test_catalog_url_shape():
    release = releases.Release(number=97, date=dt.date(2019, 11, 4))
    assert release.catalog_url().endswith("release-catalog/archive/RefSeq-release97.catalog.gz")
    assert release.catalog_url(archived=False).endswith(
        "release-catalog/RefSeq-release97.catalog.gz"
    )


NOTES_HEADER = """\
********************************************************************************
RefSeq-release97.txt       ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-notes/

\t\tNCBI Reference Sequence (RefSeq) Database

\t\t\tRelease 97
\t\t\tNovember 4, 2019

\t\tDistribution Release Notes
"""


def test_parse_notes_date():
    assert releases.parse_notes_date(NOTES_HEADER, 97) == dt.date(2019, 11, 4)


def test_parse_notes_date_rejects_wrong_release():
    with pytest.raises(releases.ReleaseError, match="expected 98"):
        releases.parse_notes_date(NOTES_HEADER, 98)


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        # Release 5's notes: a stray space before the comma.
        ("Release 5\n\t\t\tMay 2 , 2004\n", dt.date(2004, 5, 2)),
        # Release 39's notes: no space after the comma.
        ("Release 39\n\t\t\tJanuary 23,2010\n", dt.date(2010, 1, 23)),
        # Release 232's notes misspell the month as "Septmeber".
        ("Release 232\n\t\t\tSeptmeber 2, 2025\n", dt.date(2025, 9, 2)),
    ],
)
def test_parse_notes_date_tolerates_ncbi_typos(header, expected):
    assert releases.parse_notes_date(header) == expected


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("January", 1),
        ("jan", 1),
        (" march ", 3),
        ("DECEMBER", 12),
        ("Septmeber", 9),
        ("Sept", 9),
    ],
)
def test_parse_month_accepts_exact_names_and_unambiguous_prefixes(name, expected):
    assert releases.parse_month(name) == expected


@pytest.mark.parametrize("name", ["Smarch", "Jnauary", "", "Ju", "x"])
def test_parse_month_rejects_anything_looser(name):
    """A fuzzy month match could mis-date a release by a season. Refuse instead."""
    with pytest.raises(releases.ReleaseError):
        releases.parse_month(name)


def test_check_monotonic_catches_a_backdated_release():
    broken = [
        releases.Release(number=1, date=dt.date(2020, 1, 1)),
        releases.Release(number=2, date=dt.date(2019, 1, 1)),
    ]
    with pytest.raises(releases.ReleaseError, match="is not after"):
        releases.check_monotonic(broken)


def test_table_round_trips():
    table = releases.load_releases()
    assert releases.parse_table(releases.format_table(table)) == table
