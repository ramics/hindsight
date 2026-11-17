"""The leakage argument, as assertions.

These tests are the point of the tool. They encode the claim the benchmark rests on:
a November 2019 database contains SARS-CoV and its bat relatives, and cannot contain
SARS-CoV-2 or RaTG13, because neither was public.

The accession-set assertions run against a fixture cut from the real release 97
catalog. The whole-catalog assertion is marked ``network`` and checks the full
385,308-row viral set from NCBI.
"""

import datetime as dt
import gzip
import re
from pathlib import Path

import pytest

from hindsight import catalog, releases
from hindsight.accessions import Molecule
from hindsight.ncbi import NcbiClient

FIXTURE = Path(__file__).parent / "data" / "catalog_fixture.tsv.gz"

# Published 2020-01-10 onwards; the RefSeq record arrived in release 99 (March 2020).
SARS_COV_2_ACCESSIONS = {"NC_045512", "MN908947"}
# Deposited 2020-01-27.
RATG13_ACCESSIONS = {"MN996532"}

PRE_EMERGENCE_RELATIVES = {
    "NC_004718.3",  # SARS coronavirus Tor2
    "NC_014470.1",  # Bat coronavirus BM48-31/BGR/2008
    "NC_025217.1",  # Bat Hp-betacoronavirus/Zhejiang2013
    "NC_009021.1",  # Rousettus bat coronavirus HKU9
}


def fixture_rows():
    with gzip.open(FIXTURE, "rt") as handle:
        return list(catalog.iter_viral_rows(handle))


def test_december_2019_cutoff_selects_november_2019():
    chosen = releases.resolve(dt.date(2019, 12, 1))
    assert (chosen.number, chosen.date) == (97, dt.date(2019, 11, 4))


def test_pre_emergence_relatives_are_present():
    accessions = {row.accession for row in fixture_rows()}
    assert accessions >= PRE_EMERGENCE_RELATIVES


def test_sars_cov_is_pinned_at_the_version_the_release_listed():
    """Version pinning is what makes the database a record of the past.

    Release 97 listed NC_004718 at version 3. Had it listed an earlier version, that
    is the one that would have to be fetched, 15 bases shorter than today's record.
    """
    rows = {row.accession: row for row in fixture_rows()}
    sars = rows["NC_004718.3"]
    assert sars.version == 3
    assert sars.taxid == 694009
    assert "Severe acute respiratory syndrome-related" in sars.organism


def test_no_sars_cov_2_and_no_ratg13_in_the_fixture():
    bases = {row.base_accession for row in fixture_rows()}
    assert not bases & SARS_COV_2_ACCESSIONS
    assert not bases & RATG13_ACCESSIONS


# Matching on organism name needs care in both directions. A bare "coronavirus 2"
# substring also matches "Human coronavirus 229E", which is a 1960s common-cold virus
# and entirely legitimate in a 2019 database; and the pre-2020 species name "Severe
# acute respiratory syndrome-related coronavirus" must not be mistaken for the
# post-2020 "Severe acute respiratory syndrome coronavirus 2". The word boundary
# distinguishes all of them.
POST_EMERGENCE_NAME_RE = re.compile(r"\bcoronavirus 2\b|\bsars-cov-2\b|\bratg13\b")


@pytest.mark.network
def test_full_release_97_viral_set_excludes_sars_cov_2():
    """Stream the real release 97 catalog and check the whole viral set.

    This downloads ~1.2 GB and takes around nine minutes. It is the assertion that the
    benchmark's central claim holds for the actual data, not only for a hand-picked
    fixture. The counts are exact because an archived catalog is immutable: if one of
    them moves, something about the archive has changed and we want to know.
    """
    client = NcbiClient()
    release = releases.resolve(dt.date(2019, 12, 1))
    stats = catalog.ManifestStats()
    bases: set[str] = set()
    nucleotide = 0
    coronavirus_genomes = 0
    post_emergence_names = set()
    for row in catalog.iter_viral_rows(catalog.stream_catalog_lines(client, release), stats):
        bases.add(row.base_accession)
        lowered = row.organism.lower()
        if POST_EMERGENCE_NAME_RE.search(lowered):
            post_emergence_names.add(row.organism)
        if row.molecule is Molecule.NUCLEOTIDE:
            nucleotide += 1
            if "coronavirus" in lowered:
                coronavirus_genomes += 1

    assert stats.rows_viral == 385_308
    assert nucleotide == 12_179
    assert coronavirus_genomes == 45

    relative_bases = {accession.split(".")[0] for accession in PRE_EMERGENCE_RELATIVES}
    assert relative_bases <= bases, "SARS-CoV and its bat relatives must all be present"

    assert not bases & SARS_COV_2_ACCESSIONS
    assert not bases & RATG13_ACCESSIONS
    assert post_emergence_names == set()
