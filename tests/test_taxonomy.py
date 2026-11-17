import datetime as dt
import zipfile

import pytest

from hindsight import taxonomy
from hindsight.ncbi import NcbiClient

SNAPSHOT_DATES = ["2014-08-01", "2018-12-01", "2019-10-01", "2019-11-01", "2019-12-01"]


def snapshots(flavour="taxdmp"):
    return [
        taxonomy.TaxdumpSnapshot(date=dt.date.fromisoformat(d), flavour=flavour)
        for d in SNAPSHOT_DATES
    ]


def test_resolve_snapshot_takes_the_last_one_at_or_before_the_release():
    chosen = taxonomy.resolve_snapshot(dt.date(2019, 11, 4), snapshots())
    assert chosen.date == dt.date(2019, 11, 1)


def test_resolve_snapshot_allows_the_release_day_itself():
    """A snapshot taken the morning of the release is contemporaneous, not ahead."""
    chosen = taxonomy.resolve_snapshot(dt.date(2019, 11, 1), snapshots())
    assert chosen.date == dt.date(2019, 11, 1)


def test_resolve_snapshot_never_reaches_forwards():
    chosen = taxonomy.resolve_snapshot(dt.date(2019, 11, 30), snapshots())
    assert chosen.date == dt.date(2019, 11, 1)


def test_resolve_snapshot_refuses_when_the_archive_starts_too_late():
    with pytest.raises(taxonomy.TaxonomyError, match="earliest"):
        taxonomy.resolve_snapshot(dt.date(2010, 1, 1), snapshots())


def test_snapshot_url_shape():
    snapshot = taxonomy.TaxdumpSnapshot(date=dt.date(2019, 11, 1), flavour="taxdmp")
    assert snapshot.filename == "taxdmp_2019-11-01.zip"
    assert snapshot.url.endswith("taxdump_archive/taxdmp_2019-11-01.zip")


def test_list_snapshots_rejects_an_unknown_flavour():
    with pytest.raises(taxonomy.TaxonomyError, match="unknown taxdump flavour"):
        taxonomy.list_snapshots(object(), "not_a_flavour")  # type: ignore[arg-type]


# -- dump parsing --------------------------------------------------------------

NODES = (
    "1\t|\t1\t|\tno rank\t|\n694009\t|\t2509511\t|\tspecies\t|\n11118\t|\t2501931\t|\tgenus\t|\n"
)
MERGED = "123\t|\t694009\t|\n"
DELNODES = "10658\t|\n999999\t|\n"


def write_dumps(tmp_path):
    (tmp_path / "nodes.dmp").write_text(NODES)
    (tmp_path / "merged.dmp").write_text(MERGED)
    (tmp_path / "delnodes.dmp").write_text(DELNODES)
    return tmp_path


def test_read_taxids(tmp_path):
    write_dumps(tmp_path)
    assert taxonomy.read_taxids(tmp_path / "nodes.dmp") == {1, 694009, 11118}


def test_read_merged(tmp_path):
    write_dumps(tmp_path)
    assert taxonomy.read_merged(tmp_path / "merged.dmp") == {123: 694009}


def test_read_deleted(tmp_path):
    write_dumps(tmp_path)
    assert taxonomy.read_deleted(tmp_path / "delnodes.dmp") == {10658, 999999}


def test_check_taxids_classifies_every_kind_of_miss(tmp_path):
    """A miss is either merged, deleted, or genuinely unexplained."""
    write_dumps(tmp_path)
    rows = [
        {"accession": "NC_004718.3", "taxid": "694009"},  # present
        {"accession": "NC_001421.2", "taxid": "10658"},  # deleted in this snapshot
        {"accession": "NC_000000.1", "taxid": "123"},  # merged away
        {"accession": "NC_111111.1", "taxid": "4242"},  # unexplained
    ]
    check = taxonomy.check_taxids(
        rows,
        tmp_path / "nodes.dmp",
        tmp_path / "merged.dmp",
        tmp_path / "delnodes.dmp",
    )
    assert check.manifest_taxids == 4
    assert check.present == 1
    assert check.missing == [123, 4242, 10658]
    assert check.resolvable_via_merged == {123: 694009}
    assert check.deleted_in_snapshot == [10658]
    assert check.as_dict()["unexplained"] == [4242]


def test_write_seqid_map(tmp_path):
    rows = [
        {"accession": "NC_004718.3", "taxid": "694009"},
        {"accession": "YP_009639955.1", "taxid": "10658"},
    ]
    path = tmp_path / "seqid2taxid.map"
    assert taxonomy.write_seqid_map(rows, path) == 2
    assert path.read_text() == "NC_004718.3\t694009\nYP_009639955.1\t10658\n"


def test_extract_snapshot(tmp_path):
    archive = tmp_path / "taxdmp.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        for name in taxonomy.REQUIRED_DUMPS:
            zf.writestr(name, "1\t|\t1\t|\n")
        zf.writestr("readme.txt", "ignored")
    written = taxonomy.extract_snapshot(archive, tmp_path / "out")
    assert sorted(p.name for p in written) == sorted(taxonomy.REQUIRED_DUMPS)
    assert not (tmp_path / "out" / "readme.txt").exists()


def test_extract_snapshot_reports_a_missing_member(tmp_path):
    archive = tmp_path / "taxdmp.zip"
    with zipfile.ZipFile(archive, "w") as zf:
        zf.writestr("nodes.dmp", "1\t|\t1\t|\n")
    with pytest.raises(taxonomy.TaxonomyError, match="missing expected members"):
        taxonomy.extract_snapshot(archive, tmp_path / "out")


@pytest.mark.network
def test_the_archive_still_covers_the_benchmark_date():
    available = taxonomy.list_snapshots(NcbiClient(), "taxdmp")
    chosen = taxonomy.resolve_snapshot(dt.date(2019, 11, 4), available)
    assert chosen.date == dt.date(2019, 11, 1)
