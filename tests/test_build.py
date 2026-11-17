"""End-to-end build against a local catalog fixture, with no network."""

import datetime as dt
import json
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from hindsight import build as build_mod
from hindsight import catalog as catalog_mod
from hindsight import cli as cli_mod
from hindsight import provenance as provenance_mod
from hindsight import releases as releases_mod
from hindsight.accessions import Molecule

FIXTURE = Path(__file__).parent / "data" / "catalog_fixture.tsv.gz"

RELEASES = [
    releases_mod.Release(number=96, date=dt.date(2019, 9, 9)),
    releases_mod.Release(number=97, date=dt.date(2019, 11, 4)),
    releases_mod.Release(number=98, date=dt.date(2020, 1, 6)),
]


class OfflineClient:
    """Serves the catalog fixture and synthesises sequence for any accession."""

    def __init__(self):
        self.efetch_calls = 0

    def head(self, url):
        return {}

    def download(self, url, dest, **kwargs):
        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(FIXTURE, dest)
        return None

    def efetch(self, db, ids, rettype="fasta", retmode="text"):
        self.efetch_calls += 1
        letter = "A" if db == "nuccore" else "M"
        return "".join(f">{acc} synthetic {db}\n{letter * 30}\n" for acc in ids)


@pytest.fixture
def offline(monkeypatch):
    monkeypatch.setattr(
        catalog_mod, "resolve_catalog_url", lambda client, release: "file://fixture"
    )
    return OfflineClient()


def run_build(offline, tmp_path, **kwargs):
    return build_mod.build(
        offline,
        dt.date(2019, 12, 1),
        tmp_path / "db",
        tmp_path / "cache",
        releases=RELEASES,
        with_taxonomy=False,
        progress=False,
        quiet=True,
        **kwargs,
    )


def test_manifest_only_fetches_nothing(offline, tmp_path):
    result = run_build(offline, tmp_path, manifest_only=True)
    assert result.release.number == 97
    assert offline.efetch_calls == 0
    assert result.provenance_path is None
    assert result.manifest_path.exists()


def test_full_build_writes_every_output(offline, tmp_path):
    result = run_build(offline, tmp_path)
    out = result.out_dir
    for name in (
        build_mod.MANIFEST_NAME,
        build_mod.NUCLEOTIDE_FASTA,
        build_mod.PROTEIN_FASTA,
        build_mod.WITHDRAWN_NAME,
        provenance_mod.PROVENANCE_NAME,
    ):
        assert (out / name).exists(), name
    assert provenance_mod.verify(out) == []


def test_provenance_records_the_time_control(offline, tmp_path):
    result = run_build(offline, tmp_path)
    record = json.loads((result.out_dir / provenance_mod.PROVENANCE_NAME).read_text())
    control = record["time_control"]
    assert control["cutoff_date"] == "2019-12-01"
    assert control["release"] == 97
    assert control["release_date"] == "2019-11-04"
    assert control["excluded_next_release"] == {"release": 98, "release_date": "2020-01-06"}
    assert record["manifest"]["nucleotide"] > 0
    assert record["manifest"]["protein"] > 0


def test_build_is_byte_identical_on_a_rerun(offline, tmp_path):
    first = run_build(offline, tmp_path / "one")
    second = run_build(offline, tmp_path / "two")
    for name in (build_mod.MANIFEST_NAME, build_mod.NUCLEOTIDE_FASTA, build_mod.PROTEIN_FASTA):
        assert (first.out_dir / name).read_bytes() == (second.out_dir / name).read_bytes()


def test_rerunning_into_the_same_cache_refetches_nothing(offline, tmp_path):
    run_build(offline, tmp_path)
    before = offline.efetch_calls
    assert before > 0
    run_build(offline, tmp_path)
    assert offline.efetch_calls == before


def test_molecule_selection(offline, tmp_path):
    result = run_build(offline, tmp_path, molecules=(Molecule.NUCLEOTIDE,))
    assert (result.out_dir / build_mod.NUCLEOTIDE_FASTA).exists()
    assert not (result.out_dir / build_mod.PROTEIN_FASTA).exists()


def test_limit_caps_the_record_count(offline, tmp_path):
    result = run_build(offline, tmp_path, limit=2)
    text = (result.out_dir / build_mod.NUCLEOTIDE_FASTA).read_text()
    assert text.count(">") == 2


def test_verify_detects_tampering(offline, tmp_path):
    result = run_build(offline, tmp_path)
    target = result.out_dir / build_mod.NUCLEOTIDE_FASTA
    target.write_text(target.read_text() + ">NC_045512.2 leaked\nACGT\n")
    problems = provenance_mod.verify(result.out_dir)
    assert len(problems) == 1
    assert build_mod.NUCLEOTIDE_FASTA in problems[0]


def test_cutoff_before_the_archive_is_an_error(offline, tmp_path):
    with pytest.raises(releases_mod.ReleaseError):
        build_mod.build(
            offline,
            dt.date(2000, 1, 1),
            tmp_path / "db",
            tmp_path / "cache",
            releases=RELEASES,
            with_taxonomy=False,
            quiet=True,
        )


# -- CLI -----------------------------------------------------------------------


def test_cli_resolve():
    result = CliRunner().invoke(cli_mod.main, ["resolve", "2019-12-01"])
    assert result.exit_code == 0, result.output
    assert "selected release:  97 (2019-11-04)" in result.output
    assert "first excluded:    release 98 (2020-01-06)" in result.output


def test_cli_releases_json():
    result = CliRunner().invoke(cli_mod.main, ["releases", "--json", "--limit", "3"])
    assert result.exit_code == 0, result.output
    rows = json.loads(result.output)
    assert len(rows) == 3
    assert set(rows[0]) == {"release", "date"}


def test_cli_rejects_a_bad_date():
    result = CliRunner().invoke(cli_mod.main, ["resolve", "December 2019"])
    assert result.exit_code != 0


def test_cli_verify(offline, tmp_path):
    result = run_build(offline, tmp_path)
    invoked = CliRunner().invoke(cli_mod.main, ["verify", str(result.out_dir)])
    assert invoked.exit_code == 0, invoked.output
    assert "release 97 (2019-11-04)" in invoked.output


def test_cli_verify_fails_on_tampering(offline, tmp_path):
    result = run_build(offline, tmp_path)
    (result.out_dir / build_mod.MANIFEST_NAME).write_text("tampered\n")
    invoked = CliRunner().invoke(cli_mod.main, ["verify", str(result.out_dir)])
    assert invoked.exit_code == 1
