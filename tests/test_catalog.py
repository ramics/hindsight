import gzip
from pathlib import Path

import pytest

from hindsight import catalog
from hindsight.accessions import Molecule

FIXTURE = Path(__file__).parent / "data" / "catalog_fixture.tsv.gz"

SARS_COV_ROW = (
    "694009\tSevere acute respiratory syndrome-related coronavirus\t"
    "NC_004718.3\tcomplete|viral\tna\t29751\n"
)


def fixture_lines():
    with gzip.open(FIXTURE, "rt") as handle:
        yield from handle


def test_parse_row():
    row = catalog.parse_row(SARS_COV_ROW)
    assert row.taxid == 694009
    assert row.accession == "NC_004718.3"
    assert row.base_accession == "NC_004718"
    assert row.version == 3
    assert row.directories == ("complete", "viral")
    assert row.length == 29751
    assert row.molecule is Molecule.NUCLEOTIDE


def test_parse_row_rejects_wrong_column_count():
    with pytest.raises(catalog.CatalogError, match="expected 6 columns"):
        catalog.parse_row("a\tb\tc\n")


def test_is_viral_matches_the_directory_not_the_organism_name():
    viral = catalog.parse_row(SARS_COV_ROW)
    human = catalog.parse_row(
        "9606\tHomo sapiens\tNC_000001.11\tcomplete|vertebrate_mammalian\tna\t248956422\n"
    )
    assert catalog.is_viral(viral)
    assert not catalog.is_viral(human)


def test_plasmid_viral_rows_are_viral():
    """A few viral records also sit in the plasmid directory; they are still viral."""
    row = catalog.parse_row(
        "10847\tEscherichia virus phiX174\tNC_001422.1\tcomplete|plasmid|viral\tna\t5386\n"
    )
    assert catalog.is_viral(row)


def test_iter_viral_rows_filters_and_counts():
    stats = catalog.ManifestStats()
    rows = list(catalog.iter_viral_rows(fixture_lines(), stats))
    accessions = {row.accession for row in rows}

    assert "NC_004718.3" in accessions, "SARS-CoV must survive the filter"
    assert "NC_000001.11" not in accessions, "human records must not"
    assert "NC_000913.3" not in accessions, "bacterial records must not"
    assert "NC_001422.1" in accessions, "plasmid|viral records must"

    assert stats.rows_malformed == 1, "the short row is counted, not fatal"
    assert stats.rows_unclassified == 1, "the unparseable accession is counted, not fatal"
    assert stats.rows_viral == len(rows) + stats.rows_unclassified


def test_write_manifest_is_deterministic_and_sorted(tmp_path):
    first = tmp_path / "a.tsv"
    second = tmp_path / "b.tsv"
    for path in (first, second):
        stats = catalog.ManifestStats()
        catalog.write_manifest(catalog.iter_viral_rows(fixture_lines(), stats), path, stats)
    assert first.read_bytes() == second.read_bytes()

    rows = catalog.read_manifest(first)
    molecules = [row["molecule"] for row in rows]
    assert molecules == sorted(molecules), "nucleotide records come before protein ones"
    nucleotide = [r["accession"] for r in rows if r["molecule"] == "nucleotide"]
    assert nucleotide == sorted(nucleotide)


def test_manifest_stats_split_by_molecule(tmp_path):
    path = tmp_path / "manifest.tsv"
    stats = catalog.ManifestStats()
    catalog.write_manifest(catalog.iter_viral_rows(fixture_lines(), stats), path, stats)
    assert stats.nucleotide > 0
    assert stats.protein > 0
    assert stats.nucleotide + stats.protein == len(catalog.read_manifest(path))
    assert stats.nucleotide_bases > stats.nucleotide
    assert stats.taxids > 0


def test_manifest_accessions_by_molecule(tmp_path):
    path = tmp_path / "manifest.tsv"
    stats = catalog.ManifestStats()
    catalog.write_manifest(catalog.iter_viral_rows(fixture_lines(), stats), path, stats)
    nucleotide = catalog.manifest_accessions(path, Molecule.NUCLEOTIDE)
    assert "NC_004718.3" in nucleotide
    assert all(not acc.startswith("YP_") for acc in nucleotide)


def test_read_manifest_rejects_a_truncated_file(tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text("accession\ttaxid\n")
    with pytest.raises(catalog.CatalogError, match="missing columns"):
        catalog.read_manifest(path)


def test_chunk_reader_reassembles_a_gzip_stream():
    raw = gzip.compress(b"one\ntwo\nthree\n")
    chunks = iter([raw[:5], raw[5:9], raw[9:]])
    reader = catalog._ChunkReader(chunks)
    with gzip.open(reader, "rt") as handle:
        assert handle.read().splitlines() == ["one", "two", "three"]
