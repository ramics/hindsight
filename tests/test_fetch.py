import json

import pytest

from hindsight import fasta, fetch
from hindsight.accessions import Molecule
from hindsight.ncbi import NcbiClient

RECORDS = {
    "NC_004718.3": "ACGT" * 10,
    "NC_014470.1": "TTTT" * 10,
    "NC_025217.1": "GGGG" * 10,
}


class FakeClient:
    """Stands in for E-utilities; records what was asked for."""

    def __init__(self, records=None, drop=(), substitute=None):
        self.records = dict(records if records is not None else RECORDS)
        self.drop = set(drop)
        self.substitute = dict(substitute or {})
        self.calls: list[list[str]] = []

    def efetch(self, db, ids, rettype="fasta", retmode="text"):
        self.calls.append(list(ids))
        out = []
        for accession in ids:
            if accession in self.drop:
                continue
            served = self.substitute.get(accession, accession)
            sequence = self.records.get(served, "ACGT")
            out.append(f">{served} test record\n{sequence}\n")
        return "".join(out)


def test_plan_shards():
    assert fetch.plan_shards(list("abcde"), 2) == [["a", "b"], ["c", "d"], ["e"]]
    assert fetch.plan_shards([], 2) == []
    with pytest.raises(ValueError):
        fetch.plan_shards(["a"], 0)


def test_fetch_molecule_writes_shards_and_assembles(tmp_path):
    client = FakeClient()
    accessions = list(RECORDS)
    results, report = fetch.fetch_molecule(
        client, Molecule.NUCLEOTIDE, accessions, tmp_path, batch_size=2
    )

    assert report.requested == 3
    assert report.returned == 3
    assert report.withdrawn == []
    assert report.shards_total == 2

    out = tmp_path / "viral.fna"
    assert fetch.assemble(results, out) == 3
    parsed = list(fasta.parse(out.read_text()))
    assert [record.accession for record in parsed] == accessions


def test_a_resumed_run_refetches_nothing(tmp_path):
    accessions = list(RECORDS)
    first = FakeClient()
    fetch.fetch_molecule(first, Molecule.NUCLEOTIDE, accessions, tmp_path, batch_size=2)
    assert first.calls

    second = FakeClient()
    results, report = fetch.fetch_molecule(
        second, Molecule.NUCLEOTIDE, accessions, tmp_path, batch_size=2
    )
    assert second.calls == [], "nothing should be refetched"
    assert report.shards_from_cache == report.shards_total == 2
    assert report.returned == 3
    assert fetch.assemble(results, tmp_path / "out.fna") == 3


def test_a_tampered_shard_is_refetched(tmp_path):
    accessions = list(RECORDS)
    fetch.fetch_molecule(FakeClient(), Molecule.NUCLEOTIDE, accessions, tmp_path, batch_size=2)
    (tmp_path / "shard-000000.fa").write_text(">NC_004718.3 tampered\nAAAA\n")

    client = FakeClient()
    _, report = fetch.fetch_molecule(
        client, Molecule.NUCLEOTIDE, accessions, tmp_path, batch_size=2
    )
    assert client.calls, "the digest mismatch must trigger a refetch"
    assert report.shards_from_cache == 1


def test_changing_the_batch_size_invalidates_the_ledger(tmp_path):
    accessions = list(RECORDS)
    fetch.fetch_molecule(FakeClient(), Molecule.NUCLEOTIDE, accessions, tmp_path, batch_size=2)
    client = FakeClient()
    _, report = fetch.fetch_molecule(
        client, Molecule.NUCLEOTIDE, accessions, tmp_path, batch_size=3
    )
    assert report.shards_from_cache == 0
    assert client.calls == [accessions]


def test_a_withdrawn_record_is_recorded_not_fatal(tmp_path):
    client = FakeClient(drop={"NC_014470.1"})
    results, report = fetch.fetch_molecule(
        client, Molecule.NUCLEOTIDE, list(RECORDS), tmp_path, batch_size=3
    )
    assert report.withdrawn == ["NC_014470.1"]
    assert report.returned == 2
    # The shard was retried one id at a time, so the other two records survived.
    assert len(client.calls) == 2
    assert fetch.assemble(results, tmp_path / "out.fna") == 2


def test_a_record_returned_at_the_wrong_version_is_rejected(tmp_path):
    """The one failure mode that would silently break the time control."""
    client = FakeClient(substitute={"NC_004718.3": "NC_004718.4"})
    results, report = fetch.fetch_molecule(
        client, Molecule.NUCLEOTIDE, ["NC_004718.3"], tmp_path, batch_size=1
    )
    assert report.returned == 0
    assert report.withdrawn == ["NC_004718.3"]
    assert {"requested": "NC_004718.3", "returned": "NC_004718.4"} in report.version_mismatches
    assert fetch.assemble(results, tmp_path / "out.fna") == 0


def test_records_within_a_shard_keep_the_requested_order(tmp_path):
    """Output order must not depend on the order NCBI happens to answer in."""

    class Shuffling(FakeClient):
        def efetch(self, db, ids, rettype="fasta", retmode="text"):
            return super().efetch(db, list(reversed(ids)), rettype, retmode)

    accessions = list(RECORDS)
    results, _ = fetch.fetch_molecule(
        Shuffling(), Molecule.NUCLEOTIDE, accessions, tmp_path, batch_size=3
    )
    out = tmp_path / "out.fna"
    fetch.assemble(results, out)
    got = [record.accession for record in fasta.parse(out.read_text())]
    assert got == accessions


def test_ledger_is_valid_json_after_a_run(tmp_path):
    fetch.fetch_molecule(FakeClient(), Molecule.NUCLEOTIDE, list(RECORDS), tmp_path, batch_size=2)
    data = json.loads((tmp_path / fetch.LEDGER_NAME).read_text())
    assert data["batch_size"] == 2
    assert set(data["shards"]) == {"0", "1"}


def test_a_corrupt_ledger_does_not_break_the_run(tmp_path):
    (tmp_path / fetch.LEDGER_NAME).write_text("{not json")
    _, report = fetch.fetch_molecule(
        FakeClient(), Molecule.NUCLEOTIDE, list(RECORDS), tmp_path, batch_size=2
    )
    assert report.returned == 3


def test_protein_and_nucleotide_use_different_entrez_databases(tmp_path):
    class Recording(FakeClient):
        def __init__(self):
            super().__init__()
            self.dbs: list[str] = []

        def efetch(self, db, ids, rettype="fasta", retmode="text"):
            self.dbs.append(db)
            return super().efetch(db, ids, rettype, retmode)

    client = Recording()
    fetch.fetch_molecule(client, Molecule.NUCLEOTIDE, ["NC_004718.3"], tmp_path / "n")
    fetch.fetch_molecule(client, Molecule.PROTEIN, ["YP_009639955.1"], tmp_path / "p")
    assert client.dbs == ["nuccore", "protein"]


@pytest.mark.network
def test_version_pinning_is_honoured():
    """NCBI must return the revision we ask for, not the current one.

    If this ever fails, hindsight cannot do what it claims and the benchmark's
    leakage control is void. NC_004718 went from 29,736 bases at version 1 to 29,751
    at version 3.
    """
    client = NcbiClient()
    lengths = {}
    for version in (1, 3):
        accession = f"NC_004718.{version}"
        records = list(fasta.parse(client.efetch("nuccore", [accession])))
        assert len(records) == 1
        assert records[0].accession == accession
        lengths[version] = records[0].length
    assert lengths[1] == 29736
    assert lengths[3] == 29751
