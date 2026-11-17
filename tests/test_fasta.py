import pytest

from hindsight import fasta

SAMPLE = """\
>NC_004718.3 SARS coronavirus Tor2, complete genome
ACGTACGTAC
GTACGTACGT

>YP_009639955.1 some protein [Virus]
MKVLA
"""


def test_parse():
    records = list(fasta.parse(SAMPLE))
    assert [r.accession for r in records] == ["NC_004718.3", "YP_009639955.1"]
    assert records[0].sequence == "ACGTACGTACGTACGTACGT"
    assert records[0].length == 20
    assert records[1].sequence == "MKVLA"


def test_parse_empty_input():
    assert list(fasta.parse("")) == []


def test_parse_rejects_sequence_before_a_defline():
    with pytest.raises(fasta.FastaError, match="before the first defline"):
        list(fasta.parse("ACGT\n>x\nACGT\n"))


@pytest.mark.parametrize(
    ("header", "accession"),
    [
        (">NC_004718.3 SARS coronavirus", "NC_004718.3"),
        ("NC_004718.3 SARS coronavirus", "NC_004718.3"),
        (">NC_004718.3", "NC_004718.3"),
    ],
)
def test_accession_from_header(header, accession):
    assert fasta.accession_from_header(header) == accession


def test_accession_from_header_rejects_a_blank_defline():
    with pytest.raises(fasta.FastaError):
        fasta.accession_from_header(">")


def test_write_wraps_at_a_fixed_width(tmp_path):
    record = fasta.Record("test record", "A" * 165)
    path = tmp_path / "out.fa"
    assert fasta.write_file(path, [record], width=70) == 1
    lines = path.read_text().splitlines()
    assert lines[0] == ">test record"
    assert [len(line) for line in lines[1:]] == [70, 70, 25]


def test_write_round_trips(tmp_path):
    records = list(fasta.parse(SAMPLE))
    path = tmp_path / "out.fa"
    fasta.write_file(path, records)
    assert list(fasta.parse(path.read_text())) == records


def test_rewrapping_is_idempotent(tmp_path):
    """Assembling twice must give the same bytes, whatever width NCBI used."""
    first = tmp_path / "a.fa"
    second = tmp_path / "b.fa"
    fasta.write_file(first, fasta.parse(SAMPLE))
    fasta.write_file(second, fasta.parse(first.read_text()))
    assert first.read_bytes() == second.read_bytes()


def test_record_sha256_is_sequence_only():
    """A changed defline must not change the sequence digest."""
    a = fasta.Record("NC_004718.3 one description", "ACGT")
    b = fasta.Record("NC_004718.3 another description", "ACGT")
    assert a.sha256 == b.sha256


def test_sha256_file(tmp_path):
    path = tmp_path / "x.fa"
    path.write_text("abc")
    assert fasta.sha256_file(path) == (
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )
