"""ENA metadata parsing: the pairing rules that a project-wide loop gets wrong."""

from __future__ import annotations

import pytest

from hindsight.samples import ena

HEADER = "\t".join(
    [
        "run_accession",
        "study_accession",
        "sample_accession",
        "library_strategy",
        "library_layout",
        "instrument_model",
        "read_count",
        "base_count",
        "collection_date",
        "fastq_ftp",
        "fastq_bytes",
        "fastq_md5",
    ]
)


def row(accession, files, *, layout="PAIRED", read_count=10, collection_date="missing"):
    urls = ";".join(f"ftp.sra.ebi.ac.uk/vol1/fastq/x/{accession}/{name}" for name, _, _ in files)
    sizes = ";".join(str(size) for _, size, _ in files)
    md5s = ";".join(md5 for _, _, md5 in files)
    return "\t".join(
        [
            accession,
            "PRJNA1",
            "SAMN1",
            "RNA-Seq",
            layout,
            "Illumina HiSeq 4000",
            str(read_count),
            "1000",
            collection_date,
            urls,
            sizes,
            md5s,
        ]
    )


def test_paired_run_exposes_its_mates():
    text = (
        HEADER
        + "\n"
        + row("SRR1", [("SRR1_1.fastq.gz", 10, "a"), ("SRR1_2.fastq.gz", 11, "b")])
        + "\n"
    )
    run = ena.parse_filereport(text)[0]
    assert run.layout == "paired"
    assert run.mates is not None
    assert [f.name for f in run.mates] == ["SRR1_1.fastq.gz", "SRR1_2.fastq.gz"]
    assert run.single is None
    assert run.total_bytes == 21
    assert run.files[0].url.startswith("https://")


def test_merged_single_file_run_has_no_mates():
    text = HEADER + "\n" + row("SRR2", [("SRR2.fastq.gz", 99, "c")], layout="PAIRED") + "\n"
    run = ena.parse_filereport(text)[0]
    assert run.layout == "single"
    assert run.mates is None
    assert run.single is not None


def test_paired_plus_orphan_is_recognised_not_mistaken_for_single():
    # 189 of the 440 runs in PRJNA699976 look like this.
    text = (
        HEADER
        + "\n"
        + row(
            "SRR3",
            [
                ("SRR3.fastq.gz", 1, "a"),
                ("SRR3_1.fastq.gz", 2, "b"),
                ("SRR3_2.fastq.gz", 3, "c"),
            ],
        )
        + "\n"
    )
    run = ena.parse_filereport(text)[0]
    assert run.layout == "paired+orphan"
    assert run.mates is not None
    assert run.single is not None


def test_library_layout_is_not_trusted_over_the_files():
    text = HEADER + "\n" + row("SRR4", [("SRR4.fastq.gz", 5, "d")], layout="PAIRED") + "\n"
    run = ena.parse_filereport(text)[0]
    assert run.library_layout == "PAIRED"
    assert run.layout == "single"


def test_mismatched_columns_are_refused_rather_than_guessed():
    broken = "\t".join(
        [
            "SRR5",
            "PRJNA1",
            "SAMN1",
            "RNA-Seq",
            "PAIRED",
            "HiSeq",
            "10",
            "1000",
            "missing",
            "ftp.x/SRR5_1.fastq.gz;ftp.x/SRR5_2.fastq.gz",
            "10",
            "a;b",
        ]
    )
    with pytest.raises(ValueError, match="refusing to guess"):
        ena.parse_filereport(HEADER + "\n" + broken + "\n")


def test_empty_report_is_empty_not_an_error():
    assert ena.parse_filereport("") == []
    assert ena.parse_filereport(HEADER + "\n") == []


def test_missing_collection_date_is_carried_through_verbatim():
    # PRJNA699976 reports 'missing' for all 440 runs; the benchmark has to say so
    # rather than infer a date.
    text = HEADER + "\n" + row("SRR6", [("SRR6_1.fastq.gz", 1, "a"), ("SRR6_2.fastq.gz", 2, "b")])
    assert ena.parse_filereport(text)[0].collection_date == "missing"


def test_md5_of_a_known_file(tmp_path):
    path = tmp_path / "x"
    path.write_bytes(b"hindsight")
    assert ena.md5_file(path) == "2b85162ea1dc28002cc8fbebf051acd7"
