"""FASTQ handling: identifier dialects, exact subsets, and refusing bad input."""

from __future__ import annotations

import gzip

import pytest

from hindsight.samples import reads


def write_fastq(path, names, *, mate_suffix="", casava=False):
    with gzip.open(path, "wt") as handle:
        for name in names:
            header = f"@{name}{mate_suffix}"
            if casava:
                header += " 1:N:0:ACGTAC"
            handle.write(f"{header}\nACGTACGTAC\n+\nIIIIIIIIII\n")
    return path


def test_normalise_id_handles_the_dialects():
    assert reads.normalise_id("@SRR1.1/1") == "SRR1.1"
    assert reads.normalise_id("@SRR1.1/2") == "SRR1.1"
    assert reads.normalise_id("@SRR1.1 1:N:0:TAG") == "SRR1.1"
    assert reads.normalise_id("SRR1.1") == "SRR1.1"
    assert reads.normalise_id("@SRR1.1/1 1:N:0:TAG") == "SRR1.1"


def test_count_and_ids_round_trip(tmp_path):
    path = write_fastq(tmp_path / "a.fastq.gz", ["r1", "r2", "r3"], mate_suffix="/1")
    assert reads.count_reads(path) == 3
    assert reads.read_ids(path) == {"r1", "r2", "r3"}


def test_truncated_fastq_is_an_error_not_a_short_count(tmp_path):
    path = tmp_path / "bad.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\n")
    with pytest.raises(ValueError, match="not a whole number|truncated"):
        reads.count_reads(path)
    with pytest.raises(ValueError, match="truncated"):
        list(reads.iter_records(path))


def test_non_fastq_input_is_refused(tmp_path):
    path = tmp_path / "notfastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("ACGT\nACGT\nACGT\nACGT\n")
    with pytest.raises(ValueError, match="expected a header"):
        list(reads.iter_records(path))


def test_write_subset_keeps_source_order_and_exact_count(tmp_path):
    source = write_fastq(tmp_path / "s.fastq.gz", [f"r{i}" for i in range(10)], casava=True)
    dest = tmp_path / "sub.fastq.gz"
    written = reads.write_subset(source, dest, {"r7", "r2", "r5"})
    assert written == 3
    assert [reads.normalise_id(h) for h, *_ in reads.iter_records(dest)] == ["r2", "r5", "r7"]


def test_write_subset_ignores_identifiers_that_are_not_there(tmp_path):
    source = write_fastq(tmp_path / "s.fastq.gz", ["r1", "r2"])
    assert reads.write_subset(source, tmp_path / "o.fastq.gz", {"r9"}) == 0


def test_write_excluding_is_the_complement(tmp_path):
    source = write_fastq(tmp_path / "s.fastq.gz", [f"r{i}" for i in range(5)])
    kept = reads.write_excluding(source, tmp_path / "bg.fastq.gz", {"r1", "r3"})
    assert kept == 3
    assert reads.read_ids(tmp_path / "bg.fastq.gz") == {"r0", "r2", "r4"}


def test_subset_and_background_partition_the_file(tmp_path):
    names = [f"r{i}" for i in range(50)]
    source = write_fastq(tmp_path / "s.fastq.gz", names)
    viral = {"r3", "r9", "r44"}
    a = reads.write_subset(source, tmp_path / "v.fastq.gz", viral)
    b = reads.write_excluding(source, tmp_path / "b.fastq.gz", viral)
    assert a + b == len(names)
    assert reads.read_ids(tmp_path / "v.fastq.gz") | reads.read_ids(tmp_path / "b.fastq.gz") == set(
        names
    )


def test_concatenated_gzip_members_read_as_one_stream(tmp_path):
    first = write_fastq(tmp_path / "a.fastq.gz", ["a1", "a2"])
    second = write_fastq(tmp_path / "b.fastq.gz", ["b1"])
    dest = tmp_path / "joined.fastq.gz"
    reads.concatenate([first, second], dest)
    assert [reads.normalise_id(h) for h, *_ in reads.iter_records(dest)] == ["a1", "a2", "b1"]
    assert reads.count_reads(dest) == 3


def test_ids_from_sam_skips_headers_and_dedupes(tmp_path):
    path = tmp_path / "ids.sam"
    path.write_text(
        "@HD\tVN:1.6\n"
        "frag1\t99\tNC_045512.2\t1\t42\t100M\n"
        "frag1\t147\tNC_045512.2\t200\t42\t100M\n"
        "frag2/1\t0\tNC_045512.2\t300\t42\t100M\n"
    )
    assert reads.ids_from_sam(path) == {"frag1", "frag2"}


def test_hash_fraction_keeps_both_mates_in_step(tmp_path):
    names = [f"r{i:05d}" for i in range(5000)]
    one = write_fastq(tmp_path / "a_1.fastq.gz", names, mate_suffix="/1")
    two = write_fastq(tmp_path / "a_2.fastq.gz", names, mate_suffix="/2")
    kept1 = reads.write_hash_fraction(one, tmp_path / "c_1.fastq.gz", seed=7, fraction=0.2)
    kept2 = reads.write_hash_fraction(two, tmp_path / "c_2.fastq.gz", seed=7, fraction=0.2)
    assert kept1 == kept2
    assert reads.read_ids(tmp_path / "c_1.fastq.gz") == reads.read_ids(tmp_path / "c_2.fastq.gz")


def test_hash_fraction_hits_about_the_target_depth(tmp_path):
    names = [f"r{i:05d}" for i in range(20000)]
    source = write_fastq(tmp_path / "a.fastq.gz", names)
    kept = reads.write_hash_fraction(source, tmp_path / "c.fastq.gz", seed=3, fraction=0.25)
    # Binomial around 5000; three sigma is about 185.
    assert abs(kept - 5000) < 200


def test_hash_fraction_is_reproducible_and_seed_dependent(tmp_path):
    names = [f"r{i:04d}" for i in range(2000)]
    source = write_fastq(tmp_path / "a.fastq.gz", names)
    reads.write_hash_fraction(source, tmp_path / "x.fastq.gz", seed=1, fraction=0.3)
    reads.write_hash_fraction(source, tmp_path / "y.fastq.gz", seed=1, fraction=0.3)
    reads.write_hash_fraction(source, tmp_path / "z.fastq.gz", seed=2, fraction=0.3)
    assert reads.read_ids(tmp_path / "x.fastq.gz") == reads.read_ids(tmp_path / "y.fastq.gz")
    assert reads.read_ids(tmp_path / "x.fastq.gz") != reads.read_ids(tmp_path / "z.fastq.gz")


def test_hash_fraction_nests_as_the_depth_rises(tmp_path):
    """A deeper cap must be a superset of a shallower one at the same seed.

    Not required by the benchmark, but it means raising the cap later extends the
    sample rather than replacing it, so results at the old depth stay comparable.
    """
    names = [f"r{i:04d}" for i in range(3000)]
    source = write_fastq(tmp_path / "a.fastq.gz", names)
    reads.write_hash_fraction(source, tmp_path / "shallow.fastq.gz", seed=9, fraction=0.1)
    reads.write_hash_fraction(source, tmp_path / "deep.fastq.gz", seed=9, fraction=0.4)
    shallow = reads.read_ids(tmp_path / "shallow.fastq.gz")
    assert shallow <= reads.read_ids(tmp_path / "deep.fastq.gz")


def test_hash_fraction_extremes_are_pass_through_and_empty(tmp_path):
    names = [f"r{i}" for i in range(100)]
    source = write_fastq(tmp_path / "a.fastq.gz", names)
    everything = reads.write_hash_fraction(
        source, tmp_path / "all.fastq.gz", seed=1, fraction=1.0
    )
    assert everything == 100
    assert reads.write_hash_fraction(source, tmp_path / "none.fastq.gz", seed=1, fraction=0.0) == 0
    assert reads.keeps_fragment("r1", seed=1, fraction=1.5) is True
    assert reads.keeps_fragment("r1", seed=1, fraction=-0.2) is False


def test_count_redacted_counts_all_n_reads(tmp_path):
    path = tmp_path / "r.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@a\nACGTACGT\n+\nIIIIIIII\n")
        handle.write("@b\n" + "N" * 151 + "\n+\n" + "!" * 151 + "\n")
        handle.write("@c\nACGTNNNN\n+\nIIIIIIII\n")
        handle.write("@d\n" + "N" * 8 + "\n+\n" + "!" * 8 + "\n")
    # Only wholly redacted reads count; a read with some Ns is still sequence.
    assert reads.count_redacted(path) == 2


def test_count_redacted_on_a_clean_file_is_zero(tmp_path):
    path = write_fastq(tmp_path / "clean.fastq.gz", ["a", "b", "c"])
    assert reads.count_redacted(path) == 0
