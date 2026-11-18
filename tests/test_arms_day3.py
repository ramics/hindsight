"""The Day 3 driver's bookkeeping: the config it reads, the record it merges, the
summary Day 5 will take its threshold from, and the pruning that makes native depth
fit on a small machine.
"""

from __future__ import annotations

import json

import pytest

from hindsight.arms import day3

CONFIG = """
[arm]
samples_config = "config/day2.toml"
work_dir = "scratch"
out_dir = "scratch/arm-a"
provenance = "results/day3-provenance.json"

[database]
cutoff = "2019-12-01"
source_db = "db/nov2019"
kraken_db = "scratch/kraken2/nov2019"
mask_low_complexity = true

[classify]
threads = 2
confidence = 0.0
minimum_hit_groups = 2
verify_composition = false
verify_level = "viral000010"

[negatives]
max_pairs = 0
prune_after_classify = true

[tools]
kraken2 = "2.1.3"
"""


def test_load_config(tmp_path):
    path = tmp_path / "day3.toml"
    path.write_text(CONFIG)
    arm = day3.load_config(path)
    assert arm.negatives_max_pairs == 0
    assert arm.threads == 2
    assert arm.confidence == 0.0
    assert arm.minimum_hit_groups == 2
    assert arm.mask_low_complexity is True
    assert arm.verify_composition is False
    assert str(arm.kraken_db) == "scratch/kraken2/nov2019"


def test_load_config_defaults_to_native_depth(tmp_path):
    """An empty config is the paper's configuration, not a capped one."""
    path = tmp_path / "day3.toml"
    path.write_text("")
    arm = day3.load_config(path)
    assert arm.negatives_max_pairs == 0
    assert arm.confidence == 0.0
    assert arm.mask_low_complexity is True


def negative(accession, hits, pairs):
    return {"accession": accession, "false_alarm_reads": hits, "pairs": pairs}


def test_false_alarm_summary():
    summary = day3.summarise_false_alarms(
        [
            negative("SRR1", 0, 10_000_000),
            negative("SRR2", 3, 10_000_000),
            negative("SRR3", 1, 20_000_000),
        ]
    )
    assert summary["samples"] == 3
    assert summary["samples_with_at_least_one_hit"] == 2
    assert summary["reads_total"] == 4
    assert summary["reads_max_in_one_sample"] == 3
    assert summary["non_host_pairs_total"] == 40_000_000
    assert summary["reads_per_million_nonhost"] == pytest.approx(0.1)
    assert list(summary["samples_with_hits"]) == ["SRR2", "SRR3"]


def test_false_alarm_summary_of_a_clean_cohort():
    """Zero is a result, and it has to survive the summary intact.

    A false-alarm rate of exactly zero across the cohort is the case Day 5 has to
    handle most carefully — it does not mean the rate is zero, it means it is below
    what twenty samples can measure — so the summary must not quietly drop it.
    """
    summary = day3.summarise_false_alarms([negative("SRR1", 0, 10_000_000)])
    assert summary["reads_total"] == 0
    assert summary["samples_with_at_least_one_hit"] == 0
    assert summary["reads_per_million_nonhost"] == 0
    assert summary["samples_with_hits"] == {}


def test_false_alarm_summary_of_nothing_is_empty():
    assert day3.summarise_false_alarms([]) == {}


def test_merge_accumulates_negatives_across_invocations(tmp_path):
    path = tmp_path / "day3-provenance.json"
    day3.merge(
        path,
        {
            "generated_at": "2026-11-18T07:00:00+00:00",
            "wall_clock_seconds": 10,
            "negatives": [negative("SRR1", 0, 1000)],
        },
    )
    merged = day3.merge(
        path,
        {
            "generated_at": "2026-11-18T08:00:00+00:00",
            "wall_clock_seconds": 20,
            "negatives": [negative("SRR2", 2, 1000)],
        },
    )
    assert [n["accession"] for n in merged["negatives"]] == ["SRR1", "SRR2"]
    assert len(merged["runs"]) == 2
    assert merged["false_alarms"]["reads_total"] == 2
    on_disk = json.loads(path.read_text())
    assert on_disk == merged


def test_merge_replaces_a_rerun_sample_rather_than_duplicating_it(tmp_path):
    path = tmp_path / "day3-provenance.json"
    day3.merge(path, {"negatives": [negative("SRR1", 0, 1000)]})
    merged = day3.merge(path, {"negatives": [negative("SRR1", 5, 2000)]})
    assert len(merged["negatives"]) == 1
    assert merged["negatives"][0]["false_alarm_reads"] == 5
    assert merged["false_alarms"]["non_host_pairs_total"] == 2000


def test_merge_keeps_a_stage_it_was_not_given(tmp_path):
    """Running one stage must not delete another stage's results from the record."""
    path = tmp_path / "day3-provenance.json"
    day3.merge(path, {"index_sample": {"name": "idx", "pairs": 10}})
    merged = day3.merge(path, {"negatives": [negative("SRR1", 1, 1000)]})
    assert merged["index_sample"]["pairs"] == 10
    assert merged["negatives"][0]["accession"] == "SRR1"


def test_pruning_removes_reads_and_keeps_the_record(tmp_path):
    sample = tmp_path / "SRR1"
    for sub, name in (
        ("raw", "SRR1_1.fastq.gz"),
        ("qc", "SRR1_qc_1.fastq.gz"),
        ("nonhost", "SRR1_nonhost_1.fastq.gz"),
        ("sra", "SRR1.sra"),
    ):
        directory = sample / sub
        directory.mkdir(parents=True)
        (directory / name).write_bytes(b"x" * 100)
    (sample / "qc" / "fastp.json").write_text("{}")
    (sample / ".host.done").write_text("{}")

    freed = day3._prune_sample(sample)
    assert freed == 400
    assert not list(sample.glob("*/*.fastq.gz"))
    assert not list(sample.glob("*/*.sra"))
    # The things that make the sample's numbers auditable stay.
    assert (sample / "qc" / "fastp.json").exists()
    assert (sample / ".host.done").exists()


def test_pruning_a_sample_with_nothing_left_is_free(tmp_path):
    assert day3._prune_sample(tmp_path / "missing") == 0
