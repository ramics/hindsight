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


def negative_with_taxa(accession, hits, endemic, pairs):
    return {
        "accession": accession,
        "false_alarm_reads": hits,
        "false_alarm_reads_excluding_endemic": hits - endemic,
        "endemic_human_coronavirus_reads": endemic,
        "pairs": pairs,
    }


def test_false_alarm_summary_reports_both_readings():
    """A read placed on OC43 is the arm working; one on a bat sarbecovirus is not.

    Which of the two counts is "the false-alarm rate" is Day 5's decision, so the
    summary has to carry both without picking.
    """
    summary = day3.summarise_false_alarms(
        [
            negative_with_taxa("SRR1", 11, 11, 14_000_000),  # all seasonal coronavirus
            negative_with_taxa("SRR2", 2, 0, 10_000_000),  # neither endemic nor known
            negative_with_taxa("SRR3", 0, 0, 10_000_000),
        ]
    )
    assert summary["reads_total"] == 13
    assert summary["samples_with_at_least_one_hit"] == 2
    assert summary["endemic_human_coronavirus_reads_total"] == 11

    novel = summary["excluding_endemic_human_coronaviruses"]
    assert novel["reads_total"] == 2
    assert novel["samples_with_at_least_one_hit"] == 1
    assert list(novel["samples_with_hits"]) == ["SRR2"]


def test_endemic_list_excludes_sars_related_coronavirus():
    """Treating "known in 2019" as "not an alarm" would discard the benchmark's signal.

    SARS-CoV was known, and the index sample's reads land on SARS-related coronavirus.
    The endemic list is the four seasonal human coronaviruses, which are a different
    argument: they are common, they really are in respiratory samples, and they are not
    what a novel-pathogen alarm is about.
    """
    from hindsight.arms import kraken

    assert 694009 not in kraken.ENDEMIC_HUMAN_CORONAVIRUSES.values()
    assert 1335626 not in kraken.ENDEMIC_HUMAN_CORONAVIRUSES.values()  # MERS
    assert 31631 in kraken.ENDEMIC_HUMAN_CORONAVIRUSES.values()  # OC43


# ------------------------------------------------- the k-mer length is a settable knob


def test_kmer_settings_default_to_kraken2s_own_and_reach_the_builder(tmp_path, monkeypatch):
    """The README says a shorter k is an experiment someone can run; this is why.

    The knob is worth a test rather than a comment because the failure is silent: a
    config that sets ``kmer_length`` and a builder that never passes it on would produce
    a default database under an experiment's name, and every number from it would look
    like evidence about k.
    """
    from hindsight.arms import day3, kraken, taxa

    config = tmp_path / "arm.toml"
    config.write_text(
        '[arm]\nsamples_config = "config/day2.toml"\n' '[database]\nsource_db = "db/nov2019"\n'
    )
    assert day3.load_config(config).kmer_length is None
    assert day3.load_config(config).minimizer_length is None

    config.write_text(
        '[arm]\nsamples_config = "config/day2.toml"\n'
        '[database]\nsource_db = "db/nov2019"\n'
        'kraken_db = "work/kraken2/k27"\nkmer_length = 27\nminimizer_length = 24\n'
    )
    arm = day3.load_config(config)
    assert (arm.kmer_length, arm.minimizer_length) == (27, 24)

    seen: dict = {}

    def fake_build(source_db, db_dir, **kwargs):
        seen.update(kwargs)
        return {"db": str(db_dir), "files": {}, "inspect": {}}

    taxonomy = taxa.Taxonomy(parents={}, ranks={}, names={})
    monkeypatch.setattr(kraken, "build_database", fake_build)
    monkeypatch.setattr(kraken, "read_taxonomy", lambda path: taxonomy)
    monkeypatch.setattr(kraken, "clade_sets", lambda _taxonomy: taxa.CladeSets(members={}))
    day3.prepare_database(arm, log=lambda _msg: None)
    assert seen["kmer_length"] == 27
    assert seen["minimizer_length"] == 24
