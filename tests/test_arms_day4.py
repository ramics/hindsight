"""Day 4's driver: what it reads from config, and how it accumulates a record.

The merge is the part worth testing hard. Arm B's series is six assemblies and the
negatives are twenty more, so a Day 4 record is never built in one invocation — a
crash, a resumed level, a re-run of one negative — and a merge that quietly dropped the
levels already done would lose hours of compute without erroring.
"""

from __future__ import annotations

import json

from hindsight.arms import day4


def _config(tmp_path, body: str = ""):
    path = tmp_path / "day4.toml"
    path.write_text(body)
    return day4.load_config(path)


def test_defaults_are_the_documented_ones(tmp_path):
    arm = _config(tmp_path)
    assert arm.sensitivity == "very-sensitive"
    assert arm.comparison_sensitivity == "default"
    assert arm.read_level is True
    assert arm.negatives_read_level is False
    assert arm.negatives_max_pairs == 0  # native depth
    assert arm.lca_top_percent == 10.0
    assert arm.min_contig_len == 200
    assert arm.samples_config.name == "day2.toml"


def test_the_repo_config_is_what_the_write_up_describes():
    """config/day4.toml is quoted in the day's write-up, so it is checked here."""
    arm = day4.load_config("config/day4.toml")
    assert arm.sensitivity == "very-sensitive"
    assert arm.comparison_sensitivity == "default"
    assert arm.evalue == 0.001
    assert arm.max_target_seqs == 25
    assert arm.read_level is True
    assert arm.negatives_max_pairs == 0
    assert arm.negatives_read_level is False
    assert arm.prune_after_search is True
    assert arm.tool_pins == {"diamond": "2.1.9", "megahit": "1.2.9", "bowtie2": "2.5.2"}


def test_a_blank_comparison_sensitivity_means_no_comparison_run(tmp_path):
    arm = _config(tmp_path, '[search]\ncomparison_sensitivity = ""\n')
    assert arm.comparison_sensitivity is None


def test_headline_clade_is_the_family(tmp_path):
    """Day 3's recommendation, carried into Arm B rather than re-decided.

    The November 2019 taxonomy puts the dated database's nearest relative to SARS-CoV-2
    outside Betacoronavirus, so a genus-level headline would under-report both arms.
    """
    assert day4.HEADLINE_CLADE == "coronaviridae"


# ----------------------------------------------------------------------------- merge


def test_merge_writes_a_fresh_record_with_one_stamp(tmp_path):
    path = tmp_path / "day4.json"
    merged = day4.merge(path, {"generated_at": "now", "stages": ["database"], "database": {}})
    assert merged["runs"] == [
        {"generated_at": "now", "stages": ["database"], "wall_clock_seconds": None, "negatives": []}
    ]
    assert json.loads(path.read_text())["runs"] == merged["runs"]


def test_merge_keeps_levels_from_an_earlier_partial_series_run(tmp_path):
    """The expensive case: one level re-run must not erase the other five."""
    path = tmp_path / "day4.json"
    day4.merge(
        path,
        {
            "generated_at": "first",
            "series": {
                "background": {"name": "background"},
                "levels": [
                    {"name": "viral050000", "contigs_in_headline_clade": 1},
                    {"name": "viral010000", "contigs_in_headline_clade": 1},
                ],
            },
        },
    )
    merged = day4.merge(
        path,
        {
            "generated_at": "second",
            "series": {"levels": [{"name": "viral001000", "contigs_in_headline_clade": 0}]},
        },
    )
    names = [lv["name"] for lv in merged["series"]["levels"]]
    assert names == ["viral050000", "viral010000", "viral001000"]
    # The background was not in the second run and must survive it.
    assert merged["series"]["background"]["name"] == "background"
    assert len(merged["runs"]) == 2


def test_merge_replaces_a_level_that_was_run_again(tmp_path):
    path = tmp_path / "day4.json"
    day4.merge(path, {"series": {"levels": [{"name": "viral001000", "contigs": 0}]}})
    merged = day4.merge(path, {"series": {"levels": [{"name": "viral001000", "contigs": 4}]}})
    assert merged["series"]["levels"] == [{"name": "viral001000", "contigs": 4}]


def test_merge_accumulates_negatives_by_accession_and_resummarises(tmp_path):
    path = tmp_path / "day4.json"
    day4.merge(
        path,
        {
            "negatives": [
                {"accession": "SRR2", "false_alarm_contigs": 0, "nonhost_pairs": 10},
            ]
        },
    )
    merged = day4.merge(
        path,
        {
            "negatives": [
                {"accession": "SRR1", "false_alarm_contigs": 2, "nonhost_pairs": 20},
                # A re-run of SRR2 with a different answer replaces it rather than
                # appearing twice.
                {"accession": "SRR2", "false_alarm_contigs": 1, "nonhost_pairs": 10},
            ]
        },
    )
    assert [n["accession"] for n in merged["negatives"]] == ["SRR1", "SRR2"]
    assert merged["false_alarms"]["samples"] == 2
    assert merged["false_alarms"]["contigs_total"] == 3
    assert merged["false_alarms"]["non_host_pairs_total"] == 30


def test_merge_does_not_blank_an_existing_stage_with_an_empty_one(tmp_path):
    path = tmp_path / "day4.json"
    day4.merge(path, {"index_sample": {"accession": "SRR10971381"}, "database": {"db": "x"}})
    merged = day4.merge(path, {"index_sample": {}, "database": {}, "negatives": []})
    assert merged["index_sample"]["accession"] == "SRR10971381"
    assert merged["database"]["db"] == "x"


def test_pruning_keeps_the_largest_figure_across_reruns(tmp_path):
    """A second prune frees 0 bytes, and must not overwrite what the first one freed."""
    assert day4._remember_pruned(tmp_path, "SRR1", 5_000_000) == 5_000_000
    assert day4._remember_pruned(tmp_path, "SRR1", 0) == 5_000_000


def test_pruning_deletes_reads_and_leaves_the_record(tmp_path):
    sample = tmp_path / "SRR1"
    for sub in ("raw", "qc", "nonhost", "sra"):
        (sample / sub).mkdir(parents=True)
    (sample / "raw" / "a_1.fastq.gz").write_bytes(b"x" * 100)
    (sample / "nonhost" / "a_nonhost_1.fastq.gz").write_bytes(b"y" * 50)
    (sample / "sra" / "a.sra").write_bytes(b"z" * 200)
    # The read-half FASTA is Arm B's own intermediate and is pruned too.
    (sample / "nonhost" / "a.reads.fa.gz").write_bytes(b"w" * 25)
    (sample / "nonhost" / ".keep.done").write_text("marker")
    freed = day4._prune_sample(sample)
    assert freed == 375
    assert not (sample / "sra" / "a.sra").exists()
    assert (sample / "nonhost" / ".keep.done").exists()
