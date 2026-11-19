"""Day 5's arithmetic, which is the only thing between two records and a published limit.

The tests worth having here are about the places where a plausible-looking
implementation would be quietly wrong: a threshold that admits one more negative than it
means to, a zero rung that makes every sample alarm, an endemic exclusion that keeps an
infection's own spillover, and a limit of detection read off a curve that is not
monotone. Each of those would produce a table that looks finished.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from hindsight.arms import day5

REAL_DAY3 = Path("results/day3-provenance.json")
REAL_DAY4 = Path("results/day4-provenance.json")


# ------------------------------------------------------------------- the thresholds


def test_threshold_is_the_smallest_one_the_budget_allows():
    # Two negatives tie at 11, so a threshold of 11 would alarm on both: with no budget
    # the answer has to be 12 rather than 11.
    assert day5.threshold_for_budget([11, 11, 8, 2, 1], 0) == 12
    # Allowing one false alarm does not help here, because the tie means dropping to 11
    # costs two alarms, not one.
    assert day5.threshold_for_budget([11, 11, 8, 2, 1], 1) == 12
    assert day5.threshold_for_budget([11, 11, 8, 2, 1], 2) == 9
    # A cohort that produced nothing cannot justify a threshold above one read.
    assert day5.threshold_for_budget([0, 0, 0], 0) == 1
    assert day5.threshold_for_budget([], 0) == 1
    with pytest.raises(ValueError):
        day5.threshold_for_budget([1], -1)


def test_ladder_threshold_only_returns_a_pre_registered_rung():
    ladder = (0.0, 50.0, 100.0, 250.0, 500.0)
    # 226 is below 250, so 250 is the first rung that clears it.
    assert day5.ladder_threshold_for_budget([226.0, None, None], 0, ladder) == 250.0
    # With a budget of one, the arm's own unfiltered rule is already tight enough.
    assert day5.ladder_threshold_for_budget([226.0, None, None], 1, ladder) == 0.0
    # A negative above the top of the ladder cannot be excluded by it, and saying so is
    # the result: it is what the "all" framing does to Arm B.
    assert day5.ladder_threshold_for_budget([1188.0], 0, ladder) is None


def test_a_sample_with_no_contig_does_not_alarm_at_the_zero_rung():
    # The trap: `(None or 0.0) >= 0.0` is true, which would alarm on all twenty
    # negatives at the arm's own unfiltered rule.
    assert day5._arm_b_alarms({"contigs": 0, "best_bitscore": None}, 0.0) is False
    assert day5._arm_b_alarms({"contigs": 1, "best_bitscore": 40.0}, 0.0) is True
    assert day5._arm_b_alarms({"contigs": 1, "best_bitscore": 40.0}, 50.0) is False


# ------------------------------------------------------------------- the framings


def _negative(accession, pairs, clade_reads, endemic, taxa, recorded=0):
    return {
        "accession": accession,
        "pairs": pairs,
        "clades": {"coronaviridae": clade_reads, "betacoronavirus": clade_reads},
        "endemic_human_coronavirus_reads": endemic,
        "coronavirus_taxa": taxa,
        "false_alarm_reads_excluding_endemic": recorded,
    }


def test_endemic_exclusion_takes_the_larger_of_the_two_reconstructions():
    # SRR13677744's shape: one camel alphacoronavirus read in the family, and an endemic
    # read counted cumulatively outside it. Subtraction says zero, enumeration says one,
    # and the permissive answer is one — it can only make the arm look worse.
    day3 = {
        "false_alarms": {"clade": "betacoronavirus"},
        "negatives": [
            _negative("SRR_A", 18_279_498, 1, 1, [{"taxid": 1699095, "reads": 1}]),
        ],
    }
    (score,) = day5.arm_a_negative_scores(day3, "coronaviridae", "excluding_endemic")
    assert score["reads"] == 1
    assert score["estimates"] == {"subtraction": 0, "enumeration": 1}
    assert score["estimates_disagree"] is True
    # The "all" framing is the raw cumulative count and carries no estimates.
    (raw,) = day5.arm_a_negative_scores(day3, "coronaviridae", "all")
    assert raw["reads"] == 1 and "estimates" not in raw


def test_endemic_reads_are_struck_out_by_taxon():
    day3 = {
        "false_alarms": {"clade": "betacoronavirus"},
        "negatives": [
            _negative(
                "SRR_B",
                10_000_000,
                57,
                46,
                [
                    {"taxid": 277944, "reads": 46},  # NL63, endemic
                    {"taxid": 1159906, "reads": 11},  # sparrow coronavirus HKU17
                ],
            )
        ],
    }
    (score,) = day5.arm_a_negative_scores(day3, "coronaviridae", "excluding_endemic")
    assert score["reads"] == 11
    assert score["reads_per_million_nonhost"] == pytest.approx(1.1, abs=0.01)


def _contig(name, length, lca_taxid, best_taxid, bitscore):
    return {
        "contig": name,
        "length": length,
        "lca_taxid": lca_taxid,
        "lca_clades": ["viruses", "coronaviridae", "betacoronavirus"],
        "best_hit": {"taxid": best_taxid, "bitscore": bitscore},
    }


def test_attribution_decides_whether_an_infections_spillover_is_a_false_alarm():
    # The HKU1 patient's shape: one contig of HKU1 sequence whose LCA resolved up to the
    # genus. By LCA it is not endemic and counts as a false alarm; by best hit it is the
    # cold the patient had.
    sample = {
        "search": {
            "contig_detail": [
                _contig("k141_7872", 2652, 290028, 290028, 832.0),  # HKU1 / HKU1
                _contig("k141_2193", 779, 694002, 290028, 494.0),  # genus LCA / HKU1
                _contig("k141_852", 767, 694002, 502102, 226.0),  # genus LCA / rat CoV
            ]
        }
    }
    by_lca = day5.arm_b_sample_evidence(
        sample, "coronaviridae", "as_run", "excluding_endemic", "lca"
    )
    by_best = day5.arm_b_sample_evidence(
        sample, "coronaviridae", "as_run", "excluding_endemic", "best_hit"
    )
    assert by_lca["contigs"] == 2 and by_lca["best_bitscore"] == 494.0
    assert by_best["contigs"] == 1 and by_best["best_bitscore"] == 226.0
    assert by_best["endemic_contigs"] == 2
    # And the unfiltered framing keeps everything.
    everything = day5.arm_b_sample_evidence(sample, "coronaviridae", "as_run", "all", "best_hit")
    assert everything["contigs"] == 3 and everything["best_bitscore"] == 832.0


def test_the_comparison_search_is_read_from_its_own_section():
    # Day 4 records the stock-mode run beside the as-run one, with the contig detail at
    # the top of its section rather than inside a search. Reading the wrong one would
    # silently score --very-sensitive while reporting stock.
    sample = {
        "search": {"contig_detail": [_contig("a", 975, 1, 1, 57.8)]},
        "sensitivity_comparison": {
            "sensitivity": "default",
            "contig_detail": [],
            "search": {"clades": {}},
        },
    }
    as_run = day5.arm_b_sample_evidence(sample, "coronaviridae", "as_run", "all", "best_hit")
    stock = day5.arm_b_sample_evidence(sample, "coronaviridae", "comparison", "all", "best_hit")
    assert as_run["contigs"] == 1 and stock["contigs"] == 0


# --------------------------------------------------------------- limits of detection


def _curve(pairs):
    return [
        {"name": name, "viral_fragments": n, "viral_fragments_per_million": n / 10, "score": s}
        for name, n, s in pairs
    ]


def test_limit_of_detection_brackets_the_lowest_firing_rung():
    curve = _curve([("a", 1000, 50), ("b", 100, 20), ("c", 10, 1)])
    lod = day5.limit_of_detection(curve, lambda lv: lv["score"] >= 10)
    assert lod["alarms_at"] == "b"
    assert lod["viral_fragments"] == 100
    assert lod["highest_silent_level"] == "c"
    assert lod["at_or_below_series_floor"] is False
    assert lod["monotone"] is True


def test_a_limit_at_the_bottom_of_the_series_says_so():
    curve = _curve([("a", 1000, 50), ("b", 100, 20), ("c", 10, 11)])
    lod = day5.limit_of_detection(curve, lambda lv: lv["score"] >= 10)
    assert lod["alarms_at"] == "c"
    assert lod["at_or_below_series_floor"] is True
    assert lod["highest_silent_level"] is None


def test_a_non_monotone_curve_is_reported_rather_than_smoothed():
    curve = _curve([("a", 1000, 1), ("b", 100, 50), ("c", 10, 1)])
    lod = day5.limit_of_detection(curve, lambda lv: lv["score"] >= 10)
    assert lod["monotone"] is False
    assert [lv["alarm"] for lv in lod["levels"]] == [False, True, False]


# ------------------------------------------------------------------------- config


def test_config_defaults_are_the_documented_decisions():
    config = day5.load_config(Path("config/day5.toml"))
    assert config.clade == "coronaviridae"
    assert config.framing == "excluding_endemic"
    assert config.attribution == "best_hit"
    assert config.arm_b_sensitivity == "comparison"
    assert config.false_alarm_budget == 0
    assert config.denominator == "post_qc_pairs"
    assert config.bitscore_ladder == (0.0, 50.0, 100.0, 250.0, 500.0, 1000.0)


def test_config_rejects_what_it_cannot_honour(tmp_path):
    path = tmp_path / "day5.toml"
    path.write_text('[scoring]\nclade = "flaviviridae"\n')
    with pytest.raises(ValueError, match="unknown clade"):
        day5.load_config(path)
    path.write_text("[scoring]\nmystery = 1\n")
    with pytest.raises(ValueError, match="unknown key"):
        day5.load_config(path)


# --------------------------------------------------- the real records, as a regression


requires_records = pytest.mark.skipif(
    not (REAL_DAY3.exists() and REAL_DAY4.exists()),
    reason="needs the committed Day 3 and Day 4 provenance records",
)


@requires_records
def test_scoring_the_real_records_runs_no_tool(tmp_path, monkeypatch):
    # The guarantee worth asserting is that no classifier, assembler, aligner or
    # fetcher is invoked. `platform.platform()` does shell out to `uname` while
    # recording the machine, so the guard names the tools rather than forbidding
    # every subprocess.
    forbidden = {"kraken2", "megahit", "diamond", "bowtie2", "minimap2", "fasterq-dump"}
    real_run, real_popen = subprocess.run, subprocess.Popen

    def guard(real):
        def wrapper(args, *rest, **kwargs):
            parts = args if isinstance(args, list | tuple) else [args]
            names = {Path(str(part)).name for part in parts}
            if names & forbidden:
                raise AssertionError(f"Day 5 scores from the records and ran {names & forbidden}")
            return real(args, *rest, **kwargs)

        return wrapper

    monkeypatch.setattr(subprocess, "run", guard(real_run))
    monkeypatch.setattr(subprocess, "Popen", guard(real_popen))
    out = tmp_path / "day5.json"
    record = day5.score(day5.Scoring(day3=REAL_DAY3, day4=REAL_DAY4, provenance=out))
    assert json.loads(out.read_text())["parameters"]["clade"] == "coronaviridae"
    assert record["parameters"]["compute"].startswith("none")


@requires_records
def test_the_headline_numbers_are_the_ones_the_write_up_quotes():
    out = Path("results/day5-provenance.json")
    record = day5.score(day5.Scoring(day3=REAL_DAY3, day4=REAL_DAY4, provenance=out))
    arm_a = record["arm_a"]["headline"]
    arm_b = record["arm_b"]["headline"]
    assert arm_a["threshold_reads"] == 12
    assert arm_a["false_alarm_samples"] == 0
    assert arm_a["limit_of_detection"]["alarms_at"] == "viral000100"
    assert arm_a["limit_of_detection"]["per_million"] == pytest.approx(7.71, abs=0.01)
    assert arm_b["threshold_bitscore"] == 250.0
    assert arm_b["false_alarm_samples"] == 0
    assert arm_b["limit_of_detection"]["at_or_below_series_floor"] is True
    # Arm A's four reads at the bottom rung are a recovery, not an alarm; Arm B's single
    # 462 bp contig is an alarm. That is the day's cross-arm finding.
    bottom = record["cross_arm"]["levels"][-1]
    assert bottom["viral_fragments"] == 10
    assert bottom["arm_a_alarm"] is False and bottom["arm_b_alarm"] is True


@requires_records
def test_the_series_background_offsets_the_dose_axis():
    day3 = json.loads(REAL_DAY3.read_text())
    day4 = json.loads(REAL_DAY4.read_text())
    residents = day5.resident_background(day3, day4, "coronaviridae", "post_qc_pairs")
    assert residents["arm_a_background_clade_reads"] == 90
    assert residents["verified_as_sars_cov_2"] == 71
    bottom = residents["rungs"][-1]
    assert bottom["labelled_fragments"] == 10
    assert bottom["effective_fragments"] == 100
    # Nine tenths of the bottom rung's dose was never diluted, which is why the limit is
    # quoted as "at or below" rather than as a point.
    assert bottom["unlabelled_share"] == pytest.approx(0.9, abs=0.01)


@requires_records
def test_every_decision_is_recorded_with_its_alternatives():
    record = json.loads(Path("results/day5-provenance.json").read_text())
    decisions = record["decisions"]
    assert decisions["clade"]["chose"] == "coronaviridae"
    assert decisions["framing"]["chose"] == "excluding_endemic"
    assert decisions["sensitivity"]["chose"] == "stock"
    # The genus-level rows are what make the clade decision checkable: Arm B never
    # alarms there, on a sample where it rebuilt the genome.
    genus = next(r for r in decisions["clade"]["by_clade"] if r["clade"] == "betacoronavirus")
    assert genus["arm_b_lod_per_million"] is None
    # And the grid carries the rows nobody chose.
    assert len(record["operating_points"]) > 40
