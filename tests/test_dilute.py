"""The dilution series' promises: deterministic, nested, honest about its size."""

from __future__ import annotations

import pytest

from hindsight.samples import dilute


def ids(n: int, prefix: str = "frag") -> list[str]:
    return [f"{prefix}{i:06d}" for i in range(n)]


def test_ranking_is_deterministic_across_input_order():
    pool = ids(500)
    forward = dilute.rank_reads(pool, seed=7)
    backward = dilute.rank_reads(list(reversed(pool)), seed=7)
    assert forward == backward


def test_ranking_depends_on_the_seed():
    pool = ids(500)
    assert dilute.rank_reads(pool, seed=1) != dilute.rank_reads(pool, seed=2)


def test_ranking_is_pinned_to_sha256_not_to_a_runtime_rng():
    # A hand-computed expectation: if this test fails, the selection rule changed and
    # every frozen series built before the change is a different series.
    assert dilute.rank_reads(["a", "b", "c"], seed=42) == ["c", "a", "b"]


def test_levels_are_nested_prefixes_of_one_ranking():
    plan = dilute.plan_series(ids(1000), levels=(100, 10, 1), seed=3)
    hundred, ten, one = (set(level.read_ids) for level in plan)
    assert one <= ten <= hundred
    assert len(one) == 1 and len(ten) == 10 and len(hundred) == 100


def test_independent_seeds_opt_out_of_nesting():
    plan = dilute.plan_series(ids(1000), levels=(100, 10), seed=3, independent_seeds=True)
    hundred, ten = (set(level.read_ids) for level in plan)
    assert not ten <= hundred
    assert plan[0].seed != plan[1].seed


def test_level_sizes_and_names():
    plan = dilute.plan_series(ids(60000), levels=dilute.DEFAULT_LEVELS)
    assert [level.viral_reads for level in plan] == [50000, 10000, 1000, 100, 10]
    assert [level.name for level in plan] == [
        "viral050000",
        "viral010000",
        "viral001000",
        "viral000100",
        "viral000010",
    ]
    assert not any(level.truncated for level in plan)


def test_a_level_larger_than_the_pool_is_clamped_and_flagged():
    plan = dilute.plan_series(ids(30), levels=(100, 10), seed=1)
    assert plan[0].viral_reads == 30
    assert plan[0].truncated is True
    assert plan[1].truncated is False


def test_selection_digest_identifies_the_chosen_set():
    first = dilute.plan_series(ids(1000), levels=(50,), seed=11)[0]
    same = dilute.plan_series(ids(1000), levels=(50,), seed=11)[0]
    other = dilute.plan_series(ids(1000), levels=(50,), seed=12)[0]
    assert first.selection_sha256 == same.selection_sha256
    assert first.selection_sha256 != other.selection_sha256


def test_duplicate_identifiers_are_refused():
    with pytest.raises(ValueError, match="not unique"):
        dilute.plan_series(["a", "a", "b"], levels=(1,))


def test_negative_level_is_refused():
    with pytest.raises(ValueError, match="negative"):
        dilute.plan_series(ids(10), levels=(-1,))


def test_zero_level_is_a_pure_background_control():
    level = dilute.plan_series(ids(10), levels=(0,))[0]
    assert level.viral_reads == 0
    assert level.read_ids == ()
    assert level.name == "viral000000"


def test_reads_per_million_uses_the_total_the_classifier_sees():
    assert dilute.reads_per_million(1, 999_999) == pytest.approx(1.0)
    assert dilute.reads_per_million(10, 9_999_990) == pytest.approx(1.0)
    with pytest.raises(ValueError):
        dilute.reads_per_million(0, 0)
