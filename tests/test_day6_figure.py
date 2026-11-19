"""The README's figure, which is the one number-bearing artefact a reader sees first.

A figure drifts in a way a table does not: the prose beside it gets corrected, the
picture does not, and nobody diffs an SVG. So the tests here are about the two failures
that would publish a wrong picture — a point plotted from something other than the
record, and a figure left behind by a re-score — rather than about how it looks.
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pytest

FIGURE = Path("docs/figures/detection-vs-dilution.svg")
REAL_DAY5 = Path("results/day5-provenance.json")


def _load_module():
    spec = importlib.util.spec_from_file_location("day6_figure", "scripts/day6_figure.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


day6_figure = _load_module()


@pytest.fixture(scope="module")
def provenance() -> dict:
    if not REAL_DAY5.exists():
        pytest.skip("Day 5 has not been scored in this checkout")
    return json.loads(REAL_DAY5.read_text())


def test_the_committed_figure_is_what_the_script_draws_from_the_record(provenance):
    # The failure this catches is a re-score that moves a threshold and leaves the
    # picture asserting the old one.
    if not FIGURE.exists():
        pytest.skip("figure not built in this checkout")
    assert FIGURE.read_text() == day6_figure.build(provenance)


def test_every_plotted_point_comes_from_the_record(provenance):
    svg = day6_figure.build(provenance)
    levels = provenance["cross_arm"]["levels"]
    a_thresh = provenance["arm_a"]["headline"]["threshold_reads"]
    b_thresh = provenance["arm_b"]["headline"]["threshold_bitscore"]

    # Five rungs per arm, in both panels: twenty markers and no more, so a stray series
    # cannot be drawn without the count changing.
    assert svg.count("<circle") == 4 * len(levels) + 1 + 2  # markers, the ringed point, legend
    assert str(int(a_thresh)) in svg and str(int(b_thresh)) in svg

    # The figure's claim is the alarm line: every rung above it alarms in the record.
    for level in levels:
        margin_a = level["arm_a_reads"] / a_thresh
        margin_b = level["arm_b_best_bitscore"] / b_thresh
        assert level["arm_a_alarm"] == (margin_a >= 1)
        assert level["arm_b_alarm"] == (margin_b >= 1)


def test_the_crossover_marker_sits_where_the_recovery_curves_actually_cross(provenance):
    # Drawn by interpolation rather than by eye, because the rung spacing is logarithmic
    # and a marker placed at the midpoint of two rungs is a different claim.
    levels = provenance["cross_arm"]["levels"]
    low, above = levels[4], levels[3]

    def recovery(level, key):
        return 100 * level[key] / level["viral_fragments"]

    a_low, a_above = recovery(low, "arm_a_reads"), recovery(above, "arm_a_reads")
    b_low, b_above = (
        recovery(low, "arm_b_reads_recovered"),
        recovery(above, "arm_b_reads_recovered"),
    )

    # Arm B is behind at the bottom rung and ahead at the one above it: that bracket is
    # the whole reason a crossover exists, and it is what the figure annotates.
    assert b_low < a_low
    assert b_above > a_above

    t = (a_low - b_low) / ((a_low - b_low) + (b_above - a_above))
    assert 0 < t < 1


def test_the_floor_band_covers_the_rungs_the_series_cannot_resolve(provenance):
    floor = provenance["resident_background"]["resolvable_floor_per_million"]
    rungs = provenance["resident_background"]["rungs"]
    shaded = [r for r in rungs if r["labelled_per_million"] < floor]
    # Exactly one rung is inside the floor, and it is the one whose dose is 90% resident
    # background. If a later series changes that, the figure's caption changes with it.
    assert len(shaded) == 1
    assert math.isclose(shaded[0]["unlabelled_share"], 0.9, abs_tol=0.01)
    assert f"{floor:.1f}" in day6_figure.build(provenance)
