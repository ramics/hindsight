"""Organism or laboratory sequence, from the shape of the coverage.

The verdict rests entirely on two derived numbers -- how much of the reference is covered,
and how many bases the coverage takes to stop at the ends of the largest window -- so those
are what is tested here. The edge width in particular is easy to get subtly wrong: the
first implementation walked outwards starting *at* the boundary, where depth is already
low, and so never found a transition at all and reported every edge as unknown.
"""

from __future__ import annotations

from hindsight.arms import taxon


def _depth(length: int, windows: list[tuple[int, int, int]]) -> list[int]:
    """A depth array with flat-topped windows: (start, end, depth), 0-based inclusive."""
    depth = [0] * length
    for start, end, value in windows:
        for position in range(start, end + 1):
            depth[position] = value
    return depth


# ----------------------------------------------------------------- the coverage profile


def test_a_bounded_window_is_reported_as_one_window_with_its_depth():
    profile = taxon._profile(_depth(8000, [(2000, 3699, 1000)]))
    assert profile["windows"] == 1
    assert profile["covered_bases"] == 1700
    assert profile["covered_percent"] == 21.25
    biggest = profile["largest_windows"][0]
    assert (biggest["start"], biggest["end"]) == (2001, 3700)  # 1-based for reporting
    assert biggest["median_depth"] == 1000


def test_short_dropouts_inside_a_window_do_not_split_it():
    """A few uncovered bases mid-window are sampling, not structure."""
    depth = _depth(8000, [(2000, 2999, 500), (3020, 3699, 500)])
    profile = taxon._profile(depth)
    assert profile["windows"] == 1
    assert profile["largest_windows"][0]["bases"] == 1700


def test_a_real_gap_does_split_a_window():
    depth = _depth(8000, [(2000, 2999, 500), (6000, 6999, 500)])
    assert taxon._profile(depth)["windows"] == 2


def test_coverage_of_nothing_does_not_divide_by_zero():
    profile = taxon._profile([0] * 100)
    assert profile["windows"] == 0
    assert profile["covered_bases"] == 0
    assert profile["largest_windows"] == []


# --------------------------------------------------------------------- the edge width


def test_a_cloned_inserts_edges_are_one_base_wide():
    depth = _depth(8000, [(2000, 3699, 1000)])
    assert taxon._edge_width(depth, 2000, 1000, leftwards=True) == 1
    assert taxon._edge_width(depth, 3699, 1000, leftwards=False) == 1


def test_a_tapering_edge_is_reported_as_wide_rather_than_as_a_step():
    """Coverage of a genome falls away gradually; that is the case this must not call."""
    depth = [0] * 8000
    for position in range(2000, 3700):
        depth[position] = 1000
    for step, position in enumerate(range(3700, 3900)):
        depth[position] = max(0, 1000 - 5 * step)  # 200 bases to fade out
    width = taxon._edge_width(depth, 3899, 1000, leftwards=False)
    assert width is not None
    assert width > taxon.EDGE_IS_A_STEP_BASES


def test_an_edge_is_unknown_rather_than_zero_when_no_transition_is_in_range():
    depth = [1000] * 8000  # covered everywhere: no edge to find
    assert taxon._edge_width(depth, 4000, 1000, leftwards=True) is None


def test_an_empty_window_has_no_measurable_edge():
    assert taxon._edge_width([0] * 100, 50, 0, leftwards=True) is None


# ----------------------------------------------------------------------- the verdict


def _record(coverage: dict, identity: float = 98.8, fragments: int = 11227) -> dict:
    return {
        "fragments_assigned": fragments,
        "taxid": 11665,
        "taxon": "Equine infectious anemia virus",
        "reference": {"accession": "NC_001450.1"},
        "identity_percent": identity,
        "reads_mapped": 20857,
        "coverage": coverage,
    }


def test_low_identity_is_called_homology_rather_than_the_organism():
    coverage = taxon._profile(_depth(8000, [(2000, 3699, 1000)]))
    text = taxon._conclude(_record(coverage, identity=71.0))
    assert "homologous to it rather than from it" in text


def test_coverage_across_the_reference_is_called_an_organism():
    coverage = taxon._profile(_depth(8000, [(100, 7900, 50)]))
    text = taxon._conclude(_record(coverage))
    assert "what an organism present in" in text


def test_one_sharp_bounded_window_is_called_a_cloned_fragment():
    coverage = taxon._profile(_depth(8000, [(2000, 3699, 1000)]))
    text = taxon._conclude(_record(coverage))
    assert "cloned fragment of the reference" in text
    assert "vector, a plasmid or a spiked control" in text


def test_a_bounded_window_with_soft_edges_stops_short_of_calling_it_a_construct():
    depth = [0] * 8000
    for position in range(2000, 3700):
        depth[position] = 1000
    for step, position in enumerate(range(3700, 3900)):
        depth[position] = max(0, 1000 - 5 * step)
    for step, position in enumerate(range(1800, 2000)):
        depth[position] = min(1000, 5 * step)
    text = taxon._conclude(_record(taxon._profile(depth)))
    assert "not clean enough to call it a construct outright" in text
    assert "cloned fragment" not in text
