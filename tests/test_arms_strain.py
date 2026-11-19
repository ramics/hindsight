"""The strain-versus-coinfection diagnosis.

Two things here are load-bearing and neither is obvious from reading the output. The
first is the cutoff gate: this module exists to reach sequence the dated database does
not hold, so the one thing standing between it and a leak is the check that a reference
predates the cutoff. The second is the verdict logic, which has to say "one organism"
when contigs tile a genome and "more than one" when they double up on it -- and a bug
either way would be invisible, because both answers are plausible-looking prose.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from hindsight.arms import strain


class FakeClient:
    """Stands in for E-utilities, with a creation date per accession."""

    def __init__(self, dates: dict[str, str], sequences: dict[str, str] | None = None):
        self.dates = dates
        self.sequences = sequences or {}
        self.efetched: list[str] = []

    def esummary(self, db: str, ids: list[str]) -> str:
        accession = ids[0]
        date = self.dates.get(accession)
        if date is None:
            return "<eSummaryResult><DocSum></DocSum></eSummaryResult>"
        return f'<DocSum><Item Name="CreateDate" Type="String">{date}</Item></DocSum>'

    def efetch(self, db: str, ids: list[str], **kwargs) -> str:
        self.efetched.append(ids[0])
        return f">{ids[0]}\n{self.sequences.get(ids[0], 'ACGT' * 10)}\n"


def _manifest(tmp_path: Path, accessions=("NC_006577.2",)) -> Path:
    db = tmp_path / "db"
    db.mkdir()
    lines = ["accession\tsomething"] + [f"{a}\tx" for a in accessions]
    (db / "manifest.tsv").write_text("\n".join(lines) + "\n")
    return db


# ------------------------------------------------------------------ the cutoff gate


def test_a_reference_in_the_dated_manifest_needs_no_date_check(tmp_path):
    client = FakeClient(dates={})
    reference = strain.resolve_reference(
        "NC_006577.2",
        source_db=_manifest(tmp_path),
        cutoff="2019-12-01",
        dest_dir=tmp_path / "ref",
        client=client,
        log=lambda _msg: None,
    )
    assert reference.source == "dated-manifest"
    assert reference.before_cutoff is True
    assert reference.published is None


def test_a_reference_outside_the_database_is_allowed_when_it_predates_the_cutoff(tmp_path):
    """The whole point: outside RefSeq by curation policy is not outside by date."""
    client = FakeClient(dates={"DQ415911.1": "2006/06/30"})
    reference = strain.resolve_reference(
        "DQ415911.1",
        source_db=_manifest(tmp_path),
        cutoff="2019-12-01",
        dest_dir=tmp_path / "ref",
        client=client,
        log=lambda _msg: None,
    )
    assert reference.source == "pre-cutoff-genbank"
    assert reference.published == "2006-06-30"
    assert reference.before_cutoff is True
    assert client.efetched == ["DQ415911.1"]


def test_a_post_cutoff_reference_is_refused_by_default(tmp_path):
    client = FakeClient(dates={"MN908947.3": "2020/01/12"})
    with pytest.raises(strain.StrainError, match="2020-01-12"):
        strain.resolve_reference(
            "MN908947.3",
            source_db=_manifest(tmp_path),
            cutoff="2019-12-01",
            dest_dir=tmp_path / "ref",
            client=client,
            log=lambda _msg: None,
        )
    assert client.efetched == []  # refused before anything was downloaded


def test_a_post_cutoff_reference_can_be_used_deliberately_and_is_recorded_as_such(tmp_path):
    client = FakeClient(dates={"MN908947.3": "2020/01/12"})
    reference = strain.resolve_reference(
        "MN908947.3",
        source_db=_manifest(tmp_path),
        cutoff="2019-12-01",
        dest_dir=tmp_path / "ref",
        client=client,
        allow_post_cutoff=True,
        log=lambda _msg: None,
    )
    assert reference.before_cutoff is False
    assert reference.published == "2020-01-12"


def test_a_reference_with_no_creation_date_is_refused_rather_than_assumed_old(tmp_path):
    client = FakeClient(dates={})
    with pytest.raises(strain.StrainError, match="no creation date"):
        strain.resolve_reference(
            "XX999999.1",
            source_db=_manifest(tmp_path),
            cutoff="2019-12-01",
            dest_dir=tmp_path / "ref",
            client=client,
            log=lambda _msg: None,
        )


# ------------------------------------------------------------------- coverage parsing


def test_coverage_comes_from_the_megahit_header(tmp_path):
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(
        ">k141_61 flag=1 multi=10.0000 len=535\nACGT\n"
        ">k141_628 flag=1 multi=102.0000 len=2509\nACGT\n"
        ">k141_x len=100\nACGT\n"  # no multi= field
    )
    coverage = strain.contig_coverage(contigs)
    assert coverage == {"k141_61": 10.0, "k141_628": 102.0}


# -------------------------------------------------------------------- the verdict


def _panel(entries, **overrides):
    panel = {
        "reference": "DQ415911.1",
        "source": "pre-cutoff-genbank",
        "published": "2006-06-30",
        "before_cutoff": True,
        "contigs_aligned": len(entries),
        "contigs_total": len(entries),
        "aligned_bases": sum(e["alignment_length"] for e in entries),
        "weighted_identity_percent": 99.36,
        "placed_by_short_read_preset": [],
        "per_contig": {e["contig"]: e for e in entries},
    }
    panel.update(overrides)
    return panel


def _entry(contig, start, end, identity):
    return {
        "contig": contig,
        "alignment_length": end - start,
        "matching_bases": end - start,
        "identity_percent": identity,
        "reference_length": 29902,
        "reference_span": [start, end],
        "preset": "asm20",
    }


def test_contigs_that_tile_a_genome_are_not_called_a_coinfection():
    panel = _panel(
        [
            _entry("a", 1000, 3000, 99.4),
            _entry("b", 2950, 5000, 99.5),  # ends meet: a 50 bp overlap
            _entry("c", 6000, 8000, 99.2),
        ]
    )
    assert strain.duplicated_regions(panel) == []


def test_two_contigs_covering_one_region_at_different_identities_are_flagged():
    panel = _panel(
        [
            _entry("a", 1000, 3000, 99.4),
            _entry("b", 1200, 2800, 88.0),  # same gene, far lower identity
        ]
    )
    clashes = strain.duplicated_regions(panel)
    assert len(clashes) == 1
    assert sorted(clashes[0]["contigs"]) == ["a", "b"]
    assert clashes[0]["overlap_bases"] == 1600
    assert clashes[0]["identity_gap_points"] == pytest.approx(11.4)


def test_a_long_overlap_at_the_same_identity_is_not_flagged():
    """Two contigs of the same sequence are a fragmented assembly, not two organisms."""
    panel = _panel([_entry("a", 1000, 3000, 99.4), _entry("b", 1200, 2800, 99.5)])
    assert strain.duplicated_regions(panel) == []


def test_a_short_overlap_at_different_identities_is_contig_ends_meeting():
    panel = _panel([_entry("a", 1000, 3000, 99.4), _entry("b", 2950, 5000, 85.0)])
    assert strain.duplicated_regions(panel) == []


def test_the_conclusion_says_one_strain_when_every_contig_is_accounted_for():
    panel = _panel([_entry("a", 1000, 3000, 99.4)])
    record = {
        "clade": "coronaviridae",
        "cutoff": "2019-12-01",
        "duplicated_regions": [],
        "unplaced_by_every_reference": [],
    }
    text = strain._conclude(record, panel, 1)
    assert "single diverged strain rather than a coinfection" in text
    assert "before the 2019-12-01 cutoff" in text


def test_the_conclusion_names_the_preset_when_a_contig_needed_the_fallback():
    panel = _panel([_entry("a", 1000, 3000, 99.4)], placed_by_short_read_preset=["a"])
    record = {
        "clade": "coronaviridae",
        "cutoff": "2019-12-01",
        "duplicated_regions": [],
        "unplaced_by_every_reference": [],
    }
    text = strain._conclude(record, panel, 1)
    assert "aligner limit" in text


def test_the_conclusion_says_more_than_one_organism_when_a_region_doubles_up():
    panel = _panel([_entry("a", 1000, 3000, 99.4), _entry("b", 1200, 2800, 88.0)])
    record = {
        "clade": "coronaviridae",
        "cutoff": "2019-12-01",
        "duplicated_regions": strain.duplicated_regions(panel),
        "unplaced_by_every_reference": [],
    }
    text = strain._conclude(record, panel, 2)
    assert "more than one organism is present" in text


def test_the_conclusion_points_at_the_leftovers_when_the_panel_cannot_place_them():
    panel = _panel([_entry("a", 1000, 3000, 99.4)], contigs_aligned=1, contigs_total=2)
    record = {
        "clade": "coronaviridae",
        "cutoff": "2019-12-01",
        "duplicated_regions": [],
        "unplaced_by_every_reference": ["b"],
    }
    text = strain._conclude(record, panel, 2)
    assert "where to look for a second organism" in text
