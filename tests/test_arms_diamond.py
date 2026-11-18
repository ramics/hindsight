"""Arm B's promises, tested without an assembler or an aligner.

Everything here is about the parts of the arm that could be wrong silently: the LCA that
turns a cluster of near-tied protein hits into one taxon, the taxon map handed to
``diamond makedb``, the way a contig's length is read back, and the composition of the
read half. MEGAHIT and DIAMOND themselves are their authors' business.

The LCA tests carry the most weight. A sarbecovirus contig in a November 2019 database
hits SARS-CoV and several bat coronaviruses at nearly the same bitscore, and the whole
difference between "this contig is Bat coronavirus BM48-31" (precise and wrong) and
"this contig is a coronavirus" (vague and true) lives in which of those hits the arm
lets into the decision.
"""

from __future__ import annotations

import pytest

from hindsight.arms import day4, diamond
from hindsight.arms import taxa as taxa_mod

# A miniature dated taxonomy. The shape matters: BM48-31 hangs off the *family* rather
# than off Betacoronavirus, which is what the November 2019 tree actually does and the
# reason Day 3 recommended thresholding at Coronaviridae.
NODES = """\
1\t|\t1\t|\tno rank\t|
10239\t|\t1\t|\tsuperkingdom\t|
2559587\t|\t10239\t|\tclade\t|
76804\t|\t2559587\t|\torder\t|
11118\t|\t76804\t|\tfamily\t|
2501931\t|\t11118\t|\tsubfamily\t|
694002\t|\t2501931\t|\tgenus\t|
2509511\t|\t694002\t|\tsubgenus\t|
694009\t|\t2509511\t|\tspecies\t|
31631\t|\t694002\t|\tspecies\t|
864596\t|\t11118\t|\tspecies\t|
2\t|\t1\t|\tsuperkingdom\t|
562\t|\t2\t|\tspecies\t|
"""

NAMES = """\
1\t|\troot\t|\t\t|\tscientific name\t|
10239\t|\tViruses\t|\t\t|\tscientific name\t|
11118\t|\tCoronaviridae\t|\t\t|\tscientific name\t|
694002\t|\tBetacoronavirus\t|\t\t|\tscientific name\t|
2509511\t|\tSarbecovirus\t|\t\t|\tscientific name\t|
694009\t|\tSevere acute respiratory syndrome-related coronavirus\t|\t\t|\tscientific name\t|
31631\t|\tHuman coronavirus OC43\t|\t\t|\tscientific name\t|
864596\t|\tBat coronavirus BM48-31/BGR/2008\t|\t\t|\tscientific name\t|
562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|
"""


@pytest.fixture
def taxonomy_dir(tmp_path):
    directory = tmp_path / "taxonomy"
    directory.mkdir()
    (directory / "nodes.dmp").write_text(NODES)
    (directory / "names.dmp").write_text(NAMES)
    return directory


@pytest.fixture
def taxonomy(taxonomy_dir):
    return taxa_mod.read_taxonomy(taxonomy_dir)


@pytest.fixture
def clades(taxonomy):
    return taxa_mod.clade_sets(taxonomy)


# ------------------------------------------------------------------------------- LCA


def test_lca_of_near_ties_climbs_to_the_shared_ancestor(taxonomy):
    """SARS-CoV and OC43 share only the genus, so that is the honest answer."""
    assert taxonomy.lowest_common_ancestor([694009, 31631]) == 694002


def test_lca_with_bm4831_lands_on_the_family_in_the_2019_tree(taxonomy):
    """The finding Day 3 surfaced, now in Arm B's arithmetic.

    BM48-31 is the nearest relative a November 2019 database holds, and in that tree it
    is not a betacoronavirus. A contig tying across it and SARS-CoV therefore resolves
    to Coronaviridae, not to Betacoronavirus — which is exactly why this arm's headline
    clade is the family.
    """
    assert taxonomy.lowest_common_ancestor([694009, 864596]) == 11118


def test_lca_of_a_single_taxon_is_that_taxon(taxonomy):
    assert taxonomy.lowest_common_ancestor([694009]) == 694009


def test_lca_drops_taxids_the_dated_tree_never_held(taxonomy):
    """A hit on a protein whose taxon 2019 does not know must not poison the answer."""
    assert taxonomy.lowest_common_ancestor([694009, 2697049]) == 694009
    assert taxonomy.lowest_common_ancestor([2697049]) is None


def test_lca_across_superkingdoms_is_the_root(taxonomy):
    assert taxonomy.lowest_common_ancestor([694009, 562]) == 1


# ------------------------------------------------------------------- contig assignment


def test_assignment_separates_the_best_hit_from_the_lca(
    tmp_path, taxonomy, clades, blast_row, write_blast_hits
):
    """The two answers disagree, and the arm has to record both.

    The best hit says "this is SARS-CoV", which is wrong: the contig is a virus nobody
    had sequenced. The LCA over the near-ties says "a coronavirus", which is what a
    human reading the table would have concluded and is true.
    """
    hits = write_blast_hits(
        tmp_path / "hits.tsv.gz",
        [
            blast_row("k141_1", "YP_009724390.1", 694009, 900.0, pident=79.0),
            blast_row("k141_1", "YP_003858584.1", 864596, 880.0, pident=78.0),
        ],
    )
    calls = diamond.assign_contigs(hits, {"k141_1": 29800}, taxonomy, clades)
    call = calls["k141_1"]
    assert call.best_taxid == 694009
    assert "sarbecovirus" in call.best_hit_clades
    assert call.lca == 11118
    assert call.lca_clades == ["viruses", "riboviria", "nidovirales", "coronaviridae"]
    assert call.hits_considered == 2
    assert call.length == 29800


def test_a_hit_outside_the_bitscore_band_does_not_widen_the_lca(
    tmp_path, taxonomy, clades, blast_row, write_blast_hits
):
    """The band is what stops one weak, distant hit dragging the answer up to Viruses."""
    hits = write_blast_hits(
        tmp_path / "hits.tsv.gz",
        [
            blast_row("k141_1", "YP_sars", 694009, 900.0),
            blast_row("k141_1", "YP_oc43", 31631, 880.0),
            # 100.0 is far below 90% of 900, so it is evidence about the database, not
            # about this contig.
            blast_row("k141_1", "YP_ecoli", 562, 100.0),
        ],
    )
    calls = diamond.assign_contigs(hits, {"k141_1": 1200}, taxonomy, clades)
    assert calls["k141_1"].lca == 694002
    assert calls["k141_1"].hits_considered == 2
    assert calls["k141_1"].hits_total == 3


def test_widening_the_band_widens_the_answer(
    tmp_path, taxonomy, clades, blast_row, write_blast_hits
):
    """``lca_top_percent`` is a real knob, and the record says what it was set to."""
    hits = write_blast_hits(
        tmp_path / "hits.tsv.gz",
        [
            blast_row("k141_1", "YP_sars", 694009, 900.0),
            blast_row("k141_1", "YP_ecoli", 562, 100.0),
        ],
    )
    narrow = diamond.assign_contigs(hits, {"k141_1": 1200}, taxonomy, clades)
    wide = diamond.assign_contigs(hits, {"k141_1": 1200}, taxonomy, clades, lca_top_percent=99.0)
    assert narrow["k141_1"].lca == 694009
    assert wide["k141_1"].lca == 1


def test_clade_contigs_can_be_selected_either_way(
    tmp_path, taxonomy, clades, blast_row, write_blast_hits
):
    hits = write_blast_hits(
        tmp_path / "hits.tsv.gz",
        [
            blast_row("long", "YP_sars", 694009, 900.0),
            blast_row("long", "YP_bm48", 864596, 890.0),
            blast_row("short", "YP_ecoli", 562, 200.0),
        ],
    )
    calls = diamond.assign_contigs(hits, {"long": 29800, "short": 300}, taxonomy, clades)
    by_lca = diamond.clade_contigs(calls, "coronaviridae", by="lca")
    by_best = diamond.clade_contigs(calls, "betacoronavirus", by="best_hit")
    assert [c.contig for c in by_lca] == ["long"]
    # By LCA the contig is not a betacoronavirus; by best hit it is. Both are recorded.
    assert [c.contig for c in by_best] == ["long"]
    assert diamond.clade_contigs(calls, "betacoronavirus", by="lca") == []


def test_summary_counts_contigs_and_keeps_the_evidence(
    tmp_path, taxonomy, clades, blast_row, write_blast_hits
):
    hits = write_blast_hits(
        tmp_path / "hits.tsv.gz",
        [
            blast_row("long", "YP_sars", 694009, 900.0, pident=79.5),
            blast_row("long", "YP_bm48", 864596, 890.0, pident=78.0),
            blast_row("other", "YP_oc43", 31631, 400.0, pident=99.0, qlen=800),
            blast_row("bug", "YP_ecoli", 562, 300.0, qlen=500),
        ],
    )
    calls = diamond.assign_contigs(
        hits, {"long": 29800, "other": 800, "bug": 500}, taxonomy, clades
    )
    summary = diamond.summarise_calls(calls, clades, taxonomy)
    assert summary["contigs_with_a_hit"] == 3
    assert summary["clades"]["coronaviridae"]["contigs"] == 2
    assert summary["clades"]["coronaviridae"]["longest_contig"] == 29800
    assert summary["clades"]["coronaviridae"]["best_identity_percent"] == 99.0
    assert summary["clades"]["betacoronavirus"]["contigs"] == 1  # OC43 only, by LCA
    detail = {entry["contig"]: entry for entry in summary["contig_detail"]}
    assert detail["long"]["lca_taxon"] == "Coronaviridae"
    assert detail["long"]["best_hit"]["identity_percent"] == 79.5
    assert "bug" not in detail


def test_a_contig_whose_lca_left_the_family_still_appears_in_the_evidence(
    tmp_path, taxonomy, clades, blast_row, write_blast_hits
):
    """Borderline calls are exactly what a reader needs to see, so they are not hidden."""
    hits = write_blast_hits(
        tmp_path / "hits.tsv.gz",
        [
            # Best hit is a coronavirus; a near-tied bacterial hit pushes the LCA to root.
            blast_row("ambiguous", "YP_sars", 694009, 400.0),
            blast_row("ambiguous", "YP_ecoli", 562, 395.0),
        ],
    )
    calls = diamond.assign_contigs(hits, {"ambiguous": 600}, taxonomy, clades)
    assert calls["ambiguous"].lca == 1
    summary = diamond.summarise_calls(calls, clades, taxonomy)
    assert summary["clades"]["coronaviridae"]["contigs"] == 0
    assert [e["contig"] for e in summary["contig_detail"]] == ["ambiguous"]


# --------------------------------------------------------------------------- taxonmap


def test_taxonmap_covers_only_the_proteins_and_widens_to_four_columns(tmp_path):
    """DIAMOND wants prot.accession2taxid's layout, and only for records it will see."""
    (tmp_path / "viral.faa").write_text(">YP_000001.1 spike\nMKT\n>YP_000002.1 nsp\nMAA\n")
    (tmp_path / "seqid2taxid.map").write_text(
        "YP_000001.1\t694009\nYP_000002.1\t864596\nNC_004718.3\t694009\nYP_bad\tnotanumber\n"
    )
    record = diamond.write_taxonmap(
        tmp_path / "seqid2taxid.map", tmp_path / "viral.faa", tmp_path / "map.tsv"
    )
    lines = (tmp_path / "map.tsv").read_text().splitlines()
    assert lines[0] == "accession\taccession.version\ttaxid\tgi"
    assert lines[1:] == ["YP_000001\tYP_000001.1\t694009\t0", "YP_000002\tYP_000002.1\t864596\t0"]
    # The nucleotide accession is in the source map and must not reach DIAMOND.
    assert "NC_004718" not in (tmp_path / "map.tsv").read_text()
    assert record["protein_records"] == 2
    assert record["records_with_a_taxid"] == 2
    assert record["records_without_a_taxid"] == 0


def test_taxonmap_reports_a_protein_with_no_taxid(tmp_path):
    """A silent gap here shows up much later as unassigned contigs, so it is counted."""
    (tmp_path / "viral.faa").write_text(">YP_1.1 a\nMKT\n>YP_2.1 b\nMAA\n")
    (tmp_path / "seqid2taxid.map").write_text("YP_1.1\t694009\n")
    record = diamond.write_taxonmap(
        tmp_path / "seqid2taxid.map", tmp_path / "viral.faa", tmp_path / "map.tsv"
    )
    assert record["protein_records"] == 2
    assert record["records_with_a_taxid"] == 1
    assert record["records_without_a_taxid"] == 1


# --------------------------------------------------------------------------- contigs


def test_contig_lengths_and_n50(tmp_path):
    fasta = tmp_path / "c.fa"
    fasta.write_text(">k141_1 flag=1\nACGT\nACGT\n>k141_2\nAC\n>k141_3\n" + "A" * 10 + "\n")
    lengths = diamond.contig_lengths(fasta)
    assert lengths == {"k141_1": 8, "k141_2": 2, "k141_3": 10}
    # Total 20, half 10: the longest contig alone reaches it.
    assert diamond._n50(list(lengths.values())) == 10


def test_contig_lengths_of_an_empty_assembly(tmp_path):
    fasta = tmp_path / "none.fa"
    fasta.write_text("")
    assert diamond.contig_lengths(fasta) == {}
    assert diamond._n50([]) == 0


def test_contig_subset_keeps_only_the_named_contigs(tmp_path, blast_row):
    fasta = tmp_path / "c.fa"
    fasta.write_text(">keep me\nACGTACGT\n>drop\nTTTT\n>keep2\nGG\n")
    dest = tmp_path / "subset.fa"
    kept = diamond._write_contig_subset(fasta, {"keep", "keep2"}, dest)
    assert kept == 10
    text = dest.read_text()
    assert ">keep me" in text and ">keep2" in text and ">drop" not in text


# ------------------------------------------------------------------- the read half


def test_read_hits_are_pooled_across_mates_into_one_fragment(
    tmp_path, taxonomy, clades, blast_row, write_blast_hits
):
    """Arm A counts a pair as one fragment, so Arm B must too, or the table lies."""
    hits = write_blast_hits(
        tmp_path / "reads.tsv.gz",
        [
            blast_row("FRAG1/1", "YP_sars", 694009, 90.0, qlen=150),
            blast_row("FRAG1/2", "YP_sars", 694009, 85.0, qlen=150),
            blast_row("FRAG2/1", "YP_oc43", 31631, 70.0, qlen=150),
        ],
    )
    record = day4._score_read_hits(hits, taxonomy, clades)
    assert record["fragments_with_a_hit"] == 2
    assert record["fragments_assigned"] == 2
    assert record["clades"]["coronaviridae"] == 2
    assert record["clades"]["sarbecovirus"] == 1


def test_read_hits_scored_against_truth_split_detection_from_false_positives(
    tmp_path, taxonomy, clades, blast_row, write_blast_hits
):
    hits = write_blast_hits(
        tmp_path / "reads.tsv.gz",
        [
            blast_row("REAL/1", "YP_sars", 694009, 90.0, qlen=150),
            blast_row("NOISE/1", "YP_oc43", 31631, 60.0, qlen=150),
        ],
    )
    record = day4._score_read_hits(hits, taxonomy, clades, truth_ids={"REAL", "UNSEEN"})
    truth = record["truth"]
    assert truth["truth_fragments"] == 2
    assert truth["detected_by_clade"]["coronaviridae"] == 1
    assert truth["false_positive_by_clade"]["coronaviridae"] == 1
    assert truth["sensitivity"] == 0.5


def test_read_level_composition_adds_the_background_to_the_viral_part(clades):
    """The shortcut the read half depends on, written out so it can be checked.

    Valid because a translated search scores each read from the read alone. Everything
    in the viral part is ground-truth viral and nothing in the background is, so the two
    halves *are* the detections and the false positives.
    """
    background = {
        "name": "background",
        "fragments_with_a_hit": 1000,
        "clades": {"coronaviridae": 3, "betacoronavirus": 1},
    }
    viral = {
        "name": "viral000010_viral",
        "fragments_with_a_hit": 8,
        "clades": {"coronaviridae": 7, "betacoronavirus": 5},
    }
    composed = day4._compose_read_level(
        background,
        viral,
        clades,
        {"post_qc_pairs": 10_000_000},
        {"viral_fragments": 10},
    )
    assert composed["fragments_with_a_hit"] == 1008
    assert composed["clades"]["coronaviridae"] == 10
    assert composed["detected"] == 7
    assert composed["false_positive"] == 3
    assert composed["sensitivity"] == 0.7
    assert composed["detected_per_million_post_qc"] == pytest.approx(0.7)


# -------------------------------------------------------------------------- scoring


def test_score_sample_alarms_only_on_a_clade_contig(
    tmp_path, taxonomy, clades, blast_row, write_blast_hits
):
    assembly = diamond.Assembly(
        name="s", contigs=tmp_path / "c.fa", contig_count=12, longest=900, n50=400
    )
    quiet = diamond.assign_contigs(
        write_blast_hits(tmp_path / "a.tsv.gz", [blast_row("c1", "YP_ecoli", 562, 300.0)]),
        {"c1": 900},
        taxonomy,
        clades,
    )
    loud = diamond.assign_contigs(
        write_blast_hits(tmp_path / "b.tsv.gz", [blast_row("c1", "YP_sars", 694009, 900.0)]),
        {"c1": 900},
        taxonomy,
        clades,
    )
    assert day4.HEADLINE_CLADE == "coronaviridae"
    assert (
        diamond.score_sample(assembly, quiet, clades, taxonomy, headline_clade="coronaviridae")[
            "alarm"
        ]
        is False
    )
    hot = diamond.score_sample(
        assembly,
        loud,
        clades,
        taxonomy,
        denominators={"post_qc_pairs": 1_000_000},
        headline_clade="coronaviridae",
    )
    assert hot["alarm"] is True
    assert hot["contigs_in_headline_clade"] == 1
    assert hot["contigs_per_million"]["post_qc_pairs"] == pytest.approx(1.0)


def test_false_alarm_summary_counts_samples_not_reads():
    """An assembly arm alarms per sample; the summary has to say so and mean it."""
    negatives = [
        {"accession": "A", "false_alarm_contigs": 0, "nonhost_pairs": 1000},
        {
            "accession": "B",
            "false_alarm_contigs": 3,
            "false_alarm_contigs_excluding_endemic": 0,
            "endemic_human_coronavirus_contigs": 3,
            "nonhost_pairs": 2000,
        },
        {
            "accession": "C",
            "false_alarm_contigs": 1,
            "false_alarm_contigs_excluding_endemic": 1,
            "nonhost_pairs": 3000,
        },
    ]
    summary = day4.summarise_false_alarms(negatives)
    assert summary["samples"] == 3
    assert summary["samples_with_at_least_one_clade_contig"] == 2
    assert summary["contigs_total"] == 4
    assert summary["non_host_pairs_total"] == 6000
    # B's three contigs are one patient's endemic coronavirus, which is the arm working.
    novel = summary["excluding_endemic_human_coronaviruses"]
    assert novel["samples_with_at_least_one_clade_contig"] == 1
    assert novel["contigs_total"] == 1
    assert summary["endemic_human_coronavirus_contigs_total"] == 3


def test_unknown_sensitivity_is_refused_before_anything_runs(tmp_path):
    with pytest.raises(diamond.DiamondError, match="unknown sensitivity"):
        diamond.blastx(
            tmp_path / "db.dmnd",
            tmp_path / "c.fa",
            tmp_path,
            name="x",
            sensitivity="maximum-effort",
        )


def test_megahit_read_count_comes_from_its_own_log(tmp_path):
    log = tmp_path / "m.log"
    log.write_text(
        "--- [STAT] 1000 reads, 150 max length\n"
        "--- something else ---\n"
        "--- [STAT] 21800000 reads, 151 max length\n"
    )
    # MEGAHIT counts single reads across both mates; the benchmark's unit is the pair.
    assert diamond._megahit_reads_in(log) == 10_900_000
    assert diamond._megahit_reads_in(tmp_path / "absent.log") == 0
