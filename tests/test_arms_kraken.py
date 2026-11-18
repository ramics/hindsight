"""Arm A's promises, tested without a classifier.

Everything here is about the parts of the arm that could be wrong silently: clade
membership read from a *dated* tree, the composition shortcut that saves four full
classifications, and the point where the ground-truth label meets the classifier's
output. The classifier itself is Kraken2's business.
"""

from __future__ import annotations

import gzip
import json

import pytest

from hindsight.arms import kraken

# A miniature taxonomy: Viruses -> Riboviria -> ... -> Betacoronavirus -> two species,
# plus an unrelated bacterium so "everything is in the clade" cannot pass.
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
2\t|\t1\t|\tsuperkingdom\t|
562\t|\t2\t|\tspecies\t|
"""

NAMES = """\
10239\t|\tViruses\t|\t\t|\tscientific name\t|
694002\t|\tBetacoronavirus\t|\t\t|\tscientific name\t|
694009\t|\tSevere acute respiratory syndrome-related coronavirus\t|\t\t|\tscientific name\t|
562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|
694009\t|\tSARS-CoV\t|\t\t|\tsynonym\t|
"""


@pytest.fixture
def taxonomy_dir(tmp_path):
    directory = tmp_path / "taxonomy"
    directory.mkdir()
    (directory / "nodes.dmp").write_text(NODES)
    (directory / "names.dmp").write_text(NAMES)
    return directory


def test_lineage_and_descendants(taxonomy_dir):
    taxonomy = kraken.read_taxonomy(taxonomy_dir)
    assert taxonomy.lineage(694009) == [694009, 2509511, 694002, 2501931, 11118, 76804, 2559587, 10239, 1]
    assert taxonomy.descendants(694002) == {694002, 2509511, 694009, 31631}
    assert 562 not in taxonomy.descendants(10239)


def test_scientific_names_only(taxonomy_dir):
    taxonomy = kraken.read_taxonomy(taxonomy_dir)
    assert taxonomy.names[694009].startswith("Severe acute respiratory")


def test_lineage_of_an_unknown_taxid_is_empty(taxonomy_dir):
    """A taxid the dated tree never held is not an error and is not a root, either."""
    taxonomy = kraken.read_taxonomy(taxonomy_dir)
    assert taxonomy.lineage(2697049) == []
    assert taxonomy.descendants(2697049) == set()


def test_clade_sets_report_what_the_dated_tree_lacks(taxonomy_dir):
    taxonomy = kraken.read_taxonomy(taxonomy_dir)
    clades = kraken.clade_sets(
        taxonomy, {"betacoronavirus": 694002, "not_yet_named": 2697049}
    )
    assert clades.members["betacoronavirus"] == {694002, 2509511, 694009, 31631}
    assert clades.missing == ["not_yet_named (taxid 2697049)"]
    assert "not_yet_named" not in clades.members


def test_clade_counting_ignores_taxa_outside_the_clade(taxonomy_dir):
    taxonomy = kraken.read_taxonomy(taxonomy_dir)
    clades = kraken.clade_sets(taxonomy, {"betacoronavirus": 694002, "viruses": 10239})
    counts = {694009: 7, 31631: 3, 562: 1000}
    assert clades.count(counts, "betacoronavirus") == 10
    assert clades.count(counts, "viruses") == 10


def test_which_clades_a_taxid_belongs_to(taxonomy_dir):
    taxonomy = kraken.read_taxonomy(taxonomy_dir)
    clades = kraken.clade_sets(taxonomy, {"betacoronavirus": 694002, "viruses": 10239})
    assert sorted(clades.which(694009)) == ["betacoronavirus", "viruses"]
    assert clades.which(562) == []


# ------------------------------------------------------------------- output handling


def write_output(path, rows):
    with gzip.open(path, "wt") as handle:
        for state, read_id, taxid in rows:
            handle.write(f"{state}\t{read_id}\t{taxid}\t151|151\tA:1 B:2\n")


def test_counts_from_output(tmp_path):
    path = tmp_path / "out.tsv.gz"
    write_output(path, [("C", "f1", 694009), ("C", "f2", 694009), ("C", "f3", 562)])
    assert kraken.counts_from_output(path) == {694009: 2, 562: 1}


def test_scoring_splits_detections_by_the_truth_label(tmp_path, taxonomy_dir):
    taxonomy = kraken.read_taxonomy(taxonomy_dir)
    clades = kraken.clade_sets(taxonomy, {"betacoronavirus": 694002})
    path = tmp_path / "out.tsv.gz"
    write_output(
        path,
        [
            ("C", "viral1", 694009),  # true positive
            ("C", "viral2", 562),  # viral read classified, but not as a coronavirus
            ("C", "bg1", 31631),  # false positive: in the clade, not really viral
            ("C", "bg2", 562),  # neither
        ],
    )
    scored = kraken.score_against_truth(
        path, {"viral1", "viral2", "viral3"}, clades, clade="betacoronavirus"
    )
    assert scored["detected"] == 1
    assert scored["false_positive"] == 1
    assert scored["truth_fragments"] == 3
    assert scored["truth_fragments_classified_anywhere"] == 2
    # Three fragments were labelled viral; one was placed in the clade.
    assert scored["sensitivity"] == pytest.approx(1 / 3)


def test_scoring_without_a_truth_set_has_no_sensitivity(tmp_path, taxonomy_dir):
    taxonomy = kraken.read_taxonomy(taxonomy_dir)
    clades = kraken.clade_sets(taxonomy, {"betacoronavirus": 694002})
    path = tmp_path / "out.tsv.gz"
    write_output(path, [("C", "bg1", 694009)])
    scored = kraken.score_against_truth(path, set(), clades, clade="betacoronavirus")
    assert scored["sensitivity"] is None
    assert scored["false_positive"] == 1


def test_read_truth_ids_rejects_a_missing_file(tmp_path):
    with pytest.raises(kraken.KrakenError):
        kraken.read_truth_ids(tmp_path / "nope.txt")


# ---------------------------------------------------------------------- composition


def classification(name, counts, pairs, classified=None):
    result = kraken.Classification(
        name=name,
        report=f"{name}.kreport",
        classified_output=f"{name}.tsv.gz",
        pairs=pairs,
        classified=classified if classified is not None else sum(counts.values()),
    )
    result.counts = dict(counts)
    return result


def test_compose_adds_disjoint_classifications():
    background = classification("background", {562: 40, 694009: 1}, pairs=1000)
    viral = classification("viral", {694009: 9}, pairs=10)
    composed = kraken.compose([background, viral])
    assert composed["pairs"] == 1010
    assert composed["counts"] == {562: 40, 694009: 10}
    assert composed["classified"] == 50
    assert composed["unclassified"] == 960


def test_verify_composition_accepts_an_exact_match():
    composed = kraken.compose(
        [
            classification("background", {562: 40, 694009: 1}, pairs=1000),
            classification("viral", {694009: 9}, pairs=10),
        ]
    )
    direct = classification("level_direct", {562: 40, 694009: 10}, pairs=1010)
    check = kraken.verify_composition(composed, direct)
    assert check["identical"]
    assert check["taxa_disagreeing"] == 0


def test_verify_composition_names_the_taxa_that_disagree():
    composed = kraken.compose([classification("background", {562: 40}, pairs=1000)])
    direct = classification("level_direct", {562: 39, 694009: 1}, pairs=1000)
    check = kraken.verify_composition(composed, direct)
    assert not check["identical"]
    assert check["taxa_disagreeing"] == 2
    assert check["differences"][562] == [40, 39]
    assert check["differences"][694009] == [0, 1]


def test_verify_composition_notices_a_pair_count_mismatch():
    """Equal taxon counts are not enough: a level has to be the whole level."""
    composed = kraken.compose([classification("background", {562: 1}, pairs=999)])
    direct = classification("level_direct", {562: 1}, pairs=1000)
    check = kraken.verify_composition(composed, direct)
    assert not check["identical"]
    assert check["taxa_disagreeing"] == 0


# ------------------------------------------------------------------------- the record


def test_classification_record_reports_every_clade(tmp_path, taxonomy_dir):
    taxonomy = kraken.read_taxonomy(taxonomy_dir)
    clades = kraken.clade_sets(taxonomy, {"viruses": 10239, "betacoronavirus": 694002})
    result = classification("x", {694009: 5, 562: 2}, pairs=100)
    record = result.record(clades)
    assert record["clades"] == {"viruses": 5, "betacoronavirus": 5}
    assert record["top_taxa"][0] == {
        "taxid": 694009,
        "reads": 5,
        "taxon": "Severe acute respiratory syndrome-related coronavirus",
    }


def test_headline_clade_is_one_of_the_recorded_clades():
    assert kraken.HEADLINE_CLADE in kraken.CLADES


def test_no_post_2019_taxon_is_in_the_clade_list():
    """The arm cannot be scored on a taxon that did not exist at the cutoff.

    SARS-CoV-2 (2697049) and its species taxon were created in 2020. Scoring the arm
    against them would be scoring it on a label the November 2019 tree cannot carry,
    and the only honest question is whether the reads landed in the betacoronaviruses.
    """
    assert 2697049 not in kraken.CLADES.values()


def test_scoring_reports_every_clade_not_only_the_headline(tmp_path, taxonomy_dir):
    """The clade that counts as a detection is not a free choice in a dated tree.

    31631 stands in for Bat coronavirus BM48-31/BGR/2008, which the November 2019
    taxonomy files under "unclassified coronaviruses" in Coronaviridae rather than
    under Betacoronavirus, where today's tree puts it. A read matching the nearest
    relative the database holds is therefore a detection at family level and nothing
    at genus level, and the scoring has to show both or Day 5 cannot see the
    difference.
    """
    nodes = (taxonomy_dir / "nodes.dmp").read_text() + "99999\t|\t11118\t|\tspecies\t|\n"
    (taxonomy_dir / "nodes.dmp").write_text(nodes)
    taxonomy = kraken.read_taxonomy(taxonomy_dir)
    clades = kraken.clade_sets(
        taxonomy, {"coronaviridae": 11118, "betacoronavirus": 694002}
    )
    path = tmp_path / "out.tsv.gz"
    write_output(
        path,
        [
            ("C", "viral1", 694009),  # detected at both levels
            ("C", "viral2", 99999),  # detected at family level only
            ("C", "bg1", 99999),  # a family-level false positive
        ],
    )
    scored = kraken.score_against_truth(
        path, {"viral1", "viral2"}, clades, clade="betacoronavirus"
    )
    assert scored["detected"] == 1
    assert scored["by_clade"]["betacoronavirus"]["detected"] == 1
    assert scored["by_clade"]["coronaviridae"]["detected"] == 2
    assert scored["by_clade"]["betacoronavirus"]["false_positive"] == 0
    assert scored["by_clade"]["coronaviridae"]["false_positive"] == 1
    assert scored["by_clade"]["coronaviridae"]["sensitivity"] == 1.0
