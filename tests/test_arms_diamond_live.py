"""Arm B against the real DIAMOND, on a database small enough to build in a second.

The rest of Arm B's tests feed the parser rows that were written by hand, which proves
the arm's arithmetic and proves nothing about whether DIAMOND emits what the arm thinks
it emits. The column order in ``BLAST_FIELDS`` is an on-disk format, the taxon map has
a layout DIAMOND accepts silently or ignores silently, and ``staxids`` is only populated
if ``makedb`` was given the taxonomy properly. Each of those fails quietly and would be
discovered as "the arm found nothing", hours into a run.

So this builds a two-protein database with a three-node taxonomy, searches a contig that
really does translate to one of those proteins, and checks the arm reads it correctly.
Skipped where diamond is not installed, which includes CI.
"""

from __future__ import annotations

import shutil

import pytest

from hindsight.arms import diamond
from hindsight.arms import taxa as taxa_mod

pytestmark = pytest.mark.skipif(shutil.which("diamond") is None, reason="diamond is not installed")

# Two invented proteins, distinct enough that a search cannot confuse them.
SPIKE = "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNG"
NUCLEOCAPSID = "MSDNGPQNQRNAPRITFGGPSDSTGSNQNGERSGARSKQRRPQGLPNNTASWFTALTQHGKEDLKFPRGQGVPIN"

# One codon per residue, so the contig translates to SPIKE exactly in frame 1.
CODONS = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTT",
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "CGT",
    "S": "TCT",
    "T": "ACT",
    "V": "GTT",
    "W": "TGG",
    "Y": "TAT",
}

NODES = """\
1\t|\t1\t|\tno rank\t|
10239\t|\t1\t|\tsuperkingdom\t|
11118\t|\t10239\t|\tfamily\t|
694009\t|\t11118\t|\tspecies\t|
864596\t|\t11118\t|\tspecies\t|
"""

NAMES = """\
1\t|\troot\t|\t\t|\tscientific name\t|
10239\t|\tViruses\t|\t\t|\tscientific name\t|
11118\t|\tCoronaviridae\t|\t\t|\tscientific name\t|
694009\t|\tSevere acute respiratory syndrome-related coronavirus\t|\t\t|\tscientific name\t|
864596\t|\tBat coronavirus BM48-31/BGR/2008\t|\t\t|\tscientific name\t|
"""


def _reverse_translate(protein: str) -> str:
    return "".join(CODONS[residue] for residue in protein)


@pytest.fixture
def source_db(tmp_path):
    """A `hindsight build` output directory, in miniature."""
    db = tmp_path / "db"
    (db / "taxonomy").mkdir(parents=True)
    (db / "viral.faa").write_text(
        f">YP_000001.1 spike protein\n{SPIKE}\n>YP_000002.1 nucleocapsid\n{NUCLEOCAPSID}\n"
    )
    (db / "seqid2taxid.map").write_text(
        "YP_000001.1\t694009\nYP_000002.1\t864596\nNC_000001.1\t694009\n"
    )
    (db / "taxonomy" / "nodes.dmp").write_text(NODES)
    (db / "taxonomy" / "names.dmp").write_text(NAMES)
    return db


def test_the_arm_reads_what_diamond_actually_writes(source_db, tmp_path):
    """End to end: makedb with the dated taxonomy, blastx, then the arm's own scoring."""
    dmnd = tmp_path / "work" / "mini.dmnd"
    record = diamond.build_database(source_db, dmnd, threads=2, log=lambda m: None)
    assert record["protein_records_in_source"] == 2
    assert record["taxonmap"]["records_with_a_taxid"] == 2

    contigs = tmp_path / "contigs.fa"
    contigs.write_text(f">k141_1\n{_reverse_translate(SPIKE)}\n")
    hits, meta = diamond.blastx(
        dmnd,
        contigs,
        tmp_path / "search",
        name="mini",
        threads=2,
        sensitivity="very-sensitive",
        log=lambda m: None,
    )
    assert meta["hit_rows"] >= 1

    parsed = list(diamond.iter_hits(hits))
    best = max(parsed, key=lambda h: h.bitscore)
    # The columns are in the order BLAST_FIELDS claims: a wrong order would put the
    # identity in `length` and this assertion would be nonsense rather than 100.
    assert best.contig == "k141_1"
    assert best.subject == "YP_000001.1"
    assert best.pident == pytest.approx(100.0)
    assert best.length == len(SPIKE)
    assert best.qlen == len(SPIKE) * 3
    # staxids is populated only if makedb really ingested the taxon map.
    assert best.taxids == [694009]

    taxonomy = taxa_mod.read_taxonomy(source_db / "taxonomy")
    clades = taxa_mod.clade_sets(taxonomy, {"viruses": 10239, "coronaviridae": 11118})
    calls = diamond.assign_contigs(hits, diamond.contig_lengths(contigs), taxonomy, clades)
    call = calls["k141_1"]
    assert call.best_taxid == 694009
    assert "coronaviridae" in call.lca_clades
    summary = diamond.summarise_calls(calls, clades, taxonomy)
    assert summary["clades"]["coronaviridae"]["contigs"] == 1
    assert summary["clades"]["coronaviridae"]["best_identity_percent"] == pytest.approx(100.0)


def test_a_contig_with_no_homology_produces_no_hits_and_no_alarm(source_db, tmp_path):
    """The negative control for the format test: real DIAMOND, genuinely nothing there."""
    dmnd = tmp_path / "work" / "mini.dmnd"
    diamond.build_database(source_db, dmnd, threads=2, log=lambda m: None)
    contigs = tmp_path / "contigs.fa"
    # A low-complexity repeat, which is what a spurious alarm would come from.
    contigs.write_text(">k141_junk\n" + "ATCGATCGTTAGCA" * 40 + "\n")
    hits, _ = diamond.blastx(
        dmnd, contigs, tmp_path / "search", name="junk", threads=2, log=lambda m: None
    )
    taxonomy = taxa_mod.read_taxonomy(source_db / "taxonomy")
    clades = taxa_mod.clade_sets(taxonomy, {"coronaviridae": 11118})
    calls = diamond.assign_contigs(hits, diamond.contig_lengths(contigs), taxonomy, clades)
    assert diamond.clade_contigs(calls, "coronaviridae") == []


def test_an_empty_assembly_is_a_result_not_a_crash(source_db, tmp_path):
    """The bottom of the dilution series: nothing assembled, so nothing to search."""
    dmnd = tmp_path / "work" / "mini.dmnd"
    diamond.build_database(source_db, dmnd, threads=2, log=lambda m: None)
    empty = tmp_path / "none.fa"
    empty.write_text("")
    hits, meta = diamond.blastx(
        dmnd, empty, tmp_path / "search", name="empty", threads=2, log=lambda m: None
    )
    assert meta["hit_rows"] == 0
    assert meta["skipped"] == "no contigs"
    assert list(diamond.iter_hits(hits)) == []
