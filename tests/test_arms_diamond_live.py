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

import dataclasses
import shutil
from pathlib import Path

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


# --------------------------------------------------------------- contig diagnosis


def _hits_path(out_dir: Path, stage: str, name: str, arm) -> Path:
    """Where ``day4.diagnose_contigs`` looks for a stage's DIAMOND output.

    Derived from the arm's own config rather than spelled out, so that moving the
    sensitivity setting — which Day 5 did, from ``--very-sensitive`` to stock — cannot
    leave these fixtures writing to a filename nothing reads. These tests are skipped
    wherever minimap2 or diamond is absent, which includes CI, so a stale name here fails
    only in a container with the real tools and fails silently everywhere else.
    """
    path = out_dir / stage / "search" / f"{name}.{arm.sensitivity}.blastx.tsv.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 is not installed")
def test_the_contig_diagnosis_separates_the_virus_from_the_false_positives(
    source_db, tmp_path, monkeypatch, blast_row, write_blast_hits
):
    """The diagnostic's own arithmetic, against real minimap2 on a known reference.

    Two contigs are handed in: one that is a 2 kb slice of the "reference" with a handful
    of substitutions, and one unrelated sequence. The first must come back as a
    reconstruction with a plausible identity and coverage, and the second must come back
    as a contig that does not align at all — which is the arm's false positive. Getting
    this backwards would turn a false positive into a headline.
    """
    import random

    from hindsight.arms import day4
    from hindsight.samples import registry

    rng = random.Random(20191201)
    genome = "".join(rng.choice("ACGT") for _ in range(6000))
    work = tmp_path / "work"
    (work / "ref" / "truth").mkdir(parents=True)
    (work / "ref" / "truth" / "NC_TEST.1.fna").write_text(f">NC_TEST.1 reference\n{genome}\n")

    # A 2 kb stretch of the reference with 10 substitutions: ~99.5% identity, 33% of it.
    piece = list(genome[1000:3000])
    for position in rng.sample(range(len(piece)), 10):
        piece[position] = rng.choice([b for b in "ACGT" if b != piece[position]])
    unrelated = "".join(rng.choice("ACGT") for _ in range(1500))

    out = tmp_path / "arm-b"
    (out / "index" / "assembly").mkdir(parents=True)
    (out / "index" / "assembly" / "SRRTEST_nonhost.contigs.fa").write_text(
        f">k_real\n{''.join(piece)}\n>k_bogus\n{unrelated}\n"
    )
    dmnd = tmp_path / "db.dmnd"
    diamond.build_database(source_db, dmnd, threads=2, log=lambda m: None)
    # Give both contigs a coronavirus call, so the diagnosis has two candidates and has
    # to tell them apart by alignment rather than by what the search said. The fixture's
    # filename carries the sensitivity setting, so it is derived from the config the
    # diagnosis itself reads: Day 5 moved that setting from --very-sensitive to stock, and
    # a name spelled out here went looking for a file the arm no longer writes.
    arm = day4.load_config(Path("config/day4.toml"))
    hits = _hits_path(out, "index", "SRRTEST_nonhost", arm)
    write_blast_hits(
        hits,
        [
            blast_row("k_real", "YP_000001.1", 694009, 900.0, pident=79.9, qlen=2000),
            blast_row("k_bogus", "YP_000002.1", 864596, 120.0, pident=28.4, qlen=1500),
        ],
    )

    arm = dataclasses.replace(arm, work_dir=work, out_dir=out, source_db=source_db, dmnd_db=dmnd)
    samples = dataclasses.replace(
        registry.load_config(Path("config/day2.toml")),
        index_run="SRRTEST",
        truth_reference="NC_TEST.1",
    )
    taxonomy = taxa_mod.read_taxonomy(source_db / "taxonomy")
    clades = taxa_mod.clade_sets(taxonomy, {"coronaviridae": 11118})

    record = day4.diagnose_contigs(
        arm, samples, taxonomy, clades, clade="coronaviridae", log=lambda m: None
    )
    assert record["contigs_examined"] == 2
    assert record["contigs_aligning"] == 1
    assert record["contigs_not_aligning"] == 1
    best = record["best_reconstruction"]
    assert best["contig"] == "k_real"
    assert 99.0 <= best["nucleotide_identity_percent"] <= 100.0
    assert 30.0 <= best["reference_covered_percent"] <= 35.0
    assert best["best_hit_identity_percent"] == 79.9
    by_contig = {row["contig"]: row for row in record["contigs"]}
    assert by_contig["k_real"]["aligns_to_the_virus"] is True
    assert by_contig["k_bogus"]["aligns_to_the_virus"] is False
    assert "genuine false positives" in record["conclusion"]


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 is not installed")
def test_the_diagnosis_says_so_when_nothing_aligns(
    source_db, tmp_path, blast_row, write_blast_hits
):
    """Every clade contig a false positive is a result the conclusion must state."""
    import random

    from hindsight.arms import day4
    from hindsight.samples import registry

    rng = random.Random(7)
    work = tmp_path / "work"
    (work / "ref" / "truth").mkdir(parents=True)
    (work / "ref" / "truth" / "NC_TEST.1.fna").write_text(
        ">NC_TEST.1 reference\n" + "".join(rng.choice("ACGT") for _ in range(6000)) + "\n"
    )
    out = tmp_path / "arm-b"
    (out / "index" / "assembly").mkdir(parents=True)
    (out / "index" / "assembly" / "SRRTEST_nonhost.contigs.fa").write_text(
        ">k_bogus\n" + "".join(rng.choice("ACGT") for _ in range(1500)) + "\n"
    )
    dmnd = tmp_path / "db.dmnd"
    diamond.build_database(source_db, dmnd, threads=2, log=lambda m: None)
    stock = day4.load_config(Path("config/day4.toml"))
    hits = _hits_path(out, "index", "SRRTEST_nonhost", stock)
    write_blast_hits(hits, [blast_row("k_bogus", "YP_000001.1", 694009, 120.0, qlen=1500)])

    arm = dataclasses.replace(
        stock,
        work_dir=work,
        out_dir=out,
        source_db=source_db,
        dmnd_db=dmnd,
    )
    samples = dataclasses.replace(
        registry.load_config(Path("config/day2.toml")),
        index_run="SRRTEST",
        truth_reference="NC_TEST.1",
    )
    taxonomy = taxa_mod.read_taxonomy(source_db / "taxonomy")
    clades = taxa_mod.clade_sets(taxonomy, {"coronaviridae": 11118})
    record = day4.diagnose_contigs(
        arm, samples, taxonomy, clades, clade="coronaviridae", log=lambda m: None
    )
    assert record["contigs_aligning"] == 0
    assert "best_reconstruction" not in record
    assert "every one is a false positive" in record["conclusion"]
