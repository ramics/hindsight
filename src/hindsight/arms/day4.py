"""Day 4: run Arm B over the whole benchmark and write down what it found.

The arm itself is in :mod:`hindsight.arms.diamond`. This is the driver — what gets
assembled and searched, in what order, and what is recorded about each run.

Five decisions are encoded here rather than in a config file, because each one is a
claim about the benchmark and belongs where it can be read.

**The arm has two halves, and they are not redundant.** The *contig* half is Wu et al.'s
method: assemble, search the contigs, look at what the longest one hits. The *read* half
searches the individual reads against the same protein database with no assembly step.
Running both is the only way to answer the question Day 5 will actually ask about this
arm — when it stops firing, did it stop because the protein search ran out of homology,
or because the assembler ran out of coverage? One arm cannot distinguish those, and they
have completely different implications for how deep a surveillance programme must
sequence.

**The read half composes; the contig half cannot.** A translated search scores each read
against the database independently, exactly as Kraken2 does, so a dilution level's read
half is the background's hits plus that level's viral reads' hits — one expensive
background run plus five trivial ones, and
:func:`hindsight.arms.diamond.iter_hits` makes the arithmetic auditable. Assembly is the
opposite: a contig is a statement about the whole set of reads, and taking 49,990 viral
reads out changes which contigs exist. So every level is materialised and assembled
whole, which is most of this day's runtime and cannot be optimised away.

**The background is assembled and searched on its own.** It is the index patient's own
non-host reads with every ground-truth viral fragment removed, so any coronavirus contig
or read hit in it is a false positive with no labelling step needed. That is this arm's
false-positive baseline on the positive sample, the counterpart to Arm A's.

**The negatives are processed at native depth, one at a time.** Day 3 established that
this fits if each sample is fetched, filtered, searched and then deleted, and that
staging rather than compute is the cost. Day 4 pays the staging cost again because Arm A
ran in a different container and its reads are gone; see the day's write-up for what a
combined runner would have saved.
"""

from __future__ import annotations

import dataclasses
import json
import time
import tomllib
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ..samples import pipeline as samples_pipeline
from ..samples import reads, registry
from . import diamond
from .taxa import (
    clade_sets,
    machine_record,
    per_million,
    read_taxonomy,
    tool_versions,
    utc_now,
)

DEFAULT_TOOL_PINS = {
    "diamond": "2.1.9",
    "megahit": "1.2.9",
    "bowtie2": "2.5.2",
}

# Arm B's headline is drawn at the family rather than the genus, on Day 3's finding:
# the November 2019 taxonomy files *Bat coronavirus BM48-31/BGR/2008* under
# "unclassified coronaviruses" in Coronaviridae rather than in Betacoronavirus, and that
# virus is the nearest relative the dated database holds. For a homology arm this
# matters more than it did for Arm A, not less — the LCA of a sarbecovirus contig's
# near-tied hits lands on whatever node is common to SARS-CoV and the bat viruses, which
# in the 2019 tree is frequently the family. Every clade is recorded per level either
# way, so Day 5 redraws this line without a re-run.
HEADLINE_CLADE = "coronaviridae"


@dataclass(frozen=True)
class ArmConfig:
    samples_config: Path
    source_db: Path
    dmnd_db: Path
    cutoff: str
    threads: int
    sensitivity: str
    comparison_sensitivity: str | None
    evalue: float
    max_target_seqs: int
    lca_top_percent: float
    min_contig_len: int
    megahit_memory_fraction: float
    read_level: bool
    read_level_sensitivity: str
    negatives_max_pairs: int
    negatives_read_level: bool
    prune_after_search: bool
    work_dir: Path
    out_dir: Path
    provenance: Path
    tool_pins: dict[str, str | None]


def load_config(path: Path) -> ArmConfig:
    data = tomllib.loads(Path(path).read_text())
    arm = data.get("arm", {})
    database = data.get("database", {})
    search = data.get("search", {})
    assembly = data.get("assembly", {})
    read_level = data.get("read_level", {})
    negatives = data.get("negatives", {})
    work_dir = Path(arm.get("work_dir", "work"))
    return ArmConfig(
        samples_config=Path(arm.get("samples_config", "config/day2.toml")),
        source_db=Path(database.get("source_db", "db/nov2019")),
        dmnd_db=Path(database.get("dmnd_db", work_dir / "diamond" / "nov2019.dmnd")),
        cutoff=str(database.get("cutoff", "2019-12-01")),
        threads=int(search.get("threads", 4)),
        sensitivity=str(search.get("sensitivity", "very-sensitive")),
        comparison_sensitivity=search.get("comparison_sensitivity", "default") or None,
        evalue=float(search.get("evalue", 1e-3)),
        max_target_seqs=int(search.get("max_target_seqs", 25)),
        lca_top_percent=float(search.get("lca_top_percent", diamond.LCA_TOP_PERCENT)),
        min_contig_len=int(assembly.get("min_contig_len", 200)),
        megahit_memory_fraction=float(assembly.get("memory_fraction", 0.8)),
        read_level=bool(read_level.get("enabled", True)),
        read_level_sensitivity=str(read_level.get("sensitivity", "very-sensitive")),
        negatives_max_pairs=int(negatives.get("max_pairs", 0)),
        negatives_read_level=bool(negatives.get("read_level", False)),
        prune_after_search=bool(negatives.get("prune_after_search", True)),
        work_dir=work_dir,
        out_dir=Path(arm.get("out_dir", work_dir / "arm-b")),
        provenance=Path(arm.get("provenance", "results/day4-provenance.json")),
        tool_pins=dict(data.get("tools", DEFAULT_TOOL_PINS)),
    )


# ------------------------------------------------------------------------- database


def prepare_database(arm: ArmConfig, *, log=print):
    record = diamond.build_database(arm.source_db, arm.dmnd_db, threads=arm.threads, log=log)
    taxonomy = read_taxonomy(arm.source_db / "taxonomy")
    clades = clade_sets(taxonomy)
    record["clades"] = {
        name: {
            "root_taxid": clades.roots[name],
            "name_in_2019_tree": taxonomy.names.get(clades.roots[name]),
            "taxa_in_clade": len(members),
        }
        for name, members in clades.members.items()
    }
    if clades.missing:
        log("clades absent from the dated taxonomy: " + ", ".join(clades.missing))
    record["clades_absent_from_dated_taxonomy"] = clades.missing
    record["headline_clade"] = HEADLINE_CLADE
    record["protein_counts_by_clade"] = _proteins_per_clade(arm.source_db, clades)
    source_provenance = arm.source_db / "provenance.json"
    if source_provenance.exists():
        source = json.loads(source_provenance.read_text())
        record["source_release"] = {
            "time_control": source.get("time_control", {}),
            "counts": source.get("counts"),
            "taxonomy": source.get("taxonomy"),
        }
    return record, taxonomy, clades


def _proteins_per_clade(source_db: Path, clades) -> dict:
    """How many proteins the dated database actually holds for each clade.

    Arm A records minimizer counts for the same reason. "SARS-related coronavirus is in
    the database" is a weaker statement than it sounds if the entry is a handful of
    proteins, and this arm's whole sensitivity rests on how much conserved protein there
    is to match against.
    """
    seqid_map = Path(source_db) / "seqid2taxid.map"
    faa = Path(source_db) / "viral.faa"
    if not seqid_map.exists() or not faa.exists():
        return {}
    protein_ids: set[str] = set()
    with faa.open() as handle:
        for line in handle:
            if line.startswith(">"):
                protein_ids.add(line[1:].split()[0])
    taxid_of: dict[str, int] = {}
    with seqid_map.open() as handle:
        for line in handle:
            accession, _, taxid = line.rstrip("\n").partition("\t")
            if accession in protein_ids:
                try:
                    taxid_of[accession] = int(taxid)
                except ValueError:
                    continue
    counts: dict[str, int] = {}
    for name, members in clades.members.items():
        counts[name] = sum(1 for taxid in taxid_of.values() if taxid in members)
    return counts


# -------------------------------------------------------------------- one whole sample


def _contig_arm(
    arm: ArmConfig,
    r1: Path,
    r2: Path,
    out_dir: Path,
    taxonomy,
    clades,
    *,
    name: str,
    truth_ids: set[str] | None = None,
    denominators: dict[str, int] | None = None,
    sensitivity: str | None = None,
    map_reads_back: bool = True,
    log=print,
) -> dict:
    """Assemble one read pair, search the contigs, and score what came back."""
    sensitivity = sensitivity or arm.sensitivity
    assembly = diamond.assemble(
        r1,
        r2,
        out_dir / "assembly",
        name=name,
        threads=arm.threads,
        memory_fraction=arm.megahit_memory_fraction,
        min_contig_len=arm.min_contig_len,
        log=log,
    )
    hits, search_meta = diamond.blastx(
        arm.dmnd_db,
        assembly.contigs,
        out_dir / "search",
        name=name,
        threads=arm.threads,
        sensitivity=sensitivity,
        evalue=arm.evalue,
        max_target_seqs=arm.max_target_seqs,
        log=log,
    )
    lengths = diamond.contig_lengths(assembly.contigs)
    calls = diamond.assign_contigs(
        hits, lengths, taxonomy, clades, lca_top_percent=arm.lca_top_percent
    )
    record = diamond.score_sample(
        assembly,
        calls,
        clades,
        taxonomy,
        denominators=denominators,
        headline_clade=HEADLINE_CLADE,
    )
    record["search"]["meta"] = search_meta
    record["sensitivity"] = sensitivity
    clade_calls = diamond.clade_contigs(calls, HEADLINE_CLADE, by="lca")
    if map_reads_back and clade_calls:
        record["reads_on_clade_contigs"] = diamond.count_reads_on_contigs(
            assembly.contigs,
            [c.contig for c in clade_calls],
            r1,
            r2,
            out_dir / "readback",
            name=name,
            threads=arm.threads,
            truth_ids=truth_ids,
            log=log,
        )
        if denominators:
            record["reads_on_clade_contigs"]["per_million"] = {
                key: per_million(record["reads_on_clade_contigs"]["reads_aligned"], value)
                for key, value in denominators.items()
            }
    log(
        f"{name}: {assembly.contig_count:,} contigs, "
        f"{record['contigs_in_headline_clade']} in {HEADLINE_CLADE}"
        + (
            f", longest {record['search']['clades'][HEADLINE_CLADE]['longest_contig']:,} bp "
            f"at {record['search']['clades'][HEADLINE_CLADE]['best_identity_percent']}% aa identity"
            if record["contigs_in_headline_clade"]
            else ""
        )
    )
    return record


def _read_arm(
    arm: ArmConfig,
    r1: Path,
    r2: Path,
    out_dir: Path,
    taxonomy,
    clades,
    *,
    name: str,
    truth_ids: set[str] | None = None,
    log=print,
) -> dict:
    """Search the reads themselves, with no assembly step.

    The reads are converted to FASTA because DIAMOND does not need the qualities and a
    FASTA of the same reads is a third of the bytes. Both mates go in as separate
    queries and are collapsed back to fragments afterwards, so the unit matches Arm A's:
    a fragment counts once however many of its mates hit.
    """
    out_dir = Path(out_dir)
    fasta = out_dir / "reads" / f"{name}.reads.fa.gz"
    search_dir = out_dir / "search"
    # Building the FASTA means reading 22 million records and writing 3 GiB, so a search
    # that is already on disk must not pay for it. `blastx` resumes from its own marker,
    # but it cannot do that until it has been handed a query file it no longer needs.
    done = search_dir / f".{name}_reads.{arm.read_level_sensitivity}.searched.json"
    counted = search_dir / f"{name}_reads.fragments.json"
    if done.exists() and counted.exists():
        fragments = json.loads(counted.read_text())["fragments"]
    else:
        fragments = _reads_to_fasta(r1, r2, fasta)
        counted.parent.mkdir(parents=True, exist_ok=True)
        counted.write_text(json.dumps({"fragments": fragments}, indent=2) + "\n")
    hits, meta = diamond.blastx(
        arm.dmnd_db,
        fasta,
        search_dir,
        name=f"{name}_reads",
        threads=arm.threads,
        sensitivity=arm.read_level_sensitivity,
        evalue=arm.evalue,
        # A read is 50 aa of query; it does not tie across 25 subjects the way a 30 kb
        # contig does, and 10 keeps the output a manageable size over 22 million reads.
        max_target_seqs=10,
        log=log,
    )
    record = _score_read_hits(hits, taxonomy, clades, truth_ids=truth_ids)
    record["name"] = name
    record["fragments_searched"] = fragments
    record["meta"] = meta
    # 1-3 GiB of FASTA per sample whose only consumer was the search that just finished.
    fasta.unlink(missing_ok=True)
    log(
        f"{name} (reads): {record['fragments_with_a_hit']:,} fragments hit the database, "
        f"{record['clades'].get(HEADLINE_CLADE, 0):,} in {HEADLINE_CLADE}"
    )
    return record


def _reads_to_fasta(r1: Path, r2: Path, dest: Path) -> int:
    """Write both mates of a pair into one FASTA, returning the fragment count."""
    import gzip

    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    seen: set[str] = set()
    with gzip.open(dest, "wt") as out:
        for mate, source in ((1, r1), (2, r2)):
            for header, sequence, _plus, _qual in reads.iter_records(source):
                fragment = reads.normalise_id(header)
                seen.add(fragment)
                out.write(f">{fragment}/{mate}\n{sequence}\n")
    return len(seen)


def _score_read_hits(
    hits_path: Path,
    taxonomy,
    clades,
    *,
    truth_ids: set[str] | None = None,
) -> dict:
    """Collapse per-mate protein hits into per-fragment clade calls.

    Each fragment's hits are pooled across its mates and an LCA is taken over the
    near-tied ones, the same rule the contig half uses. Pooling the mates is the right
    unit: Arm A classifies a pair as one fragment, and a benchmark table that counted
    one arm's mates and the other arm's fragments would be comparing two different
    things.
    """
    from collections import Counter

    best: dict[str, float] = {}
    rows: dict[str, list[tuple[float, list[int]]]] = {}
    for hit in diamond.iter_hits(hits_path):
        fragment = hit.contig.rsplit("/", 1)[0]
        rows.setdefault(fragment, []).append((hit.bitscore, hit.taxids))
        if hit.bitscore > best.get(fragment, 0.0):
            best[fragment] = hit.bitscore

    assignments: dict[str, int] = {}
    for fragment, entries in rows.items():
        cutoff = best[fragment] * (1.0 - diamond.LCA_TOP_PERCENT / 100.0)
        taxids = [taxid for score, ids in entries if score >= cutoff for taxid in ids]
        lca = taxonomy.lowest_common_ancestor(taxids)
        if lca is not None:
            assignments[fragment] = lca

    counts: Counter[int] = Counter(assignments.values())
    record: dict = {
        "fragments_with_a_hit": len(rows),
        "fragments_assigned": len(assignments),
        "clades": {
            name: sum(1 for taxid in assignments.values() if taxid in members)
            for name, members in clades.members.items()
        },
        "top_lca_taxa": [
            {"taxid": taxid, "fragments": n, "taxon": taxonomy.names.get(taxid)}
            for taxid, n in counts.most_common(15)
        ],
    }
    if truth_ids is not None:
        detected = {}
        false_positive = {}
        for name, members in clades.members.items():
            in_clade = {f for f, taxid in assignments.items() if taxid in members}
            detected[name] = len(in_clade & truth_ids)
            false_positive[name] = len(in_clade - truth_ids)
        record["truth"] = {
            "truth_fragments": len(truth_ids),
            "detected_by_clade": detected,
            "false_positive_by_clade": false_positive,
            "detected": detected.get(HEADLINE_CLADE, 0),
            "false_positive": false_positive.get(HEADLINE_CLADE, 0),
            "sensitivity": (
                detected.get(HEADLINE_CLADE, 0) / len(truth_ids) if truth_ids else None
            ),
        }
    return record


# ---------------------------------------------------------------------- index sample


def run_index(
    arm: ArmConfig,
    samples_config: registry.Config,
    day2: dict,
    taxonomy,
    clades,
    *,
    log=print,
) -> dict:
    accession = samples_config.index_run
    sample_dir = arm.work_dir / "samples" / accession
    r1 = sample_dir / "nonhost" / f"{accession}_nonhost_1.fastq.gz"
    r2 = sample_dir / "nonhost" / f"{accession}_nonhost_2.fastq.gz"
    for path in (r1, r2):
        if not path.exists():
            raise diamond.DiamondError(f"{path} is missing; re-run the Day 2 pipeline")
    truth_ids = _truth_ids(sample_dir)
    denominators = (day2.get("dilution_series") or {}).get("denominators", {})

    record = _contig_arm(
        arm,
        r1,
        r2,
        arm.out_dir / "index",
        taxonomy,
        clades,
        name=f"{accession}_nonhost",
        truth_ids=truth_ids,
        denominators=denominators,
        log=log,
    )
    record["accession"] = accession
    record["role"] = "index"
    record["truth_fragments_in_nonhost"] = len(truth_ids)

    # The same sample at DIAMOND's stock sensitivity. The arm runs at --very-sensitive
    # and says why; this is the evidence for what that choice is worth, measured rather
    # than asserted.
    if arm.comparison_sensitivity and arm.comparison_sensitivity != arm.sensitivity:
        log(f"{accession}: comparison run at --{arm.comparison_sensitivity}")
        record["sensitivity_comparison"] = _contig_arm(
            arm,
            r1,
            r2,
            arm.out_dir / "index",
            taxonomy,
            clades,
            name=f"{accession}_nonhost",
            truth_ids=truth_ids,
            denominators=denominators,
            sensitivity=arm.comparison_sensitivity,
            map_reads_back=False,
            log=log,
        )
    if arm.read_level:
        record["read_level"] = _index_read_level(
            arm, samples_config, day2, r1, r2, truth_ids, taxonomy, clades, log=log
        )
    return record


def _index_read_level(
    arm: ArmConfig,
    samples_config: registry.Config,
    day2: dict,
    r1: Path,
    r2: Path,
    truth_ids: set[str],
    taxonomy,
    clades,
    *,
    log=print,
) -> dict:
    """The index sample's read half, composed rather than searched whole.

    Searching the 10.9 million non-host pairs costs hours, and the series stage has to
    search the background — the same reads with the 62,117 viral fragments taken out —
    anyway. Day 2 built that background by exclusion, so the two sets are an exact
    partition of the non-host reads, and the index sample's read half is the background's
    hits plus the viral fragments' hits. The viral part is 61,219 fragments, which is
    seconds.

    The saving is one full read-level search, and the arithmetic is the same one Arm A
    used for its series. The composition is not assumed: ``read_level_partition`` in the
    record states both halves' fragment counts and the non-host total they must sum to,
    so a reader can check it.
    """
    series_out = arm.out_dir / "series"
    parts = arm.work_dir / "series" / "parts"
    bg1, bg2 = parts / "background_1.fastq.gz", parts / "background_2.fastq.gz"
    for path in (bg1, bg2):
        if not path.exists():
            raise diamond.DiamondError(f"{path} is missing; re-run the Day 2 pipeline")
    background = _read_arm(arm, bg1, bg2, series_out, taxonomy, clades, name="background", log=log)

    # Every truth fragment that survived host filtering, pulled out of the non-host set
    # rather than taken from the label: 898 of the 62,117 align to the human genome and
    # are not in these files, and the arm must score what it was actually given.
    viral_dir = arm.out_dir / "index" / "viral-fragments"
    viral_dir.mkdir(parents=True, exist_ok=True)
    v1 = viral_dir / "index_viral_1.fastq.gz"
    v2 = viral_dir / "index_viral_2.fastq.gz"
    if not (v1.exists() and v2.exists()):
        for source, dest in ((r1, v1), (r2, v2)):
            reads.write_subset(source, dest, truth_ids)
    viral = _read_arm(
        arm,
        v1,
        v2,
        arm.out_dir / "index",
        taxonomy,
        clades,
        name="index_viral",
        log=log,
    )
    for stale in (v1, v2):
        stale.unlink(missing_ok=True)

    sample = next(
        (s for s in day2.get("samples", []) if s["accession"] == samples_config.index_run),
        {},
    )
    nonhost_pairs = sample.get("nonhost_pairs")
    composed = _compose_read_level(
        background,
        viral,
        clades,
        (day2.get("dilution_series") or {}).get("denominators", {}),
        {"viral_fragments": viral["fragments_searched"]},
    )
    composed["read_level_partition"] = {
        "background_fragments": background["fragments_searched"],
        "viral_fragments": viral["fragments_searched"],
        "sum": (background["fragments_searched"] or 0) + (viral["fragments_searched"] or 0),
        "nonhost_pairs_recorded_by_day2": nonhost_pairs,
        "partition_is_exact": (
            (background["fragments_searched"] or 0) + (viral["fragments_searched"] or 0)
            == nonhost_pairs
        ),
        "truth_fragments_in_qc_set": len(truth_ids),
        "truth_fragments_lost_to_the_host_filter": (
            len(truth_ids) - (viral["fragments_searched"] or 0)
        ),
    }
    composed["background"] = background
    composed["viral_part"] = viral
    if not composed["read_level_partition"]["partition_is_exact"]:
        log(
            "WARNING: background + viral fragments do not sum to Day 2's non-host count; "
            "the read half's composition is not exact and must not be quoted"
        )
    return composed


def _truth_ids(sample_dir: Path) -> set[str]:
    path = Path(sample_dir) / "truth" / "viral_fragment_ids.txt"
    if not path.exists():
        raise diamond.DiamondError(f"no ground-truth label file at {path}")
    with path.open() as handle:
        return {line.strip() for line in handle if line.strip()}


# --------------------------------------------------------------------------- series


def run_series(
    arm: ArmConfig,
    day2: dict,
    taxonomy,
    clades,
    *,
    levels: Sequence[str] | None = None,
    log=print,
) -> dict:
    """Assemble and search the background and each dilution level, whole.

    Every level is materialised from its frozen parts, assembled, searched and then
    deleted, so peak disk is one level rather than five. The read half of the arm reuses
    the background's search for every level, which is valid for the same reason Arm A's
    composition is: a per-read search has no cross-read dependence.
    """
    series = day2.get("dilution_series")
    if not series:
        raise diamond.DiamondError("the Day 2 provenance record holds no dilution series")
    out_dir = arm.out_dir / "series"
    denominators = series.get("denominators", {})
    parts = arm.work_dir / "series" / "parts"

    record: dict = {"levels": []}

    # The background on its own: this arm's false-positive baseline on the positive
    # sample, since every ground-truth viral fragment has been removed from it.
    bg1, bg2 = parts / "background_1.fastq.gz", parts / "background_2.fastq.gz"
    for path in (bg1, bg2):
        if not path.exists():
            raise diamond.DiamondError(f"{path} is missing; re-run the Day 2 pipeline")
    background = _contig_arm(arm, bg1, bg2, out_dir, taxonomy, clades, name="background", log=log)
    background["role"] = "background (viral fragments removed; any clade hit is a false positive)"
    record["background"] = background

    background_reads: dict | None = None
    if arm.read_level:
        background_reads = _read_arm(
            arm, bg1, bg2, out_dir, taxonomy, clades, name="background", log=log
        )
        record["background_read_level"] = background_reads

    wanted = set(levels) if levels else None
    materialised_dir = arm.work_dir / "series" / "materialised"
    for level in series["levels"]:
        name = level["name"]
        if wanted and name not in wanted:
            continue
        log(f"{name}: materialising {level['total_fragments']:,} fragments")
        r1, r2 = samples_pipeline.materialise_level(series, name, materialised_dir, verify=True)
        # This level's own ground truth, which needs no alignment and no label file: a
        # fragment in a level is viral exactly when it came from the viral part. So the
        # reads that map back onto a clade contig split cleanly into detections and
        # false positives, in the same unit Arm A reports.
        level_truth = reads.read_ids(_viral_part(series, name, 1))
        try:
            entry = _contig_arm(
                arm,
                r1,
                r2,
                out_dir,
                taxonomy,
                clades,
                name=name,
                truth_ids=level_truth,
                denominators=denominators,
                log=log,
            )
        finally:
            for path in (r1, r2):
                Path(path).unlink(missing_ok=True)
        entry["name"] = name
        entry["viral_fragments"] = level["viral_fragments"]
        entry["total_fragments"] = level["total_fragments"]
        entry["per_million"] = level.get("per_million", {})
        entry["viral_fragments_per_million"] = level.get("viral_fragments_per_million")

        if arm.read_level and background_reads is not None:
            viral_reads = _read_level_viral_part(
                arm, series, name, out_dir, taxonomy, clades, log=log
            )
            entry["read_level"] = _compose_read_level(
                background_reads, viral_reads, clades, denominators, level
            )
        record["levels"].append(entry)
    return record


def _viral_part(series: dict, level_name: str, mate: int) -> Path:
    """The path to one mate of a level's viral part, as Day 2 froze it.

    A level is stored as parts rather than as one file, which is what lets both halves
    of this arm treat the viral reads separately from the background: the read half
    searches them on their own, and the contig half uses their identifiers as the
    level's ground truth.
    """
    level = next((lv for lv in series["levels"] if lv["name"] == level_name), None)
    if level is None:
        raise diamond.DiamondError(f"no level {level_name} in the series")
    key = "read_1" if mate == 1 else "read_2"
    candidates = [p for p in level["parts"][key] if "background" not in Path(p["path"]).name]
    if len(candidates) != 1:
        raise diamond.DiamondError(
            f"{level_name}: expected one viral part for mate {mate}, got {len(candidates)}"
        )
    path = Path(candidates[0]["path"])
    if not path.exists():
        raise diamond.DiamondError(f"{path} is missing; re-run the Day 2 pipeline")
    return path


def _read_level_viral_part(
    arm: ArmConfig, series: dict, level_name: str, out_dir: Path, taxonomy, clades, *, log=print
) -> dict:
    """Search just this level's viral reads, which is what composition needs."""
    return _read_arm(
        arm,
        _viral_part(series, level_name, 1),
        _viral_part(series, level_name, 2),
        out_dir,
        taxonomy,
        clades,
        name=f"{level_name}_viral",
        log=log,
    )


def _compose_read_level(
    background: dict, viral: dict, clades, denominators: dict, level: dict
) -> dict:
    """Add a level's viral read hits to the shared background's.

    Valid because a translated search scores each read from the read alone. Both halves
    are kept separately in the record as well, so the sum is checkable rather than
    asserted: every fragment in the viral part is ground-truth viral and every fragment
    in the background is not, so detection and false positives need no labelling step.
    """
    detected = {name: viral["clades"].get(name, 0) for name in clades.members}
    false_positive = {name: background["clades"].get(name, 0) for name in clades.members}
    return {
        "composed_from": {
            "background": background["name"],
            "viral_part": viral["name"],
        },
        "fragments_with_a_hit": background["fragments_with_a_hit"] + viral["fragments_with_a_hit"],
        "clades": {name: detected[name] + false_positive[name] for name in clades.members},
        "detected_by_clade": detected,
        "false_positive_by_clade": false_positive,
        "detected": detected.get(HEADLINE_CLADE, 0),
        "false_positive": false_positive.get(HEADLINE_CLADE, 0),
        "sensitivity": (
            detected.get(HEADLINE_CLADE, 0) / level["viral_fragments"]
            if level.get("viral_fragments")
            else None
        ),
        "detected_per_million_post_qc": per_million(
            detected.get(HEADLINE_CLADE, 0), denominators.get("post_qc_pairs")
        ),
    }


# ------------------------------------------------------------------------- negatives


def run_negatives(
    arm: ArmConfig,
    samples_config: registry.Config,
    taxonomy,
    clades,
    *,
    only: Sequence[str] | None = None,
    log=print,
) -> list[dict]:
    """Fetch, filter, assemble, search and discard each negative in turn, at native depth.

    One sample is on disk at a time, in the pinned order, so an interrupted cohort
    resumes where it stopped and a sample already searched is skipped without
    re-fetching it.
    """
    specs = [s for s in registry.samples_from_config(samples_config) if not s.is_index]
    if only:
        wanted = set(only)
        specs = [s for s in specs if s.accession in wanted]
    deep_config = dataclasses.replace(samples_config, negatives_max_pairs=arm.negatives_max_pairs)
    out_dir = arm.out_dir / "negatives"
    records: list[dict] = []
    for position, spec in enumerate(specs, start=1):
        marker = out_dir / "search" / f".{spec.accession}_nonhost.{arm.sensitivity}.searched.json"
        log(f"[{position}/{len(specs)}] {spec.accession}")
        if marker.exists():
            log(f"{spec.accession}: already searched; not re-fetching")
        else:
            started = time.monotonic()
            sample = samples_pipeline.process_sample(spec, deep_config, truth_index=None, log=log)
            _write_sample_note(out_dir, spec.accession, sample, time.monotonic() - started)
        record = _search_one_negative(arm, spec.accession, taxonomy, clades, log=log)
        records.append(record)
        if arm.prune_after_search:
            freed = _prune_sample(arm.work_dir / "samples" / spec.accession)
            record["pruned_bytes"] = _remember_pruned(out_dir, spec.accession, freed)
            log(f"{spec.accession}: pruned {freed / 2**30:.2f} GiB")
    return records


def _search_one_negative(arm: ArmConfig, accession: str, taxonomy, clades, *, log=print) -> dict:
    sample_dir = arm.work_dir / "samples" / accession
    r1 = sample_dir / "nonhost" / f"{accession}_nonhost_1.fastq.gz"
    r2 = sample_dir / "nonhost" / f"{accession}_nonhost_2.fastq.gz"
    out_dir = arm.out_dir / "negatives"
    record = _contig_arm(
        arm,
        r1,
        r2,
        out_dir,
        taxonomy,
        clades,
        name=f"{accession}_nonhost",
        log=log,
    )
    record["accession"] = accession
    record["role"] = "negative"
    nonhost_pairs = record["assembly"].get("read_pairs_in") or 0
    reads_note = out_dir / f"{accession}.reads.json"
    if reads_note.exists():
        note = json.loads(reads_note.read_text())
        record["reads"] = {
            key: note.get(key)
            for key in (
                "raw_pairs",
                "qc_pairs",
                "redacted_pairs",
                "capped_pairs",
                "cap_fraction",
                "nonhost_pairs",
                "prepare_seconds",
                "source",
            )
        }
        nonhost_pairs = note.get("nonhost_pairs") or nonhost_pairs
    record["nonhost_pairs"] = nonhost_pairs
    pruned_note = out_dir / f"{accession}.pruned.json"
    if pruned_note.exists():
        record["pruned_bytes"] = json.loads(pruned_note.read_text()).get("freed_bytes", 0)

    contigs = record["contigs_in_headline_clade"]
    record["false_alarm_contigs"] = contigs
    # The endemic human coronaviruses, separated out for the same reason Arm A separates
    # them: a contig whose LCA is HKU1 in a respiratory patient is the arm working.
    detail = record["search"].get("contig_detail", [])
    record["endemic_human_coronavirus_contigs"] = sum(
        1 for entry in detail if entry.get("lca_taxid") in clades.endemic
    )
    record["false_alarm_contigs_excluding_endemic"] = sum(
        1
        for entry in detail
        if entry.get("lca_taxid") not in clades.endemic
        and HEADLINE_CLADE in (entry.get("lca_clades") or [])
    )
    if arm.negatives_read_level and arm.read_level:
        record["read_level"] = _read_arm(
            arm, r1, r2, out_dir, taxonomy, clades, name=f"{accession}_nonhost", log=log
        )
    log(
        f"{accession}: {contigs} contigs in {HEADLINE_CLADE} of "
        f"{record['assembly']['contig_count']:,} assembled"
    )
    return record


def _write_sample_note(out_dir: Path, accession: str, sample, prepare_seconds: float) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{accession}.reads.json").write_text(
        json.dumps(
            {
                **samples_pipeline._sample_record(sample),
                "prepare_seconds": round(prepare_seconds, 1),
            },
            indent=2,
        )
        + "\n"
    )


def _remember_pruned(out_dir: Path, accession: str, freed: int) -> int:
    """Keep the largest number of bytes this sample's pruning has ever freed."""
    marker = Path(out_dir) / f"{accession}.pruned.json"
    previous = 0
    if marker.exists():
        previous = json.loads(marker.read_text()).get("freed_bytes", 0)
    total = max(previous, freed)
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps({"freed_bytes": total}, indent=2) + "\n")
    return total


def _prune_sample(sample_dir: Path) -> int:
    """Delete every read file of a finished negative, keeping markers and reports."""
    freed = 0
    patterns = ("*.fastq.gz", "*.fastq", "*.sra", "*.sra.*", "*.fa.gz")
    for sub in ("raw", "qc", "capped", "nonhost", "tmp", "sra"):
        directory = Path(sample_dir) / sub
        if not directory.exists():
            continue
        for pattern in patterns:
            for path in directory.glob(pattern):
                if path.is_file():
                    freed += path.stat().st_size
                    path.unlink()
    return freed


def summarise_false_alarms(negatives: list[dict]) -> dict:
    """The number Day 5's threshold is drawn from, with its denominators stated.

    An assembly arm alarms per *sample*, not per read: the thing an analyst sees is a
    contig, and one coronavirus contig in a sample is one alarm whether it is backed by
    80 reads or 8,000. So the headline here is the count of samples with at least one
    clade contig, and the read-level figures are recorded beside it for Day 5's table.
    """
    if not negatives:
        return {}
    contigs = [n.get("false_alarm_contigs", 0) for n in negatives]
    novel = [n.get("false_alarm_contigs_excluding_endemic", 0) for n in negatives]
    nonhost = [n.get("nonhost_pairs", 0) or 0 for n in negatives]
    with_hits = [n for n, c in zip(negatives, contigs, strict=False) if c]
    with_novel = [n for n, c in zip(negatives, novel, strict=False) if c]
    return {
        "samples": len(negatives),
        "non_host_pairs_total": sum(nonhost),
        "clade": HEADLINE_CLADE,
        "unit": "contigs; an assembly arm alarms per sample",
        "samples_with_at_least_one_clade_contig": len(with_hits),
        "contigs_total": sum(contigs),
        "contigs_max_in_one_sample": max(contigs) if contigs else 0,
        "samples_with_hits": {
            n["accession"]: n.get("false_alarm_contigs", 0)
            for n in sorted(with_hits, key=lambda n: -n.get("false_alarm_contigs", 0))
        },
        "excluding_endemic_human_coronaviruses": {
            "contigs_total": sum(novel),
            "samples_with_at_least_one_clade_contig": len(with_novel),
            "contigs_max_in_one_sample": max(novel) if novel else 0,
            "samples_with_hits": {
                n["accession"]: n.get("false_alarm_contigs_excluding_endemic", 0)
                for n in sorted(
                    with_novel,
                    key=lambda n: -n.get("false_alarm_contigs_excluding_endemic", 0),
                )
            },
        },
        "endemic_human_coronavirus_contigs_total": sum(
            n.get("endemic_human_coronavirus_contigs", 0) for n in negatives
        ),
    }


# ------------------------------------------------------------------- the whole day


def run_day4(
    arm: ArmConfig,
    *,
    stages: Sequence[str] = ("database", "index", "series", "negatives"),
    only_negatives: Sequence[str] | None = None,
    only_levels: Sequence[str] | None = None,
    log=print,
) -> dict:
    started = time.monotonic()
    samples_config = registry.load_config(arm.samples_config)
    day2_path = arm.work_dir / "day2-provenance.json"
    if not day2_path.exists():
        day2_path = Path("results/day2-provenance.json")
    day2 = json.loads(day2_path.read_text())

    record: dict = {
        "tool": {"name": "hindsight arm-b", "version": _version()},
        "generated_at": utc_now(),
        "stages": list(stages),
        "machine": machine_record(arm.work_dir),
        "tools": tool_versions(arm.tool_pins),
        "parameters": {
            "cutoff": arm.cutoff,
            "source_db": str(arm.source_db),
            "dmnd_db": str(arm.dmnd_db),
            "threads": arm.threads,
            "sensitivity": arm.sensitivity,
            "comparison_sensitivity": arm.comparison_sensitivity,
            "evalue": arm.evalue,
            "max_target_seqs": arm.max_target_seqs,
            "lca_top_percent": arm.lca_top_percent,
            "min_contig_len": arm.min_contig_len,
            "read_level": arm.read_level,
            "read_level_sensitivity": arm.read_level_sensitivity,
            "negatives_max_pairs": arm.negatives_max_pairs,
            "negatives_depth": (
                "native" if arm.negatives_max_pairs == 0 else f"{arm.negatives_max_pairs} pairs"
            ),
            "negatives_read_level": arm.negatives_read_level,
            "headline_clade": HEADLINE_CLADE,
            "day2_provenance": str(day2_path),
        },
    }

    database_record, taxonomy, clades = prepare_database(arm, log=log)
    record["database"] = database_record

    if "index" in stages:
        record["index_sample"] = run_index(arm, samples_config, day2, taxonomy, clades, log=log)
    if "series" in stages:
        record["series"] = run_series(arm, day2, taxonomy, clades, levels=only_levels, log=log)
    if "negatives" in stages:
        record["negatives"] = run_negatives(
            arm, samples_config, taxonomy, clades, only=only_negatives, log=log
        )
        record["false_alarms"] = summarise_false_alarms(record["negatives"])

    record["wall_clock_seconds"] = round(time.monotonic() - started, 1)
    return merge(arm.provenance, record)


def merge(path: Path, fresh: dict) -> dict:
    """Merge a run into the record on disk, the way Days 2 and 3 do."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        merged = dict(fresh)
        merged["runs"] = [_stamp(fresh)]
        path.write_text(json.dumps(merged, indent=2) + "\n")
        return merged

    existing = json.loads(path.read_text())
    merged = dict(existing)
    for key, value in fresh.items():
        if key == "negatives":
            by_accession = {n["accession"]: n for n in existing.get("negatives", [])}
            for negative in value:
                by_accession[negative["accession"]] = negative
            merged["negatives"] = sorted(by_accession.values(), key=lambda n: n["accession"])
        elif key == "series" and value and existing.get("series"):
            # A partial series run (one level) must not drop the levels already done.
            levels = {lv["name"]: lv for lv in existing["series"].get("levels", [])}
            for level in value.get("levels", []):
                levels[level["name"]] = level
            combined = dict(existing["series"])
            combined.update({k: v for k, v in value.items() if k != "levels"})
            combined["levels"] = [
                levels[name] for name in _level_order(existing, value) if name in levels
            ]
            merged["series"] = combined
        elif key in ("index_sample", "series", "database") and not value:
            continue
        else:
            merged[key] = value
    if merged.get("negatives"):
        merged["false_alarms"] = summarise_false_alarms(merged["negatives"])
    runs = list(existing.get("runs", []))
    runs.append(_stamp(fresh))
    merged["runs"] = runs
    path.write_text(json.dumps(merged, indent=2) + "\n")
    return merged


def _level_order(existing: dict, fresh: dict) -> list[str]:
    order: list[str] = []
    for source in (existing.get("series") or {}, fresh):
        for level in source.get("levels", []) or []:
            if level["name"] not in order:
                order.append(level["name"])
    return order


def _stamp(record: dict) -> dict:
    return {
        "generated_at": record.get("generated_at"),
        "stages": record.get("stages"),
        "wall_clock_seconds": record.get("wall_clock_seconds"),
        "negatives": [n["accession"] for n in record.get("negatives", [])],
    }


def _version() -> str:
    from .. import __version__

    return __version__
