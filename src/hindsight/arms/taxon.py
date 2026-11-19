"""Ask whether a classifier's call is an organism or a laboratory sequence.

A read count is not evidence of an infection. Arm A's largest viral call in one
pre-COVID negative was 11,227 reads of *Equine infectious anemia virus*, a horse
lentivirus, in a human respiratory library -- four times the count of the HKU1
infection the same sample genuinely had. The usual explanations for that shape of
result are host-derived repeats or bare homology to reverse transcriptase, and both
would have made it a false positive to count against the arm.

Neither was true. The reads are 99% identical to the EIAV genome and confined to a
single 1,679 bp window inside *pol*, with zero coverage of gag, env or either LTR and
edges that go from nothing to 500-fold depth in two bases. That is not a virus; it is a
cloned fragment, and the thing to do with it is name it rather than argue about it.

**Where the reads sit on the reference is the discriminator, and a count cannot see it.**
An organism replicating in a patient gives coverage across its genome, uneven but
everywhere. A construct gives a bounded window at high identity and high depth, with
sharp edges. So this module takes the reads a classifier assigned to one taxon, aligns
them to a reference for that taxon, and reports the coverage profile and the shape of
its edges alongside the verdict those imply.

Everything here is a diagnostic, run after an arm has finished, and nothing it computes
feeds back into an arm's score. References go through the same cutoff gate the strain
panel uses, so a diagnosis cannot quietly introduce post-cutoff sequence.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from ..samples import reads as reads_mod
from . import kraken, strain

# A covered stretch is treated as one window if the unseen stretch between it and the
# next is shorter than this. Short dropouts inside a covered region are sampling, not
# structure.
WINDOW_JOIN_BASES = 50

# Below this, "covers most of the reference" is not a fair description, and a bounded
# window is the better reading.
WHOLE_GENOME_COVERAGE_PERCENT = 70.0

# An edge is measured rather than asserted: the number of bases depth takes to go from
# most of the window's own median to almost none of it. A cloned insert ends where the
# insert ends, so that transition is as short as the read length allows; coverage of a
# real genome tapers over hundreds of bases or does not end at all. The width is reported
# whatever it is, so a reader can draw the line somewhere else than here.
EDGE_HIGH_FRACTION = 0.5
EDGE_LOW_FRACTION = 0.05
EDGE_SEARCH_BASES = 400
# Within a read length of the boundary, a transition is as sharp as this assay can show.
EDGE_IS_A_STEP_BASES = 20


class TaxonError(RuntimeError):
    """The classifier output, the reads or the reference could not be used."""


def diagnose_taxon(
    arm,
    samples_config,
    clades,
    *,
    taxid: int,
    reference: str,
    accession: str | None = None,
    include_descendants: bool = True,
    allow_post_cutoff: bool = False,
    client=None,
    log=print,
) -> dict:
    """Say whether the reads a classifier put in ``taxid`` look like an organism."""
    from ..ncbi import NcbiClient

    client = client or NcbiClient()
    out_dir = arm.out_dir / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)
    accession = accession or samples_config.index_run
    name = f"{accession}_nonhost"
    stage = "index" if accession == samples_config.index_run else "negatives"
    classified = arm.out_dir / stage / f"{name}.classified.tsv.gz"
    if not classified.exists():
        raise TaxonError(f"{classified} is missing; classify {accession} first")

    taxonomy = kraken.read_taxonomy(arm.source_db / "taxonomy")
    wanted_taxids = taxonomy.descendants(taxid) if include_descendants else {taxid}
    fragments = {
        fragment
        for fragment, assigned in kraken.iter_output(classified)
        if assigned in wanted_taxids
    }
    record: dict = {
        "sample": accession,
        "taxid": taxid,
        "taxon": taxonomy.names.get(taxid),
        "clades": clades.which(taxid),
        "descendant_taxa_included": len(wanted_taxids) if include_descendants else 1,
        "fragments_assigned": len(fragments),
        "cutoff": arm.cutoff,
        "note": ("a diagnostic, run after the arm; nothing here feeds back into the arm's score"),
    }
    if not fragments:
        record["conclusion"] = f"no read in {accession} was assigned to {taxid}"
        return record

    sample_dir = arm.work_dir / "samples" / accession
    mates = []
    for mate in (1, 2):
        source = sample_dir / "nonhost" / f"{accession}_nonhost_{mate}.fastq.gz"
        if not source.exists():
            raise TaxonError(
                f"{source} is missing; this sample's reads were pruned, so re-stage it "
                f"before diagnosing a taxon in it"
            )
        dest = out_dir / f"{accession}.taxon{taxid}_{mate}.fastq.gz"
        mates.append((dest, reads_mod.write_subset(source, dest, fragments)))
    record["fragments_extracted"] = [count for _dest, count in mates]

    resolved = strain.resolve_reference(
        reference,
        source_db=arm.source_db,
        cutoff=arm.cutoff,
        dest_dir=arm.work_dir / "ref" / "panel",
        client=client,
        allow_post_cutoff=allow_post_cutoff,
        log=log,
    )
    record["reference"] = {
        "accession": resolved.accession,
        "source": resolved.source,
        "published": resolved.published,
        "before_cutoff": resolved.before_cutoff,
    }

    log(f"aligning {len(fragments):,} {record['taxon'] or taxid} fragments to {reference}")
    depth, stats = _align_and_profile(
        resolved.path,
        mates[0][0],
        mates[1][0],
        out_dir / f"{accession}.taxon{taxid}",
        threads=arm.threads,
    )
    record.update(stats)
    record["coverage"] = _profile(depth)
    record["conclusion"] = _conclude(record)
    log(record["conclusion"])
    return record


def _align_and_profile(
    reference: Path, r1: Path, r2: Path, prefix: Path, *, threads: int
) -> tuple[list[int], dict]:
    """bowtie2 --local, then per-position depth over the whole reference.

    ``--local`` and bowtie2 for the same reason Day 3's background diagnostic uses them:
    the host filter and the truth label were both local alignments, and a diagnostic that
    changed aligner would be answering a slightly different question.
    """
    index = prefix.with_suffix(".bt2idx")
    log_path = prefix.with_suffix(".bowtie2.log")
    with log_path.open("wb") as errlog:
        built = subprocess.run(  # noqa: S603
            ["bowtie2-build", "--threads", str(threads), str(reference), str(index)],
            stdout=errlog,
            stderr=errlog,
            check=False,
        )
        if built.returncode != 0:
            raise TaxonError(f"bowtie2-build failed ({built.returncode}); see {log_path}")
        bam = prefix.with_suffix(".bam")
        bowtie = subprocess.Popen(  # noqa: S603
            [
                "bowtie2",
                "--local",
                "-p",
                str(threads),
                "-x",
                str(index),
                "-1",
                str(r1),
                "-2",
                str(r2),
            ],
            stdout=subprocess.PIPE,
            stderr=errlog,
        )
        sort = subprocess.Popen(  # noqa: S603
            ["samtools", "sort", "-@", str(max(1, threads - 1)), "-o", str(bam), "-"],
            stdin=bowtie.stdout,
            stdout=errlog,
            stderr=errlog,
        )
        bowtie.stdout.close()
        codes = {"samtools sort": sort.wait(), "bowtie2": bowtie.wait()}
    if any(codes.values()):
        raise TaxonError(f"diagnostic alignment failed: {codes} (see {log_path})")

    stats_text = subprocess.run(  # noqa: S603
        ["samtools", "stats", str(bam)], capture_output=True, text=True, check=True
    ).stdout
    stats: dict = {}
    for line in stats_text.splitlines():
        if line.startswith("SN\treads mapped:"):
            stats["reads_mapped"] = int(line.split("\t")[2])
        elif line.startswith("SN\tsequences:"):
            stats["reads_examined"] = int(line.split("\t")[2])
        elif line.startswith("SN\terror rate:"):
            rate = float(line.split("\t")[2].split("#")[0].strip())
            stats["identity_percent"] = round(100 * (1 - rate), 2)
    if stats.get("reads_examined"):
        stats["mapped_fraction"] = round(stats["reads_mapped"] / stats["reads_examined"], 4)
    # With nothing mapped, samtools reports an error rate of 0, which this would otherwise
    # publish as "100% identity" beside "0 reads aligned" — the most misleading pair of
    # numbers the diagnostic can produce, and the exact case that matters most, since a
    # taxon whose reads do not align at all is the one being called spurious.
    if not stats.get("reads_mapped"):
        stats["identity_percent"] = None

    lengths = subprocess.run(  # noqa: S603
        ["samtools", "view", "-H", str(bam)], capture_output=True, text=True, check=True
    ).stdout
    ref_len = 0
    for line in lengths.splitlines():
        if line.startswith("@SQ"):
            for field in line.split("\t"):
                if field.startswith("LN:"):
                    ref_len = max(ref_len, int(field[3:]))
    depth = [0] * ref_len
    subprocess.run(["samtools", "index", str(bam)], check=True)  # noqa: S603
    text = subprocess.run(  # noqa: S603
        ["samtools", "depth", "-a", str(bam)], capture_output=True, text=True, check=True
    ).stdout
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3:
            position = int(fields[1])
            if 1 <= position <= ref_len:
                depth[position - 1] = int(fields[2])
    return depth, stats


def _profile(depth: list[int]) -> dict:
    """Covered windows, their depths, and whether their edges are steps or slopes."""
    ref_len = len(depth)
    covered = sum(1 for d in depth if d > 0)
    windows: list[list[int]] = []
    for position, value in enumerate(depth):
        if value <= 0:
            continue
        if windows and position - windows[-1][1] <= WINDOW_JOIN_BASES:
            windows[-1][1] = position
        else:
            windows.append([position, position])
    described = []
    for start, end in sorted(windows, key=lambda w: -(w[1] - w[0]))[:10]:
        inside = [d for d in depth[start : end + 1] if d > 0]
        inside.sort()
        median = inside[len(inside) // 2] if inside else 0
        left = _edge_width(depth, start, median, leftwards=True)
        right = _edge_width(depth, end, median, leftwards=False)
        described.append(
            {
                "start": start + 1,
                "end": end + 1,
                "bases": end - start + 1,
                "median_depth": median,
                "max_depth": max(depth[start : end + 1]) if end >= start else 0,
                "left_edge_width_bases": left,
                "right_edge_width_bases": right,
                "edges_are_steps": (
                    left is not None
                    and right is not None
                    and left <= EDGE_IS_A_STEP_BASES
                    and right <= EDGE_IS_A_STEP_BASES
                ),
            }
        )
    return {
        "reference_length": ref_len,
        "covered_bases": covered,
        "covered_percent": round(100 * covered / ref_len, 2) if ref_len else None,
        "mean_depth": round(sum(depth) / ref_len, 1) if ref_len else None,
        "max_depth": max(depth) if depth else 0,
        "windows": len(windows),
        "largest_windows": described,
    }


def _edge_width(depth: list[int], boundary: int, median: int, *, leftwards: bool) -> int | None:
    """Bases taken to fall from most of ``median`` to almost none of it at one boundary.

    Walking outwards from the boundary rather than inwards, because the question is how
    abruptly coverage stops. ``None`` when the transition is not found inside the search
    window, which happens when the window abuts another covered stretch.
    """
    if median <= 0:
        return None
    high = EDGE_HIGH_FRACTION * median
    low = EDGE_LOW_FRACTION * median
    step = -1 if leftwards else 1
    # Walk outwards from well inside the window, so that the high-coverage side of the
    # transition is seen before the low one. Starting at the boundary itself would begin
    # on ground that is already low and find no transition at all.
    inside = [boundary - step * offset for offset in range(EDGE_SEARCH_BASES, 0, -1)]
    outside = [boundary + step * offset for offset in range(EDGE_SEARCH_BASES)]
    last_high: int | None = None
    for index, position in enumerate(inside + outside):
        if not 0 <= position < len(depth):
            continue
        value = depth[position]
        if value >= high:
            last_high = index
        elif value <= low and last_high is not None:
            return index - last_high
    return None


def _conclude(record: dict) -> str:
    coverage = record["coverage"]
    taxon = record["taxon"] or record["taxid"]
    reference = record["reference"]["accession"]
    identity = record.get("identity_percent")
    mapped = record.get("reads_mapped", 0)
    if not mapped:
        return (
            f"{record['fragments_assigned']:,} fragments were assigned to {taxon}, and "
            f"**not one of them aligns to {reference}**. The call rests on shared k-mers "
            f"rather than on the organism's sequence being present: whatever these reads "
            f"are, they are not this."
        )
    head = (
        f"{record['fragments_assigned']:,} fragments were assigned to {taxon}; "
        f"{mapped:,} reads align to {reference} at "
        f"{identity}% identity, covering {coverage['covered_percent']}% of it "
        f"({coverage['covered_bases']:,} of {coverage['reference_length']:,} bases) "
        f"at a mean depth of {coverage['mean_depth']}."
    )
    if identity is not None and identity < 90:
        return (
            head + " The identity is too low for this to be the organism itself: these "
            "reads are homologous to it rather than from it."
        )
    if coverage["covered_percent"] and coverage["covered_percent"] >= WHOLE_GENOME_COVERAGE_PERCENT:
        return (
            head + " Coverage spans the reference, which is what an organism present in "
            "the sample looks like."
        )
    biggest = coverage["largest_windows"][0] if coverage["largest_windows"] else None
    if not biggest:
        return head + " Nothing aligned well enough to describe."
    share = round(100 * biggest["bases"] / max(1, coverage["covered_bases"]), 1)
    edges = biggest["edges_are_steps"]
    shape = (
        f" {share}% of the covered bases are one {biggest['bases']:,} bp window "
        f"({biggest['start']:,}-{biggest['end']:,}) at median depth "
        f"{biggest['median_depth']}"
    )
    if edges:
        return (
            head
            + shape
            + f", and its coverage stops within {biggest['left_edge_width_bases']} and "
            f"{biggest['right_edge_width_bases']} bases at its two ends. A genome in a "
            "patient is not covered that way: this is a cloned fragment of the reference "
            "-- a vector, a plasmid or a spiked control -- and not the organism."
        )
    return (
        head + shape + ". Coverage is confined to part of the reference rather than spanning "
        "it, so this is a fragment of the organism's sequence rather than evidence the "
        "organism was present; the window's edges are not clean enough to call it a "
        "construct outright."
    )
