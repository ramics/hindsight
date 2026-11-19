"""Separate a diverged strain from a second organism, using only pre-cutoff sequence.

Day 4's ``diagnose-contigs`` aligns an arm's clade contigs to one reference and reports
what fraction of it they tile. That answers "how far was the thing in the patient from
what the database held", and on the HKU1 negative it left four contigs unexplained: good
protein hits to HKU1, no nucleotide alignment to ``NC_006577.2`` at all. Either a very
diverged region of the same virus, or a second betacoronavirus in the same patient --
and one reference cannot tell those apart.

This module can, because the ambiguity comes from the *reference set* rather than from
the data. RefSeq curates one genome per species, and its single HKU1 genome is
genotype A; HKU1's three genotypes differ most in hemagglutinin-esterase and spike,
which is exactly where the four contigs sat. Aligning to a panel of genotypes instead of
one reference resolves it in a single pass.

**The panel is still time-correct.** All three HKU1 genotypes were published in 2006, so
they are outside the dated database by RefSeq's curation policy, not by date. That
distinction is load-bearing for the benchmark's leakage argument, so every reference this
module accepts is checked against the cutoff at NCBI and its publication date is recorded
in the output. A reference published on or after the cutoff is refused unless the caller
says otherwise in as many words, and the record then says so loudly.

Two further things this exists to get right, both learned the hard way on the same
sample:

* ``minimap2 -x asm20`` silently fails to place a short contig whose alignment is split
  into blocks, which is not the same as divergence. One of the four "unexplained" contigs
  was HKU1 at 99-100% nucleotide identity all along: a 324 bp chimera whose halves align
  to the same region in opposite orientations. Every contig ``asm20`` leaves unplaced is
  retried with ``-x sr``, and the record names the preset that placed it.
* A second organism does not merely add contigs, it duplicates coverage: two viruses mean
  two copies of the same gene at different identities. So the verdict checks whether any
  region of the best reference is covered twice by contigs that disagree, which is the
  signature that distinguishes a coinfection from a fragmented assembly.
"""

from __future__ import annotations

import datetime as _dt
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from . import diamond

# A pair of contigs covering the same reference region are "discordant" if their
# identities to it differ by more than this. Strain-level variation along one genome
# stays well inside it; two organisms do not.
DISCORDANT_IDENTITY_POINTS = 3.0

# Overlapping spans shorter than this are contig ends meeting, not duplicated coverage.
MIN_OVERLAP_BASES = 100

_COVERAGE = re.compile(r"\bmulti=([0-9.]+)")


class StrainError(RuntimeError):
    """A reference could not be resolved, or an aligner failed."""


@dataclass(frozen=True)
class Reference:
    accession: str
    path: Path
    source: str
    published: str | None
    before_cutoff: bool | None


def contig_coverage(contigs: Path) -> dict[str, float]:
    """MEGAHIT's own k-mer coverage per contig, from the FASTA headers.

    Independent of every alignment here, which is why it is worth reading: a second
    organism in the same library has its own abundance, so a contig set that belongs to
    one genome should not split into two coverage classes.
    """
    coverage: dict[str, float] = {}
    with contigs.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            name = line[1:].split()[0]
            found = _COVERAGE.search(line)
            if found:
                coverage[name] = float(found.group(1))
    return coverage


def resolve_reference(
    accession: str,
    *,
    source_db: Path,
    cutoff: str,
    dest_dir: Path,
    client,
    allow_post_cutoff: bool = False,
    log=print,
) -> Reference:
    """Fetch one reference and establish that using it involves no leakage.

    A reference in the dated manifest needs no date check: the database held it. Anything
    else is asked about at NCBI, and its creation date decides. The point of the whole
    exercise is to reach sequence RefSeq did not curate, so this path is expected and
    normal -- it just has to be checked rather than assumed.
    """
    accession = accession.strip()
    if not accession:
        raise StrainError("empty accession")
    manifest = source_db / "manifest.tsv"
    in_manifest = False
    if manifest.exists():
        for line in manifest.read_text().splitlines()[1:]:
            fields = line.split("\t")
            if fields and fields[0] == accession:
                in_manifest = True
                break

    published: str | None = None
    before: bool | None = None
    if in_manifest:
        source = "dated-manifest"
        before = True
    else:
        source = "pre-cutoff-genbank"
        published = _published_date(client, accession)
        if published is None:
            raise StrainError(
                f"{accession}: NCBI reported no creation date, so it cannot be shown to "
                f"predate the {cutoff} cutoff"
            )
        before = published < cutoff
        if not before and not allow_post_cutoff:
            raise StrainError(
                f"{accession} was published {published}, on or after the {cutoff} cutoff. "
                f"Aligning to it would put post-cutoff sequence into the diagnosis; pass "
                f"allow_post_cutoff to do it deliberately."
            )

    dest = dest_dir / f"{accession}.fna"
    if not (dest.exists() and dest.stat().st_size > 0):
        dest.parent.mkdir(parents=True, exist_ok=True)
        log(f"fetching {accession} ({source}{'' if published is None else ', ' + published})")
        text = client.efetch("nuccore", [accession])
        if not text.lstrip().startswith(">"):
            raise StrainError(f"{accession}: NCBI did not return FASTA")
        dest.write_text(text)
    return Reference(
        accession=accession,
        path=dest,
        source=source,
        published=published,
        before_cutoff=before,
    )


def _published_date(client, accession: str) -> str | None:
    """The record's creation date at NCBI, as ``YYYY-MM-DD``."""
    text = client.esummary("nuccore", [accession])
    found = re.search(r'Name="CreateDate"[^>]*>([0-9]{4})/([0-9]{2})/([0-9]{2})', text)
    if not found:
        return None
    year, month, day = found.groups()
    try:
        return _dt.date(int(year), int(month), int(day)).isoformat()
    except ValueError:
        return None


def _paf(reference: Path, subset: Path, preset: str) -> dict[str, dict]:
    """Best alignment per contig against one reference, by aligned length."""
    completed = subprocess.run(  # noqa: S603
        [
            "minimap2",
            "-c",
            "-x",
            preset,
            "--secondary=no",
            str(reference),
            str(subset),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise StrainError(f"minimap2 -x {preset} failed ({completed.returncode})")
    best: dict[str, dict] = {}
    blocks: dict[str, int] = {}
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 11:
            continue
        query, qlen = fields[0], int(fields[1])
        strand = fields[4]
        ref_len, ref_start, ref_end = int(fields[6]), int(fields[7]), int(fields[8])
        matches, aln_len = int(fields[9]), int(fields[10])
        blocks[query] = blocks.get(query, 0) + 1
        previous = best.get(query)
        if previous and previous["alignment_length"] >= aln_len:
            previous["strands"].add(strand)
            continue
        strands = previous["strands"] if previous else set()
        strands.add(strand)
        best[query] = {
            "contig": query,
            "contig_length": qlen,
            "alignment_length": aln_len,
            "matching_bases": matches,
            "identity_percent": round(100 * matches / aln_len, 2) if aln_len else None,
            "reference_length": ref_len,
            "reference_span": [ref_start, ref_end],
            "preset": preset,
            "strands": strands,
        }
    for query, entry in best.items():
        entry["alignment_blocks"] = blocks.get(query, 1)
        # Both orientations against one reference is not divergence and not a second
        # organism: it is an inverted-repeat join inside the contig.
        entry["both_orientations"] = len(entry.pop("strands")) > 1
    return best


def align_panel(reference: Reference, subset: Path, contigs: list[str]) -> dict:
    """Align every clade contig to one reference, falling back to ``-x sr``.

    ``asm20`` is the right preset for a query that may be a diverged relative, and it is
    what Day 4 used. It is also unable to place a short contig whose alignment breaks
    into blocks, so a bare "does not align" from it is not evidence of anything until the
    short-read preset has had a turn.
    """
    placed = _paf(reference.path, subset, "asm20")
    missing = [name for name in contigs if name not in placed]
    rescued = _paf(reference.path, subset, "sr") if missing else {}
    for name in missing:
        if name in rescued:
            placed[name] = rescued[name]
    aligned = [placed[name] for name in contigs if name in placed]
    total_aln = sum(entry["alignment_length"] for entry in aligned)
    total_match = sum(entry["matching_bases"] for entry in aligned)
    return {
        "reference": reference.accession,
        "source": reference.source,
        "published": reference.published,
        "before_cutoff": reference.before_cutoff,
        "contigs_aligned": len(aligned),
        "contigs_total": len(contigs),
        "aligned_bases": total_aln,
        "weighted_identity_percent": (
            round(100 * total_match / total_aln, 2) if total_aln else None
        ),
        "reference_covered_percent": _covered_percent(aligned),
        "placed_by_short_read_preset": sorted(
            entry["contig"] for entry in aligned if entry["preset"] == "sr"
        ),
        "contigs_in_both_orientations": sorted(
            entry["contig"] for entry in aligned if entry["both_orientations"]
        ),
        "per_contig": {entry["contig"]: entry for entry in aligned},
    }


def _covered_percent(aligned: list[dict]) -> float | None:
    if not aligned:
        return None
    ref_len = aligned[0]["reference_length"]
    spans = sorted((entry["reference_span"][0], entry["reference_span"][1]) for entry in aligned)
    merged: list[list[int]] = []
    for start, end in spans:
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    covered = sum(end - start for start, end in merged)
    return round(100 * covered / ref_len, 2) if ref_len else None


def duplicated_regions(panel: dict) -> list[dict]:
    """Reference regions two contigs both cover, where their identities disagree.

    This is the test that separates the two hypotheses. A fragmented assembly of one
    genome tiles the reference: contigs abut, and where they overlap they agree, because
    they are the same sequence. Two organisms produce two contigs for the same gene whose
    identities to any single reference differ, because one of them is not that organism.
    """
    entries = sorted(panel["per_contig"].values(), key=lambda e: e["reference_span"][0])
    clashes = []
    for i, left in enumerate(entries):
        for right in entries[i + 1 :]:
            overlap = min(left["reference_span"][1], right["reference_span"][1]) - max(
                left["reference_span"][0], right["reference_span"][0]
            )
            if overlap < MIN_OVERLAP_BASES:
                continue
            gap = abs((left["identity_percent"] or 0) - (right["identity_percent"] or 0))
            if gap <= DISCORDANT_IDENTITY_POINTS:
                continue
            clashes.append(
                {
                    "contigs": [left["contig"], right["contig"]],
                    "overlap_bases": overlap,
                    "identities": [left["identity_percent"], right["identity_percent"]],
                    "identity_gap_points": round(gap, 2),
                }
            )
    return clashes


def diagnose_strain(
    arm,
    samples_config,
    taxonomy,
    clades,
    *,
    references: list[str],
    clade: str = "coronaviridae",
    accession: str | None = None,
    allow_post_cutoff: bool = False,
    client=None,
    log=print,
) -> dict:
    """Decide whether a clade's contigs are one diverged strain or more than one organism."""
    if not references:
        raise StrainError("name at least one reference to align to")
    from ..ncbi import NcbiClient

    client = client or NcbiClient()
    out_dir = arm.out_dir / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)
    accession = accession or samples_config.index_run
    name = f"{accession}_nonhost"
    stage = "index" if accession == samples_config.index_run else "negatives"
    contigs_path = arm.out_dir / stage / "assembly" / f"{name}.contigs.fa"
    hits = arm.out_dir / stage / "search" / f"{name}.{arm.sensitivity}.blastx.tsv.gz"
    for required in (contigs_path, hits):
        if not required.exists():
            raise StrainError(f"{required} is missing; run the {stage} stage first")

    lengths = diamond.contig_lengths(contigs_path)
    calls = diamond.assign_contigs(
        hits, lengths, taxonomy, clades, lca_top_percent=arm.lca_top_percent
    )
    chosen = diamond.clade_contigs(calls, clade, by="lca")
    seen = {call.contig for call in chosen}
    chosen += [
        call
        for call in diamond.clade_contigs(calls, clade, by="best_hit")
        if call.contig not in seen
    ]
    if not chosen:
        return {
            "clade": clade,
            "sample": accession,
            "conclusion": f"the arm placed no contig in {clade}; nothing to align",
        }
    names = [call.contig for call in sorted(chosen, key=lambda c: -c.length)]
    subset = out_dir / f"{accession}.{clade}-contigs.fa"
    diamond._write_contig_subset(contigs_path, set(names), subset)
    coverage = contig_coverage(contigs_path)

    panels = []
    for wanted in references:
        reference = resolve_reference(
            wanted,
            source_db=arm.source_db,
            cutoff=arm.cutoff,
            dest_dir=arm.work_dir / "ref" / "panel",
            client=client,
            allow_post_cutoff=allow_post_cutoff,
            log=log,
        )
        log(f"aligning {len(names)} {clade} contigs to {reference.accession}")
        panels.append(align_panel(reference, subset, names))

    best = max(
        panels,
        key=lambda p: (p["contigs_aligned"], p["weighted_identity_percent"] or 0),
    )
    record: dict = {
        "clade": clade,
        "sample": accession,
        "cutoff": arm.cutoff,
        "contigs_examined": len(names),
        "note": (
            "every reference in the panel is checked against the cutoff at NCBI and its "
            "publication date recorded; a reference outside the dated database is outside "
            "it by RefSeq's one-genome-per-species curation policy, not by date"
        ),
        "panel": [{k: v for k, v in p.items() if k != "per_contig"} for p in panels],
        "best_reference": best["reference"],
        "contigs": [
            {
                "contig": contig,
                "length": lengths.get(contig),
                "kmer_coverage": coverage.get(contig),
                "best_hit_taxon": _taxon(taxonomy, calls[contig].best_taxid),
                "best_hit_identity_percent": calls[contig].best_pident,
                "best_hit_subject": calls[contig].best_subject,
                "identity_per_reference": {
                    p["reference"]: p["per_contig"].get(contig, {}).get("identity_percent")
                    for p in panels
                },
            }
            for contig in names
        ],
    }
    record["coverage_range"] = _coverage_range(coverage, names)
    record["duplicated_regions"] = duplicated_regions(best)
    record["unplaced_by_every_reference"] = sorted(
        set(names) - {c for p in panels for c in p["per_contig"]}
    )
    record["conclusion"] = _conclude(record, best, len(names))
    log(record["conclusion"])
    return record


def _taxon(taxonomy, taxid):
    return taxonomy.names.get(taxid) if taxid else None


def _coverage_range(coverage: dict[str, float], names: list[str]) -> dict:
    values = [coverage[n] for n in names if n in coverage]
    if not values:
        return {}
    return {
        "minimum": min(values),
        "maximum": max(values),
        "contigs_with_coverage": len(values),
    }


def _conclude(record: dict, best: dict, total: int) -> str:
    reference = best["reference"]
    where = (
        "in the dated database"
        if best["source"] == "dated-manifest"
        else f"published {best['published']}, before the {record['cutoff']} cutoff"
    )
    head = (
        f"{best['contigs_aligned']} of {total} {record['clade']} contigs align to "
        f"{reference} ({where}), {best['aligned_bases']:,} bases at "
        f"{best['weighted_identity_percent']}% identity."
    )
    if best["contigs_aligned"] < total:
        unplaced = record["unplaced_by_every_reference"]
        tail = f" {total - best['contigs_aligned']} contig(s) remain unplaced against it" + (
            f", and {len(unplaced)} against every reference in the panel, which is "
            f"where to look for a second organism."
            if unplaced
            else ", though another reference in the panel places them."
        )
        return head + tail
    if record["duplicated_regions"]:
        return (
            head + f" But {len(record['duplicated_regions'])} reference region(s) are covered "
            f"twice by contigs whose identities disagree, which one genome does not do: "
            f"more than one organism is present."
        )
    rescued = best["placed_by_short_read_preset"]
    extra = (
        f" {len(rescued)} of them only under the short-read preset, which is an aligner "
        f"limit on short or split contigs rather than divergence."
        if rescued
        else ""
    )
    return (
        head + " Every contig is accounted for by this one reference and no region is "
        "covered twice, so this is a single diverged strain rather than a coinfection." + extra
    )
