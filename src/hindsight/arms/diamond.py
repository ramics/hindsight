"""Arm B: assemble, then search translated protein against a time-correct protein set.

This arm automates what Wu et al. actually did. They did not run a classifier. They
took the reads left after host filtering, assembled them de novo with MEGAHIT, and
matched the longest contig against prior coronavirus sequences — at 79-88% nucleotide
identity, which is to say by homology rather than by exact match. The whole reason this
arm exists is that Arm A's question ("would a classifier have called it?") and this one
("would an assembly-plus-homology-search have called it?") are different questions, and
the second is the one the discovery actually answered.

Four things in here are decisions rather than mechanics.

**The search is translated, and that is the point.** SARS-CoV-2 is 12-21% divergent in
nucleotide from anything in the November 2019 database, which is enough to defeat exact
35-mers over much of the genome. Protein is far better conserved: the ORF1ab
polyprotein of a sarbecovirus is recognisable across that distance where the nucleotide
sequence is not. So the arm searches contigs against the dated *protein* set —
``viral.faa``, 373,129 records at the versions release 97 held — with
``diamond blastx``.

**DIAMOND runs at ``--very-sensitive``, unlike Arm A which runs at Kraken2's
defaults**, and that asymmetry is deliberate rather than sloppy. Kraken2's defaults are
what a 2019 pipeline would have run. DIAMOND's default mode is tuned for high-identity
matches between closely related proteins, and running a *homology* arm in a mode that
cannot see homology would measure the mode and not the method. The honest way to hold
both is to record the alternative: the index sample and the series are also searched in
default mode, so Day 5 can see exactly what the sensitivity setting is worth and the
write-up can quote both.

**A dilution level has to be materialised and assembled whole.** Arm A composes a level
from its parts because Kraken2 assigns each fragment independently of the rest of the
file. Assembly is the opposite: a contig is a statement about the *set* of reads, and
removing 49,990 viral reads changes which contigs exist, their length and their
coverage. There is nothing to compose. Five levels therefore cost five full assemblies,
which is most of this arm's runtime and is not optimisable away.

**The arm's output is contigs; the benchmark's unit is reads.** A contig is the thing
an analyst looks at, so the alarm is recorded per contig — how many, how long, how
similar to what. But Day 5 needs one table across both arms, in reads per million, so
the clade-assigned contigs are made into a bowtie2 index and the non-host reads are
mapped back to them. A read that aligns to a contig the arm placed in Coronaviridae is
this arm's equivalent of a read Kraken2 placed in Coronaviridae.
"""

from __future__ import annotations

import gzip
import json
import re
import shutil
import subprocess
import time
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ..samples import reads
from .taxa import CladeSets, Taxonomy, per_million

# DIAMOND's tabular output, in the order the arm writes and reads it. Recorded as a
# constant because the column order is part of the on-disk format: a re-run with a
# different --outfmt would silently change what column 5 means.
BLAST_FIELDS = (
    "qseqid",
    "sseqid",
    "staxids",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
)

# How wide a band of hits the LCA considers, as a fraction of the best bitscore for that
# contig. 10% is the convention MEGAN established and most published pipelines use. It
# matters here more than usual: a sarbecovirus contig in a 2019 database hits SARS-CoV
# and several bat coronaviruses at similar scores, and taking only the single best hit
# would make the arm's answer depend on which of a cluster of near-ties won.
LCA_TOP_PERCENT = 10.0

SENSITIVITY_FLAGS: dict[str, list[str]] = {
    "default": [],
    "fast": ["--fast"],
    "mid-sensitive": ["--mid-sensitive"],
    "sensitive": ["--sensitive"],
    "more-sensitive": ["--more-sensitive"],
    "very-sensitive": ["--very-sensitive"],
    "ultra-sensitive": ["--ultra-sensitive"],
}


class DiamondError(RuntimeError):
    pass


# --------------------------------------------------------------------------- database


def write_taxonmap(seqid_map: Path, protein_fasta: Path, dest: Path) -> dict:
    """Turn hindsight's ``seqid2taxid.map`` into the file ``diamond makedb`` wants.

    DIAMOND expects NCBI's ``prot.accession2taxid`` layout: a header line, then
    ``accession``, ``accession.version``, ``taxid``, ``gi``. hindsight writes a
    two-column ``accession<TAB>taxid`` map covering both molecules, so this narrows it
    to the accessions actually present in ``viral.faa`` and widens it to four columns.

    Narrowing matters: the map holds 385,308 entries across nucleotide and protein, and
    handing DIAMOND nucleotide accessions it will never see is how a taxonomy mapping
    ends up looking complete when it is not. The count of protein records that got a
    taxid is returned so a gap shows up in the record rather than as unassigned contigs
    later.
    """
    seqid_map, protein_fasta, dest = Path(seqid_map), Path(protein_fasta), Path(dest)
    wanted: set[str] = set()
    with protein_fasta.open() as handle:
        for line in handle:
            if line.startswith(">"):
                wanted.add(line[1:].split()[0])
    written = 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    with seqid_map.open() as src, dest.open("w") as out:
        out.write("accession\taccession.version\ttaxid\tgi\n")
        for line in src:
            accession, _, taxid = line.rstrip("\n").partition("\t")
            if accession not in wanted:
                continue
            try:
                value = int(taxid)
            except ValueError:
                continue
            bare = accession.split(".")[0]
            out.write(f"{bare}\t{accession}\t{value}\t0\n")
            written += 1
    return {
        "path": str(dest),
        "protein_records": len(wanted),
        "records_with_a_taxid": written,
        "records_without_a_taxid": len(wanted) - written,
    }


def build_database(
    source_db: Path,
    dmnd_path: Path,
    *,
    threads: int = 4,
    log=print,
) -> dict:
    """Make a DIAMOND database from a ``hindsight build`` output directory.

    Only the taxon *map* is handed to DIAMOND, not ``nodes.dmp`` and ``names.dmp``,
    which is a deliberate narrowing. DIAMOND can store a tree inside the database and do
    its own LCA, but this arm computes the LCA itself from the dated dumps (see
    :func:`assign_contigs`) because the 2019 tree's shape is the thing under test — Day
    3's finding was that *Bat coronavirus BM48-31/BGR/2008* sits outside Betacoronavirus
    in November 2019 and inside it today. Two copies of the tree can only ever disagree,
    and the one the scorer reads is the one that decides the answer.

    Dropping them also removes a failure mode: ``makedb`` validates every rank string in
    ``nodes.dmp`` against its own list and aborts on one it does not recognise, so a
    dated taxdump carrying an unfashionable rank would fail the build for no benefit.
    ``--taxonmap`` alone is what populates ``staxids`` in the output, which is all the
    arm reads.

    Resumable the way Arm A's build is: a finished ``.dmnd`` is left alone.
    """
    source_db, dmnd_path = Path(source_db), Path(dmnd_path)
    faa = source_db / "viral.faa"
    seqid_map = source_db / "seqid2taxid.map"
    taxonomy_src = source_db / "taxonomy"
    for required in (faa, seqid_map, taxonomy_src / "nodes.dmp", taxonomy_src / "names.dmp"):
        if not required.exists():
            raise DiamondError(
                f"{required} is missing; build the protein database first with "
                "`hindsight build --date 2019-12-01 --out <db> --molecules protein`"
            )

    record: dict = {"db": str(dmnd_path), "source_db": str(source_db)}
    with faa.open() as handle:
        record["protein_records_in_source"] = sum(1 for line in handle if line.startswith(">"))

    taxonmap = dmnd_path.parent / "prot.accession2taxid"
    if not taxonmap.exists():
        log(f"writing the dated protein taxon map to {taxonmap}")
    record["taxonmap"] = (
        write_taxonmap(seqid_map, faa, taxonmap)
        if not taxonmap.exists()
        else {"path": str(taxonmap), "reused": True}
    )

    if dmnd_path.exists():
        log(f"diamond database already built at {dmnd_path}")
        record["reused"] = True
    else:
        record["reused"] = False
        dmnd_path.parent.mkdir(parents=True, exist_ok=True)
        log(f"diamond makedb -> {dmnd_path} on {threads} threads")
        start = time.monotonic()
        log_path = dmnd_path.with_suffix(".makedb.log")
        command = [
            "diamond",
            "makedb",
            "--in",
            str(faa),
            "--db",
            str(dmnd_path),
            "--threads",
            str(threads),
            "--taxonmap",
            str(taxonmap),
        ]
        with log_path.open("ab") as errlog:
            completed = subprocess.run(  # noqa: S603
                command, stdout=errlog, stderr=errlog, check=False
            )
        if completed.returncode != 0:
            tail = log_path.read_text(errors="replace")[-2000:]
            raise DiamondError(f"diamond makedb failed ({completed.returncode}):\n{tail}")
        record["build_seconds"] = round(time.monotonic() - start, 1)
        record["log"] = str(log_path)

    actual = dmnd_path if dmnd_path.exists() else dmnd_path.with_suffix(".dmnd")
    if actual.exists():
        record["bytes"] = actual.stat().st_size
    return record


# --------------------------------------------------------------------------- assembly


@dataclass
class Assembly:
    """One MEGAHIT run over one read pair."""

    name: str
    contigs: Path
    contig_count: int = 0
    total_bases: int = 0
    longest: int = 0
    n50: int = 0
    seconds: float = 0.0
    reads_in: int = 0

    def record(self) -> dict:
        return {
            "name": self.name,
            "contigs": str(self.contigs),
            "contig_count": self.contig_count,
            "total_bases": self.total_bases,
            "longest_contig": self.longest,
            "n50": self.n50,
            "seconds": round(self.seconds, 1),
            "read_pairs_in": self.reads_in,
        }


def contig_lengths(path: Path) -> dict[str, int]:
    """Contig name -> length, from a FASTA, without holding the sequence."""
    lengths: dict[str, int] = {}
    name: str | None = None
    total = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:  # type: ignore[operator]
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = total
                name = line[1:].split()[0]
                total = 0
            else:
                total += len(line.strip())
    if name is not None:
        lengths[name] = total
    return lengths


def _n50(lengths: Sequence[int]) -> int:
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2
    running = 0
    for length in ordered:
        running += length
        if running >= half:
            return length
    return ordered[-1]


def assemble(
    r1: Path,
    r2: Path,
    out_dir: Path,
    *,
    name: str,
    threads: int = 4,
    memory_fraction: float = 0.8,
    min_contig_len: int = 200,
    extra: Sequence[str] = (),
    resume: bool = True,
    log=print,
) -> Assembly:
    """Assemble one paired FASTQ with MEGAHIT and keep only the contigs.

    MEGAHIT's defaults are kept — the k-mer list included — because this arm is meant to
    be the method Wu et al. used, and they used it as it comes. The intermediate
    directory is deleted on success: it is several times the size of the contigs, holds
    nothing the arm reads again, and at native depth across twenty negatives it is the
    difference between fitting on this machine and not.

    MEGAHIT refuses to write into an existing output directory, so a resumed run is
    detected from the contigs file rather than by letting MEGAHIT decide.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    contigs = out_dir / f"{name}.contigs.fa"
    marker = out_dir / f".{name}.assembled.json"
    megahit_dir = out_dir / f"{name}.megahit"

    if resume and contigs.exists() and marker.exists():
        log(f"{name}: already assembled")
        payload = json.loads(marker.read_text())
        return Assembly(
            name=name,
            contigs=contigs,
            contig_count=payload.get("contig_count", 0),
            total_bases=payload.get("total_bases", 0),
            longest=payload.get("longest_contig", 0),
            n50=payload.get("n50", 0),
            seconds=payload.get("seconds", 0.0),
            reads_in=payload.get("read_pairs_in", 0),
        )

    if megahit_dir.exists():
        shutil.rmtree(megahit_dir)
    log(f"{name}: megahit on {threads} threads")
    start = time.monotonic()
    log_path = out_dir / f"{name}.megahit.log"
    command = [
        "megahit",
        "-1",
        str(r1),
        "-2",
        str(r2),
        "-o",
        str(megahit_dir),
        "-t",
        str(threads),
        "-m",
        str(memory_fraction),
        "--min-contig-len",
        str(min_contig_len),
        *extra,
    ]
    with log_path.open("ab") as errlog:
        completed = subprocess.run(  # noqa: S603
            command, stdout=errlog, stderr=errlog, check=False
        )
    produced = megahit_dir / "final.contigs.fa"
    if completed.returncode != 0:
        tail = log_path.read_text(errors="replace")[-2000:]
        # An assembly that finds nothing is not an error, but MEGAHIT exits non-zero
        # when it cannot build a single contig, which happens at the bottom of the
        # dilution series and is a *result* rather than a failure. Distinguish them by
        # what MEGAHIT said rather than by the exit code alone.
        if "Too few contigs" in tail or "empty" in tail.lower():
            log(f"{name}: megahit produced no contigs")
            contigs.write_text("")
        else:
            raise DiamondError(f"{name}: megahit failed ({completed.returncode}):\n{tail}")
    else:
        shutil.copyfile(produced, contigs)

    lengths = list(contig_lengths(contigs).values())
    result = Assembly(
        name=name,
        contigs=contigs,
        contig_count=len(lengths),
        total_bases=sum(lengths),
        longest=max(lengths) if lengths else 0,
        n50=_n50(lengths),
        seconds=time.monotonic() - start,
        reads_in=_megahit_reads_in(log_path),
    )
    marker.write_text(json.dumps(result.record(), indent=2) + "\n")
    if megahit_dir.exists():
        shutil.rmtree(megahit_dir)
    log(
        f"{name}: {result.contig_count:,} contigs, longest {result.longest:,} bp, "
        f"N50 {result.n50:,} ({result.seconds / 60:.1f} min)"
    )
    return result


def _megahit_reads_in(log_path: Path) -> int:
    """How many read pairs MEGAHIT actually took in, from its own log.

    Recorded because it is the denominator the assembly worked from, and because it is
    a cheap cross-check on the reads side: if it disagrees with the pipeline's non-host
    count, one of the two is wrong about the file.
    """
    if not log_path.exists():
        return 0
    text = log_path.read_text(errors="replace")
    counts = re.findall(r"([0-9]+) reads, [0-9]+ max length", text)
    if not counts:
        return 0
    # MEGAHIT reports this once per input file and counts single reads across both
    # mates, so the last line is the total and a pair is two reads.
    return int(counts[-1]) // 2


# ----------------------------------------------------------------------------- search


def blastx(
    dmnd_path: Path,
    contigs: Path,
    out_dir: Path,
    *,
    name: str,
    threads: int = 4,
    sensitivity: str = "very-sensitive",
    evalue: float = 1e-3,
    max_target_seqs: int = 25,
    block_size: float | None = None,
    resume: bool = True,
    log=print,
) -> tuple[Path, dict]:
    """Search contigs against the dated protein database, keeping the tabular hits.

    ``max_target_seqs`` is 25 rather than DIAMOND's default 25-per-query-frame
    behaviour being left implicit: the LCA needs the *cluster* of near-tied hits, not
    the single best, and 25 is comfortably more than the number of sarbecovirus and bat
    coronavirus proteins in a 2019 database that a real coronavirus contig ties across.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = out_dir / f"{name}.{sensitivity}.blastx.tsv.gz"
    marker = out_dir / f".{name}.{sensitivity}.searched.json"
    if sensitivity not in SENSITIVITY_FLAGS:
        raise DiamondError(f"unknown sensitivity {sensitivity!r}; have {sorted(SENSITIVITY_FLAGS)}")

    if resume and hits.exists() and marker.exists():
        log(f"{name}: already searched ({sensitivity})")
        return hits, json.loads(marker.read_text())

    if not Path(contigs).exists() or Path(contigs).stat().st_size == 0:
        log(f"{name}: no contigs to search")
        with gzip.open(hits, "wt"):
            pass
        payload = {
            "name": name,
            "sensitivity": sensitivity,
            "contigs": 0,
            "hit_rows": 0,
            "seconds": 0.0,
            "skipped": "no contigs",
        }
        marker.write_text(json.dumps(payload, indent=2) + "\n")
        return hits, payload

    log(f"{name}: diamond blastx ({sensitivity}) on {threads} threads")
    start = time.monotonic()
    log_path = out_dir / f"{name}.{sensitivity}.diamond.log"
    raw = out_dir / f"{name}.{sensitivity}.blastx.tsv"
    command = [
        "diamond",
        "blastx",
        "--db",
        str(dmnd_path),
        "--query",
        str(contigs),
        "--out",
        str(raw),
        "--outfmt",
        "6",
        *BLAST_FIELDS,
        "--threads",
        str(threads),
        "--evalue",
        str(evalue),
        "--max-target-seqs",
        str(max_target_seqs),
        *SENSITIVITY_FLAGS[sensitivity],
    ]
    if block_size:
        command += ["--block-size", str(block_size)]
    with log_path.open("ab") as errlog:
        completed = subprocess.run(  # noqa: S603
            command, stdout=errlog, stderr=errlog, check=False
        )
    if completed.returncode != 0:
        raw.unlink(missing_ok=True)
        tail = log_path.read_text(errors="replace")[-2000:]
        raise DiamondError(f"{name}: diamond blastx failed ({completed.returncode}):\n{tail}")

    rows = 0
    with raw.open("rb") as src, gzip.open(hits, "wb") as sink:
        for line in src:
            rows += 1
            sink.write(line)
    raw.unlink(missing_ok=True)

    payload = {
        "name": name,
        "sensitivity": sensitivity,
        "evalue": evalue,
        "max_target_seqs": max_target_seqs,
        "fields": list(BLAST_FIELDS),
        "hit_rows": rows,
        "seconds": round(time.monotonic() - start, 1),
        "hits": str(hits),
        "log": str(log_path),
    }
    marker.write_text(json.dumps(payload, indent=2) + "\n")
    log(f"{name}: {rows:,} hit rows ({payload['seconds'] / 60:.1f} min)")
    return hits, payload


@dataclass
class Hit:
    contig: str
    subject: str
    taxids: list[int]
    pident: float
    length: int
    evalue: float
    bitscore: float
    qlen: int


def iter_hits(path: Path):
    """Yield :class:`Hit` rows from a kept tabular search output."""
    index = {field: i for i, field in enumerate(BLAST_FIELDS)}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:  # type: ignore[operator]
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(BLAST_FIELDS):
                continue
            raw_taxids = fields[index["staxids"]]
            taxids = []
            for piece in raw_taxids.replace(" ", "").split(";"):
                if piece.isdigit():
                    taxids.append(int(piece))
            yield Hit(
                contig=fields[index["qseqid"]],
                subject=fields[index["sseqid"]],
                taxids=taxids,
                pident=float(fields[index["pident"]]),
                length=int(fields[index["length"]]),
                evalue=float(fields[index["evalue"]]),
                bitscore=float(fields[index["bitscore"]]),
                qlen=int(fields[index["qlen"]]),
            )


# ------------------------------------------------------------------ contig assignment


@dataclass
class ContigCall:
    """What the arm decided one contig is, and the evidence it decided from."""

    contig: str
    length: int
    lca: int | None
    best_taxid: int | None
    best_subject: str = ""
    best_pident: float = 0.0
    best_bitscore: float = 0.0
    best_evalue: float = 0.0
    best_aln_length: int = 0
    hits_considered: int = 0
    hits_total: int = 0
    lca_clades: list[str] = field(default_factory=list)
    best_hit_clades: list[str] = field(default_factory=list)

    def record(self, taxonomy: Taxonomy) -> dict:
        return {
            "contig": self.contig,
            "length": self.length,
            "lca_taxid": self.lca,
            "lca_taxon": taxonomy.names.get(self.lca) if self.lca else None,
            "lca_rank": taxonomy.ranks.get(self.lca) if self.lca else None,
            "lca_clades": self.lca_clades,
            "best_hit": {
                "subject": self.best_subject,
                "taxid": self.best_taxid,
                "taxon": taxonomy.names.get(self.best_taxid) if self.best_taxid else None,
                "identity_percent": self.best_pident,
                "alignment_length_aa": self.best_aln_length,
                "bitscore": self.best_bitscore,
                "evalue": self.best_evalue,
                "clades": self.best_hit_clades,
            },
            "hits_in_lca_band": self.hits_considered,
            "hits_total": self.hits_total,
        }


def assign_contigs(
    hits_path: Path,
    lengths: dict[str, int],
    taxonomy: Taxonomy,
    clades: CladeSets,
    *,
    lca_top_percent: float = LCA_TOP_PERCENT,
) -> dict[str, ContigCall]:
    """Decide a taxon per contig: best hit, and an LCA over the near-tied hits.

    Both are kept because they answer different questions and can disagree in a way
    that matters here. A sarbecovirus contig's best hit in a 2019 database is some
    particular bat coronavirus protein, which is a precise-looking answer that is
    precisely wrong — the contig is not that virus. The LCA over everything within
    ``lca_top_percent`` of the best bitscore says "a sarbecovirus", which is both
    vaguer and true, and is what a human reading the BLAST table would have concluded.

    Day 5 can score either. The arm's own headline uses the LCA.
    """
    grouped: dict[str, list[Hit]] = {}
    for hit in iter_hits(hits_path):
        grouped.setdefault(hit.contig, []).append(hit)

    calls: dict[str, ContigCall] = {}
    for contig, contig_hits in grouped.items():
        contig_hits.sort(key=lambda h: (-h.bitscore, h.evalue, h.subject))
        best = contig_hits[0]
        cutoff = best.bitscore * (1.0 - lca_top_percent / 100.0)
        band = [h for h in contig_hits if h.bitscore >= cutoff]
        band_taxids = [taxid for h in band for taxid in h.taxids]
        lca = taxonomy.lowest_common_ancestor(band_taxids)
        best_taxid = best.taxids[0] if best.taxids else None
        calls[contig] = ContigCall(
            contig=contig,
            length=lengths.get(contig, best.qlen),
            lca=lca,
            best_taxid=best_taxid,
            best_subject=best.subject,
            best_pident=best.pident,
            best_bitscore=best.bitscore,
            best_evalue=best.evalue,
            best_aln_length=best.length,
            hits_considered=len(band),
            hits_total=len(contig_hits),
            lca_clades=clades.which(lca) if lca else [],
            best_hit_clades=clades.which(best_taxid) if best_taxid else [],
        )
    return calls


def clade_contigs(calls: dict[str, ContigCall], clade: str, *, by: str = "lca") -> list[ContigCall]:
    """Every contig the arm placed in ``clade``, longest first."""
    key = "lca_clades" if by == "lca" else "best_hit_clades"
    chosen = [call for call in calls.values() if clade in getattr(call, key)]
    return sorted(chosen, key=lambda c: (-c.length, c.contig))


def summarise_calls(
    calls: dict[str, ContigCall],
    clades: CladeSets,
    taxonomy: Taxonomy,
    *,
    detail_clade: str = "coronaviridae",
    detail_limit: int = 50,
) -> dict:
    """Per-clade contig and base counts, by LCA and by best hit, plus the detail table."""
    summary: dict = {
        "contigs_with_a_hit": len(calls),
        "clades": {},
        "clades_by_best_hit": {},
    }
    for name in clades.members:
        by_lca = clade_contigs(calls, name, by="lca")
        by_best = clade_contigs(calls, name, by="best_hit")
        summary["clades"][name] = {
            "contigs": len(by_lca),
            "bases": sum(c.length for c in by_lca),
            "longest_contig": by_lca[0].length if by_lca else 0,
            "best_identity_percent": max((c.best_pident for c in by_lca), default=None),
        }
        summary["clades_by_best_hit"][name] = {
            "contigs": len(by_best),
            "bases": sum(c.length for c in by_best),
            "longest_contig": by_best[0].length if by_best else 0,
        }
    detail = clade_contigs(calls, detail_clade, by="lca")
    # Union with the best-hit set, so a contig whose LCA climbed out of the family but
    # whose best hit is a coronavirus still appears in the evidence table. That gap is
    # exactly where a borderline call lives, and it should not be invisible.
    seen = {c.contig for c in detail}
    for call in clade_contigs(calls, detail_clade, by="best_hit"):
        if call.contig not in seen:
            detail.append(call)
    summary["detail_clade"] = detail_clade
    summary["contigs_in_detail_clade"] = len(detail)
    summary["contig_detail"] = [
        call.record(taxonomy)
        for call in sorted(detail, key=lambda c: (-c.length, c.contig))[:detail_limit]
    ]
    lca_counts = Counter(call.lca for call in calls.values() if call.lca)
    summary["top_lca_taxa"] = [
        {"taxid": taxid, "contigs": n, "taxon": taxonomy.names.get(taxid)}
        for taxid, n in lca_counts.most_common(15)
    ]
    return summary


# ------------------------------------------------------- reads back onto the contigs


def count_reads_on_contigs(
    contigs: Path,
    names: Sequence[str],
    r1: Path,
    r2: Path,
    out_dir: Path,
    *,
    name: str,
    threads: int = 4,
    truth_ids: set[str] | None = None,
    log=print,
) -> dict:
    """Count the non-host reads whose best home in the assembly is a clade contig.

    This is what makes Arm B's numbers comparable with Arm A's. Arm A counts reads
    placed in a clade; this counts reads that align to a contig the arm placed in that
    clade. Both are "reads this arm would have shown an analyst".

    **Two passes, and the second one is the point.** Mapping quality means "unique among
    the sequences you gave the aligner", not "this read came from here". An index holding
    only the clade contigs leaves every read from the rest of the library nowhere else to
    go, so reads land on whatever thirty bases they happen to share and are scored
    confidently for it. Measured: one negative's 23 coronavirus contigs took 1,435,005
    read pairs against a 23-contig index and 2,885 against the whole 7,438-contig
    assembly, a factor of 498, where Arm A's count for that sample was 2,399.

    So the subset index is used only as a candidate filter, which is sound because adding
    references to an index can take reads away from a contig but never add them, and the
    candidates are then re-adjudicated against the whole assembly. That costs one cheap
    pass over the library plus one expensive pass over a few hundred thousand reads,
    rather than the whole library against the whole assembly.

    Local alignment, matching the host filter, and the full mapping-quality distribution
    is recorded rather than a single threshold applied: Day 5 should be able to redraw
    the line without a re-run. ``truth_ids`` splits the adjudicated reads into real viral
    fragments and the rest, which is this arm's sensitivity and its false-positive count
    on one sample.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    record: dict = {"name": name, "contigs_used": len(names)}
    if not names:
        record["note"] = "no contigs in the clade; nothing to map"
        record["reads_aligned"] = 0
        return record

    # This is the one step in the arm that needs the reads themselves, and the negatives
    # are pruned as each finishes, so resumability and prune-as-you-go collide here and
    # nowhere else: a resumed run finds the assembly and the blastx output cached and the
    # reads gone. Only a negative with contigs in the clade reaches this far, which is
    # why it survived nine clean samples before showing itself. Persisting the record
    # makes the two features compatible.
    marker = out_dir / f".{name}.readback.json"
    if marker.exists():
        cached = json.loads(marker.read_text())
        log(f"{name}: reusing the cached read mapping ({cached.get('reads_aligned', 0)} reads)")
        return cached
    if not r1.exists() or not r2.exists():
        accession = name.removesuffix("_nonhost")
        raise DiamondError(
            f"{name}: the reads have been pruned and there is no cached read mapping at "
            f"{marker}. Re-stage the sample, then run this stage again:\n"
            f"  hindsight samples run --only {accession}"
        )

    subset = out_dir / f"{name}.clade-contigs.fa"
    kept = _write_contig_subset(contigs, set(names), subset)
    record["contig_bases"] = kept
    index = out_dir / f"{name}.clade-contigs"
    build_log = out_dir / f"{name}.bowtie2-build.log"
    build_log.unlink(missing_ok=True)
    _bowtie2_build(subset, index, threads=threads, log_path=build_log, name=name)

    log(f"{name}: pass 1, candidate reads against {len(names)} clade contigs")
    align_log = out_dir / f"{name}.bowtie2.log"
    start = time.monotonic()
    candidates = out_dir / f"{name}.candidates.fa"
    found = _align_and_collect(
        ["bowtie2", "--local", "--no-unal", "-p", str(threads), "-x", str(index)]
        + ["-1", str(r1), "-2", str(r2)],
        align_log,
        name=name,
        keep=None,
        write_fasta=candidates,
    )
    record["candidate_reads"] = len(found["reads"])
    record["candidate_mate_alignments"] = sum(found["mapqs"].values())
    record["pass1_alignment_summary"] = (
        align_log.read_text(errors="replace").strip().splitlines()[-6:]
    )
    log(f"{name}: {len(found['reads']):,} candidate reads; adjudicating on the whole assembly")

    whole = out_dir / f"{name}.whole-assembly"
    _bowtie2_build(contigs, whole, threads=threads, log_path=build_log, name=name)
    adjudicate_log = out_dir / f"{name}.bowtie2-adjudicate.log"
    final = _align_and_collect(
        # -f because the candidates arrive as FASTA: the first pass saved the aligner's
        # own SEQ field, which carries no quality string.
        ["bowtie2", "--local", "--no-unal", "-f", "-p", str(threads), "-x", str(whole)]
        + ["-U", str(candidates)],
        adjudicate_log,
        name=name,
        keep=set(names),
        write_fasta=None,
    )
    aligned = final["reads"]
    mapqs = final["mapqs"]

    record["seconds"] = round(time.monotonic() - start, 1)
    record["reads_aligned"] = len(aligned)
    record["mapq_distribution"] = dict(sorted(mapqs.items()))
    record["mate_alignments"] = sum(mapqs.values())
    record["alignment_summary"] = (
        adjudicate_log.read_text(errors="replace").strip().splitlines()[-6:]
    )
    # Kept because the ratio is the measurement that justifies the second pass, and
    # because a future change to either index would show up here first.
    record["subset_index_inflation"] = (
        round(record["candidate_reads"] / len(aligned), 1) if aligned else None
    )
    if truth_ids is not None:
        hit_truth = aligned & truth_ids
        record["truth_fragments"] = len(truth_ids)
        record["detected"] = len(hit_truth)
        record["false_positive"] = len(aligned) - len(hit_truth)
        record["sensitivity"] = len(hit_truth) / len(truth_ids) if truth_ids else None
    marker.write_text(json.dumps(record, indent=2) + "\n")
    # Both indexes and the candidate FASTA are rebuildable from the contigs and the name
    # list, and the whole-assembly index is the largest thing this step makes.
    for stale in out_dir.glob(f"{name}.clade-contigs*.bt2"):
        stale.unlink(missing_ok=True)
    for stale in out_dir.glob(f"{name}.whole-assembly*.bt2"):
        stale.unlink(missing_ok=True)
    candidates.unlink(missing_ok=True)
    return record


def _bowtie2_build(fasta: Path, index: Path, *, threads: int, log_path: Path, name: str) -> None:
    with log_path.open("ab") as errlog:
        completed = subprocess.run(  # noqa: S603
            ["bowtie2-build", "--threads", str(threads), str(fasta), str(index)],
            stdout=errlog,
            stderr=errlog,
            check=False,
        )
    if completed.returncode != 0:
        raise DiamondError(f"{name}: bowtie2-build on {fasta.name} failed ({completed.returncode})")


def _align_and_collect(
    argv: Sequence[str],
    log_path: Path,
    *,
    name: str,
    keep: set[str] | None,
    write_fasta: Path | None,
) -> dict:
    """Run bowtie2, counting primary alignments and optionally saving the reads.

    ``keep`` restricts the count to alignments landing on those references; ``None``
    counts every alignment. ``write_fasta`` saves each counted mate's sequence, which is
    how the candidate set reaches the second pass without a second trip through the
    FASTQ. The sequence is written as the aligner reported it, which is reverse
    complemented for a mate on the minus strand — immaterial, since the next pass aligns
    it to the same assembly either way.
    """
    mapqs: Counter[int] = Counter()
    seen: set[str] = set()
    out = write_fasta.open("w") if write_fasta is not None else None
    try:
        with log_path.open("wb") as errlog:
            bowtie = subprocess.Popen(  # noqa: S603
                list(argv),
                stdout=subprocess.PIPE,
                stderr=errlog,
                text=True,
                bufsize=1 << 20,
            )
            assert bowtie.stdout is not None
            for line in bowtie.stdout:
                if line.startswith("@"):
                    continue
                fields = line.split("\t", 11)
                if len(fields) < 10 or fields[2] == "*":
                    continue
                flag = int(fields[1])
                if flag & 0x100 or flag & 0x800:  # secondary or supplementary
                    continue
                if keep is not None and fields[2] not in keep:
                    continue
                mapqs[int(fields[4])] += 1
                read_id = reads.normalise_id(fields[0])
                seen.add(read_id)
                if out is not None and fields[9] != "*":
                    # Mates of a pair share a normalised id, so the mate number keeps the
                    # two records distinct in the candidate FASTA.
                    mate = 2 if flag & 0x80 else 1
                    out.write(f">{read_id}/{mate}\n{fields[9]}\n")
            bowtie.stdout.close()
            rc = bowtie.wait()
    finally:
        if out is not None:
            out.close()
    if rc != 0:
        tail = log_path.read_text(errors="replace")[-1500:]
        raise DiamondError(f"{name}: bowtie2 failed ({rc}):\n{tail}")
    return {"mapqs": mapqs, "reads": seen}


def _write_contig_subset(contigs: Path, names: set[str], dest: Path) -> int:
    kept = 0
    writing = False
    with Path(contigs).open() as src, Path(dest).open("w") as out:
        for line in src:
            if line.startswith(">"):
                writing = line[1:].split()[0] in names
            if writing:
                out.write(line)
                if not line.startswith(">"):
                    kept += len(line.strip())
    return kept


# ------------------------------------------------------------------------- the score


# Bitscore cuts the evidence profile is reported against. Not thresholds — the arm
# applies none, and `alarm` below is the unfiltered rule. They are reporting bands, fixed
# in advance at round numbers rather than chosen to fit a curve, so that Day 5 can read
# the series and the negatives against the same ladder and place one cut itself.
BITSCORE_BANDS = (50.0, 100.0, 250.0, 500.0, 1000.0)


def evidence_profile(
    calls: dict[str, ContigCall],
    clade: str,
    *,
    bands: Sequence[float] = BITSCORE_BANDS,
) -> dict:
    """How strong the clade evidence is, not merely whether there is any.

    A translated search against a database of distant relatives always finds *something*:
    on this cohort every assembly, the spike-free background included, yields a handful of
    short contigs whose best protein match sits near 30% identity at a bitscore in the
    forties. They are real remote homology and the arm is right to report them, but a rule
    that alarms on "any contig in the family" therefore alarms on a clean sample, and the
    bare contig count cannot tell a 29 kb genome from that floor.

    Worse, the count is not even monotone in dose. At high coverage the genome assembles
    into one contig; at intermediate coverage it shatters into twenty short ones, so the
    count *rises* as the signal falls. Bases and best bitscore both behave; the count does
    not, and this profile records all three so the write-up cannot lean on the wrong one.
    """
    chosen = clade_contigs(calls, clade, by="lca")
    scores = sorted((call.best_bitscore or 0.0) for call in chosen)
    profile = {
        "contigs": len(chosen),
        "bases": sum(call.length for call in chosen),
        "longest_contig": chosen[0].length if chosen else 0,
        "best_bitscore": max(scores, default=None),
        "bands": [],
    }
    for cut in bands:
        above = [call for call in chosen if (call.best_bitscore or 0.0) >= cut]
        profile["bands"].append(
            {
                "min_bitscore": cut,
                "contigs": len(above),
                "bases": sum(call.length for call in above),
                "longest_contig": max((call.length for call in above), default=0),
            }
        )
    return profile


def score_sample(
    assembly: Assembly,
    calls: dict[str, ContigCall],
    clades: CladeSets,
    taxonomy: Taxonomy,
    *,
    denominators: dict[str, int] | None = None,
    headline_clade: str = "coronaviridae",
) -> dict:
    """Pull one sample's assembly and search together into the numbers Day 5 reads."""
    summary = summarise_calls(calls, clades, taxonomy)
    record = {
        "assembly": assembly.record(),
        "search": summary,
        "headline_clade": headline_clade,
        # The unfiltered rule, recorded as the arm actually behaves: any contig the LCA
        # placed in the headline clade. `evidence` is what distinguishes this firing from
        # the background floor, and `docs/day4-arm-b.md` reports both.
        "alarm": summary["clades"].get(headline_clade, {}).get("contigs", 0) > 0,
        "contigs_in_headline_clade": summary["clades"].get(headline_clade, {}).get("contigs", 0),
        "evidence": evidence_profile(calls, headline_clade),
    }
    if denominators:
        record["contigs_per_million"] = {
            key: per_million(record["contigs_in_headline_clade"], value)
            for key, value in denominators.items()
        }
    return record
