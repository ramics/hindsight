"""The Day 2 pipeline: accession in, frozen benchmark samples out.

One command runs all of it, every stage is resumable, and every number the later arms
depend on lands in a JSON record beside the data. The stages, in order:

1. **fetch** — ENA FASTQ, resumed until whole, MD5 verified.
2. **qc** — fastp: adapters, quality, length. The JSON report is kept.
3. **truth** (index sample only) — align to the modern SARS-CoV-2 genome to *label*
   which fragments are viral. This is the one place the pipeline is allowed to know
   something from 2020, and it is used only to decide what the ground truth is. The
   detection arms on Days 3 and 4 never see this label and never see this reference.
4. **host** — align to GRCh38 and keep the fragments where neither mate mapped.
   Host filtering is where a careless pipeline loses viral reads, so the viral label
   is computed *before* this step and the overlap is reported: how many known-viral
   fragments the host filter would have thrown away.
5. **background** (index sample only) — the non-host fragments with the viral ones
   removed. This is the matrix the dilution series is held in: the index sample's own
   background, not a synthetic one, so the series varies viral abundance and nothing
   else.
6. **dilute** (index sample only) — background plus a seeded subsample of the viral
   fragments, at each level.

The negatives go through 1, 2 and 4 only. They have no viral label and no dilution.
"""

from __future__ import annotations

import datetime as dt
import json
import platform
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import requests

from .. import __version__
from . import dilute as dilute_mod
from . import ena, reads, registry, sra, tools

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


@dataclass
class Timing:
    stage: str
    seconds: float
    note: str = ""


@dataclass
class SampleResult:
    accession: str
    role: str
    layout: str = ""
    source: str = ""
    raw_bytes: int = 0
    raw_pairs: int = 0
    qc_pairs: int = 0
    redacted_pairs: int = 0
    non_redacted_pairs: int = 0
    capped_pairs: int = 0
    pruned_bytes: int = 0
    cap_fraction: float = 1.0
    host_pairs: int = 0
    nonhost_pairs: int = 0
    viral_pairs: int = 0
    viral_lost_to_host: int = 0
    background_pairs: int = 0
    viral_fraction_raw: float = 0.0
    viral_fraction_qc: float = 0.0
    files: dict[str, dict] = field(default_factory=dict)
    timings: list[Timing] = field(default_factory=list)
    ena_metadata: dict = field(default_factory=dict)


class Stage:
    """A tiny resumability helper: a stage is done when its marker exists."""

    def __init__(self, directory: Path, name: str):
        self.marker = Path(directory) / f".{name}.done"

    @property
    def done(self) -> bool:
        return self.marker.exists()

    def finish(self, payload: dict | None = None) -> None:
        self.marker.parent.mkdir(parents=True, exist_ok=True)
        self.marker.write_text(json.dumps(payload or {}, indent=2) + "\n")

    def load(self) -> dict:
        return json.loads(self.marker.read_text()) if self.marker.exists() else {}


def fetch_reference(accession: str, dest: Path, *, session: requests.Session | None = None) -> Path:
    """Fetch one nucleotide record as FASTA, for ground-truth labelling only."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    sess = session or requests.Session()
    response = sess.get(
        EUTILS,
        params={"db": "nuccore", "id": accession, "rettype": "fasta", "retmode": "text"},
        timeout=120,
    )
    response.raise_for_status()
    text = response.text
    # E-utilities signals a rate-limit rejection with HTTP 200 and a JSON body, which
    # would otherwise be written into the FASTA as a sequence.
    if not text.startswith(">"):
        raise RuntimeError(f"efetch did not return FASTA for {accession}: {text[:200]!r}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(text)
    return dest


def build_index(fasta: Path, prefix: Path, *, threads: int, log: Path) -> Path:
    prefix = Path(prefix)
    if (prefix.parent / f"{prefix.name}.1.bt2").exists():
        return prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    tools.run(
        ["bowtie2-build", "--threads", str(threads), str(fasta), str(prefix)],
        log=log,
    )
    return prefix


def _pipe_host_filter(
    *,
    index: Path,
    r1: Path,
    r2: Path,
    out1: Path,
    out2: Path,
    threads: int,
    extra: list[str],
    log: Path,
) -> None:
    """bowtie2 | samtools: keep fragments where *neither* mate mapped.

    ``--un-conc-gz`` is the obvious shortcut and the wrong one: it writes pairs that
    failed to align *concordantly*, which includes pairs with one mate mapped
    perfectly well to the human genome. Flag 12 (both mates unmapped) is the honest
    test, and 256 drops secondary alignments.
    """
    out1.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    sam_threads = max(1, threads // 2)
    with open(log, "ab") as errlog:
        bowtie = subprocess.Popen(  # noqa: S603
            [
                "bowtie2",
                "-p",
                str(threads),
                "-x",
                str(index),
                "-1",
                str(r1),
                "-2",
                str(r2),
                *extra,
            ],
            stdout=subprocess.PIPE,
            stderr=errlog,
        )
        view = subprocess.Popen(  # noqa: S603
            ["samtools", "view", "-@", str(sam_threads), "-u", "-f", "12", "-F", "256", "-"],
            stdin=bowtie.stdout,
            stdout=subprocess.PIPE,
            stderr=errlog,
        )
        bowtie.stdout.close()
        fastq = subprocess.Popen(  # noqa: S603
            [
                "samtools",
                "fastq",
                "-@",
                str(sam_threads),
                "-n",
                "-c",
                "6",
                "-1",
                str(out1),
                "-2",
                str(out2),
                "-0",
                "/dev/null",
                "-s",
                "/dev/null",
                "-",
            ],
            stdin=view.stdout,
            stdout=subprocess.DEVNULL,
            stderr=errlog,
        )
        view.stdout.close()
        codes = {
            "samtools fastq": fastq.wait(),
            "samtools view": view.wait(),
            "bowtie2": bowtie.wait(),
        }
    bad = {name: code for name, code in codes.items() if code != 0}
    if bad:
        raise RuntimeError(f"host filtering failed: {bad} (see {log})")


def _truth_label(
    *,
    index: Path,
    r1: Path,
    r2: Path,
    threads: int,
    min_mapq: int,
    extra: list[str],
    log: Path,
    ids_out: Path,
) -> set[str]:
    """Align to the truth reference and write the identifiers that mapped.

    A fragment counts as viral if either mate aligns at mapping quality ``min_mapq``
    or better, excluding secondary and supplementary records. Against a 30 kb
    reference, mapping quality is doing little work beyond discarding the handful of
    low-complexity repeats, but it is recorded so the threshold is visible.
    """
    log.parent.mkdir(parents=True, exist_ok=True)
    ids_out.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "ab") as errlog, open(ids_out, "wb") as out:
        bowtie = subprocess.Popen(  # noqa: S603
            [
                "bowtie2",
                "-p",
                str(threads),
                "-x",
                str(index),
                "-1",
                str(r1),
                "-2",
                str(r2),
                "--no-unal",
                *extra,
            ],
            stdout=subprocess.PIPE,
            stderr=errlog,
        )
        view = subprocess.Popen(  # noqa: S603
            [
                "samtools",
                "view",
                "-@",
                "1",
                "-F",
                "2308",  # unmapped, secondary, supplementary
                "-q",
                str(min_mapq),
                "-",
            ],
            stdin=bowtie.stdout,
            stdout=subprocess.PIPE,
            stderr=errlog,
        )
        bowtie.stdout.close()
        cut = subprocess.Popen(["cut", "-f1"], stdin=view.stdout, stdout=subprocess.PIPE)  # noqa: S603, S607
        view.stdout.close()
        sort = subprocess.Popen(  # noqa: S603, S607
            ["sort", "-u"],
            stdin=cut.stdout,
            stdout=out,
            env={"LC_ALL": "C", "PATH": "/usr/bin:/bin"},
        )
        cut.stdout.close()
        codes = {
            "sort": sort.wait(),
            "cut": cut.wait(),
            "samtools view": view.wait(),
            "bowtie2": bowtie.wait(),
        }
    bad = {name: code for name, code in codes.items() if code != 0}
    if bad:
        raise RuntimeError(f"truth labelling failed: {bad} (see {log})")
    return {reads.normalise_id(line) for line in ids_out.read_text().splitlines() if line.strip()}


def process_sample(
    spec: registry.SampleSpec,
    config: registry.Config,
    *,
    truth_index: Path | None,
    log=print,
) -> SampleResult:
    """Run every stage that applies to one sample, skipping what is already done."""
    sample_dir = config.work_dir / "samples" / spec.accession
    raw_dir = sample_dir / "raw"
    qc_dir = sample_dir / "qc"
    host_dir = sample_dir / "nonhost"
    logs = sample_dir / "logs"
    result = SampleResult(accession=spec.accession, role=spec.role)

    # ---- fetch -----------------------------------------------------------------
    runs = ena.filereport(spec.accession)
    if not runs:
        raise RuntimeError(f"ENA knows nothing about {spec.accession}")
    run = runs[0]
    mates = run.mates
    if mates is None:
        raise RuntimeError(f"{spec.accession}: no _1/_2 pair in the deposit ({run.layout})")
    result.layout = run.layout
    result.ena_metadata = {
        "study": run.study,
        "sample": run.sample,
        "library_strategy": run.library_strategy,
        "library_layout": run.library_layout,
        "instrument_model": run.instrument_model,
        "ena_read_count": run.read_count,
        "ena_base_count": run.base_count,
        "collection_date": run.collection_date,
        "deposit_layout": run.layout,
        "files": [{"name": f.name, "bytes": f.size, "md5": f.md5} for f in run.files],
    }

    stage = Stage(sample_dir, "fetch")
    started = time.monotonic()
    r1, r2, source = _fetch_pair(
        spec.accession,
        mates,
        raw_dir,
        source=config.fetch_source,
        threads=config.threads,
        tmp_dir=config.work_dir / "tmp",
        already_done=stage.done,
        log=log,
    )
    result.source = source
    result.raw_bytes = r1.stat().st_size + r2.stat().st_size
    if not stage.done:
        stage.finish({"bytes": result.raw_bytes, "source": source})
        result.timings.append(Timing("fetch", time.monotonic() - started, source))
    else:
        result.source = stage.load().get("source", source)

    # ---- qc --------------------------------------------------------------------
    qc1 = qc_dir / f"{spec.accession}_qc_1.fastq.gz"
    qc2 = qc_dir / f"{spec.accession}_qc_2.fastq.gz"
    stage = Stage(sample_dir, "qc")
    if not stage.done:
        started = time.monotonic()
        qc_dir.mkdir(parents=True, exist_ok=True)
        tools.run(
            [
                "fastp",
                "-i",
                str(r1),
                "-I",
                str(r2),
                "-o",
                str(qc1),
                "-O",
                str(qc2),
                "-w",
                str(min(config.threads, 16)),
                "--json",
                str(qc_dir / "fastp.json"),
                "--html",
                str(qc_dir / "fastp.html"),
                *config.fastp_args,
            ],
            log=logs / "fastp.log",
        )
        stage.finish()
        result.timings.append(Timing("qc", time.monotonic() - started))
    fastp_report = json.loads((qc_dir / "fastp.json").read_text())
    result.raw_pairs = fastp_report["summary"]["before_filtering"]["total_reads"] // 2
    result.qc_pairs = fastp_report["summary"]["after_filtering"]["total_reads"] // 2

    # ---- archive redaction -----------------------------------------------------
    # Counted from the raw FASTQ rather than inferred: NCBI distributes the reads it
    # redacted as full-length runs of N, and in the index sample half the run is
    # exactly that. Every per-million figure downstream needs this number to say
    # which denominator it used.
    redaction_stage = Stage(sample_dir, "redaction")
    if redaction_stage.done:
        result.redacted_pairs = redaction_stage.load()["redacted_pairs"]
    else:
        started = time.monotonic()
        result.redacted_pairs = reads.count_redacted(r1)
        redaction_stage.finish({"redacted_pairs": result.redacted_pairs})
        result.timings.append(Timing("redaction", time.monotonic() - started))
    result.non_redacted_pairs = result.raw_pairs - result.redacted_pairs

    # ---- truth labelling (index only) -----------------------------------------
    viral_ids: set[str] = set()
    ids_path = sample_dir / "truth" / "viral_fragment_ids.txt"
    if spec.is_index:
        stage = Stage(sample_dir, "truth")
        if not stage.done:
            if truth_index is None:
                raise RuntimeError("the index sample needs a truth reference index")
            started = time.monotonic()
            viral_ids = _truth_label(
                index=truth_index,
                r1=qc1,
                r2=qc2,
                threads=config.threads,
                min_mapq=config.min_truth_mapq,
                extra=list(config.bowtie2_truth_args),
                log=logs / "bowtie2-truth.log",
                ids_out=ids_path,
            )
            stage.finish({"viral_fragments": len(viral_ids)})
            result.timings.append(Timing("truth", time.monotonic() - started))
        else:
            viral_ids = {
                reads.normalise_id(line)
                for line in ids_path.read_text().splitlines()
                if line.strip()
            }
        result.viral_pairs = len(viral_ids)
        result.viral_fraction_raw = (
            result.viral_pairs / result.raw_pairs if result.raw_pairs else 0.0
        )
        result.viral_fraction_qc = result.viral_pairs / result.qc_pairs if result.qc_pairs else 0.0

    # ---- depth cap (negatives only) -------------------------------------------
    # The negatives exist to put a number on how often a method cries wolf, and that
    # number is per million reads. Twenty negatives at their native depth are 18 GiB
    # of FASTQ and 120 million pairs of alignment; capped at a recorded depth they are
    # a few GB, they finish in minutes, and every sample contributes the same weight
    # to the false-alarm estimate instead of the deepest one dominating it. The cap
    # costs sensitivity to genuinely rare background hits, which the write-up has to
    # say, and `max_pairs = 0` turns it off for a machine with room.
    filter_from_1, filter_from_2 = qc1, qc2
    if not spec.is_index and config.negatives_max_pairs:
        cap_dir = sample_dir / "capped"
        cap1 = cap_dir / f"{spec.accession}_capped_1.fastq.gz"
        cap2 = cap_dir / f"{spec.accession}_capped_2.fastq.gz"
        stage = Stage(sample_dir, "cap")
        if not stage.done:
            started = time.monotonic()
            fraction = min(1.0, config.negatives_max_pairs / max(result.qc_pairs, 1))
            kept1 = reads.write_hash_fraction(
                qc1, cap1, seed=config.negatives_seed, fraction=fraction
            )
            kept2 = reads.write_hash_fraction(
                qc2, cap2, seed=config.negatives_seed, fraction=fraction
            )
            if kept1 != kept2:
                raise RuntimeError(f"{spec.accession}: capped mates disagree, {kept1} vs {kept2}")
            stage.finish({"capped_pairs": kept1, "fraction": fraction})
            result.timings.append(Timing("cap", time.monotonic() - started))
        capped = stage.load()
        result.capped_pairs = capped["capped_pairs"]
        result.cap_fraction = capped["fraction"]
        filter_from_1, filter_from_2 = cap1, cap2

    # ---- host filtering --------------------------------------------------------
    nh1 = host_dir / f"{spec.accession}_nonhost_1.fastq.gz"
    nh2 = host_dir / f"{spec.accession}_nonhost_2.fastq.gz"
    stage = Stage(sample_dir, "host")
    if not stage.done:
        started = time.monotonic()
        _pipe_host_filter(
            index=config.host_index,
            r1=filter_from_1,
            r2=filter_from_2,
            out1=nh1,
            out2=nh2,
            threads=config.threads,
            extra=list(config.bowtie2_host_args),
            log=logs / "bowtie2-host.log",
        )
        stage.finish()
        result.timings.append(Timing("host", time.monotonic() - started))
    counts_stage = Stage(sample_dir, "counts")
    if counts_stage.done:
        counts = counts_stage.load()
        result.nonhost_pairs = counts["nonhost_pairs"]
    else:
        result.nonhost_pairs = reads.count_reads(nh1)
        counts_stage.finish({"nonhost_pairs": result.nonhost_pairs})
    filtered_input = result.capped_pairs or result.qc_pairs
    result.host_pairs = filtered_input - result.nonhost_pairs

    # ---- background and dilution (index only) ---------------------------------
    if spec.is_index:
        nonhost_ids = reads.read_ids(nh1)
        kept_viral = viral_ids & nonhost_ids
        result.viral_lost_to_host = len(viral_ids) - len(kept_viral)
        background_ids_count = len(nonhost_ids) - len(kept_viral)
        result.background_pairs = background_ids_count
        del nonhost_ids

    result.files = {
        "raw_1": _describe(r1),
        "raw_2": _describe(r2),
        "qc_1": _describe(qc1),
        "qc_2": _describe(qc2),
        "nonhost_1": _describe(nh1),
        "nonhost_2": _describe(nh2),
    }
    if spec.is_index and ids_path.exists():
        result.files["viral_fragment_ids"] = _describe(ids_path)

    # ---- pruning ---------------------------------------------------------------
    # The negatives' raw and QC'd FASTQ are 18 GiB between them and nothing
    # downstream reads either: Days 3 and 4 consume the host-filtered files. The
    # digests above are taken first, so the provenance record still describes what was
    # deleted, and the container plus the recorded parameters regenerate it.
    if config.prune_intermediates and not spec.is_index:
        freed = 0
        for path in (r1, r2, qc1, qc2):
            if path.exists():
                freed += path.stat().st_size
                path.unlink()
        for extra in (sample_dir / "capped").glob("*.fastq.gz"):
            freed += extra.stat().st_size
            extra.unlink()
        if freed:
            log(f"{spec.accession}: pruned {freed / 2**30:.2f} GiB of intermediates")
            result.pruned_bytes = freed
    return result


def _fetch_pair(
    accession: str,
    mates: tuple[ena.RunFile, ena.RunFile],
    raw_dir: Path,
    *,
    source: str,
    threads: int,
    tmp_dir: Path,
    already_done: bool,
    log=print,
) -> tuple[Path, Path, str]:
    """Get this run's FASTQ pair, from whichever source suits the run and the disk.

    ``auto`` prefers SRA's Open Data mirror, which is an order of magnitude faster,
    and falls back to ENA's ready-made FASTQ when converting the container locally
    would not fit in the free space — which is exactly the case for the index sample
    on a small instance. The choice is returned so it lands in the provenance record:
    the two sources produce different read *names*, and a reader comparing two builds
    needs to know which they have.
    """
    gz1 = raw_dir / f"{accession}_1.fastq.gz"
    gz2 = raw_dir / f"{accession}_2.fastq.gz"
    if already_done and gz1.exists() and gz2.exists():
        return gz1, gz2, source

    if source in ("auto", "sra-odp"):
        try:
            obj = sra.locate(accession)
            free = shutil.disk_usage(raw_dir.parent if raw_dir.exists() else ".").free
            whole = sra.conversion_fits(obj.size, free)
            log(
                f"{accession}: {obj.service}, {obj.size / 2**30:.2f} GiB container, "
                f"{free / 2**30:.1f} GiB free"
            )
            container = sra.download(obj, raw_dir.parent / "sra", log=log)
            if whole:
                r1, r2 = sra.fasterq_dump_paired(
                    container,
                    raw_dir,
                    accession=accession,
                    threads=threads,
                    tmp_dir=tmp_dir,
                    log_path=raw_dir.parent / "logs" / "fasterq-dump.log",
                )
                container.unlink(missing_ok=True)
                return r1, r2, "sra-odp"
            # Not enough room to convert the run whole, so convert it in spot ranges:
            # peak disk becomes one chunk instead of one run, at the cost of running
            # the slower single-threaded dumper once per chunk.
            log(f"{accession}: converting in chunks; the whole run would not fit")
            r1, r2 = sra.fastq_dump_chunked(
                container,
                raw_dir,
                accession=accession,
                threads=threads,
                log=log,
            )
            container.unlink(missing_ok=True)
            return r1, r2, "sra-odp-chunked"
        except Exception as exc:  # noqa: BLE001 - any failure here is a fallback, not fatal
            log(f"{accession}: SRA mirror unavailable ({exc}); falling back to ENA FASTQ")

    paths = [ena.download(mate, raw_dir, verify=True, log=log) for mate in mates]
    return paths[0], paths[1], "ena"


def _describe(path: Path) -> dict:
    path = Path(path)
    if not path.exists():
        return {}
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": reads.sha256_file(path),
    }


def build_series(
    config: registry.Config,
    *,
    index_accession: str,
    viral_ids: set[str],
    denominators: dict[str, int] | None = None,
    log=print,
) -> dict:
    """Write the frozen dilution series and return its record.

    Each level is the background (non-host, viral removed) plus the level's viral
    fragments, both mates. Levels are written by concatenating gzip members, so the
    background is compressed once rather than once per level.
    """
    sample_dir = config.work_dir / "samples" / index_accession
    out_dir = config.work_dir / "series"
    out_dir.mkdir(parents=True, exist_ok=True)
    nh1 = sample_dir / "nonhost" / f"{index_accession}_nonhost_1.fastq.gz"
    nh2 = sample_dir / "nonhost" / f"{index_accession}_nonhost_2.fastq.gz"
    qc1 = sample_dir / "qc" / f"{index_accession}_qc_1.fastq.gz"
    qc2 = sample_dir / "qc" / f"{index_accession}_qc_2.fastq.gz"

    parts = out_dir / "parts"
    parts.mkdir(parents=True, exist_ok=True)
    bg1, bg2 = parts / "background_1.fastq.gz", parts / "background_2.fastq.gz"
    stage = Stage(out_dir, "background")
    if not stage.done:
        log("writing the background (non-host fragments, viral removed)")
        written1 = reads.write_excluding(nh1, bg1, viral_ids)
        written2 = reads.write_excluding(nh2, bg2, viral_ids)
        if written1 != written2:
            raise RuntimeError(f"background mates disagree: {written1} vs {written2}")
        stage.finish({"background_pairs": written1})
    background_pairs = stage.load()["background_pairs"]

    plan = dilute_mod.plan_series(
        sorted(viral_ids), levels=config.levels, seed=config.dilution_seed
    )
    levels_record = []
    for level in plan:
        keep = set(level.read_ids)
        v1 = parts / f"{level.name}_viral_1.fastq.gz"
        v2 = parts / f"{level.name}_viral_2.fastq.gz"
        lvl_stage = Stage(out_dir, level.name)
        if not lvl_stage.done:
            log(f"{level.name}: extracting {level.viral_reads} viral fragments")
            got1 = reads.write_subset(qc1, v1, keep)
            got2 = reads.write_subset(qc2, v2, keep)
            if not (got1 == got2 == level.viral_reads):
                raise RuntimeError(
                    f"{level.name}: wanted {level.viral_reads} fragments, "
                    f"found {got1}/{got2} in the QC'd FASTQs"
                )
            out1 = out_dir / f"{index_accession}_{level.name}_1.fastq.gz"
            out2 = out_dir / f"{index_accession}_{level.name}_2.fastq.gz"
            reads.concatenate([bg1, v1], out1)
            reads.concatenate([bg2, v2], out2)
            lvl_stage.finish({"viral_reads": level.viral_reads})
        out1 = out_dir / f"{index_accession}_{level.name}_1.fastq.gz"
        out2 = out_dir / f"{index_accession}_{level.name}_2.fastq.gz"
        levels_record.append(
            {
                "name": level.name,
                "requested_viral_fragments": level.requested,
                "viral_fragments": level.viral_reads,
                "truncated": level.truncated,
                "seed": level.seed,
                "selection_sha256": level.selection_sha256,
                "background_fragments": background_pairs,
                "total_fragments": background_pairs + level.viral_reads,
                "viral_fragments_per_million": dilute_mod.reads_per_million(
                    level.viral_reads, background_pairs
                ),
                # Three denominators, because three are defensible and they differ by
                # a factor of two. The headline figure is the first: it is the only
                # one the detection arms see, so it is the one that converts into a
                # required sequencing depth. The raw figure is the one comparable in
                # kind to Wu et al.'s 0.19% — though not directly, since theirs was
                # computed before NCBI redacted half the run.
                "per_million": {
                    key: level.viral_reads / value * 1_000_000
                    for key, value in (denominators or {}).items()
                    if value
                },
                "files": {"read_1": _describe(out1), "read_2": _describe(out2)},
            }
        )
    return {
        "index_run": index_accession,
        "seed": config.dilution_seed,
        "selection_rule": (
            "fragments ranked by sha256(f'{seed}:{fragment_id}'), ascending; "
            "levels are nested prefixes of that ranking"
        ),
        "nested": True,
        "background_fragments": background_pairs,
        "viral_fragments_available": len(viral_ids),
        "denominators": denominators or {},
        "headline_denominator": "in_file",
        "levels": levels_record,
    }


def run_day2(
    config: registry.Config,
    *,
    only: list[str] | None = None,
    skip_series: bool = False,
    log=print,
) -> dict:
    """Everything: fetch, QC, label, filter, freeze. Idempotent and resumable."""
    started_wall = time.monotonic()
    started_at = dt.datetime.now(dt.UTC)
    config.work_dir.mkdir(parents=True, exist_ok=True)

    tool_versions = tools.probe_all(
        {
            "curl": config.tool_pins.get("curl"),
            "fastp": config.tool_pins.get("fastp"),
            "bowtie2": config.tool_pins.get("bowtie2"),
            "samtools": config.tool_pins.get("samtools"),
        }
    )
    drift = [t for t in tool_versions if not t.matches_pin]
    for t in tool_versions:
        log(f"{t.name} {t.version}" + ("" if t.matches_pin else f" (config pins {t.pinned})"))

    specs = registry.samples_from_config(config)
    if only:
        wanted = set(only)
        specs = [s for s in specs if s.accession in wanted]

    truth_index = None
    if any(s.is_index for s in specs):
        ref_dir = config.work_dir / "ref" / "truth"
        fasta = fetch_reference(config.truth_reference, ref_dir / f"{config.truth_reference}.fna")
        truth_index = build_index(
            fasta,
            ref_dir / config.truth_reference,
            threads=config.threads,
            log=config.work_dir / "logs" / "bowtie2-build.log",
        )

    results: list[SampleResult] = []
    for spec in specs:
        log(f"=== {spec.accession} ({spec.role}) ===")
        results.append(process_sample(spec, config, truth_index=truth_index, log=log))

    series = None
    index_result = next((r for r in results if r.role == "index"), None)
    if index_result and not skip_series:
        ids_path = (
            config.work_dir
            / "samples"
            / index_result.accession
            / "truth"
            / "viral_fragment_ids.txt"
        )
        viral_ids = {
            reads.normalise_id(line) for line in ids_path.read_text().splitlines() if line.strip()
        }
        series = build_series(
            config,
            index_accession=index_result.accession,
            viral_ids=viral_ids,
            denominators={
                # "in_file" is filled per level by build_series from the level's own
                # total; these are the sample-wide alternatives.
                "raw_spots": index_result.raw_pairs,
                "non_redacted_spots": index_result.non_redacted_pairs,
                "post_qc_pairs": index_result.qc_pairs,
            },
            log=log,
        )

    record = {
        "tool": {
            "name": "hindsight samples",
            "version": __version__,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
        "generated_at": started_at.isoformat(timespec="seconds"),
        "command": list(sys.argv),
        "wall_clock_seconds": round(time.monotonic() - started_wall, 1),
        "machine": machine_record(config.work_dir),
        "tools": [asdict(t) for t in tool_versions],
        "tool_pin_drift": [t.name for t in drift],
        "parameters": {
            "fastp_args": list(config.fastp_args),
            "bowtie2_host_args": list(config.bowtie2_host_args),
            "bowtie2_truth_args": list(config.bowtie2_truth_args),
            "min_truth_mapq": config.min_truth_mapq,
            "threads": config.threads,
            "host_index": str(config.host_index),
            "truth_reference": config.truth_reference,
            "dilution_levels": list(config.levels),
            "dilution_seed": config.dilution_seed,
            "negatives_seed": config.negatives_seed,
            "negatives_count": config.negatives_count,
            "negatives_max_pairs": config.negatives_max_pairs,
            "prune_intermediates": config.prune_intermediates,
        },
        "samples": [_sample_record(r) for r in results],
        "dilution_series": series,
    }
    out = config.work_dir / "day2-provenance.json"
    out.write_text(json.dumps(record, indent=2, default=str) + "\n")
    log(f"wrote {out}")
    return record


def _sample_record(result: SampleResult) -> dict:
    data = asdict(result)
    data["timings"] = [asdict(t) for t in result.timings]
    return data


def machine_record(work_dir: Path) -> dict:
    """What this ran on, because "$40 of spot" is a claim that needs a machine.

    Captured from the host rather than configured, so a reproduction on different
    hardware records the difference instead of inheriting the original's numbers.
    """
    info: dict = {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_count": None,
        "memory_gib": None,
        "instance_type": None,
    }
    try:
        import os

        info["cpu_count"] = os.cpu_count()
    except Exception:  # pragma: no cover - defensive
        pass
    try:
        meminfo = Path("/proc/meminfo").read_text()
        for line in meminfo.splitlines():
            if line.startswith("MemTotal:"):
                info["memory_gib"] = round(int(line.split()[1]) / 1024 / 1024, 1)
                break
    except Exception:
        pass
    try:
        model = Path("/proc/cpuinfo").read_text()
        for line in model.splitlines():
            if line.lower().startswith("model name"):
                info["cpu_model"] = line.split(":", 1)[1].strip()
                break
    except Exception:
        pass
    # On EC2 the instance type is worth having verbatim; IMDSv2 needs a token and is
    # absent everywhere else, so a failure here is normal and silent.
    try:
        token = requests.put(
            "http://169.254.169.254/latest/api/token",
            headers={"X-aws-ec2-metadata-token-ttl-seconds": "60"},
            timeout=1,
        )
        if token.ok:
            meta = requests.get(
                "http://169.254.169.254/latest/meta-data/instance-type",
                headers={"X-aws-ec2-metadata-token": token.text},
                timeout=1,
            )
            if meta.ok:
                info["instance_type"] = meta.text.strip()
    except Exception:
        pass
    try:
        usage = shutil.disk_usage(work_dir)
        info["work_dir_disk"] = {
            "total_gib": round(usage.total / 2**30, 1),
            "used_gib": round(usage.used / 2**30, 1),
            "free_gib": round(usage.free / 2**30, 1),
        }
    except Exception:
        pass
    return info
