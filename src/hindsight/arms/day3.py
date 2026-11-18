"""Day 3: run Arm A over the whole benchmark and write down what it found.

The arm itself is in :mod:`hindsight.arms.kraken`. This is the driver — what gets
classified, in what order, and what is recorded about each run.

Three decisions are encoded here rather than in a config file, because each one is a
claim about the benchmark and belongs where it can be read:

**The negatives are classified at native depth.** Day 2 capped them at 6 million pairs
because 17 GiB of host-filtered FASTQ does not fit on the machine this runs on. The
cap degrades exactly the quantity Day 5 needs — a false-alarm rate is about rare
events, and a cap removes rare events — so Day 3 removes it, and pays for it by
coupling the two days: each negative is fetched, filtered, classified and then
*deleted*, one at a time. What persists is the classifier's report and its classified
reads, which are megabytes. Peak disk is therefore one negative's worth, not twenty.

**A dilution level is composed, not materialised.** See
:func:`hindsight.arms.kraken.compose`. The composition is checked against a directly
classified level before any of it is believed.

**The index sample is classified undiluted as well.** The dilution series tops out at
50,000 of its 62,117 viral fragments, so the series alone never answers "what does the
arm do on the sample as it was". That run is the benchmark's upper bound, and the
number the write-up will quote first.
"""

from __future__ import annotations

import dataclasses
import json
import subprocess
import time
import tomllib
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ..samples import pipeline as samples_pipeline
from ..samples import reads, registry
from . import kraken

DEFAULT_TOOL_PINS = {"kraken2": "2.1.3"}


@dataclass(frozen=True)
class ArmConfig:
    samples_config: Path
    source_db: Path
    kraken_db: Path
    cutoff: str
    threads: int
    confidence: float
    minimum_hit_groups: int
    mask_low_complexity: bool
    negatives_max_pairs: int
    verify_level: str
    verify_composition: bool
    prune_after_classify: bool
    work_dir: Path
    out_dir: Path
    provenance: Path
    tool_pins: dict[str, str | None]


def load_config(path: Path) -> ArmConfig:
    data = tomllib.loads(Path(path).read_text())
    arm = data.get("arm", {})
    classify = data.get("classify", {})
    database = data.get("database", {})
    negatives = data.get("negatives", {})
    work_dir = Path(arm.get("work_dir", "work"))
    return ArmConfig(
        samples_config=Path(arm.get("samples_config", "config/day2.toml")),
        source_db=Path(database.get("source_db", "db/nov2019")),
        kraken_db=Path(database.get("kraken_db", work_dir / "kraken2" / "nov2019")),
        cutoff=str(database.get("cutoff", "2019-12-01")),
        threads=int(classify.get("threads", 4)),
        confidence=float(classify.get("confidence", 0.0)),
        minimum_hit_groups=int(classify.get("minimum_hit_groups", 2)),
        mask_low_complexity=bool(database.get("mask_low_complexity", True)),
        negatives_max_pairs=int(negatives.get("max_pairs", 0)),
        verify_level=str(classify.get("verify_level", "viral000010")),
        verify_composition=bool(classify.get("verify_composition", True)),
        prune_after_classify=bool(negatives.get("prune_after_classify", True)),
        work_dir=work_dir,
        out_dir=Path(arm.get("out_dir", work_dir / "arm-a")),
        provenance=Path(arm.get("provenance", "results/day3-provenance.json")),
        tool_pins=dict(data.get("tools", DEFAULT_TOOL_PINS)),
    )


# ------------------------------------------------------------------------- database


def prepare_database(arm: ArmConfig, *, log=print) -> tuple[dict, kraken.CladeSets]:
    record = kraken.build_database(
        arm.source_db,
        arm.kraken_db,
        threads=arm.threads,
        mask=arm.mask_low_complexity,
        log=log,
    )
    taxonomy = kraken.read_taxonomy(arm.source_db / "taxonomy")
    clades = kraken.clade_sets(taxonomy)
    record["clades"] = {
        name: {
            "root_taxid": clades.roots[name],
            "name_in_2019_tree": taxonomy.names.get(clades.roots[name]),
            "taxa_in_clade": len(members),
        }
        for name, members in clades.members.items()
    }
    if clades.missing:
        # A clade the dated tree does not contain is a real finding, not a config
        # error: it means the group had not been named yet, and a report written
        # against that tree can never place a read there.
        log("clades absent from the dated taxonomy: " + ", ".join(clades.missing))
    record["clades_absent_from_dated_taxonomy"] = clades.missing
    record["headline_clade"] = kraken.HEADLINE_CLADE
    source_provenance = arm.source_db / "provenance.json"
    if source_provenance.exists():
        source = json.loads(source_provenance.read_text())
        record["source_release"] = {
            "cutoff": source.get("cutoff"),
            "release": (source.get("release") or {}),
            "first_excluded_release": source.get("first_excluded_release"),
            "counts": source.get("counts"),
            "taxonomy": source.get("taxonomy"),
        }
    return record, clades


# --------------------------------------------------------------------------- series


def classify_series(
    arm: ArmConfig,
    day2: dict,
    clades: kraken.CladeSets,
    *,
    log=print,
) -> dict:
    """Classify the background once and each level's viral part once, then compose."""
    series = day2.get("dilution_series")
    if not series:
        raise kraken.KrakenError("the Day 2 provenance record holds no dilution series")
    out_dir = arm.out_dir / "series"

    background = _classify_parts(arm, series, "background", out_dir, log=log)
    bg_record = background.record(clades)
    # Every background fragment is non-viral by construction — the series' background
    # is the non-host set with the ground-truth viral fragments removed — so a clade
    # hit here is a false positive with no labelling step needed.
    bg_record["role"] = "background (viral fragments removed; any clade hit is a false positive)"

    levels = []
    for level in series["levels"]:
        name = level["name"]
        viral = _classify_parts(arm, series, name, out_dir, log=log, viral_only=True)
        viral_record = viral.record(clades)
        composed = kraken.compose([background, viral])
        clade_counts = {
            clade: clades.count({int(k): v for k, v in composed["counts"].items()}, clade)
            for clade in clades.members
        }
        # Every fragment in a viral part is ground-truth viral by construction, and
        # every fragment in the background is not, so detection and false positives
        # need no labelling step here — they are the two files' own clade counts.
        detected_by_clade = {
            name: clades.count(viral.counts, name) for name in clades.members
        }
        false_positive_by_clade = {
            name: clades.count(background.counts, name) for name in clades.members
        }
        detected = detected_by_clade[kraken.HEADLINE_CLADE]
        false_positive = false_positive_by_clade[kraken.HEADLINE_CLADE]
        levels.append(
            {
                "name": name,
                "viral_fragments": level["viral_fragments"],
                "total_fragments": level["total_fragments"],
                "per_million": level.get("per_million", {}),
                "viral_fragments_per_million": level.get("viral_fragments_per_million"),
                "viral_part": viral_record,
                "composed": {
                    "pairs": composed["pairs"],
                    "classified": composed["classified"],
                    "unclassified": composed["unclassified"],
                    "clades": clade_counts,
                },
                "detected_viral_reads": detected,
                "false_positive_reads": false_positive,
                "detected_by_clade": detected_by_clade,
                "false_positive_by_clade": false_positive_by_clade,
                "sensitivity": (
                    detected / level["viral_fragments"] if level["viral_fragments"] else None
                ),
                "detected_per_million_post_qc": _per_million(
                    detected, series.get("denominators", {}).get("post_qc_pairs")
                ),
            }
        )
        log(
            f"{name}: {detected:,} of {level['viral_fragments']:,} viral fragments placed in "
            f"{kraken.HEADLINE_CLADE}, {false_positive:,} background false positives"
        )

    record = {"background": bg_record, "levels": levels}
    if arm.verify_composition:
        record["composition_check"] = _verify(arm, series, clades, background, log=log)
    return record


def _per_million(count: int, denominator: int | None) -> float | None:
    if not denominator:
        return None
    return count * 1_000_000 / denominator


def _classify_parts(
    arm: ArmConfig,
    series: dict,
    level_name: str,
    out_dir: Path,
    *,
    log=print,
    viral_only: bool = False,
) -> kraken.Classification:
    """Classify one part pair of the series: the background, or one level's viral reads.

    The parts are the unit the series is stored in, so nothing is concatenated or
    rewritten here — the files on disk are the input.
    """
    if level_name == "background":
        r1 = arm.work_dir / "series" / "parts" / "background_1.fastq.gz"
        r2 = arm.work_dir / "series" / "parts" / "background_2.fastq.gz"
        name = "background"
    else:
        level = next((lv for lv in series["levels"] if lv["name"] == level_name), None)
        if level is None:
            raise kraken.KrakenError(f"no level {level_name} in the series")
        viral_parts = {
            mate: [p for p in level["parts"][mate] if "background" not in Path(p["path"]).name]
            for mate in ("read_1", "read_2")
        }
        if not viral_only:
            raise kraken.KrakenError("levels are classified as background + viral parts")
        if len(viral_parts["read_1"]) != 1 or len(viral_parts["read_2"]) != 1:
            raise kraken.KrakenError(
                f"{level_name}: expected one viral part per mate, got "
                f"{len(viral_parts['read_1'])}/{len(viral_parts['read_2'])}"
            )
        r1 = Path(viral_parts["read_1"][0]["path"])
        r2 = Path(viral_parts["read_2"][0]["path"])
        name = f"{level_name}_viral"
    for path in (r1, r2):
        if not path.exists():
            raise kraken.KrakenError(f"{path} is missing; re-run the Day 2 pipeline")
    return kraken.classify(
        arm.kraken_db,
        r1,
        r2,
        out_dir,
        name=name,
        threads=arm.threads,
        confidence=arm.confidence,
        minimum_hit_groups=arm.minimum_hit_groups,
        log=log,
    )


def _verify(
    arm: ArmConfig,
    series: dict,
    clades: kraken.CladeSets,
    background: kraken.Classification,
    *,
    log=print,
) -> dict:
    """Materialise one level, classify it whole, and check the composition matches.

    This is the expensive check in the day — a full 10.9 M-pair classification whose
    only purpose is to test an assumption — and it is worth it. Everything the series
    reports is composed, so if Kraken2's per-read independence ever stopped holding,
    every dilution number in the paper would be wrong in a way nothing else here would
    catch.
    """
    level_name = arm.verify_level
    log(f"composition check: materialising {level_name} and classifying it whole")
    materialised = arm.work_dir / "series" / "materialised"
    r1, r2 = samples_pipeline.materialise_level(series, level_name, materialised, verify=True)
    try:
        direct = kraken.classify(
            arm.kraken_db,
            r1,
            r2,
            arm.out_dir / "verify",
            name=f"{level_name}_direct",
            threads=arm.threads,
            confidence=arm.confidence,
            minimum_hit_groups=arm.minimum_hit_groups,
            log=log,
        )
        viral = _classify_parts(arm, series, level_name, arm.out_dir / "series", viral_only=True)
        composed = kraken.compose([background, viral])
        result = kraken.verify_composition(composed, direct)
        result["clades_direct"] = {
            name: clades.count(direct.counts, name) for name in clades.members
        }
        log(
            "composition check: "
            + (
                "identical"
                if result["identical"]
                else f"MISMATCH on {result['taxa_disagreeing']} taxa"
            )
        )
        return result
    finally:
        # 2 GB of materialised FASTQ whose only purpose was this check.
        for path in (r1, r2):
            Path(path).unlink(missing_ok=True)


# ---------------------------------------------------------------------- index sample


def classify_index(
    arm: ArmConfig,
    samples_config: registry.Config,
    day2: dict,
    clades: kraken.CladeSets,
    *,
    log=print,
) -> dict:
    accession = samples_config.index_run
    sample_dir = arm.work_dir / "samples" / accession
    r1 = sample_dir / "nonhost" / f"{accession}_nonhost_1.fastq.gz"
    r2 = sample_dir / "nonhost" / f"{accession}_nonhost_2.fastq.gz"
    for path in (r1, r2):
        if not path.exists():
            raise kraken.KrakenError(f"{path} is missing; re-run the Day 2 pipeline")
    result = kraken.classify(
        arm.kraken_db,
        r1,
        r2,
        arm.out_dir / "index",
        name=f"{accession}_nonhost",
        threads=arm.threads,
        confidence=arm.confidence,
        minimum_hit_groups=arm.minimum_hit_groups,
        log=log,
    )
    record = result.record(clades)
    truth_ids = kraken.read_truth_ids(sample_dir / "truth" / "viral_fragment_ids.txt")
    sample = next(
        (s for s in day2.get("samples", []) if s["accession"] == accession),
        {},
    )
    # The label is applied to the classifier's output, never to its input.
    nonhost_viral = sample.get("viral_pairs")
    scored = kraken.score_against_truth(result.classified_output, truth_ids, clades)
    scored["truth_fragments_in_nonhost"] = nonhost_viral
    if nonhost_viral:
        scored["sensitivity"] = scored["detected"] / nonhost_viral
    record["truth"] = scored
    denominators = (day2.get("dilution_series") or {}).get("denominators", {})
    record["detected_per_million"] = {
        key: _per_million(scored["detected"], value) for key, value in denominators.items()
    }
    log(
        f"{accession}: {scored['detected']:,} of {nonhost_viral:,} ground-truth viral "
        f"fragments placed in {kraken.HEADLINE_CLADE} "
        f"({100 * (scored['sensitivity'] or 0):.1f}%), "
        f"{scored['false_positive']:,} false positives"
    )
    return record


# ----------------------------------------------------------- the background's hits


def diagnose_background_hits(
    arm: ArmConfig,
    samples_config: registry.Config,
    clades: kraken.CladeSets,
    *,
    taxid: int = 694009,
    log=print,
) -> dict:
    """Find out what the index sample's own background is actually made of.

    The series' background is the non-host set with every ground-truth viral fragment
    removed, so a read the arm places in the sarbecoviruses should not be there. Two
    explanations fit, and they have opposite consequences: either the classifier is
    wrong about those reads, in which case the arm has a false-positive rate on its own
    sample, or the *label* is wrong about them, in which case the sample has more viral
    reads than the benchmark credits it with.

    So this asks the question the label was asked, with the threshold removed: take
    those reads, align them to the same modern reference, and report the alignment rate
    and the mapping qualities. A read that aligns but at mapping quality below the
    label's cut-off was never a classifier error; it was a fragment the label declined
    to count.

    This is a diagnostic, not part of the arm. It is the one place in Day 3 that
    touches the post-2019 reference, and it touches it after the arm has finished and
    only to interpret the arm's output.
    """
    out_dir = arm.out_dir / "diagnostics"
    out_dir.mkdir(parents=True, exist_ok=True)
    classified = arm.out_dir / "series" / "background.classified.tsv.gz"
    if not classified.exists():
        raise kraken.KrakenError(f"{classified} is missing; classify the series first")

    wanted = {
        fragment for fragment, assigned in kraken.iter_output(classified) if assigned == taxid
    }
    taxonomy = kraken.read_taxonomy(arm.source_db / "taxonomy")
    record: dict = {
        "taxid": taxid,
        "taxon": taxonomy.names.get(taxid),
        "clades": clades.which(taxid),
        "fragments": len(wanted),
        "truth_reference": samples_config.truth_reference,
        "label_min_mapq": samples_config.min_truth_mapq,
    }
    if not wanted:
        record["conclusion"] = "no background fragments were assigned to this taxon"
        return record
    (out_dir / f"background_{taxid}_ids.txt").write_text("\n".join(sorted(wanted)) + "\n")

    parts = arm.work_dir / "series" / "parts"
    extracted = []
    for mate in (1, 2):
        dest = out_dir / f"background_{taxid}_{mate}.fastq.gz"
        count = reads.write_subset(parts / f"background_{mate}.fastq.gz", dest, wanted)
        extracted.append(count)
    record["fragments_extracted"] = extracted

    index = arm.work_dir / "ref" / "truth" / samples_config.truth_reference
    sam = out_dir / f"background_{taxid}.sam"
    log_path = out_dir / f"background_{taxid}.bowtie2.log"
    with open(sam, "wb") as out, open(log_path, "wb") as errlog:
        bowtie = subprocess.Popen(  # noqa: S603
            [
                "bowtie2",
                "--local",
                "-p",
                str(arm.threads),
                "-x",
                str(index),
                "-1",
                str(out_dir / f"background_{taxid}_1.fastq.gz"),
                "-2",
                str(out_dir / f"background_{taxid}_2.fastq.gz"),
            ],
            stdout=subprocess.PIPE,
            stderr=errlog,
        )
        view = subprocess.Popen(  # noqa: S603
            ["samtools", "view", "-F", "2308", "-"],
            stdin=bowtie.stdout,
            stdout=out,
            stderr=errlog,
        )
        bowtie.stdout.close()
        codes = {"samtools": view.wait(), "bowtie2": bowtie.wait()}
    if any(codes.values()):
        raise kraken.KrakenError(f"diagnostic alignment failed: {codes} (see {log_path})")

    mapqs: Counter[int] = Counter()
    aligned_fragments: set[str] = set()
    with sam.open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if len(fields) < 5 or fields[2] == "*":
                continue
            mapqs[int(fields[4])] += 1
            aligned_fragments.add(reads.normalise_id(fields[0]))
    record["alignment_summary"] = log_path.read_text(errors="replace").strip().splitlines()
    record["fragments_aligning_to_the_modern_reference"] = len(aligned_fragments)
    record["mapq_distribution"] = dict(sorted(mapqs.items()))
    record["max_mapq"] = max(mapqs) if mapqs else None
    record["all_below_the_label_threshold"] = (
        bool(mapqs) and max(mapqs) < samples_config.min_truth_mapq
    )
    if record["all_below_the_label_threshold"]:
        record["conclusion"] = (
            f"{len(aligned_fragments)} of {len(wanted)} fragments align to "
            f"{samples_config.truth_reference}, every mate below the label's "
            f"MAPQ >= {samples_config.min_truth_mapq} cut-off. These are viral reads the "
            "ground-truth label declined to count, not classifier false positives, so "
            "the sample's viral content is slightly larger than the benchmark credits "
            "and the arm's apparent false-positive rate on its own background is not one."
        )
    else:
        record["conclusion"] = (
            f"{len(aligned_fragments)} of {len(wanted)} fragments align, with mapping "
            f"quality up to {record['max_mapq']} — at or above the label's cut-off, so "
            "the label and the classifier disagree for another reason. Investigate "
            "before Day 5 uses this background for anything."
        )
    log(record["conclusion"])
    for path in (sam,):
        path.unlink(missing_ok=True)
    return record


# ------------------------------------------------------------------------- negatives


def classify_negatives(
    arm: ArmConfig,
    samples_config: registry.Config,
    clades: kraken.CladeSets,
    *,
    only: Sequence[str] | None = None,
    log=print,
) -> list[dict]:
    """Fetch, filter, classify and discard each negative in turn, at native depth.

    One sample is on disk at a time. The order is the pinned order, so an interrupted
    cohort resumes where it stopped, and a sample already classified is skipped
    without re-fetching it.
    """
    specs = [s for s in registry.samples_from_config(samples_config) if not s.is_index]
    if only:
        wanted = set(only)
        specs = [s for s in specs if s.accession in wanted]
    deep_config = dataclasses.replace(samples_config, negatives_max_pairs=arm.negatives_max_pairs)
    records: list[dict] = []
    for position, spec in enumerate(specs, start=1):
        out_dir = arm.out_dir / "negatives"
        marker = out_dir / f".{spec.accession}_nonhost.done"
        log(f"[{position}/{len(specs)}] {spec.accession}")
        if marker.exists():
            log(f"{spec.accession}: already classified; not re-fetching")
        else:
            started = time.monotonic()
            sample = samples_pipeline.process_sample(spec, deep_config, truth_index=None, log=log)
            prepare_seconds = time.monotonic() - started
            _write_sample_note(out_dir, spec.accession, sample, prepare_seconds)
        record = _classify_one_negative(arm, spec.accession, clades, log=log)
        records.append(record)
        if arm.prune_after_classify:
            record["pruned_bytes"] = _prune_sample(arm.work_dir / "samples" / spec.accession)
            log(f"{spec.accession}: pruned {record['pruned_bytes'] / 2**30:.2f} GiB")
    return records


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


def _classify_one_negative(
    arm: ArmConfig, accession: str, clades: kraken.CladeSets, *, log=print
) -> dict:
    sample_dir = arm.work_dir / "samples" / accession
    r1 = sample_dir / "nonhost" / f"{accession}_nonhost_1.fastq.gz"
    r2 = sample_dir / "nonhost" / f"{accession}_nonhost_2.fastq.gz"
    out_dir = arm.out_dir / "negatives"
    result = kraken.classify(
        arm.kraken_db,
        r1,
        r2,
        out_dir,
        name=f"{accession}_nonhost",
        threads=arm.threads,
        confidence=arm.confidence,
        minimum_hit_groups=arm.minimum_hit_groups,
        log=log,
    )
    record = result.record(clades)
    record["accession"] = accession
    record["role"] = "negative"
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
    hits = clades.count(result.counts, kraken.HEADLINE_CLADE)
    record["false_alarm_reads"] = hits
    record["false_alarm_per_million_nonhost"] = _per_million(hits, result.pairs)
    log(f"{accession}: {hits:,} reads in {kraken.HEADLINE_CLADE} of {result.pairs:,} non-host")
    return record


def _prune_sample(sample_dir: Path) -> int:
    """Delete every read file of a finished negative, keeping markers and reports.

    This is what makes native depth fit: the classifier's output is kept, the reads it
    read are not. The sample's digests were recorded before deletion and the pinned
    parameters regenerate it.
    """
    freed = 0
    patterns = ("*.fastq.gz", "*.fastq", "*.sra", "*.sra.*")
    # "sra" holds the downloaded container, which is 1-3 GiB a sample and the single
    # largest thing left behind. Day 2 could leave it: it pruned FASTQ and the cap
    # meant the container was the only copy worth keeping. At native depth, twenty
    # containers is the disk this day does not have.
    for sub in ("raw", "qc", "capped", "nonhost", "tmp", "sra"):
        directory = sample_dir / sub
        if not directory.exists():
            continue
        for pattern in patterns:
            for path in directory.glob(pattern):
                if path.is_file():
                    freed += path.stat().st_size
                    path.unlink()
    return freed


# ------------------------------------------------------------------------ the whole day


def run_day3(
    arm: ArmConfig,
    *,
    stages: Sequence[str] = ("database", "index", "series", "negatives"),
    only_negatives: Sequence[str] | None = None,
    log=print,
) -> dict:
    started = time.monotonic()
    samples_config = registry.load_config(arm.samples_config)
    day2_path = arm.work_dir / "day2-provenance.json"
    if not day2_path.exists():
        day2_path = Path("results/day2-provenance.json")
    day2 = json.loads(day2_path.read_text())

    record: dict = {
        "tool": {"name": "hindsight arm-a", "version": _version()},
        "generated_at": kraken.utc_now(),
        "stages": list(stages),
        "machine": kraken.machine_record(arm.work_dir),
        "tools": kraken.tool_versions(arm.tool_pins),
        "parameters": {
            "cutoff": arm.cutoff,
            "source_db": str(arm.source_db),
            "kraken_db": str(arm.kraken_db),
            "mask_low_complexity": arm.mask_low_complexity,
            "threads": arm.threads,
            "confidence": arm.confidence,
            "minimum_hit_groups": arm.minimum_hit_groups,
            "negatives_max_pairs": arm.negatives_max_pairs,
            "negatives_depth": (
                "native" if arm.negatives_max_pairs == 0 else f"{arm.negatives_max_pairs} pairs"
            ),
            "headline_clade": kraken.HEADLINE_CLADE,
            "day2_provenance": str(day2_path),
        },
    }

    # The database is prepared whatever stages were asked for: every other stage needs
    # the clade sets that come with it, and a built database is reused rather than
    # rebuilt, so there is nothing to skip.
    database_record, clades = prepare_database(arm, log=log)
    record["database"] = database_record

    if "index" in stages:
        record["index_sample"] = classify_index(arm, samples_config, day2, clades, log=log)
    if "series" in stages:
        record["series"] = classify_series(arm, day2, clades, log=log)
        # What the background's own clade hits are, which decides whether they are the
        # arm's false positives or the label's false negatives.
        record["background_diagnosis"] = diagnose_background_hits(
            arm, samples_config, clades, log=log
        )
    if "negatives" in stages:
        record["negatives"] = classify_negatives(
            arm, samples_config, clades, only=only_negatives, log=log
        )
        record["false_alarms"] = summarise_false_alarms(record["negatives"])

    record["wall_clock_seconds"] = round(time.monotonic() - started, 1)
    merged = merge(arm.provenance, record)
    return merged


def summarise_false_alarms(negatives: list[dict]) -> dict:
    """The number Day 5's threshold is drawn from, with its denominators stated."""
    if not negatives:
        return {}
    hits = [n.get("false_alarm_reads", 0) for n in negatives]
    nonhost = [n.get("pairs", 0) for n in negatives]
    with_hits = [n for n, h in zip(negatives, hits, strict=False) if h]
    return {
        "samples": len(negatives),
        "non_host_pairs_total": sum(nonhost),
        "clade": kraken.HEADLINE_CLADE,
        "samples_with_at_least_one_hit": len(with_hits),
        "reads_total": sum(hits),
        "reads_max_in_one_sample": max(hits) if hits else 0,
        "reads_median_per_sample": sorted(hits)[len(hits) // 2] if hits else 0,
        "reads_per_million_nonhost": (
            sum(hits) * 1_000_000 / sum(nonhost) if sum(nonhost) else None
        ),
        "samples_with_hits": {
            n["accession"]: n.get("false_alarm_reads", 0)
            for n in sorted(with_hits, key=lambda n: -n.get("false_alarm_reads", 0))
        },
    }


def merge(path: Path, fresh: dict) -> dict:
    """Merge a run into the record on disk, the way Day 2 does.

    A cohort of twenty negatives is not built in one invocation — a crash, a resume, a
    rerun of one sample — so the record accumulates and keeps a stamp per invocation
    rather than being overwritten by whichever run finished last.
    """
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
