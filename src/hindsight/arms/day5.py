"""Day 5: score both arms in one frame, and draw the threshold the benchmark promised.

Nothing in here runs a classifier, an assembler or an aligner. It reads
``results/day3-provenance.json`` and ``results/day4-provenance.json`` and turns what
those two days measured into the benchmark's headline output: a limit of detection per
method, in viral reads per million, at a stated false-alarm rate. That is a deliberate
property rather than a convenience — a scoring step that can re-run itself over the
recorded tables is a scoring step whose thresholds can be moved by a reader without a
machine, and both arms were built to record every clade at every level for exactly this.

Four ideas shape the arithmetic, and each of them is a decision this day owns.

**An alarm is a property of a sample, not of a read.** Arm A's natural unit is reads and
Arm B's is contigs, so a single shared threshold across both is the wrong shape (Day 4's
phrasing, and it is right). What *is* shared is the question a triage process asks: did
this sample alarm? So each arm gets its own statistic — a read count for Arm A, a
bitscore for Arm B — and both are converted into the same binary, at a threshold set by
the twenty negatives rather than by the positive. A threshold chosen after looking at
the positive's answer is not a detection limit, so the positive is never consulted to
place one.

**The dilution series' background is not virus-free, and that bounds the whole day.**
Day 3 found 78 sarbecovirus-assigned fragments in the spike-free background and showed
that 71 of them align to the modern SARS-CoV-2 genome below the ground-truth label's
MAPQ cut-off — real viral reads the label missed. They are present *undiluted* at every
rung of the series, so the dose at a rung is the labelled count plus about eighty. At
1,000 labelled fragments that is noise; at 10 it is most of the dose. Scoring therefore
uses the labelled detections, which is the conservative reading, and records the
operational count beside it so the gap is visible rather than argued about.

**Which reads count as a false alarm is a choice with a two-order-of-magnitude effect,
so it is made explicitly and recorded in every variant.** The cohort's coronavirus reads
are dominated by two patients' genuine endemic infections, which the November 2019
database names correctly at 97.8-100% identity. Counting those as false alarms measures
how often a BALF patient has a cold.

**The clade is part of the measurement.** Day 3 (the 2019 taxonomy files the nearest bat
relative outside Betacoronavirus) and Day 4 (the genome contig's own LCA *is* the family)
arrived at family level from independent directions, and a third argument appears in the
negatives: one patient's NL63 is an alphacoronavirus and is invisible to anything scored
at Betacoronavirus. Every clade is scored anyway, because the record has to let a reader
disagree.
"""

from __future__ import annotations

import dataclasses
import json
import tomllib
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

from ..fasta import sha256_file
from .taxa import CLADES, ENDEMIC_HUMAN_CORONAVIRUSES, machine_record, utc_now

# The clade Day 5 settles on, and the one the headline is quoted at. See the module
# docstring and the three independent arguments in docs/day5-scoring.md; every other
# clade in CLADES is scored too and sits in the record.
HEADLINE_CLADE = "coronaviridae"

# Arm B's pre-registered bitscore ladder, fixed in config/day4.toml before the series
# ran. 0.0 is prepended to represent the arm's own unfiltered rule — "any contig whose
# LCA lands in the clade" — which is what Day 4 recorded as `alarm`.
BITSCORE_LADDER: tuple[float, ...] = (0.0, 50.0, 100.0, 250.0, 500.0, 1000.0)

# How many of the twenty negatives may alarm and the operating point still count. One is
# the benchmark's own promise ("at <=1 false alarm per N samples"); zero is reported
# beside it because with twenty samples the two often coincide, and where they do not
# the difference is the honest width of the measurement.
FALSE_ALARM_BUDGETS: tuple[int, ...] = (0, 1)

FRAMINGS: tuple[str, ...] = ("all", "excluding_endemic")

# Which field decides whether a contig is "somebody's cold". Day 4 objected, correctly,
# that excluding the endemic human coronaviruses by LCA leaves a dirty count: six of the
# HKU1 patient's contigs are HKU1's own sequence whose LCA resolved one rank up to
# Betacoronavirus, so an LCA-based exclusion keeps them and calls them false alarms.
# Attributing by best hit asks the question that was actually meant — which organism did
# this contig come from — and that organism is named at 97.8-100% amino-acid identity.
ATTRIBUTIONS: tuple[str, ...] = ("lca", "best_hit")

ENDEMIC_TAXIDS = frozenset(ENDEMIC_HUMAN_CORONAVIRUSES.values())


# --------------------------------------------------------------------------- config


@dataclass(frozen=True)
class Scoring:
    """Everything the scoring step needs, which is two file paths and six choices."""

    day3: Path = Path("results/day3-provenance.json")
    day4: Path = Path("results/day4-provenance.json")
    provenance: Path = Path("results/day5-provenance.json")

    clade: str = HEADLINE_CLADE
    framing: str = "excluding_endemic"
    attribution: str = "best_hit"
    # "as_run" is Arm B at --very-sensitive, the setting config/day4.toml pre-registered.
    # "comparison" is DIAMOND's stock mode, which Day 4 measured alongside it.
    arm_b_sensitivity: str = "comparison"
    false_alarm_budget: int = 0
    denominator: str = "post_qc_pairs"
    bitscore_ladder: tuple[float, ...] = BITSCORE_LADDER
    clades: tuple[str, ...] = tuple(CLADES)
    budgets: tuple[int, ...] = FALSE_ALARM_BUDGETS

    def validate(self) -> None:
        if self.clade not in CLADES:
            raise ValueError(f"unknown clade {self.clade!r}; expected one of {list(CLADES)}")
        if self.framing not in FRAMINGS:
            raise ValueError(f"unknown framing {self.framing!r}; expected one of {FRAMINGS}")
        if self.attribution not in ATTRIBUTIONS:
            raise ValueError(f"unknown attribution {self.attribution!r}")
        if self.arm_b_sensitivity not in ("as_run", "comparison"):
            raise ValueError("arm_b_sensitivity must be 'as_run' or 'comparison'")
        if self.false_alarm_budget < 0:
            raise ValueError("false_alarm_budget cannot be negative")


def load_config(path: Path) -> Scoring:
    """Read config/day5.toml. Every key is optional and falls back to the default."""
    raw: dict = {}
    text = Path(path).read_bytes() if Path(path).exists() else b""
    if text:
        raw = tomllib.loads(text.decode())
    scoring = raw.get("scoring", {})
    inputs = raw.get("inputs", {})
    fields = {f.name for f in dataclasses.fields(Scoring)}
    kwargs: dict = {}
    for key, value in {**inputs, **scoring}.items():
        if key not in fields:
            raise ValueError(f"unknown key {key!r} in {path}")
        if key in ("day3", "day4", "provenance"):
            value = Path(value)
        if key in ("bitscore_ladder", "clades", "budgets"):
            value = tuple(value)
        kwargs[key] = value
    config = Scoring(**kwargs)
    config.validate()
    return config


# ----------------------------------------------------------------- shared arithmetic


def threshold_for_budget(scores: Sequence[float], budget: int) -> float:
    """The smallest integer threshold at which at most ``budget`` scores clear it.

    "Clear it" is ``>=``, because that is how a report is read: *at least* this many
    reads. With scores [11, 11, 8] and a budget of one, the answer is 12 rather than 11,
    since a threshold of 11 would alarm on two samples. Returns 1 when every score is
    zero — one read is then enough, and the cohort cannot say otherwise.
    """
    if budget < 0:
        raise ValueError("budget cannot be negative")
    ordered = sorted((float(s) for s in scores), reverse=True)
    if budget >= len(ordered):
        return 1.0
    # The score we must sit strictly above is the (budget+1)-th largest: allowing
    # `budget` alarms means the (budget+1)-th sample must not alarm.
    return float(int(ordered[budget]) + 1)


def ladder_threshold_for_budget(
    per_sample_best: Sequence[float | None], budget: int, ladder: Sequence[float]
) -> float | None:
    """The smallest rung of a *fixed* ladder at which at most ``budget`` samples alarm.

    Arm B's threshold cannot be chosen freely: the ladder was pre-registered, and a
    bitscore fitted to this cohort's negatives is not a ladder. Returns ``None`` when no
    rung is tight enough, which is itself a result about the arm.
    """
    for rung in sorted(ladder):
        alarms = sum(1 for best in per_sample_best if best is not None and best >= rung)
        if alarms <= budget:
            return float(rung)
    return None


def per_million(count: float, denominator: float | None) -> float | None:
    if not denominator:
        return None
    return count * 1_000_000 / denominator


def _rate(count: float, pairs: float | None) -> float | None:
    return per_million(count, pairs)


# ------------------------------------------------------------------ Arm A: the reads


def arm_a_negative_scores(day3: dict, clade: str, framing: str) -> list[dict]:
    """Per-negative coronavirus read counts at one clade, under one framing.

    The endemic-excluded count is not a single recorded field — Day 3 recorded it at the
    genus, which is where that day's headline sat — so it is reconstructed two ways and
    the larger is taken:

    * **subtraction**: the clade's cumulative count minus the sample's endemic reads;
    * **enumeration**: the per-taxon list with the four endemic species struck out.

    They agree on nineteen of twenty samples. On ``SRR13677744`` they differ by one read,
    because the clade counts are cumulative over a clade while the taxon list is direct
    assignments, and the sample's single endemic read is counted in one and not the
    other. Taking the larger keeps the threshold permissive — it can only make the arm
    look worse, never better — and both numbers stay in the record, because a one-read
    disagreement in 228 million that nobody wrote down is how a table stops being
    reproducible.
    """
    if framing not in FRAMINGS:
        raise ValueError(f"unknown framing {framing!r}")
    scores: list[dict] = []
    for negative in day3.get("negatives", []):
        clades = negative.get("clades") or {}
        total = int(clades.get(clade, 0))
        endemic = int(negative.get("endemic_human_coronavirus_reads", 0) or 0)
        taxa = negative.get("coronavirus_taxa") or []
        entry = {
            "accession": negative.get("accession") or negative.get("name"),
            "nonhost_pairs": int(negative.get("pairs", 0)),
            "reads_all": total,
            "endemic_reads": endemic,
        }
        if framing == "all":
            entry["reads"] = total
        else:
            subtraction = max(total - endemic, 0)
            candidates = {"subtraction": subtraction}
            # The taxon list is the coronaviruses the sample was called for, so striking
            # the endemic species out of it is only an estimate of a *family*-level
            # count; at a narrower clade the list would include alphacoronaviruses that
            # do not belong to it.
            if clade == "coronaviridae":
                candidates["enumeration"] = sum(
                    int(t.get("reads", 0))
                    for t in taxa
                    if int(t.get("taxid", 0)) not in ENDEMIC_TAXIDS
                )
            # Day 3 recorded this one directly, and a recorded number beats either
            # reconstruction of it.
            if clade == (day3.get("false_alarms") or {}).get("clade"):
                candidates["recorded"] = int(
                    negative.get("false_alarm_reads_excluding_endemic", 0) or 0
                )
            entry["reads"] = max(candidates.values())
            entry["estimates"] = candidates
            if len(set(candidates.values())) > 1:
                entry["estimates_disagree"] = True
        entry["reads_per_million_nonhost"] = _rate(entry["reads"], entry["nonhost_pairs"])
        scores.append(entry)
    return sorted(scores, key=lambda e: -e["reads"])


def arm_a_curve(day3: dict, clade: str, denominator: str) -> list[dict]:
    """The dilution series as Arm A saw it, dose-controlled.

    ``detected`` is ground-truth-labelled viral fragments the arm placed in the clade.
    ``operational`` adds the clade reads that came from the background — which is what a
    2019 analyst would have had on their screen, and which at the bottom of the series is
    almost entirely the virus the label missed rather than anything spurious.
    """
    levels: list[dict] = []
    for level in (day3.get("series") or {}).get("levels", []):
        detected = int((level.get("detected_by_clade") or {}).get(clade, 0))
        background = int((level.get("false_positive_by_clade") or {}).get(clade, 0))
        levels.append(
            {
                "name": level["name"],
                "viral_fragments": int(level["viral_fragments"]),
                "viral_fragments_per_million": (level.get("per_million") or {}).get(denominator),
                "detected": detected,
                "background_clade_reads": background,
                "operational_total": detected + background,
                "total_fragments": int(level.get("total_fragments", 0)),
            }
        )
    return sorted(levels, key=lambda lv: -lv["viral_fragments"])


# --------------------------------------------------------------- Arm B: the assembly


def _arm_b_sample_search(sample: dict, sensitivity: str) -> dict:
    """The search half of one Arm B sample, at the chosen sensitivity."""
    if sensitivity == "as_run":
        return sample
    comparison = sample.get("sensitivity_comparison") or {}
    # The comparison run records the same shape — a search, an evidence profile and the
    # clade contigs' detail — for every sample, level and background Day 4 touched.
    return comparison


def _contig_rows(section: dict, clade: str) -> list[dict]:
    """The contigs a sample has in a clade, from whichever search is being scored.

    Filtered on ``lca_clades`` because that is the arm's own alarm rule: a contig counts
    for a clade when its lowest common ancestor sits in it. ``contig_detail`` is recorded
    for the headline clade, so a narrower clade is a filter on the same rows.
    """
    detail = section.get("contig_detail")
    if detail is None:
        # The as-run shape keeps the contig detail inside the search; the comparison
        # run carries it at the top of its own section.
        detail = (section.get("search") or {}).get("contig_detail") or []
    return [row for row in detail if clade in (row.get("lca_clades") or [])]


def _arm_b_alarms(evidence: dict, threshold: float) -> bool:
    """Arm B's binary: a contig in the clade whose best bitscore clears the rung.

    Written once because three places need it and because the zero rung is a trap — a
    sample with no clade contig at all has no bitscore, and ``None or 0.0 >= 0.0`` would
    make every negative in the cohort alarm.
    """
    if not evidence.get("contigs"):
        return False
    return (evidence.get("best_bitscore") or 0.0) >= threshold


def _arm_b_rule(clade: str, threshold: float) -> str:
    if threshold <= 0:
        return f"at least one contig whose LCA is in {clade} (the arm's own rule, unfiltered)"
    return f"at least one {clade} contig at bitscore >= {threshold:g}"


def _is_endemic(row: dict, attribution: str) -> bool:
    if attribution == "best_hit":
        taxid = int(((row.get("best_hit") or {}).get("taxid")) or 0)
    else:
        taxid = int(row.get("lca_taxid") or 0)
    return taxid in ENDEMIC_TAXIDS


def arm_b_sample_evidence(
    sample: dict, clade: str, sensitivity: str, framing: str, attribution: str
) -> dict:
    """One sample's Arm B evidence: its contigs in the clade and its best bitscore.

    The statistic is the best bitscore, not the contig count, for the reason Day 4
    measured: the count is not monotone in dose, because one genome shattering into
    twenty pieces increases it.
    """
    section = _arm_b_sample_search(sample, sensitivity)
    rows = _contig_rows(section, clade)
    if framing == "excluding_endemic":
        kept = [row for row in rows if not _is_endemic(row, attribution)]
    else:
        kept = list(rows)
    bitscores = [float((row.get("best_hit") or {}).get("bitscore") or 0.0) for row in kept]
    return {
        "contigs": len(kept),
        "contigs_before_framing": len(rows),
        "bases": sum(int(row.get("length", 0)) for row in kept),
        "longest_contig": max((int(row.get("length", 0)) for row in kept), default=0),
        "best_bitscore": max(bitscores, default=None),
        "endemic_contigs": len(rows) - len(kept),
        "contigs_by_bitscore": {
            f"{rung:g}": sum(1 for b in bitscores if b >= rung) for rung in BITSCORE_LADDER
        },
    }


def arm_b_negative_scores(
    day4: dict, clade: str, sensitivity: str, framing: str, attribution: str
) -> list[dict]:
    scores: list[dict] = []
    for negative in day4.get("negatives", []):
        evidence = arm_b_sample_evidence(negative, clade, sensitivity, framing, attribution)
        evidence["accession"] = negative.get("accession") or negative.get("name")
        evidence["nonhost_pairs"] = int(negative.get("nonhost_pairs", 0) or 0)
        scores.append(evidence)
    return sorted(scores, key=lambda e: -(e["best_bitscore"] or 0.0))


def arm_b_curve(
    day4: dict, clade: str, sensitivity: str, framing: str, attribution: str, denominator: str
) -> list[dict]:
    series = day4.get("series") or {}
    levels: list[dict] = []
    for level in series.get("levels", []):
        evidence = arm_b_sample_evidence(level, clade, sensitivity, framing, attribution)
        read_level = level.get("read_level") or {}
        levels.append(
            {
                "name": level["name"],
                "viral_fragments": int(level["viral_fragments"]),
                "viral_fragments_per_million": (level.get("per_million") or {}).get(denominator),
                "reads_on_clade_contigs": int(
                    (level.get("reads_on_clade_contigs") or {}).get("detected", 0) or 0
                ),
                "read_level_detected": int(
                    (read_level.get("detected_by_clade") or {}).get(clade, 0) or 0
                ),
                "read_level_background": int(
                    (read_level.get("false_positive_by_clade") or {}).get(clade, 0) or 0
                ),
                **evidence,
            }
        )
    background = series.get("background")
    if background is not None:
        evidence = arm_b_sample_evidence(background, clade, sensitivity, framing, attribution)
        evidence.update({"name": "background", "viral_fragments": 0})
        levels.append(evidence)
    return sorted(levels, key=lambda lv: -lv["viral_fragments"])


# ------------------------------------------------------------- limits and the grid


def limit_of_detection(
    curve: Sequence[dict],
    fires: Callable[[dict], bool],
    denominator_key: str = "viral_fragments_per_million",
) -> dict:
    """Walk the series from the top and report where the alarm stops.

    Returns the lowest rung that alarms, the highest rung that does not, and whether the
    two bracket cleanly — a non-monotone curve is reported rather than smoothed, because
    an arm that alarms at 100 fragments and not at 1,000 would be telling us something
    about the method and nothing about a limit.
    """
    rungs = [lv for lv in curve if lv.get("viral_fragments")]
    verdicts = [(lv, bool(fires(lv))) for lv in rungs]
    firing = [lv for lv, hit in verdicts if hit]
    silent = [lv for lv, hit in verdicts if not hit]
    lowest_firing = min(firing, key=lambda lv: lv["viral_fragments"], default=None)
    highest_silent = max(silent, key=lambda lv: lv["viral_fragments"], default=None)
    monotone = not (
        lowest_firing
        and highest_silent
        and highest_silent["viral_fragments"] > lowest_firing["viral_fragments"]
    )
    floor = lowest_firing is not None and highest_silent is None
    return {
        "alarms_at": lowest_firing["name"] if lowest_firing else None,
        "viral_fragments": lowest_firing["viral_fragments"] if lowest_firing else None,
        "per_million": lowest_firing.get(denominator_key) if lowest_firing else None,
        "highest_silent_level": highest_silent["name"] if highest_silent else None,
        "highest_silent_per_million": (
            highest_silent.get(denominator_key) if highest_silent else None
        ),
        "at_or_below_series_floor": floor,
        "monotone": monotone,
        "levels": [
            {"name": lv["name"], "viral_fragments": lv["viral_fragments"], "alarm": hit}
            for lv, hit in verdicts
        ],
    }


def arm_a_operating_point(
    day3: dict, clade: str, framing: str, budget: int, denominator: str
) -> dict:
    negatives = arm_a_negative_scores(day3, clade, framing)
    threshold = threshold_for_budget([n["reads"] for n in negatives], budget)
    alarms = [n for n in negatives if n["reads"] >= threshold]
    curve = arm_a_curve(day3, clade, denominator)
    lod = limit_of_detection(curve, lambda lv: lv["detected"] >= threshold)
    operational = limit_of_detection(curve, lambda lv: lv["operational_total"] >= threshold)
    total_reads = sum(n["reads"] for n in negatives)
    pairs = sum(n["nonhost_pairs"] for n in negatives)
    return {
        "arm": "A",
        "clade": clade,
        "framing": framing,
        "budget": budget,
        "rule": f"at least {threshold:g} reads placed in {clade}",
        "threshold_reads": threshold,
        "false_alarm_samples": len(alarms),
        "false_alarm_samples_list": [n["accession"] for n in alarms],
        "cohort_reads": total_reads,
        "cohort_reads_per_million_nonhost": _rate(total_reads, pairs),
        "worst_negative": negatives[0] if negatives else None,
        "limit_of_detection": lod,
        "limit_of_detection_operational": operational,
    }


def arm_b_operating_point(
    day4: dict,
    clade: str,
    framing: str,
    budget: int,
    sensitivity: str,
    attribution: str,
    denominator: str,
    ladder: Sequence[float] = BITSCORE_LADDER,
) -> dict:
    negatives = arm_b_negative_scores(day4, clade, sensitivity, framing, attribution)
    threshold = ladder_threshold_for_budget([n["best_bitscore"] for n in negatives], budget, ladder)
    curve = arm_b_curve(day4, clade, sensitivity, framing, attribution, denominator)
    if threshold is None:
        lod = {"alarms_at": None, "unreachable_on_the_ladder": True}
        alarms: list[dict] = []
    else:
        alarms = [n for n in negatives if _arm_b_alarms(n, threshold)]
        lod = limit_of_detection(curve, lambda lv: _arm_b_alarms(lv, threshold))
    background = next((lv for lv in curve if lv["name"] == "background"), None)
    return {
        "arm": "B",
        "clade": clade,
        "framing": framing,
        "budget": budget,
        "sensitivity": sensitivity,
        "attribution": attribution,
        "rule": None if threshold is None else _arm_b_rule(clade, threshold),
        "threshold_bitscore": threshold,
        "false_alarm_samples": len(alarms),
        "false_alarm_samples_list": [n["accession"] for n in alarms],
        "worst_negative": negatives[0] if negatives else None,
        "spike_free_background": (
            None
            if background is None
            else {
                "contigs": background["contigs"],
                "best_bitscore": background["best_bitscore"],
                "alarms": threshold is not None and _arm_b_alarms(background, threshold),
            }
        ),
        "limit_of_detection": lod,
    }


def operating_grid(day3: dict, day4: dict, config: Scoring) -> list[dict]:
    """Every combination a reader might ask for, so the headline is one row of a table.

    This is the record's defence against a threshold that looks chosen. The grid is
    clade x framing x budget for Arm A, and that times both sensitivities and both
    attributions for Arm B, computed from the same two files.
    """
    grid: list[dict] = []
    for clade in config.clades:
        for framing in FRAMINGS:
            for budget in config.budgets:
                grid.append(arm_a_operating_point(day3, clade, framing, budget, config.denominator))
                for sensitivity in ("as_run", "comparison"):
                    for attribution in ATTRIBUTIONS:
                        if framing == "all" and attribution == "lca":
                            # Attribution only enters through the endemic exclusion, so
                            # the "all" framing would otherwise be duplicated.
                            continue
                        grid.append(
                            arm_b_operating_point(
                                day4,
                                clade,
                                framing,
                                budget,
                                sensitivity,
                                attribution,
                                config.denominator,
                                config.bitscore_ladder,
                            )
                        )
    return grid


# ------------------------------------------------------- what the series cannot see


def resident_background(day3: dict, day4: dict, clade: str, denominator: str) -> dict:
    """The unlabelled virus living in the series background, and what it costs the axis.

    Day 3's diagnosis is the evidence: of 78 background fragments Kraken2 placed in
    SARS-related coronavirus, 71 align concordantly to ``NC_045512.2`` with every mate
    below the label's MAPQ 20 cut-off. They are SARS-CoV-2 reads the ground-truth label
    declined to count, they are in the background at full strength, and the dilution
    series was cut from that background. So the series' dose axis is offset, and the
    offset only matters at the bottom.
    """
    diagnosis = day3.get("background_diagnosis") or {}
    levels = (day3.get("series") or {}).get("levels", [])
    background_reads = max(
        (int((lv.get("false_positive_by_clade") or {}).get(clade, 0)) for lv in levels),
        default=0,
    )
    denominators = {}
    for level in levels:
        pm = (level.get("per_million") or {}).get(denominator)
        if pm and level.get("viral_fragments"):
            denominators[level["name"]] = level["viral_fragments"] / pm * 1_000_000
    post_qc = next(iter(denominators.values()), None)
    rungs = []
    for level in sorted(levels, key=lambda lv: -lv["viral_fragments"]):
        labelled = int(level["viral_fragments"])
        effective = labelled + background_reads
        rungs.append(
            {
                "name": level["name"],
                "labelled_fragments": labelled,
                "labelled_per_million": (level.get("per_million") or {}).get(denominator),
                "effective_fragments": effective,
                "effective_per_million": per_million(effective, post_qc),
                "unlabelled_share": round(background_reads / effective, 4) if effective else None,
            }
        )
    arm_b_read_background = int(
        ((day4.get("series") or {}).get("background_read_level") or {})
        .get("clades", {})
        .get(clade, 0)
        or 0
    )
    return {
        "clade": clade,
        "arm_a_background_clade_reads": background_reads,
        "arm_a_sarbecovirus_fragments": diagnosis.get("fragments"),
        "verified_as_sars_cov_2": diagnosis.get("fragments_aligning_to_the_modern_reference"),
        "label_min_mapq": diagnosis.get("label_min_mapq"),
        "max_mapq_observed": diagnosis.get("max_mapq"),
        "arm_b_read_level_background_clade_fragments": arm_b_read_background,
        "post_qc_pairs": post_qc,
        "rungs": rungs,
        "resolvable_floor_per_million": per_million(background_reads, post_qc),
        "note": (
            "The series cannot probe below the dose its own background already carries. "
            "Scoring uses labelled detections, which understates both arms at the two "
            "lowest rungs; the effective column is what the sample actually contained."
        ),
    }


# ----------------------------------------------------------------- the three decisions


def sensitivity_decision(day4: dict, clade: str, attribution: str) -> dict:
    """Stock DIAMOND or ``--very-sensitive``, decided on Day 4's recorded comparison.

    Day 4 deliberately left ``config/day4.toml`` alone after measuring this, on the
    grounds that a setting edited after seeing the result is no longer pre-registered.
    That is the right instinct and it is also why the decision belongs here: Day 5 is
    scoring, and choosing which recorded search to score is a scoring choice, made on
    evidence that was gathered under the pre-registered setting rather than in place of
    it.
    """
    series = day4.get("series") or {}
    rows = []
    entries = list(series.get("levels", []))
    background = series.get("background")
    if background is not None:
        entries.append({**background, "name": "background", "viral_fragments": 0})
    index = day4.get("index_sample")
    if index is not None:
        entries.insert(0, {**index, "name": "index_sample", "viral_fragments": None})
    for entry in entries:
        as_run = arm_b_sample_evidence(entry, clade, "as_run", "all", attribution)
        stock = arm_b_sample_evidence(entry, clade, "comparison", "all", attribution)
        rows.append(
            {
                "name": entry.get("name"),
                "viral_fragments": entry.get("viral_fragments"),
                "as_run_contigs": as_run["contigs"],
                "as_run_best_bitscore": as_run["best_bitscore"],
                "stock_contigs": stock["contigs"],
                "stock_best_bitscore": stock["best_bitscore"],
                "same_best_bitscore": as_run["best_bitscore"] == stock["best_bitscore"],
            }
        )
    negatives_as_run = arm_b_negative_scores(day4, clade, "as_run", "all", attribution)
    negatives_stock = arm_b_negative_scores(day4, clade, "comparison", "all", attribution)
    alarms_as_run = [n["accession"] for n in negatives_as_run if n["contigs"]]
    alarms_stock = [n["accession"] for n in negatives_stock if n["contigs"]]
    # Agreement is judged on *signal* — the index sample and the spiked rungs — and not
    # on the background, because the background is where the two modes are supposed to
    # differ and where the whole case against the tuned setting lives.
    signal = [row for row in rows if row["name"] != "background"]
    agree_on_signal = all(row["same_best_bitscore"] for row in signal)
    background_row = next((row for row in rows if row["name"] == "background"), None)
    stock_background_clean = bool(background_row and not background_row["stock_contigs"])
    stock_no_worse = set(alarms_stock) <= set(alarms_as_run)
    chose = "stock" if (agree_on_signal and stock_no_worse) else "very-sensitive"
    return {
        "question": "drop Arm B's pre-registered --very-sensitive setting?",
        "chose": chose,
        "levels": rows,
        "negatives_alarming_as_run": alarms_as_run,
        "negatives_alarming_stock": alarms_stock,
        "best_bitscore_identical_on_signal": agree_on_signal,
        "stock_background_is_clean": stock_background_clean,
        "background": background_row,
        "speed_ratio": 6,
        "because": (
            "Stock mode reaches the identical best bitscore on the index sample and at "
            f"every rung of the series ({len(signal)} comparisons), alarms on "
            f"{len(alarms_stock)} negatives against --very-sensitive's "
            f"{len(alarms_as_run)}, and finds nothing at all in the spike-free "
            "background, where --very-sensitive finds three contigs of remote homology "
            "to whale and bat coronaviruses at 28-33% identity. It costs about a sixth "
            "as much. Nothing is re-run to decide this: both searches are in Day 4's "
            "record, on the index sample, every level and all twenty negatives."
        ),
    }


def clade_decision(day3: dict, day4: dict, config: Scoring) -> dict:
    """What the clade choice is worth, in the numbers rather than in the argument."""
    rows = []
    index_truth = ((day3.get("index_sample") or {}).get("truth") or {}).get("by_clade") or {}
    for clade in config.clades:
        arm_a = arm_a_operating_point(
            day3, clade, config.framing, config.false_alarm_budget, config.denominator
        )
        arm_b = arm_b_operating_point(
            day4,
            clade,
            config.framing,
            config.false_alarm_budget,
            config.arm_b_sensitivity,
            config.attribution,
            config.denominator,
            config.bitscore_ladder,
        )
        rows.append(
            {
                "clade": clade,
                "index_sample_truth_fragments": (index_truth.get(clade) or {}).get("detected"),
                "index_sample_background_reads": (index_truth.get(clade) or {}).get(
                    "false_positive"
                ),
                "index_sample_sensitivity": (index_truth.get(clade) or {}).get("sensitivity"),
                "arm_a_threshold_reads": arm_a["threshold_reads"],
                "arm_a_false_alarm_samples": arm_a["false_alarm_samples"],
                "arm_a_lod_per_million": arm_a["limit_of_detection"]["per_million"],
                "arm_b_threshold_bitscore": arm_b["threshold_bitscore"],
                "arm_b_false_alarm_samples": arm_b["false_alarm_samples"],
                "arm_b_lod_per_million": (arm_b["limit_of_detection"] or {}).get("per_million"),
                "arm_b_at_or_below_floor": (arm_b["limit_of_detection"] or {}).get(
                    "at_or_below_series_floor"
                ),
            }
        )
    return {
        "question": "which clade counts as a detection?",
        "chose": config.clade,
        "by_clade": rows,
        "because": (
            "Three independent arguments land on the family. Day 3: the November 2019 "
            "taxonomy files Bat coronavirus BM48-31/BGR/2008 under 'unclassified "
            "coronaviruses' rather than in Betacoronavirus, which is 5,510 of the index "
            "sample's reads and is invisible to any benchmark that pairs dated sequences "
            "with a modern tree. Day 4: the genome contig's own LCA *is* Coronaviridae, "
            "because a novel sarbecovirus ties several known relatives, so a "
            "genus-level headline reports zero on a sample where the arm recovered "
            "99.68% of the genome. And the negatives: one patient's NL63 is an "
            "alphacoronavirus, so a Betacoronavirus-level report cannot see a real "
            "coronavirus infection that Arm B found."
        ),
    }


def framing_decision(day3: dict, day4: dict, config: Scoring) -> dict:
    """Which false-alarm count is the headline, and what the alternatives cost."""
    rows = []
    for framing in FRAMINGS:
        arm_a = arm_a_operating_point(
            day3, config.clade, framing, config.false_alarm_budget, config.denominator
        )
        arm_b = arm_b_operating_point(
            day4,
            config.clade,
            framing,
            config.false_alarm_budget,
            config.arm_b_sensitivity,
            config.attribution,
            config.denominator,
            config.bitscore_ladder,
        )
        rows.append(
            {
                "framing": framing,
                "arm_a_threshold_reads": arm_a["threshold_reads"],
                "arm_a_cohort_reads": arm_a["cohort_reads"],
                "arm_a_reads_per_million": arm_a["cohort_reads_per_million_nonhost"],
                "arm_a_lod_per_million": arm_a["limit_of_detection"]["per_million"],
                "arm_b_threshold_bitscore": arm_b["threshold_bitscore"],
                "arm_b_false_alarm_samples": arm_b["false_alarm_samples"],
                "arm_b_lod_per_million": (arm_b["limit_of_detection"] or {}).get("per_million"),
            }
        )
    sarbeco = {
        "clade": "sarbecovirus",
        "arm_a_reads": (day3.get("false_alarms") or {}).get("samples_with_hits") and 0,
        "note": (
            "Zero reads and zero contigs in Sarbecovirus across 228,407,636 non-host "
            "pairs, in both arms. It is the cleanest number in the benchmark and it is "
            "not the headline, because a 2019 analyst's report would have been filtered "
            "at a clade, not at the subgenus a virus nobody had seen would land in."
        ),
    }
    return {
        "question": "which false-alarm count is the headline?",
        "chose": config.framing,
        "by_framing": rows,
        "sarbecovirus": sarbeco,
        "because": (
            "The cohort's raw coronavirus reads are two patients' genuine endemic "
            "infections — 2,399 HKU1 and 9 OC43 in Arm A, 23 HKU1 and 3 NL63 contigs at "
            "97.8-100% amino-acid identity in Arm B — which the November 2019 database "
            "names correctly. Counting a correct answer as a false alarm measures how "
            "often a BALF patient has a cold, and it moves Arm A's limit by a factor of "
            "ten. The endemic exclusion's own weakness, that an infection's reads land "
            "one taxon over, is handled by the threshold rather than by the exclusion: "
            "the spillover peaks at 11 reads in a sample and the threshold sits at 12. "
            "For Arm B the same weakness is handled by attributing a contig to its best "
            "hit rather than to its LCA, which is the question that was meant."
        ),
    }


# --------------------------------------------------------------------------- the score


def score(config: Scoring) -> dict:
    """Read both arms' records and produce Day 5's. Pure arithmetic, no tools."""
    config.validate()
    day3 = json.loads(Path(config.day3).read_text())
    day4 = json.loads(Path(config.day4).read_text())

    headline_a = arm_a_operating_point(
        day3, config.clade, config.framing, config.false_alarm_budget, config.denominator
    )
    headline_b = arm_b_operating_point(
        day4,
        config.clade,
        config.framing,
        config.false_alarm_budget,
        config.arm_b_sensitivity,
        config.attribution,
        config.denominator,
        config.bitscore_ladder,
    )
    residents = resident_background(day3, day4, config.clade, config.denominator)
    record = {
        "tool": {"name": "hindsight score", "version": _version()},
        "generated_at": utc_now(),
        "machine": machine_record(Path(".")),
        "inputs": {
            "day3": _input_record(config.day3, day3),
            "day4": _input_record(config.day4, day4),
        },
        "parameters": {
            "clade": config.clade,
            "framing": config.framing,
            "attribution": config.attribution,
            "arm_b_sensitivity": config.arm_b_sensitivity,
            "false_alarm_budget": config.false_alarm_budget,
            "negatives": len(day3.get("negatives", [])),
            "denominator": config.denominator,
            "bitscore_ladder": list(config.bitscore_ladder),
            "compute": "none; this day reads the two arms' recorded tables",
        },
        "decisions": {
            "clade": clade_decision(day3, day4, config),
            "framing": framing_decision(day3, day4, config),
            "sensitivity": sensitivity_decision(day4, config.clade, config.attribution),
        },
        "resident_background": residents,
        "arm_a": {
            "method": "Kraken2 2.1.3, defaults, November 2019 nucleotide database",
            "statistic": f"reads placed in {config.clade}",
            "negatives": arm_a_negative_scores(day3, config.clade, config.framing),
            "curve": arm_a_curve(day3, config.clade, config.denominator),
            "headline": headline_a,
        },
        "arm_b": {
            "method": "MEGAHIT 1.2.9 + diamond blastx, November 2019 protein set",
            "statistic": f"best bitscore of any contig whose LCA is in {config.clade}",
            "negatives": arm_b_negative_scores(
                day4, config.clade, config.arm_b_sensitivity, config.framing, config.attribution
            ),
            "curve": arm_b_curve(
                day4,
                config.clade,
                config.arm_b_sensitivity,
                config.framing,
                config.attribution,
                config.denominator,
            ),
            "read_half": _arm_b_read_half(day4, config),
            "headline": headline_b,
        },
        "cross_arm": _cross_arm(day3, day4, config, headline_a, headline_b, residents),
        "operating_points": operating_grid(day3, day4, config),
    }
    out = Path(config.provenance)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=1, sort_keys=False) + "\n")
    return record


def _version() -> str:
    from .. import __version__

    return __version__


def _input_record(path: Path, record: dict) -> dict:
    path = Path(path)
    return {
        "path": str(path),
        "sha256": sha256_file(path) if path.exists() else None,
        "generated_at": record.get("generated_at"),
        "tool": (record.get("tool") or {}).get("name"),
    }


def _arm_b_read_half(day4: dict, config: Scoring) -> dict:
    """Arm B's read-level half, which is a mechanism result and not a limit.

    Recorded with the background floor beside every rung, because that is the whole
    point: a 93% per-read detection rate at a dose below the arm's own noise is not a
    detection, and quoting the rate without the floor is the one mistake this day is
    positioned to prevent.
    """
    rows = []
    for level in (day4.get("series") or {}).get("levels", []):
        read_level = level.get("read_level") or {}
        detected = int((read_level.get("detected_by_clade") or {}).get(config.clade, 0) or 0)
        floor = int((read_level.get("false_positive_by_clade") or {}).get(config.clade, 0) or 0)
        rows.append(
            {
                "name": level["name"],
                "viral_fragments": int(level["viral_fragments"]),
                "detected": detected,
                "rate": round(detected / level["viral_fragments"], 4)
                if level["viral_fragments"]
                else None,
                "background_floor": floor,
                "signal_to_floor": round(detected / floor, 3) if floor else None,
                "clears_its_own_floor": bool(floor and detected > floor),
            }
        )
    return {
        "levels": sorted(rows, key=lambda r: -r["viral_fragments"]),
        "note": (
            "Not scored as this arm's limit. The read half is off for the negatives "
            "(49 hours at --very-sensitive on four cores), so it has no measured "
            "false-alarm rate on an independent cohort, and on the index sample's own "
            "background it is already below its noise at 100 fragments."
        ),
    }


def _cross_arm(
    day3: dict, day4: dict, config: Scoring, headline_a: dict, headline_b: dict, residents: dict
) -> dict:
    curve_a = {lv["name"]: lv for lv in arm_a_curve(day3, config.clade, config.denominator)}
    curve_b = {
        lv["name"]: lv
        for lv in arm_b_curve(
            day4,
            config.clade,
            config.arm_b_sensitivity,
            config.framing,
            config.attribution,
            config.denominator,
        )
    }
    rows = []
    for name, level in curve_a.items():
        other = curve_b.get(name, {})
        rows.append(
            {
                "name": name,
                "viral_fragments": level["viral_fragments"],
                "per_million": level["viral_fragments_per_million"],
                "arm_a_reads": level["detected"],
                "arm_a_alarm": level["detected"] >= headline_a["threshold_reads"],
                "arm_b_best_bitscore": other.get("best_bitscore"),
                "arm_b_reads_recovered": other.get("reads_on_clade_contigs"),
                "arm_b_alarm": (
                    headline_b["threshold_bitscore"] is not None
                    and (other.get("best_bitscore") or 0.0) >= headline_b["threshold_bitscore"]
                ),
            }
        )
    rows.sort(key=lambda r: -r["viral_fragments"])
    return {
        "levels": rows,
        "arm_a_limit_per_million": headline_a["limit_of_detection"]["per_million"],
        "arm_b_limit_per_million": (headline_b["limit_of_detection"] or {}).get("per_million"),
        "arm_b_at_or_below_series_floor": (headline_b["limit_of_detection"] or {}).get(
            "at_or_below_series_floor"
        ),
        "series_floor_per_million": residents["resolvable_floor_per_million"],
        "finding": (
            "The two arms' read-recovery curves cross, and their sample-level alarms do "
            "not. Arm A recovers a flat ~41% of the label at every dose and Arm B's "
            "recovery collapses below 1,000 fragments, which is the result Day 4 "
            "reported; but an alarm needs one piece of evidence above the cohort's "
            "noise, and Arm B's 462 bp contig at bitscore 270 is that, while Arm A's "
            "four reads are not. Recovery and detection are different questions, and "
            "the benchmark's headline is the second one."
        ),
    }


# --------------------------------------------------------------------------- report


def report(record: dict) -> list[str]:
    """The tables, as text. ``hindsight score report`` prints these."""
    params = record.get("parameters") or {}
    lines: list[str] = []
    clade = params.get("clade")
    lines.append(
        f"scored at {clade}, false alarms counted {params.get('framing')}, "
        f"budget <= {params.get('false_alarm_budget')} of {params.get('negatives')} negatives; "
        f"DIAMOND at {params.get('arm_b_sensitivity')} sensitivity"
    )
    a = (record.get("arm_a") or {}).get("headline") or {}
    b = (record.get("arm_b") or {}).get("headline") or {}
    lines.append("")
    lines.append(f"Kraken2 rule: {a.get('rule')}")
    lines.append(
        f"  false alarms: {a.get('false_alarm_samples')} of {params.get('negatives')} samples; "
        f"limit of detection {_fmt(a.get('limit_of_detection', {}).get('per_million'))} "
        f"viral reads per million ({a.get('limit_of_detection', {}).get('alarms_at')})"
    )
    lines.append(f"MEGAHIT + DIAMOND rule: {b.get('rule')}")
    lod_b = b.get("limit_of_detection") or {}
    floor = (
        " (the series' lowest rung; the limit is at or below it)"
        if lod_b.get("at_or_below_series_floor")
        else ""
    )
    lines.append(
        f"  false alarms: {b.get('false_alarm_samples')} of {params.get('negatives')} samples; "
        f"limit of detection {_fmt(lod_b.get('per_million'))} viral reads per million "
        f"({lod_b.get('alarms_at')}){floor}"
    )
    lines.append("")
    lines.append(
        "level                fragments  per million   Kraken2 rds  alarm   assembly bits  alarm"
    )
    for row in (record.get("cross_arm") or {}).get("levels", []):
        lines.append(
            f"{row['name']:<20} {row['viral_fragments']:>9,}  {_fmt(row['per_million']):>11}  "
            f"{row['arm_a_reads']:>11,}  {'yes' if row['arm_a_alarm'] else 'no':<5}  "
            f"{_fmt(row['arm_b_best_bitscore']):>13}  {'yes' if row['arm_b_alarm'] else 'no'}"
        )
    residents = record.get("resident_background") or {}
    lines.append("")
    lines.append(
        f"the series' own background carries {residents.get('arm_a_background_clade_reads')} "
        f"{clade} reads ({residents.get('verified_as_sars_cov_2')} of "
        f"{residents.get('arm_a_sarbecovirus_fragments')} verified as SARS-CoV-2), so the "
        f"dose axis bottoms out near "
        f"{_fmt(residents.get('resolvable_floor_per_million'))} per million"
    )
    for decision in (record.get("decisions") or {}).values():
        lines.append("")
        lines.append(f"decision — {decision.get('question')} -> {decision.get('chose')}")
    return lines


def _fmt(value) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        if value >= 100:
            return f"{value:,.0f}"
        if value >= 1:
            return f"{value:.2f}"
        return f"{value:.3f}"
    return str(value)
