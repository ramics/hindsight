"""Which samples the benchmark uses, and how the negatives were chosen.

The positive is not a choice: SRR10971381 is *the* index sample, the BALF
metatranscriptome Wu et al. assembled the first SARS-CoV-2 genome from.

The negatives are a choice, and the choice has to be written down. PRJNA699976 holds
440 runs from the Shi et al. 2022 pre-COVID Wuhan respiratory study — 408 patient
BALF metatranscriptomes collected May 2016 to December 2017 at the same hospital, plus
28 healthy-volunteer and 4 blank controls. The benchmark uses a 20-run subset of the
408 patients, because a false-alarm rate estimated from 20 samples is already enough
to place a threshold for a one-week spike, and 408 runs is 487 GB.

Picking 20 "representative" runs by eye is how a negative set ends up quietly
selected for being easy. Instead the subset is a seeded ranking over the sorted run
accessions, the same ``sha256``-ranking the dilution series uses, restricted to runs
whose deposit exposes a ``_1``/``_2`` pair so the layout matches the positive. The
seed and the resulting accession list are both recorded, so the subset is auditable
and regenerable without this machine.

Two properties of PRJNA699976 are load-bearing enough to repeat in the methods:

* **The deposit is host-depleted.** Shi et al. submitted non-human reads only. So the
  negatives cannot be compared with the positive on raw-read composition, and a
  false-alarm rate per raw read is not available for them. It is available per
  non-host read, which is the denominator the detection arms actually work in.
* **There are no per-run collection dates.** ENA reports ``collection_date`` as
  ``missing`` for all 440 runs, and the supplementary tables do not map patient to
  run. The 2016-17 window comes from the paper's text, at cohort level, and cannot be
  verified or stratified per run from the archive.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

from .dilute import rank_reads


@dataclass(frozen=True)
class SampleSpec:
    """One run the benchmark will process."""

    accession: str
    role: str  # "index" or "negative"
    label: str = ""
    notes: str = ""

    @property
    def is_index(self) -> bool:
        return self.role == "index"


@dataclass(frozen=True)
class Config:
    """The whole Day 2 configuration, as read from a TOML file."""

    index_run: str
    negatives_project: str
    negatives_count: int
    negatives_seed: int
    negatives_manifest: Path | None
    negatives_pinned: tuple[str, ...]
    negatives_max_pairs: int
    truth_reference: str
    levels: tuple[int, ...]
    dilution_seed: int
    threads: int
    fetch_source: str
    prune_intermediates: bool
    host_index: Path
    work_dir: Path
    tool_pins: dict[str, str | None] = field(default_factory=dict)
    fastp_args: tuple[str, ...] = ()
    bowtie2_host_args: tuple[str, ...] = ()
    bowtie2_truth_args: tuple[str, ...] = ()
    min_truth_mapq: int = 20


def load_config(path: Path) -> Config:
    with open(path, "rb") as handle:
        raw = tomllib.load(handle)
    samples = raw.get("samples", {})
    negatives = raw.get("negatives", {})
    dilution = raw.get("dilution", {})
    run = raw.get("run", {})
    manifest = negatives.get("manifest")
    return Config(
        index_run=samples["index_run"],
        negatives_project=negatives.get("project", "PRJNA699976"),
        negatives_count=int(negatives.get("count", 20)),
        negatives_seed=int(negatives.get("seed", 20191201)),
        negatives_manifest=(path.parent / manifest) if manifest else None,
        negatives_pinned=tuple(negatives.get("pinned_runs", ())),
        negatives_max_pairs=int(negatives.get("max_pairs", 0)),
        truth_reference=samples.get("truth_reference", "NC_045512.2"),
        levels=tuple(int(x) for x in dilution.get("levels", (50000, 10000, 1000, 100, 10))),
        dilution_seed=int(dilution.get("seed", 20191201)),
        threads=int(run.get("threads", 4)),
        fetch_source=str(raw.get("fetch", {}).get("source", "auto")),
        prune_intermediates=bool(run.get("prune_intermediates", False)),
        host_index=Path(run["host_index"]).expanduser(),
        work_dir=Path(run.get("work_dir", "work")).expanduser(),
        tool_pins=dict(raw.get("tools", {})),
        fastp_args=tuple(raw.get("fastp", {}).get("args", ())),
        bowtie2_host_args=tuple(raw.get("bowtie2", {}).get("host_args", ())),
        bowtie2_truth_args=tuple(raw.get("bowtie2", {}).get("truth_args", ())),
        min_truth_mapq=int(raw.get("bowtie2", {}).get("min_truth_mapq", 20)),
    )


def read_negatives_manifest(path: Path) -> list[dict[str, str]]:
    """Read the annotated ENA filereport for the negatives project.

    Columns are whatever the filereport was asked for, plus the BioSample attributes
    (``biological_replicate``, ``health_state``) that distinguish the 408 patients
    from the healthy and blank controls, which the filereport alone cannot.
    """
    rows: list[dict[str, str]] = []
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"), strict=False)))
    return rows


def patient_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """The 408 patient runs: disease cohort, excluding volunteers and blanks."""
    return [r for r in rows if r.get("biological_replicate", "").startswith("Patient_")]


def paired_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Runs whose deposit exposes a ``_1``/``_2`` pair.

    189 of the 440 runs do; the rest were deposited as a single merged file. Taking
    only paired runs keeps one read layout through the pipeline, and is recorded as a
    selection criterion rather than left as an accident of which runs happened to
    work.
    """
    out = []
    for row in rows:
        names = [u.rsplit("/", 1)[-1] for u in row.get("fastq_ftp", "").split(";") if u]
        acc = row.get("run_accession", "")
        if f"{acc}_1.fastq.gz" in names and f"{acc}_2.fastq.gz" in names:
            out.append(row)
    return out


def choose_negatives(
    rows: list[dict[str, str]],
    *,
    count: int,
    seed: int,
) -> list[str]:
    """Pick ``count`` negative runs by seeded ranking over the eligible accessions."""
    eligible = sorted(r["run_accession"] for r in paired_rows(patient_rows(rows)))
    if len(eligible) < count:
        raise ValueError(
            f"only {len(eligible)} eligible paired patient runs, asked for {count}"
        )
    return sorted(rank_reads(eligible, seed=seed)[:count])


def samples_from_config(config: Config) -> list[SampleSpec]:
    """The full sample list: the index run plus the chosen negatives.

    ``pinned_runs`` in the config short-circuits the selection. Once a selection has
    been made and frozen, later runs should reproduce *that* list rather than re-derive
    it, so that adding a column to the manifest cannot silently change the cohort.
    """
    specs = [
        SampleSpec(
            accession=config.index_run,
            role="index",
            label="Wu et al. index BALF (SARS-CoV-2 positive)",
        )
    ]
    if config.negatives_pinned:
        chosen = list(config.negatives_pinned)
        how = "pinned in config"
    elif config.negatives_manifest and config.negatives_manifest.exists():
        rows = read_negatives_manifest(config.negatives_manifest)
        chosen = choose_negatives(
            rows, count=config.negatives_count, seed=config.negatives_seed
        )
        how = f"seeded ranking (seed {config.negatives_seed})"
    else:
        raise ValueError(
            "no negatives: provide [negatives].pinned_runs or a manifest to select from"
        )
    specs += [
        SampleSpec(
            accession=acc,
            role="negative",
            label="pre-COVID Wuhan BALF negative",
            notes=how,
        )
        for acc in chosen
    ]
    return specs
