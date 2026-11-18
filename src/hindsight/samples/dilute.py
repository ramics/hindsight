"""The dilution series: which reads are kept at each level, and why that is stable.

A limit of detection needs a ladder, and the ladder has to be the same ladder every
time anyone rebuilds it. Two properties matter more than they first appear:

**Determinism that does not depend on a language runtime.** ``random.sample`` is
reproducible for a given CPython, which is not the same as reproducible. Here the
selection is a *ranking*: each read identifier is scored by
``sha256(f"{seed}:{read_id}")`` and the lowest scores win. Any implementation of
SHA-256 in any language reproduces the set from the seed and the identifier list
alone, which is what a methods section should be able to promise.

**Nesting.** Because the levels are prefixes of one ranking, the 10-read level is a
subset of the 100-read level, which is a subset of the 1,000, and so on. A detection
curve built on nested levels measures sensitivity falling away as reads are removed.
Built on independent draws it also measures the luck of each draw, and with 10 reads
at the bottom that noise is the whole signal. ``independent_seeds=True`` is available
for anyone who wants to quantify that draw-to-draw variance deliberately, and the
default stays nested.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

#: Viral read-pair counts the benchmark reports a detection result for. The top of the
#: ladder is below the index sample's own viral count, so every level is a genuine
#: subsample rather than an extrapolation.
DEFAULT_LEVELS: tuple[int, ...] = (50_000, 10_000, 1_000, 100, 10)

#: Seed for the ranking. Fixed in the repository rather than generated per run: a seed
#: nobody can read off the page is not a recorded seed.
DEFAULT_SEED = 20191201


def _score(seed: int, read_id: str) -> str:
    return hashlib.sha256(f"{seed}:{read_id}".encode()).hexdigest()


def rank_reads(read_ids: Iterable[str], *, seed: int = DEFAULT_SEED) -> list[str]:
    """Return the identifiers in the seed's ranking order.

    Ties on the digest are broken by the identifier itself, so the order does not
    depend on the order the identifiers arrived in.
    """
    return sorted(read_ids, key=lambda rid: (_score(seed, rid), rid))


@dataclass(frozen=True)
class Level:
    """One rung: how many viral read pairs, which ones, and how to check them."""

    viral_reads: int
    requested: int
    seed: int
    read_ids: tuple[str, ...] = field(repr=False)
    truncated: bool = False

    @property
    def name(self) -> str:
        return level_name(self.requested)

    @property
    def selection_sha256(self) -> str:
        """Digest of the selected identifiers, in ranking order, newline separated.

        Two builds that agree on this digest selected exactly the same reads, which is
        a cheaper thing to compare than two multi-gigabyte FASTQ files.
        """
        joined = "\n".join(self.read_ids) + "\n" if self.read_ids else ""
        return hashlib.sha256(joined.encode()).hexdigest()


def level_name(requested: int) -> str:
    """A short, sortable, filesystem-safe name for a level."""
    if requested == 0:
        return "viral000000"
    return f"viral{requested:06d}"


def plan_series(
    read_ids: Sequence[str],
    *,
    levels: Sequence[int] = DEFAULT_LEVELS,
    seed: int = DEFAULT_SEED,
    independent_seeds: bool = False,
) -> list[Level]:
    """Decide, for each requested level, which viral reads it holds.

    ``read_ids`` is the set of viral read identifiers found in the index sample. A
    level asking for more reads than exist is clamped to what exists and marked
    ``truncated`` — silently serving a thinner level than the name claims would put a
    wrong x-coordinate on the detection curve.
    """
    unique = sorted(set(read_ids))
    if len(unique) != len(read_ids):
        duplicates = len(read_ids) - len(unique)
        raise ValueError(f"viral read identifiers are not unique ({duplicates} repeated)")

    ranked = rank_reads(unique, seed=seed)
    plan: list[Level] = []
    for requested in levels:
        if requested < 0:
            raise ValueError(f"level must not be negative: {requested}")
        level_seed = seed + requested if independent_seeds else seed
        pool = rank_reads(unique, seed=level_seed) if independent_seeds else ranked
        take = min(requested, len(pool))
        plan.append(
            Level(
                viral_reads=take,
                requested=requested,
                seed=level_seed,
                read_ids=tuple(pool[:take]),
                truncated=take < requested,
            )
        )
    return plan


def reads_per_million(viral_reads: int, background_reads: int) -> float:
    """Viral read pairs per million read pairs in the file, the benchmark's x-axis.

    Expressed against the total the file actually holds (background plus viral), since
    that is what a classifier sees. The provenance record also carries the figure
    against the raw pre-QC total, because that is the number a sequencing run would be
    reported with.
    """
    total = viral_reads + background_reads
    if total <= 0:
        raise ValueError("a sample with no reads has no viral fraction")
    return viral_reads / total * 1_000_000
