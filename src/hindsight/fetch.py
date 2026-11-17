"""Resumable, version-checked retrieval of the records a manifest names.

A build fetches hundreds of thousands of records over a rate-limited API, so it will
be interrupted. The design follows from that:

* Work is split into fixed shards of accessions, in manifest order.
* Each finished shard is written to the cache and recorded in a ledger with its
  SHA-256. A resumed run re-verifies the digest and skips the shard.
* Every record returned is checked against the accession *and version* requested.
  A record that comes back at a different version is a leak, not a near-miss, so it
  is rejected rather than accepted with a warning.
* Accessions NCBI no longer serves are recorded as withdrawn, with the batch retried
  one accession at a time so a single dead record cannot mask a whole shard.
"""

from __future__ import annotations

import json
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from . import fasta
from .accessions import Molecule
from .ncbi import NcbiClient, NcbiError

DEFAULT_BATCH_SIZE = {Molecule.NUCLEOTIDE: 100, Molecule.PROTEIN: 400}
LEDGER_NAME = "fetch_state.json"


class FetchError(RuntimeError):
    pass


@dataclass
class ShardResult:
    index: int
    path: Path
    accessions: list[str]
    returned: list[str]
    withdrawn: list[str]
    sha256: str
    from_cache: bool = False


@dataclass
class FetchReport:
    molecule: Molecule
    requested: int = 0
    returned: int = 0
    withdrawn: list[str] = field(default_factory=list)
    version_mismatches: list[dict[str, str]] = field(default_factory=list)
    shards_total: int = 0
    shards_from_cache: int = 0

    def as_dict(self) -> dict:
        return {
            "molecule": self.molecule.value,
            "requested": self.requested,
            "returned": self.returned,
            "withdrawn_count": len(self.withdrawn),
            "withdrawn": sorted(self.withdrawn),
            "version_mismatches": self.version_mismatches,
            "shards_total": self.shards_total,
            "shards_from_cache": self.shards_from_cache,
        }


def plan_shards(accessions: Sequence[str], batch_size: int) -> list[list[str]]:
    if batch_size < 1:
        raise ValueError("batch_size must be at least 1")
    return [list(accessions[i : i + batch_size]) for i in range(0, len(accessions), batch_size)]


class _Ledger:
    """Shard bookkeeping, flushed after every shard so a kill -9 loses one shard."""

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self.data: dict = {"batch_size": None, "shards": {}}
        if self.path.exists():
            try:
                self.data = json.loads(self.path.read_text())
            except json.JSONDecodeError:
                self.data = {"batch_size": None, "shards": {}}
        self.data.setdefault("shards", {})

    def matches_plan(self, batch_size: int) -> bool:
        return self.data.get("batch_size") in (None, batch_size)

    def reset(self, batch_size: int) -> None:
        self.data = {"batch_size": batch_size, "shards": {}}

    def get(self, index: int) -> dict | None:
        return self.data["shards"].get(str(index))

    def put(self, result: ShardResult) -> None:
        self.data["shards"][str(result.index)] = {
            "accessions": result.accessions,
            "returned": result.returned,
            "withdrawn": result.withdrawn,
            "sha256": result.sha256,
            "file": result.path.name,
        }
        self.flush()

    def flush(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(self.data, indent=2, sort_keys=True) + "\n")
        tmp.replace(self.path)


def _verify_batch(
    requested: Sequence[str], records: list[fasta.Record]
) -> tuple[list[fasta.Record], list[str], list[dict[str, str]]]:
    """Match returned records to the exact accession.version requested."""
    wanted = {acc: i for i, acc in enumerate(requested)}
    wanted_bases = {acc.split(".", 1)[0]: acc for acc in requested}

    accepted: dict[str, fasta.Record] = {}
    mismatches: list[dict[str, str]] = []
    for record in records:
        accession = record.accession
        if accession in wanted:
            accepted[accession] = record
            continue
        base = accession.split(".", 1)[0]
        if base in wanted_bases:
            mismatches.append({"requested": wanted_bases[base], "returned": accession})
        else:
            mismatches.append({"requested": "", "returned": accession})

    ordered = [accepted[acc] for acc in requested if acc in accepted]
    withdrawn = [acc for acc in requested if acc not in accepted]
    return ordered, withdrawn, mismatches


def _fetch_batch(client: NcbiClient, db: str, accessions: Sequence[str]) -> list[fasta.Record]:
    text = client.efetch(db=db, ids=list(accessions))
    return list(fasta.parse(text))


def fetch_shard(
    client: NcbiClient,
    molecule: Molecule,
    accessions: Sequence[str],
    index: int,
    shard_dir: Path,
    report: FetchReport,
) -> ShardResult:
    db = molecule.entrez_db
    records = _fetch_batch(client, db, accessions)
    ordered, withdrawn, mismatches = _verify_batch(accessions, records)

    if withdrawn:
        # Could be withdrawn records, or one bad id poisoning the batch. Ask one at a
        # time so we can tell the difference and keep the rest of the shard.
        still_missing: list[str] = []
        recovered: dict[str, fasta.Record] = {}
        for accession in withdrawn:
            try:
                single = _fetch_batch(client, db, [accession])
            except NcbiError:
                still_missing.append(accession)
                continue
            got, missing_again, more_mismatches = _verify_batch([accession], single)
            mismatches.extend(more_mismatches)
            if got:
                recovered[accession] = got[0]
            else:
                still_missing.extend(missing_again)
        if recovered:
            merged = {record.accession: record for record in ordered} | recovered
            ordered = [merged[acc] for acc in accessions if acc in merged]
        withdrawn = still_missing

    path = Path(shard_dir) / f"shard-{index:06d}.fa"
    fasta.write_file(path, ordered)

    report.version_mismatches.extend(mismatches)
    report.withdrawn.extend(withdrawn)
    report.returned += len(ordered)

    return ShardResult(
        index=index,
        path=path,
        accessions=list(accessions),
        returned=[record.accession for record in ordered],
        withdrawn=withdrawn,
        sha256=fasta.sha256_file(path),
    )


def fetch_molecule(
    client: NcbiClient,
    molecule: Molecule,
    accessions: Sequence[str],
    shard_dir: Path,
    *,
    batch_size: int | None = None,
    progress: bool = False,
) -> tuple[list[ShardResult], FetchReport]:
    """Fetch every accession for one molecule type, resuming completed shards."""
    shard_dir = Path(shard_dir)
    shard_dir.mkdir(parents=True, exist_ok=True)
    size = batch_size or DEFAULT_BATCH_SIZE[molecule]
    shards = plan_shards(accessions, size)

    ledger = _Ledger(shard_dir / LEDGER_NAME)
    if not ledger.matches_plan(size):
        ledger.reset(size)
    ledger.data["batch_size"] = size

    report = FetchReport(molecule=molecule, requested=len(accessions), shards_total=len(shards))
    results: list[ShardResult] = []

    iterator: Iterator[tuple[int, list[str]]] = enumerate(shards)
    if progress:
        from tqdm import tqdm

        iterator = tqdm(list(enumerate(shards)), desc=f"fetch {molecule.value}", unit="shard")

    for index, batch in iterator:
        cached = ledger.get(index)
        path = shard_dir / f"shard-{index:06d}.fa"
        if (
            cached
            and cached.get("accessions") == batch
            and path.exists()
            and cached.get("sha256") == fasta.sha256_file(path)
        ):
            result = ShardResult(
                index=index,
                path=path,
                accessions=batch,
                returned=list(cached.get("returned", [])),
                withdrawn=list(cached.get("withdrawn", [])),
                sha256=cached["sha256"],
                from_cache=True,
            )
            report.returned += len(result.returned)
            report.withdrawn.extend(result.withdrawn)
            report.shards_from_cache += 1
            results.append(result)
            continue

        result = fetch_shard(client, molecule, batch, index, shard_dir, report)
        ledger.put(result)
        results.append(result)

    return results, report


def assemble(
    results: Sequence[ShardResult],
    out_path: Path,
    *,
    width: int = fasta.LINE_WIDTH,
) -> int:
    """Concatenate shards into one FASTA, rewrapped to a fixed width.

    Shards are already in manifest order, and records within a shard are stored in the
    order they were requested, so streaming them straight through gives a
    byte-identical file on every run without holding the database in memory.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with open(out_path, "w") as handle:
        for result in sorted(results, key=lambda r: r.index):
            text = result.path.read_text()
            written += fasta.write(handle, fasta.parse(text), width=width)
    return written
