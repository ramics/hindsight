"""The provenance record that makes a build auditable.

Everything needed to answer "what exactly is in this database, and why is it entitled
to call itself November 2019" goes in one JSON file beside the sequences: the cutoff
asked for, the release chosen and its date, the next release that was deliberately
excluded, the taxonomy snapshot, counts, withdrawn records, and a SHA-256 of every
output file.
"""

from __future__ import annotations

import datetime as dt
import json
import platform
import sys
from pathlib import Path

from . import __version__
from .fasta import sha256_file

PROVENANCE_NAME = "provenance.json"


def checksum_outputs(out_dir: Path, names: list[str]) -> dict[str, dict]:
    out_dir = Path(out_dir)
    result: dict[str, dict] = {}
    for name in names:
        path = out_dir / name
        if not path.exists():
            continue
        result[name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    return result


def build_record(
    *,
    cutoff: dt.date,
    release,
    next_release,
    catalog_url: str,
    catalog_bytes: int | None,
    manifest_stats: dict,
    fetch_reports: list[dict],
    taxonomy: dict | None,
    outputs: dict[str, dict],
    command: list[str] | None = None,
) -> dict:
    return {
        "tool": {
            "name": "hindsight",
            "version": __version__,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
        "generated_at": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
        "command": command if command is not None else list(sys.argv),
        "time_control": {
            "cutoff_date": cutoff.isoformat(),
            "rule": "last RefSeq release issued strictly before the cutoff date",
            "release": release.number,
            "release_date": release.date.isoformat(),
            "excluded_next_release": (
                {"release": next_release.number, "release_date": next_release.date.isoformat()}
                if next_release
                else None
            ),
            "catalog_url": catalog_url,
            "catalog_bytes": catalog_bytes,
        },
        "manifest": manifest_stats,
        "fetch": fetch_reports,
        "taxonomy": taxonomy,
        "outputs": outputs,
    }


def write(out_dir: Path, record: dict) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / PROVENANCE_NAME
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return path


def read(out_dir: Path) -> dict:
    return json.loads((Path(out_dir) / PROVENANCE_NAME).read_text())


def verify(out_dir: Path) -> list[str]:
    """Re-check every recorded output digest. Returns a list of problems."""
    record = read(out_dir)
    problems: list[str] = []
    for name, meta in sorted(record.get("outputs", {}).items()):
        path = Path(out_dir) / name
        if not path.exists():
            problems.append(f"{name}: missing")
            continue
        actual = sha256_file(path)
        if actual != meta.get("sha256"):
            problems.append(f"{name}: sha256 {actual} != recorded {meta.get('sha256')}")
    return problems
