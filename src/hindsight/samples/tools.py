"""Running the external tools, and writing down exactly which ones ran.

"Reproducible with pinned versions" is a claim about the record, not the intention.
Every external tool this pipeline invokes is looked up once, its version string
captured verbatim from the tool itself, and compared against the version the config
pins. A mismatch is reported rather than ignored: fastp 0.23.4 and fastp 0.26 do not
trim identically, and a detection limit that moves because a trimmer changed is a
result about the trimmer.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

VERSION_COMMANDS: dict[str, list[str]] = {
    "fastp": ["--version"],
    "bowtie2": ["--version"],
    "samtools": ["--version"],
    "pigz": ["--version"],
    "curl": ["--version"],
    "diamond": ["--version"],
    "megahit": ["--version"],
    "minimap2": ["--version"],
}

VERSION_PATTERNS: dict[str, str] = {
    "fastp": r"fastp\s+([0-9][0-9.]*)",
    "bowtie2": r"version\s+([0-9][0-9.]*)",
    "samtools": r"samtools\s+([0-9][0-9.]*)",
    "pigz": r"pigz\s+([0-9][0-9.]*)",
    "curl": r"curl\s+([0-9][0-9.]*)",
    # "diamond version 2.1.9" — no leading "v", unlike megahit's "MEGAHIT v1.2.9".
    "diamond": r"diamond version\s+([0-9][0-9.]*)",
    "megahit": r"MEGAHIT\s+v([0-9][0-9.]*)",
    # minimap2 prints a bare "2.26-r1175" and nothing else.
    "minimap2": r"^([0-9][0-9.]*)",
}


@dataclass(frozen=True)
class ToolVersion:
    name: str
    path: str
    version: str
    raw: str
    pinned: str | None = None

    @property
    def matches_pin(self) -> bool:
        return self.pinned is None or self.version == self.pinned


class MissingTool(RuntimeError):
    pass


def probe(name: str, *, pinned: str | None = None) -> ToolVersion:
    path = shutil.which(name)
    if path is None:
        raise MissingTool(f"{name} is not on PATH")
    proc = subprocess.run(  # noqa: S603
        [path, *VERSION_COMMANDS.get(name, ["--version"])],
        capture_output=True,
        check=False,
    )
    # Bytes, not text: curl's --version banner carries the locale's own encoding and
    # a strict UTF-8 decode of it brought the whole pipeline down at the first step.
    raw = (proc.stdout + proc.stderr).decode("utf-8", errors="replace").strip()
    first = raw.splitlines()[0] if raw else ""
    pattern = VERSION_PATTERNS.get(name, r"([0-9][0-9.]*)")
    match = re.search(pattern, raw)
    return ToolVersion(
        name=name,
        path=path,
        version=match.group(1) if match else "unknown",
        raw=first,
        pinned=pinned,
    )


def probe_all(pins: dict[str, str | None]) -> list[ToolVersion]:
    return [probe(name, pinned=pin) for name, pin in pins.items()]


def run(
    argv: list[str],
    *,
    log: Path | None = None,
    stdout: Path | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """Run a command, tee its chatter to a log file, and fail loudly.

    Tool logs are outputs, not noise: bowtie2's alignment summary is where the host
    fraction comes from, and fastp's JSON is where the QC loss comes from.
    """
    out_handle = open(stdout, "wb") if stdout else subprocess.PIPE  # noqa: SIM115
    try:
        proc = subprocess.run(  # noqa: S603
            argv,
            stdout=out_handle,
            stderr=subprocess.PIPE,
            check=False,
        )
    finally:
        if stdout:
            out_handle.close()
    if log is not None:
        log.parent.mkdir(parents=True, exist_ok=True)
        with open(log, "ab") as handle:
            handle.write(f"$ {' '.join(argv)}\n".encode())
            handle.write(proc.stderr or b"")
            handle.write(b"\n")
    if check and proc.returncode != 0:
        tail = (proc.stderr or b"").decode(errors="replace").strip().splitlines()[-5:]
        raise RuntimeError(f"{argv[0]} exited {proc.returncode}: " + " / ".join(tail))
    return proc
