"""Talking to ENA: what a run contains, and getting the bytes down intact.

ENA rather than SRA because ENA serves submitted FASTQ over plain HTTPS with a
published MD5 per file, which means a download can be verified rather than trusted,
and resumed rather than restarted. The SRA toolkit would work too and needs a local
configuration step that is one more thing to get wrong in a reproduction.

Two things here are defensive rather than decorative:

* **Resume with verification.** Long transfers through a proxy or a flaky link get
  cut, and a truncated gzip is not always an obvious error — it can decompress to a
  plausible prefix. Every file is resumed until it reaches its advertised size and
  then checked against the advertised MD5.
* **Pairing is per run, not per project.** In PRJNA699976, some runs expose ``_1`` and
  ``_2`` plus an orphan single-end file, and others expose one merged file. A loop
  that assumes ``_1``/``_2`` quietly skips half the project.
"""

from __future__ import annotations

import hashlib
import subprocess
import time
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import requests

PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

FILEREPORT_FIELDS = (
    "run_accession",
    "experiment_accession",
    "study_accession",
    "sample_accession",
    "library_strategy",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "read_count",
    "base_count",
    "collection_date",
    "fastq_ftp",
    "fastq_bytes",
    "fastq_md5",
)


@dataclass(frozen=True)
class RunFile:
    url: str
    size: int
    md5: str

    @property
    def name(self) -> str:
        return self.url.rsplit("/", 1)[-1]


@dataclass(frozen=True)
class Run:
    """One sequencing run as ENA describes it, with its files paired up."""

    accession: str
    study: str
    sample: str
    library_strategy: str
    library_layout: str
    instrument_model: str
    read_count: int
    base_count: int
    collection_date: str
    files: tuple[RunFile, ...]

    @property
    def mates(self) -> tuple[RunFile, RunFile] | None:
        """The ``_1``/``_2`` pair, if this run exposes one."""
        by_suffix = {f.name: f for f in self.files}
        one = by_suffix.get(f"{self.accession}_1.fastq.gz")
        two = by_suffix.get(f"{self.accession}_2.fastq.gz")
        return (one, two) if one and two else None

    @property
    def single(self) -> RunFile | None:
        """The merged/unpaired file, if this run exposes one."""
        for f in self.files:
            if f.name == f"{self.accession}.fastq.gz":
                return f
        return None

    @property
    def layout(self) -> str:
        """``paired``, ``single``, or ``paired+orphan`` — what the files actually are.

        ENA's own ``library_layout`` describes the experiment, not the deposit. A run
        sequenced as PAIRED can still be deposited as one merged file, and the
        pipeline has to follow the files.
        """
        if self.mates and self.single:
            return "paired+orphan"
        if self.mates:
            return "paired"
        if self.single:
            return "single"
        return "unknown"

    @property
    def total_bytes(self) -> int:
        return sum(f.size for f in self.files)


def _parse_row(header: Sequence[str], row: Sequence[str]) -> Run:
    rec = dict(zip(header, row, strict=False))
    urls = [u for u in rec.get("fastq_ftp", "").split(";") if u]
    sizes = [s for s in rec.get("fastq_bytes", "").split(";") if s]
    md5s = [m for m in rec.get("fastq_md5", "").split(";") if m]
    if not (len(urls) == len(sizes) == len(md5s)):
        raise ValueError(
            f"{rec.get('run_accession')}: ENA listed {len(urls)} files, {len(sizes)} "
            f"sizes and {len(md5s)} checksums; refusing to guess the pairing"
        )
    files = tuple(
        RunFile(url=f"https://{u}" if not u.startswith("http") else u, size=int(s), md5=m)
        for u, s, m in zip(urls, sizes, md5s, strict=True)
    )
    return Run(
        accession=rec["run_accession"],
        study=rec.get("study_accession", ""),
        sample=rec.get("sample_accession", ""),
        library_strategy=rec.get("library_strategy", ""),
        library_layout=rec.get("library_layout", ""),
        instrument_model=rec.get("instrument_model", ""),
        read_count=int(rec.get("read_count") or 0),
        base_count=int(rec.get("base_count") or 0),
        collection_date=rec.get("collection_date", ""),
        files=files,
    )


def parse_filereport(text: str) -> list[Run]:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []
    header = lines[0].split("\t")
    return [_parse_row(header, ln.split("\t")) for ln in lines[1:]]


def filereport(
    accession: str,
    *,
    session: requests.Session | None = None,
    attempts: int = 5,
    log=None,
) -> list[Run]:
    """Ask ENA what a run, experiment, sample or project contains.

    Retried, because this is a one-line metadata lookup standing in front of hours of
    work: a cohort run asks ENA once per sample, and an unretried read timeout on the
    eleventh of twenty negatives takes the whole run down after the ten before it have
    already been fetched, assembled and searched. The file downloads have always had
    their own retry loop; this call was the gap.
    """
    sess = session or requests.Session()
    last: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            response = sess.get(
                PORTAL,
                params={
                    "accession": accession,
                    "result": "read_run",
                    "fields": ",".join(FILEREPORT_FIELDS),
                    "format": "tsv",
                },
                timeout=120,
            )
            response.raise_for_status()
            return parse_filereport(response.text)
        except (requests.RequestException, OSError) as exc:
            last = exc
            if attempt == attempts:
                break
            delay = min(60.0, 5.0 * 2 ** (attempt - 1))
            if log is not None:
                log(
                    f"{accession}: ENA metadata attempt {attempt} failed ({exc}); "
                    f"retrying in {delay:.0f}s"
                )
            time.sleep(delay)
    raise RuntimeError(f"{accession}: ENA metadata unavailable after {attempts} attempts") from last


def md5_file(path: Path, *, chunk: int = 1 << 20) -> str:
    digest = hashlib.md5()  # noqa: S324 - ENA publishes MD5, not a security claim
    with open(path, "rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def download(
    file: RunFile,
    dest_dir: Path,
    *,
    attempts: int = 60,
    verify: bool = True,
    log=None,
) -> Path:
    """Fetch one file, resuming until it is whole, then verify it.

    Returns the path. A file already present at the right size with the right MD5 is
    left alone, which is what makes the whole pipeline resumable: interrupt it, run the
    same command again, and it picks up where it stopped.
    """
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / file.name
    note = log or (lambda msg: None)

    for attempt in range(1, attempts + 1):
        have = path.stat().st_size if path.exists() else 0
        if have == file.size:
            break
        if have > file.size:
            note(f"{file.name}: local copy is larger than ENA's; starting over")
            path.unlink()
            have = 0
        note(f"{file.name}: {have}/{file.size} bytes, attempt {attempt}")
        # curl rather than requests: resume, retry and sparse-write behaviour that is
        # already correct, and it is the same command a reader can run by hand.
        subprocess.run(  # noqa: S603, S607
            [
                "curl",
                "--silent",
                "--show-error",
                "--location",
                "--continue-at",
                "-",
                "--max-time",
                "900",
                "--output",
                str(path),
                file.url,
            ],
            check=False,
        )
        if (path.stat().st_size if path.exists() else 0) == have and attempt < attempts:
            time.sleep(min(30, 2 * attempt))
    else:
        raise RuntimeError(f"{file.name}: gave up after {attempts} attempts")

    size = path.stat().st_size
    if size != file.size:
        raise RuntimeError(f"{file.name}: {size} bytes on disk, ENA says {file.size}")
    if verify:
        got = md5_file(path)
        if got != file.md5:
            raise RuntimeError(f"{file.name}: MD5 {got}, ENA says {file.md5}")
        note(f"{file.name}: {size} bytes, MD5 verified")
    return path
