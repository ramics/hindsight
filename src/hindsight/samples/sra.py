"""The fast path: SRA's Open Data mirror on S3, converted locally to FASTQ.

ENA's FASTQ is the convenient source — ready-made, MD5-published, no toolkit — and on
a good link it is the right one. It is not always a good link. From the machine this
benchmark was first built on, ``ftp.sra.ebi.ac.uk`` delivered about 1 MB/s with the
transfer dropped every few minutes, which puts the index sample at half an hour and
the twenty negatives at five hours. The same data from SRA's Open Data mirror
(``sra-pub-run-odp.s3.amazonaws.com``, no credentials, no sign-up) arrived at 24 MB/s.

The cost of the fast path is a conversion step: the mirror serves the ``.sra``
container, so ``fasterq-dump`` has to turn it back into FASTQ, and that step wants
scratch space in proportion to the run. Measured on the index sample: converting
28.3 M spots needs roughly 19 GB for the uncompressed pair, plus temporary space of
the same order again. ``--stdout`` does not avoid it — piping the interleaved stream
into ``pigz`` looks like it should stream, but ``fasterq-dump`` buffers through its
temporary directory first, which filled 20 GB in ten seconds and then died with
"storage exhausted" having emitted nothing.

So the source is chosen per run, by size:

* a run whose conversion will not fit in the free space goes through ENA's ready-made
  FASTQ, slowly but with no scratch at all;
* everything else comes from the mirror and is converted locally, which for the
  20 negatives is minutes rather than hours.

``fasterq-dump`` also refuses to start when it estimates the disk is too small, and
its estimate is conservative, so the limits are set explicitly rather than left to the
default guess.

Both sources are kept, selectable per run, and which one served a file is recorded.
They are not identical: ENA FASTQ carries the submitter's original read names, while
``fasterq-dump`` emits ``<accession>.<spot>`` names. Fragment identifiers are only
ever compared within one sample, so this does not matter to the benchmark — but it
means a dilution series built from one source cannot be diffed read-name-wise against
one built from the other, and the provenance record says which was used.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import requests

SDL = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"

#: Public, anonymous, and deliberately virtual-host style. The path-style endpoint
#: (``s3.us-east-1.amazonaws.com/sra-pub-run-odp/...``) invites any ambient AWS
#: credentials into the request and is then refused with AccessDenied, which looks
#: like the data is gone rather than like the request was signed by mistake.
ODP_URL = "https://sra-pub-run-odp.s3.amazonaws.com/sra/{acc}/{acc}"


@dataclass(frozen=True)
class RunObject:
    accession: str
    url: str
    size: int
    md5: str | None
    service: str


def locate(accession: str, *, session: requests.Session | None = None) -> RunObject:
    """Ask NCBI's Storage Delivery Locator where a run lives and what it weighs."""
    sess = session or requests.Session()
    response = sess.get(SDL, params={"acc": accession}, timeout=120)
    response.raise_for_status()
    payload = response.json()
    results = payload.get("result") or []
    if not results:
        raise RuntimeError(f"SDL returned no result for {accession}: {payload}")
    bundle = results[0]
    if bundle.get("status") not in (200, "200"):
        raise RuntimeError(f"SDL refused {accession}: {bundle.get('status')} {bundle.get('msg')}")
    for entry in bundle.get("files", []):
        if entry.get("type") != "sra":
            continue
        for location in entry.get("locations", []):
            if location.get("service") == "s3":
                return RunObject(
                    accession=accession,
                    url=ODP_URL.format(acc=accession),
                    size=int(entry["size"]),
                    md5=entry.get("md5"),
                    service="sra-pub-run-odp",
                )
    raise RuntimeError(f"no S3-hosted .sra object for {accession}")


def md5_file(path: Path, *, chunk: int = 1 << 20) -> str:
    import hashlib

    digest = hashlib.md5()  # noqa: S324 - NCBI publishes MD5; not a security claim
    with open(path, "rb") as handle:
        while block := handle.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def download(obj: RunObject, dest_dir: Path, *, attempts: int = 20, log=None) -> Path:
    """Fetch the ``.sra`` container, resuming, and verify size and MD5."""
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / f"{obj.accession}.sra"
    note = log or (lambda msg: None)
    for attempt in range(1, attempts + 1):
        have = path.stat().st_size if path.exists() else 0
        if have == obj.size:
            break
        if have > obj.size:
            path.unlink()
            have = 0
        note(f"{obj.accession}.sra: {have}/{obj.size} bytes, attempt {attempt}")
        subprocess.run(  # noqa: S603, S607
            [
                "curl",
                "--silent",
                "--show-error",
                "--location",
                "--continue-at",
                "-",
                "--max-time",
                "1800",
                "--output",
                str(path),
                obj.url,
            ],
            check=False,
        )
        if (path.stat().st_size if path.exists() else 0) == have and attempt < attempts:
            time.sleep(min(30, 2 * attempt))
    else:
        raise RuntimeError(f"{obj.accession}.sra: gave up after {attempts} attempts")
    size = path.stat().st_size
    if size != obj.size:
        raise RuntimeError(f"{obj.accession}.sra: {size} bytes, SDL says {obj.size}")
    if obj.md5:
        got = md5_file(path)
        if got != obj.md5:
            raise RuntimeError(f"{obj.accession}.sra: MD5 {got}, SDL says {obj.md5}")
        note(f"{obj.accession}.sra: {size} bytes, MD5 verified")
    return path


#: Rule of thumb from the measurements above: converting a ``.sra`` container needs
#: free space of roughly 14 times its own size — about 7x for the uncompressed FASTQ
#: pair and as much again for ``fasterq-dump``'s temporary files — before the
#: uncompressed pair is replaced by its gzipped form.
CONVERSION_SPACE_FACTOR = 14


def conversion_fits(
    sra_bytes: int, free_bytes: int, *, factor: int = CONVERSION_SPACE_FACTOR
) -> bool:
    """Whether a local conversion of this run fits in the free space.

    Checked before downloading rather than after: finding out that a 2.7 GB container
    cannot be converted *after* fetching it wastes the fetch, and finding out halfway
    through the conversion leaves a filled disk and a half-written FASTQ.
    """
    return free_bytes >= sra_bytes * factor


def fasterq_dump_paired(
    sra_path: Path,
    out_dir: Path,
    *,
    accession: str,
    threads: int = 4,
    tmp_dir: Path | None = None,
    compress_level: int = 6,
    log_path: Path | None = None,
) -> tuple[Path, Path]:
    """Convert a ``.sra`` container to a gzipped FASTQ pair.

    Written as two steps — convert, then compress each mate in turn and delete the
    uncompressed file as soon as it has a compressed counterpart — so that the peak
    disk use is the uncompressed pair plus one compressed mate, rather than the pair
    plus both.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    plain1 = out_dir / f"{accession}_1.fastq"
    plain2 = out_dir / f"{accession}_2.fastq"
    r1 = out_dir / f"{accession}_1.fastq.gz"
    r2 = out_dir / f"{accession}_2.fastq.gz"
    # Checked before the tools are: a resumed run that already has its FASTQ should
    # not need sra-tools installed at all.
    if r1.exists() and r2.exists():
        return r1, r2

    # Exactly one of the pair can exist here: a download or a conversion that stopped
    # partway. pigz refuses to overwrite, so leaving it in place fails the conversion
    # at its last step -- after the dump -- and the failure reads as the mirror being
    # unavailable. Half a FASTQ is worth nothing anyway.
    for stale in (r1, r2, plain1, plain2):
        stale.unlink(missing_ok=True)

    tmp_dir = Path(tmp_dir or out_dir.parent / "tmp")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    for tool in ("fasterq-dump", "pigz"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"{tool} is not on PATH; the SRA path needs sra-tools and pigz")

    # Held open across the subprocess call rather than wrapped in `with`, so both
    # the dumper's stdout and stderr land in one log in order. The directory may not
    # exist yet on a sample whose first stage this is, and a missing log directory
    # must not be reported as "the SRA mirror is unavailable".
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
    log_handle = open(log_path, "ab") if log_path else subprocess.DEVNULL  # noqa: SIM115
    try:
        proc = subprocess.run(  # noqa: S603
            [
                "fasterq-dump",
                "--split-3",
                "--skip-technical",
                "-e",
                str(threads),
                "-O",
                str(out_dir),
                "-t",
                str(tmp_dir),
                # The default check is conservative and refuses machines that do fit;
                # conversion_fits() above is the check this pipeline relies on.
                "--size-check",
                "off",
                "--disk-limit",
                "500G",
                "--disk-limit-tmp",
                "500G",
                str(sra_path),
            ],
            stdout=log_handle if log_path else subprocess.DEVNULL,
            stderr=log_handle if log_path else subprocess.DEVNULL,
            check=False,
        )
    finally:
        if log_path:
            log_handle.close()
    if proc.returncode != 0:
        raise RuntimeError(f"fasterq-dump exited {proc.returncode} on {accession}")
    if not (plain1.exists() and plain2.exists()):
        raise RuntimeError(
            f"{accession}: fasterq-dump wrote no _1/_2 pair "
            f"({[p.name for p in out_dir.glob(f'{accession}*')]})"
        )
    for plain, gz in ((plain1, r1), (plain2, r2)):
        subprocess.run(  # noqa: S603, S607
            ["pigz", "-p", str(threads), f"-{compress_level}", str(plain)],
            check=True,
        )
        if not gz.exists():
            raise RuntimeError(f"{accession}: pigz produced no {gz.name}")
    return r1, r2


def fastq_dump_chunked(
    sra_path: Path,
    out_dir: Path,
    *,
    accession: str,
    chunk_spots: int = 2_000_000,
    threads: int = 4,
    compress_level: int = 6,
    log=None,
) -> tuple[Path, Path]:
    """Convert a container in spot ranges, so peak disk use is one chunk, not one run.

    This is the path for a run too big to convert whole — the index sample needs about
    19 GB for its uncompressed pair, which a cheap instance does not have spare. The
    legacy ``fastq-dump`` takes a spot range (``-N``/``-X``) and uses no temporary
    space, so each chunk is dumped uncompressed, compressed with ``pigz``, appended to
    the output, and deleted. Peak extra disk is about 1.5 GB per million spots in the
    chunk.

    Appending compressed chunks is legitimate rather than a trick: the gzip format is
    concatenable, so the result decompresses to exactly the concatenation of its
    parts. ``gzip -t`` validates it and every FASTQ reader reads it.

    Avoid ``fastq-dump --gzip`` instead of this: its internal single-threaded
    compressor ran at about 35 MB/min here, where uncompressed dumping plus ``pigz``
    runs at roughly 100 times that.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    note = log or (lambda msg: None)
    r1 = out_dir / f"{accession}_1.fastq.gz"
    r2 = out_dir / f"{accession}_2.fastq.gz"
    if r1.exists() and r2.exists():
        return r1, r2
    for tool in ("fastq-dump", "pigz"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"{tool} is not on PATH; the SRA path needs sra-tools and pigz")

    part1 = out_dir / f"{accession}_1.fastq.gz.part"
    part2 = out_dir / f"{accession}_2.fastq.gz.part"
    marker = out_dir / f".{accession}.chunks.json"
    state = json.loads(marker.read_text()) if marker.exists() else {"next_spot": 1, "spots": 0}
    if state["next_spot"] > 1 and not (part1.exists() and part2.exists()):
        # The parts went away but the bookkeeping did not; starting over is the only
        # safe reading of that.
        state = {"next_spot": 1, "spots": 0}
    elif state["next_spot"] > 1:
        # A process killed midway through appending a chunk leaves a partial gzip
        # member on the end of a part file, which the next append would follow with a
        # complete copy of the same chunk: duplicated reads behind a truncated member,
        # in a file that is otherwise 2 GB of plausible FASTQ. The marker therefore
        # records the byte length each part had when the last chunk was *finished*,
        # and a resume truncates back to it.
        for part, key in ((part1, "part1_bytes"), (part2, "part2_bytes")):
            recorded = state.get(key)
            if recorded is None:
                note(f"{accession}: no recorded part length; converting from the start")
                state = {"next_spot": 1, "spots": 0}
                part1.unlink(missing_ok=True)
                part2.unlink(missing_ok=True)
                break
            if part.stat().st_size > recorded:
                note(
                    f"{accession}: discarding {part.stat().st_size - recorded} bytes of "
                    f"a half-written chunk from {part.name}"
                )
                with open(part, "r+b") as handle:
                    handle.truncate(recorded)

    scratch = out_dir / f".{accession}.chunk"
    while True:
        start = state["next_spot"]
        stop = start + chunk_spots - 1
        if scratch.exists():
            shutil.rmtree(scratch)
        scratch.mkdir(parents=True)
        subprocess.run(  # noqa: S603, S607
            [
                "fastq-dump",
                "--split-3",
                "--skip-technical",
                "-N",
                str(start),
                "-X",
                str(stop),
                "-O",
                str(scratch),
                str(sra_path),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        chunk1 = scratch / f"{sra_path.stem}_1.fastq"
        chunk2 = scratch / f"{sra_path.stem}_2.fastq"
        if not chunk1.exists() or chunk1.stat().st_size == 0:
            shutil.rmtree(scratch)
            break
        for chunk, part in ((chunk1, part1), (chunk2, part2)):
            with open(part, "ab") as out:
                subprocess.run(  # noqa: S603, S607
                    ["pigz", "-p", str(threads), f"-{compress_level}", "-c", str(chunk)],
                    stdout=out,
                    check=True,
                )
        with open(chunk1, "rb") as handle:
            lines = sum(1 for _ in handle)
        state = {
            "next_spot": stop + 1,
            "spots": state["spots"] + lines // 4,
            "part1_bytes": part1.stat().st_size,
            "part2_bytes": part2.stat().st_size,
        }
        marker.write_text(json.dumps(state, indent=2) + "\n")
        note(f"{accession}: converted {state['spots']} spots")
        shutil.rmtree(scratch)

    if not part1.exists():
        raise RuntimeError(f"{accession}: no spots converted from {sra_path}")
    part1.rename(r1)
    part2.rename(r2)
    marker.write_text(json.dumps({**state, "complete": True}, indent=2) + "\n")
    return r1, r2


def sdl_record(obj: RunObject) -> dict:
    return {
        "accession": obj.accession,
        "service": obj.service,
        "url": obj.url,
        "bytes": obj.size,
        "md5": obj.md5,
    }


def write_marker(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")
