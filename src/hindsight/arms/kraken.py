"""Arm A: Kraken2 against a time-correct nucleotide database.

The arm is deliberately plain. It is the thing a surveillance pipeline in late 2019
would have done — exact k-mer matching against the curated viral reference set, at the
classifier's own default thresholds — and the benchmark's whole point is that the
reference set is the one that existed then. No tuning, no custom thresholds, no
post-hoc filtering. If a default-settings run would not have flagged the index sample,
that is the result.

Three things in here are less obvious than they look.

**The database is built from hindsight's output, not Kraken2's downloader.**
``kraken2-build --download-library viral`` fetches *today's* viral RefSeq, which is the
leak this whole project exists to prevent. Instead the November 2019 ``viral.fna`` goes
into ``library/``, the dated ``seqid2taxid.map`` goes in beside it — Kraken2 uses that
file verbatim when it is present, skipping its own accession lookups — and the
contemporaneous ``nodes.dmp``/``names.dmp`` become the taxonomy. The resulting database
is a 2019 database in every component, including the shape of the tree the reports are
written against.

**Low-complexity masking is kept**, because ``kraken2-build`` masks by default and an
unmasked library inflates exactly the number Day 5 is trying to measure: a poly-A tail
matching a poly-A tail in some reference genome is a false alarm that says nothing
about either sequence.

**A level of the dilution series is classified by composition, not by materialising
it.** Kraken2 assigns each fragment independently of every other fragment in the input,
so the classification of ``background + viral`` is exactly the classification of
``background`` concatenated with the classification of ``viral``. The background is
10.9 M pairs and the same for all five levels; the viral parts are 50,000 pairs at the
top level and ten at the bottom. Composing therefore costs one background run plus five
trivial ones instead of five full runs. ``verify_composition`` checks the claim against
a materialised level rather than asserting it, because "reads are independent" is the
kind of thing that is true until a tool has a batching bug.
"""

from __future__ import annotations

import datetime as dt
import gzip
import json
import os
import platform
import shutil
import subprocess
import time
from collections import Counter
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from ..samples import reads, tools

# Kraken2 ships its helpers in a libexec directory that is not on PATH; kraken2-build
# puts its own directory there before calling them, so a direct call has to find them.
KRAKEN2_LIBEXEC_CANDIDATES = (
    Path("/usr/lib/kraken2"),
    Path("/usr/libexec/kraken2"),
    Path("/usr/local/libexec/kraken2"),
)

# The clades a betacoronavirus detection would land in, from the widest useful net to
# the narrowest. Taxids are stable integers, but their *presence* in a 2019 tree is
# not something to assume, so every one is checked against the dated nodes.dmp and a
# missing one is reported rather than silently scoring zero.
#
# No SARS-CoV-2 taxon appears here. It did not exist in the November 2019 taxonomy, and
# the question is not "did Kraken2 say SARS-CoV-2" — it could not have — but "did
# Kraken2 put these reads in the betacoronaviruses, where a human analyst would have
# looked twice".
CLADES: dict[str, int] = {
    "viruses": 10239,
    "riboviria": 2559587,
    "nidovirales": 76804,
    "coronaviridae": 11118,
    "orthocoronavirinae": 2501931,
    "betacoronavirus": 694002,
    "sarbecovirus": 2509511,
    "sars_related_cov": 694009,
}

# The narrowest clade whose membership counts as a detection for the headline number.
# A read placed anywhere in Betacoronavirus is a read that would have made an analyst
# look; a read placed in Coronaviridae but no deeper is reported too, and Day 5 can
# draw the line where it likes from the per-clade counts.
HEADLINE_CLADE = "betacoronavirus"


class KrakenError(RuntimeError):
    pass


def libexec_dir() -> Path:
    for candidate in KRAKEN2_LIBEXEC_CANDIDATES:
        if (candidate / "build_db").exists():
            return candidate
    raise KrakenError(
        "cannot find Kraken2's helper directory (build_db, k2mask); looked in "
        + ", ".join(str(c) for c in KRAKEN2_LIBEXEC_CANDIDATES)
    )


# --------------------------------------------------------------------------- taxonomy


@dataclass(frozen=True)
class Taxonomy:
    """Just enough of a dated taxdump to decide clade membership."""

    parents: dict[int, int]
    ranks: dict[int, str]
    names: dict[int, str]

    def lineage(self, taxid: int) -> list[int]:
        """Taxid to root. Returns [] for a taxid the dated tree does not hold."""
        if taxid not in self.parents:
            return []
        chain = [taxid]
        seen = {taxid}
        while True:
            parent = self.parents.get(chain[-1])
            if parent is None or parent in seen or parent == chain[-1]:
                break
            chain.append(parent)
            seen.add(parent)
        return chain

    def descendants(self, root: int) -> set[int]:
        """Every taxid at or below ``root``, by walking children once."""
        children: dict[int, list[int]] = {}
        for node, parent in self.parents.items():
            if node != parent:
                children.setdefault(parent, []).append(node)
        out: set[int] = set()
        if root not in self.parents:
            return out
        stack = [root]
        while stack:
            node = stack.pop()
            if node in out:
                continue
            out.add(node)
            stack.extend(children.get(node, ()))
        return out


def read_taxonomy(taxonomy_dir: Path) -> Taxonomy:
    taxonomy_dir = Path(taxonomy_dir)
    nodes = taxonomy_dir / "nodes.dmp"
    names_path = taxonomy_dir / "names.dmp"
    if not nodes.exists():
        raise KrakenError(f"no nodes.dmp in {taxonomy_dir}")
    parents: dict[int, int] = {}
    ranks: dict[int, str] = {}
    with nodes.open() as handle:
        for line in handle:
            fields = [f.strip() for f in line.split("|")]
            if len(fields) < 3:
                continue
            taxid, parent, rank = int(fields[0]), int(fields[1]), fields[2]
            parents[taxid] = parent
            ranks[taxid] = rank
    names: dict[int, str] = {}
    if names_path.exists():
        with names_path.open() as handle:
            for line in handle:
                fields = [f.strip() for f in line.split("|")]
                if len(fields) < 4 or fields[3] != "scientific name":
                    continue
                names[int(fields[0])] = fields[1]
    return Taxonomy(parents=parents, ranks=ranks, names=names)


@dataclass
class CladeSets:
    """Clade name -> every taxid in that clade, in the dated tree."""

    members: dict[str, set[int]]
    missing: list[str] = field(default_factory=list)
    roots: dict[str, int] = field(default_factory=dict)
    # The dated tree's names, carried along so a record can say "Enterovirus D68"
    # rather than "1 taxid 42789". The name is the 2019 one, which is the point: a
    # report written against a dated tree should be read in that tree's words.
    names: dict[int, str] = field(default_factory=dict)

    def which(self, taxid: int) -> list[str]:
        return [name for name, members in self.members.items() if taxid in members]

    def count(self, counts: dict[int, int], clade: str) -> int:
        members = self.members.get(clade, set())
        return sum(n for taxid, n in counts.items() if taxid in members)


def clade_sets(taxonomy: Taxonomy, clades: dict[str, int] | None = None) -> CladeSets:
    wanted = dict(clades or CLADES)
    members: dict[str, set[int]] = {}
    missing: list[str] = []
    for name, root in wanted.items():
        found = taxonomy.descendants(root)
        if not found:
            missing.append(f"{name} (taxid {root})")
            continue
        members[name] = found
    return CladeSets(members=members, missing=missing, roots=wanted, names=taxonomy.names)


# --------------------------------------------------------------------------- database


def _library_fasta(db_dir: Path) -> Path:
    return db_dir / "library" / "viral.fna"


def build_database(
    source_db: Path,
    db_dir: Path,
    *,
    threads: int = 4,
    mask: bool = True,
    kmer_length: int | None = None,
    minimizer_length: int | None = None,
    log=print,
) -> dict:
    """Assemble a Kraken2 database from a ``hindsight build`` output directory.

    Resumable in the only way that matters: a finished ``hash.k2d`` is left alone.
    """
    source_db, db_dir = Path(source_db), Path(db_dir)
    fna = source_db / "viral.fna"
    seqid_map = source_db / "seqid2taxid.map"
    taxonomy_src = source_db / "taxonomy"
    for required in (fna, seqid_map, taxonomy_src / "nodes.dmp", taxonomy_src / "names.dmp"):
        if not required.exists():
            raise KrakenError(f"{required} is missing; build the database with taxonomy first")

    record: dict = {"db": str(db_dir), "source_db": str(source_db), "masked": bool(mask)}
    if (db_dir / "hash.k2d").exists():
        log(f"kraken2 database already built at {db_dir}")
        record["reused"] = True
    else:
        record["reused"] = False
        db_dir.mkdir(parents=True, exist_ok=True)
        (db_dir / "library").mkdir(exist_ok=True)
        (db_dir / "taxonomy").mkdir(exist_ok=True)
        for dump in ("nodes.dmp", "names.dmp"):
            shutil.copyfile(taxonomy_src / dump, db_dir / "taxonomy" / dump)
        # Kraken2 reads this file verbatim when it exists and skips building its own
        # map from accession2taxid, which is the whole reason a dated map can be used.
        shutil.copyfile(seqid_map, db_dir / "seqid2taxid.map")

        target = _library_fasta(db_dir)
        if not target.exists():
            shutil.copyfile(fna, target)
        if mask and not (target.with_suffix(target.suffix + ".masked")).exists():
            log("masking low-complexity sequence in the library (k2mask)")
            _mask_in_place(target, threads=threads)

        start = time.monotonic()
        log(f"building the kraken2 database in {db_dir} on {threads} threads")
        env = dict(os.environ)
        env["PATH"] = f"{libexec_dir()}{os.pathsep}{env.get('PATH', '')}"
        command = [
            "kraken2-build",
            "--build",
            "--db",
            str(db_dir),
            "--threads",
            str(threads),
        ]
        if kmer_length:
            command += ["--kmer-len", str(kmer_length)]
        if minimizer_length:
            command += ["--minimizer-len", str(minimizer_length)]
        log_path = db_dir / "build.log"
        with log_path.open("ab") as errlog:
            completed = subprocess.run(  # noqa: S603
                command, env=env, stdout=errlog, stderr=errlog, check=False
            )
        if completed.returncode != 0:
            tail = log_path.read_text(errors="replace")[-2000:]
            raise KrakenError(f"kraken2-build failed ({completed.returncode}):\n{tail}")
        record["build_seconds"] = round(time.monotonic() - start, 1)

    record["files"] = {}
    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        path = db_dir / name
        if path.exists():
            record["files"][name] = {"bytes": path.stat().st_size}
    record["inspect"] = inspect_database(db_dir)
    record["library"] = _library_audit(db_dir, taxonomy_src)
    return record


def _library_audit(db_dir: Path, taxonomy_src: Path) -> dict:
    """Account for every library record, including the ones Kraken2 did not index.

    Kraken2 skips a sequence whose taxid is not in the taxonomy it was given, and says
    so only as a count in its build log. A silently smaller database is the kind of
    thing that turns up later as an unexplained miss, so the difference is reconciled
    here against the dated tree: how many records went in, how many came out, and which
    accessions account for the gap.
    """
    from .. import taxonomy as taxonomy_mod

    fna = _library_fasta(db_dir)
    audit: dict = {}
    if fna.exists():
        with fna.open() as handle:
            audit["records_in_library"] = sum(1 for line in handle if line.startswith(">"))
    log_path = db_dir / "build.log"
    if log_path.exists():
        for line in log_path.read_text(errors="replace").splitlines():
            if line.startswith("Completed processing of"):
                parts = line.split()
                audit["sequences_indexed"] = int(parts[3])
                audit["base_pairs_indexed"] = int(parts[5])

    known = taxonomy_mod.read_taxids(taxonomy_src / "nodes.dmp")
    unmapped: list[dict] = []
    seqid_map = db_dir / "seqid2taxid.map"
    nucleotide_ids: set[str] = set()
    if fna.exists():
        with fna.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    nucleotide_ids.add(line[1:].split()[0])
    if seqid_map.exists():
        with seqid_map.open() as handle:
            for line in handle:
                accession, _, taxid = line.partition("\t")
                if accession not in nucleotide_ids:
                    continue
                try:
                    value = int(taxid)
                except ValueError:
                    continue
                if value not in known:
                    unmapped.append({"accession": accession, "taxid": value})
    audit["records_whose_taxid_is_absent_from_the_dated_tree"] = unmapped
    indexed = audit.get("sequences_indexed")
    in_library = audit.get("records_in_library")
    if indexed is not None and in_library is not None:
        audit["records_not_indexed"] = in_library - indexed
        audit["gap_explained_by_absent_taxids"] = (in_library - indexed) == len(unmapped)
    return audit


def _mask_in_place(fasta: Path, *, threads: int) -> None:
    """What ``mask_low_complexity.sh`` does, without needing it on PATH.

    k2mask lowercases the masked stretches and Kraken2's own script then rewrites them
    to ``x``; the sed is part of the behaviour, not a cosmetic step, because a
    lowercase base is still a base to the k-mer extractor.
    """
    masker = libexec_dir() / "k2mask"
    tmp = fasta.with_suffix(fasta.suffix + ".tmp")
    with tmp.open("wb") as out:
        masker_proc = subprocess.Popen(  # noqa: S603
            [str(masker), "-in", str(fasta), "-outfmt", "fasta", "-threads", str(threads)],
            stdout=subprocess.PIPE,
        )
        sed = subprocess.Popen(  # noqa: S603
            ["sed", "-e", "/^>/!s/[a-z]/x/g"], stdin=masker_proc.stdout, stdout=out
        )
        masker_proc.stdout.close()
        sed_rc = sed.wait()
        masker_rc = masker_proc.wait()
    if masker_rc != 0 or sed_rc != 0:
        tmp.unlink(missing_ok=True)
        raise KrakenError(f"masking failed (k2mask {masker_rc}, sed {sed_rc})")
    tmp.replace(fasta)
    fasta.with_suffix(fasta.suffix + ".masked").write_text("")


def inspect_database(db_dir: Path) -> dict:
    """Summarise a built database, and record what it holds for the clades that matter.

    ``kraken2-inspect --skip-counts`` prints only the header in 2.1.3, so the taxon
    count comes from the full listing — which is also the only place the database says
    how much distinguishing material it actually holds for a clade. "SARS-related
    coronavirus is in the database" is weaker than it sounds if the entry is one short
    genome; the minimizer count is the honest version of that statement, and it is
    recorded for every clade the arm scores.
    """
    listing = db_dir / "inspect.txt"
    try:
        with listing.open("w") as out:
            completed = subprocess.run(  # noqa: S603
                ["kraken2-inspect", "--db", str(db_dir)],
                stdout=out,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
                timeout=3600,
            )
    except (OSError, subprocess.TimeoutExpired) as exc:  # pragma: no cover - environment
        return {"error": str(exc)}
    if completed.returncode != 0:
        return {"error": (completed.stderr or "").strip()[-500:]}

    header: dict[str, str] = {}
    taxa = 0
    minimizers: dict[int, int] = {}
    with listing.open() as handle:
        for line in handle:
            if line.startswith("#"):
                if ":" in line:
                    key, _, value = line.lstrip("# ").partition(":")
                    header[key.strip()] = value.strip()
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            taxa += 1
            try:
                minimizers[int(fields[4])] = int(fields[1])
            except ValueError:
                continue
    return {
        "taxa_in_db": taxa,
        "header": header,
        "clade_minimizers": {
            name: minimizers.get(taxid) for name, taxid in CLADES.items() if taxid in minimizers
        },
        "listing": str(listing),
    }


# -------------------------------------------------------------------------- classify


@dataclass
class Classification:
    """One Kraken2 run over one read pair."""

    name: str
    report: Path
    classified_output: Path
    pairs: int = 0
    classified: int = 0
    unclassified: int = 0
    seconds: float = 0.0
    counts: dict[int, int] = field(default_factory=dict)

    def record(self, clades: CladeSets) -> dict:
        return {
            "name": self.name,
            "pairs": self.pairs,
            "classified": self.classified,
            "unclassified": self.unclassified,
            "seconds": round(self.seconds, 1),
            "report": str(self.report),
            "classified_output": str(self.classified_output),
            "clades": {name: clades.count(self.counts, name) for name in clades.members},
            "top_taxa": _top_taxa(self.counts, 15, clades.names),
        }


def _top_taxa(
    counts: dict[int, int], limit: int, names: dict[int, str] | None = None
) -> list[dict]:
    """What a sample was called, most-assigned first.

    Worth recording for every sample and not only the index one: a negative control is
    negative for SARS-CoV-2, not free of viruses, and the only way to tell a
    coronavirus false alarm from somebody's unrelated infection is to look at what the
    reads were actually called.
    """
    names = names or {}
    return [
        {"taxid": taxid, "reads": n, "taxon": names.get(taxid)}
        for taxid, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:limit]
    ]


def classify(
    db_dir: Path,
    r1: Path,
    r2: Path,
    out_dir: Path,
    *,
    name: str,
    threads: int = 4,
    confidence: float = 0.0,
    minimum_hit_groups: int = 2,
    extra: Sequence[str] = (),
    resume: bool = True,
    log=print,
) -> Classification:
    """Classify one paired FASTQ and keep the report plus the *classified* reads only.

    Unclassified lines are dropped on the way out. Against a viral-only database
    almost every read in a clinical metatranscriptome is unclassified, so the full
    per-read output is gigabytes of the word "U" and the report already carries the
    count. What is kept is every classified read with its k-mer hit map, which is what
    Day 5 needs to recompute a stricter confidence threshold without re-running
    anything: raising ``--confidence`` can only ever unclassify a read, never classify
    a new one, so the permissive run's output is a superset of every stricter one's.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    report = out_dir / f"{name}.kreport"
    classified_output = out_dir / f"{name}.classified.tsv.gz"
    marker = out_dir / f".{name}.done"

    if resume and marker.exists() and report.exists() and classified_output.exists():
        log(f"{name}: already classified")
        payload = json.loads(marker.read_text())
        result = Classification(
            name=name,
            report=report,
            classified_output=classified_output,
            pairs=payload.get("pairs", 0),
            classified=payload.get("classified", 0),
            unclassified=payload.get("unclassified", 0),
            seconds=payload.get("seconds", 0.0),
        )
        result.counts = counts_from_output(classified_output)
        return result

    for stale in (report, classified_output):
        stale.unlink(missing_ok=True)

    command = [
        "kraken2",
        "--db",
        str(db_dir),
        "--threads",
        str(threads),
        "--paired",
        "--gzip-compressed",
        "--confidence",
        str(confidence),
        "--minimum-hit-groups",
        str(minimum_hit_groups),
        "--report",
        str(report),
        # No --output: Kraken2 writes per-read output to stdout by default, and
        # `--output -` does not mean stdout — it *suppresses* the per-read output
        # entirely, which costs a run to discover because the report still appears and
        # looks right. The report is a separate file and is unaffected.
        *extra,
        str(r1),
        str(r2),
    ]
    log(f"{name}: kraken2 on {threads} threads")
    start = time.monotonic()
    err_path = out_dir / f"{name}.kraken2.log"
    with err_path.open("ab") as errlog, gzip.open(classified_output, "wt") as sink:
        proc = subprocess.Popen(  # noqa: S603
            command, stdout=subprocess.PIPE, stderr=errlog, text=True, bufsize=1 << 20
        )
        counts: Counter[int] = Counter()
        classified = 0
        total = 0
        assert proc.stdout is not None
        for line in proc.stdout:
            total += 1
            if line[0] != "C":
                continue
            classified += 1
            sink.write(line)
            # C \t id \t taxid \t length \t hitmap
            fields = line.split("\t", 3)
            if len(fields) > 2:
                try:
                    counts[int(fields[2])] += 1
                except ValueError:
                    # Kraken2 writes "name (taxid 123)" when --use-names is on; this
                    # arm does not use it, but refuse to guess rather than miscount.
                    raise KrakenError(f"{name}: cannot read a taxid from {fields[2]!r}") from None
        proc.stdout.close()
        rc = proc.wait()
    if rc != 0:
        tail = err_path.read_text(errors="replace")[-2000:]
        raise KrakenError(f"{name}: kraken2 failed ({rc}):\n{tail}")

    result = Classification(
        name=name,
        report=report,
        classified_output=classified_output,
        pairs=total,
        classified=classified,
        unclassified=total - classified,
        seconds=time.monotonic() - start,
        counts=dict(counts),
    )
    marker.write_text(
        json.dumps(
            {
                "pairs": result.pairs,
                "classified": result.classified,
                "unclassified": result.unclassified,
                "seconds": round(result.seconds, 1),
            },
            indent=2,
        )
        + "\n"
    )
    log(
        f"{name}: {result.pairs:,} pairs, {result.classified:,} classified "
        f"({result.seconds / 60:.1f} min)"
    )
    return result


def iter_output(path: Path) -> Iterator[tuple[str, int]]:
    """Yield (fragment id, taxid) from a kept classified-output file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:  # type: ignore[operator]
        for line in handle:
            fields = line.split("\t", 3)
            if len(fields) < 3:
                continue
            yield fields[1], int(fields[2])


def counts_from_output(path: Path) -> dict[int, int]:
    counts: Counter[int] = Counter()
    for _, taxid in iter_output(path):
        counts[taxid] += 1
    return dict(counts)


# ------------------------------------------------------------------------ composition


def compose(parts: Iterable[Classification]) -> dict:
    """Add up classifications of disjoint inputs.

    Valid because Kraken2 assigns each fragment from the fragment alone: no
    normalisation, no abundance estimation, no dependence on what else is in the file.
    ``verify_composition`` is the check that this is still true of the installed
    version.
    """
    total = Counter()
    pairs = classified = 0
    for part in parts:
        total.update(part.counts)
        pairs += part.pairs
        classified += part.classified
    return {
        "pairs": pairs,
        "classified": classified,
        "unclassified": pairs - classified,
        "counts": dict(total),
    }


def verify_composition(composed: dict, direct: Classification) -> dict:
    """Compare a composed level against a directly classified one, taxon by taxon."""
    composed_counts = {int(k): v for k, v in composed["counts"].items()}
    direct_counts = dict(direct.counts)
    differences = {
        taxid: [composed_counts.get(taxid, 0), direct_counts.get(taxid, 0)]
        for taxid in set(composed_counts) | set(direct_counts)
        if composed_counts.get(taxid, 0) != direct_counts.get(taxid, 0)
    }
    return {
        "level": direct.name,
        "pairs_composed": composed["pairs"],
        "pairs_direct": direct.pairs,
        "classified_composed": composed["classified"],
        "classified_direct": direct.classified,
        "taxa_disagreeing": len(differences),
        "differences": dict(sorted(differences.items())[:20]),
        "identical": not differences and composed["pairs"] == direct.pairs,
    }


# ----------------------------------------------------------------------------- truth


def score_against_truth(
    classified_output: Path,
    viral_ids: set[str],
    clades: CladeSets,
    *,
    clade: str = HEADLINE_CLADE,
) -> dict:
    """Split a sample's classified reads by whether they are really viral.

    The ground-truth label never reaches the classifier; it is applied here, afterwards,
    to the classifier's own output. ``detected`` is a read the arm placed in the clade
    *and* that really is SARS-CoV-2; ``false_positive`` is one the arm placed in the
    clade that is not.

    Every clade is scored, not only the headline one, because which clade counts as a
    detection is not a free choice in a time-correct benchmark. In the November 2019
    taxonomy *Bat coronavirus BM48-31/BGR/2008* — the nearest relative the database
    holds — sits under "unclassified coronaviruses" in Coronaviridae, not under
    Betacoronavirus, where today's tree puts it. Reads that match it are therefore
    detections at family level and invisible at genus level, and the gap between the
    two columns is that fact. Day 5 draws its own line, from this table.
    """
    detected: Counter[str] = Counter()
    false_positive: Counter[str] = Counter()
    classified_viral = 0
    for fragment_id, taxid in iter_output(classified_output):
        truly_viral = fragment_id in viral_ids
        if truly_viral:
            classified_viral += 1
        for name, members in clades.members.items():
            if taxid in members:
                if truly_viral:
                    detected[name] += 1
                else:
                    false_positive[name] += 1
    return {
        "clade": clade,
        "truth_fragments": len(viral_ids),
        "truth_fragments_classified_anywhere": classified_viral,
        "detected": detected.get(clade, 0),
        "false_positive": false_positive.get(clade, 0),
        "sensitivity": (detected.get(clade, 0) / len(viral_ids)) if viral_ids else None,
        "by_clade": {
            name: {
                "detected": detected.get(name, 0),
                "false_positive": false_positive.get(name, 0),
                "sensitivity": (detected.get(name, 0) / len(viral_ids)) if viral_ids else None,
            }
            for name in clades.members
        },
    }


def read_truth_ids(path: Path) -> set[str]:
    path = Path(path)
    if not path.exists():
        raise KrakenError(f"no ground-truth label file at {path}")
    with path.open() as handle:
        return {line.strip() for line in handle if line.strip()}


# ------------------------------------------------------------------------ provenance


def machine_record(work_dir: Path) -> dict:
    usage = shutil.disk_usage(Path(work_dir))
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "memory_gib": _memory_gib(),
        "instance_type": os.environ.get("HINDSIGHT_INSTANCE_TYPE"),
        "work_dir_disk": {
            "total_gib": round(usage.total / 2**30, 1),
            "used_gib": round(usage.used / 2**30, 1),
            "free_gib": round(usage.free / 2**30, 1),
        },
    }


def _memory_gib() -> float | None:
    try:
        with open("/proc/meminfo") as handle:
            for line in handle:
                if line.startswith("MemTotal:"):
                    return round(int(line.split()[1]) / 2**20, 1)
    except OSError:  # pragma: no cover - environment
        pass
    return None


def tool_versions(pins: dict[str, str | None]) -> list[dict]:
    out = []
    for name, pin in pins.items():
        try:
            probed = tools.probe(name, pinned=pin)
        except tools.MissingTool as exc:
            out.append({"name": name, "error": str(exc), "pinned": pin})
            continue
        out.append(
            {
                "name": probed.name,
                "path": probed.path,
                "version": probed.version,
                "raw": probed.raw,
                "pinned": probed.pinned,
            }
        )
    return out


def describe_file(path: Path) -> dict:
    path = Path(path)
    if not path.exists():
        return {}
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": reads.sha256_file(path),
    }


def utc_now() -> str:
    return dt.datetime.now(dt.UTC).replace(microsecond=0).isoformat()
