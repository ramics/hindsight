"""What both detection arms need: a dated taxonomy, the clades, and the provenance.

This started inside :mod:`hindsight.arms.kraken` and moved out when Arm B needed the
same machinery. Nothing in here knows what a classifier or an aligner is — it reads a
taxdump, decides clade membership in *that* tree, and records what the machine and the
tools were. Both arms score against the same clade sets as a result, which is what makes
Day 5 able to put their numbers in one table.

The one idea worth keeping in mind: **a taxon is an answer, and answers have dates.**
Clade membership is computed from the November 2019 ``nodes.dmp``, not from today's
tree, because the arms' output is only meaningful in the tree the database was built
against. *Bat coronavirus BM48-31/BGR/2008* is the example that bites — in 2019 it sits
under "unclassified coronaviruses" in Coronaviridae, and today it sits in Sarbecovirus.
Scoring a 2019 report against a 2026 tree would silently move 5,000 reads.
"""

from __future__ import annotations

import datetime as dt
import os
import platform
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from ..samples import reads, tools

# The clades a betacoronavirus detection would land in, from the widest useful net to
# the narrowest. Taxids are stable integers, but their *presence* in a 2019 tree is
# not something to assume, so every one is checked against the dated nodes.dmp and a
# missing one is reported rather than silently scoring zero.
#
# No SARS-CoV-2 taxon appears here. It did not exist in the November 2019 taxonomy, and
# the question is not "did the arm say SARS-CoV-2" — it could not have — but "did the
# arm put these reads in the betacoronaviruses, where a human analyst would have looked
# twice".
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

# The four endemic human coronaviruses, each already a named species in the November
# 2019 tree. They are not a list of "things to ignore": a read placed on one of them is
# the arm working, and these viruses really are in respiratory samples. They are
# separated out because the benchmark's question is about a *novel* pathogen, and a
# reader has to be able to see both numbers — coronavirus reads, and coronavirus reads
# that are not somebody's cold. Which one is the false-alarm rate is Day 5's call.
#
# SARS-CoV and MERS are deliberately absent. Both were known in 2019, but neither is
# endemic, and reads landing on SARS-related coronavirus are exactly what the index
# sample produces: treating that taxon as "known, therefore not an alarm" would discard
# the benchmark's own positive signal.
ENDEMIC_HUMAN_CORONAVIRUSES: dict[str, int] = {
    "hcov_229e": 11137,
    "hcov_nl63": 277944,
    "hcov_oc43": 31631,
    "hcov_hku1": 290028,
    # The species OC43 and HKU1 sit under in the 2019 tree; a read assigned to the
    # species rather than the strain lands here.
    "betacoronavirus_1": 694003,
}

# The narrowest clade whose membership counts as a detection for the headline number.
# Day 3 found that the family level is the defensible one in a dated tree and said so;
# this stays at the genus because both arms record every clade and Day 5 draws the line.
HEADLINE_CLADE = "betacoronavirus"


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

    def lowest_common_ancestor(self, taxids: list[int]) -> int | None:
        """The deepest taxon that every one of ``taxids`` sits at or below.

        Arm B needs this and Arm A does not: Kraken2 does its own LCA internally and
        hands back one taxid per read, whereas a translated search hands back a list of
        protein hits and the arm has to decide what they jointly mean.

        Taxids absent from the dated tree are dropped rather than failing the call: a
        hit on a protein whose taxon the 2019 taxonomy does not hold is information
        about the archive, not about this contig, and :func:`clade_sets` already
        reports that case.
        """
        lineages = [self.lineage(taxid) for taxid in taxids]
        lineages = [chain for chain in lineages if chain]
        if not lineages:
            return None
        # Root-first, so the common prefix is the shared lineage.
        shared = list(reversed(lineages[0]))
        for chain in lineages[1:]:
            other = list(reversed(chain))
            keep = 0
            for mine, theirs in zip(shared, other, strict=False):
                if mine != theirs:
                    break
                keep += 1
            shared = shared[:keep]
            if not shared:
                return None
        return shared[-1] if shared else None


def read_taxonomy(taxonomy_dir: Path) -> Taxonomy:
    taxonomy_dir = Path(taxonomy_dir)
    nodes = taxonomy_dir / "nodes.dmp"
    names_path = taxonomy_dir / "names.dmp"
    if not nodes.exists():
        raise FileNotFoundError(f"no nodes.dmp in {taxonomy_dir}")
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
    endemic: set[int] = field(default_factory=set)
    # The dated tree's names, carried along so a record can say "Enterovirus D68"
    # rather than "1 taxid 42789". The name is the 2019 one, which is the point: a
    # report written against a dated tree should be read in that tree's words.
    names: dict[int, str] = field(default_factory=dict)

    def which(self, taxid: int) -> list[str]:
        return [name for name, members in self.members.items() if taxid in members]

    def count(self, counts: dict[int, int], clade: str) -> int:
        members = self.members.get(clade, set())
        return sum(n for taxid, n in counts.items() if taxid in members)


def endemic_taxids(taxonomy: Taxonomy) -> set[int]:
    """Every taxid at or below one of the endemic human coronaviruses."""
    out: set[int] = set()
    for root in ENDEMIC_HUMAN_CORONAVIRUSES.values():
        out |= taxonomy.descendants(root)
    return out


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
    return CladeSets(
        members=members,
        missing=missing,
        roots=wanted,
        names=taxonomy.names,
        endemic=endemic_taxids(taxonomy),
    )


def top_taxa(counts: dict[int, int], limit: int, names: dict[int, str] | None = None) -> list[dict]:
    """What a sample was called, most-assigned first.

    Worth recording for every sample and not only the index one: a negative control is
    negative for SARS-CoV-2, not free of viruses, and the only way to tell a coronavirus
    false alarm from somebody's unrelated infection is to look at what the reads were
    actually called.
    """
    names = names or {}
    return [
        {"taxid": taxid, "reads": n, "taxon": names.get(taxid)}
        for taxid, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[:limit]
    ]


def per_million(count: int, denominator: int | None) -> float | None:
    if not denominator:
        return None
    return count * 1_000_000 / denominator


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
