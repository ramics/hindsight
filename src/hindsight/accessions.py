"""Accession parsing and molecule classification.

RefSeq accession prefixes are documented at
https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/
The viral release directory contains genomic records (NC_, AC_) and protein records
(YP_, NP_, AP_, WP_, XP_); we need the split because they are fetched from different
Entrez databases.
"""

from __future__ import annotations

import re
from enum import Enum

# RefSeq accessions come in two shapes: the familiar two-letter prefix with digits
# (NC_004718), and the WGS/assembly forms that carry a GenBank accession after the
# underscore (NZ_CP012345, NZ_ABCD01000001). Both must parse, and both must carry a
# version suffix.
ACCESSION_RE = re.compile(r"^([A-Z]{2}_[A-Z0-9]{4,}|[A-Z]{1,4}[0-9]{5,})\.([0-9]+)$")

NUCLEOTIDE_PREFIXES = frozenset({"NC", "AC", "NG", "NT", "NW", "NZ", "NM", "NR", "XM", "XR"})
PROTEIN_PREFIXES = frozenset({"NP", "XP", "YP", "AP", "WP", "ZP"})


class Molecule(str, Enum):
    NUCLEOTIDE = "nucleotide"
    PROTEIN = "protein"

    @property
    def entrez_db(self) -> str:
        return "nuccore" if self is Molecule.NUCLEOTIDE else "protein"


class AccessionError(ValueError):
    """Raised for an accession that is not a versioned RefSeq accession."""


def split_version(accession: str) -> tuple[str, int]:
    """Split ``NC_004718.3`` into ``("NC_004718", 3)``.

    The catalog always carries the version; an unversioned accession is an error
    rather than something to paper over, because the version *is* the time control.
    """
    match = ACCESSION_RE.match(accession.strip())
    if not match:
        raise AccessionError(f"not a versioned RefSeq accession: {accession!r}")
    return match.group(1), int(match.group(2))


def prefix(accession: str) -> str:
    base, _ = split_version(accession)
    head = base.split("_", 1)[0] if "_" in base else base[:2]
    return head.upper()


def classify(accession: str) -> Molecule:
    """Return the molecule type implied by a RefSeq accession prefix."""
    head = prefix(accession)
    if head in NUCLEOTIDE_PREFIXES:
        return Molecule.NUCLEOTIDE
    if head in PROTEIN_PREFIXES:
        return Molecule.PROTEIN
    raise AccessionError(f"unrecognised RefSeq accession prefix {head!r} in {accession!r}")
