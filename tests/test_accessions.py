import pytest

from hindsight.accessions import (
    AccessionError,
    Molecule,
    classify,
    prefix,
    split_version,
)


@pytest.mark.parametrize(
    ("accession", "base", "version"),
    [
        ("NC_004718.3", "NC_004718", 3),
        ("YP_009639955.1", "YP_009639955", 1),
        ("AC_000001.1", "AC_000001", 1),
        ("NC_045512.2", "NC_045512", 2),
        ("  NC_004718.3  ", "NC_004718", 3),
    ],
)
def test_split_version(accession, base, version):
    assert split_version(accession) == (base, version)


@pytest.mark.parametrize("accession", ["NC_004718", "", "nonsense", "NC_004718.x", "4718.1"])
def test_split_version_requires_a_version(accession):
    """An unversioned accession is an error: the version is the time control."""
    with pytest.raises(AccessionError):
        split_version(accession)


@pytest.mark.parametrize(
    ("accession", "molecule"),
    [
        ("NC_004718.3", Molecule.NUCLEOTIDE),
        ("AC_000001.1", Molecule.NUCLEOTIDE),
        ("NZ_CP012345.1", Molecule.NUCLEOTIDE),
        ("YP_009639955.1", Molecule.PROTEIN),
        ("NP_619791.1", Molecule.PROTEIN),
        ("AP_000001.1", Molecule.PROTEIN),
        ("WP_000012345.1", Molecule.PROTEIN),
    ],
)
def test_classify(accession, molecule):
    assert classify(accession) is molecule


def test_classify_rejects_an_unknown_prefix():
    with pytest.raises(AccessionError, match="unrecognised"):
        classify("QQ_000001.1")


def test_entrez_db_mapping():
    assert Molecule.NUCLEOTIDE.entrez_db == "nuccore"
    assert Molecule.PROTEIN.entrez_db == "protein"


def test_prefix():
    assert prefix("NC_004718.3") == "NC"
    assert prefix("YP_009639955.1") == "YP"
