"""hindsight: build a time-correct viral reference database for a given cutoff date.

The database is reconstructed from NCBI's archived RefSeq release catalogs, which
record every accession *at the version it held* in that release. Sequences are then
fetched at those exact versions, so the result is the record set as it stood, not a
present-day approximation of it.
"""

__version__ = "0.1.0"

USER_AGENT = f"hindsight/{__version__} (https://github.com/ramics/hindsight)"
