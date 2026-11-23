# What hindsight does not give you

Stated plainly, because a benchmark whose database has unexamined properties is not yet a
benchmark. Read this before writing a methods section.

**Release granularity, not day granularity.** Releases are roughly two months apart
(97 is 4 November 2019, 98 is 6 January 2020 — a 63-day gap). A 1 December 2019 cutoff
therefore yields the database as of 4 November, not 30 November. The error is one-sided
and that is deliberate: hindsight never includes something published after the cutoff,
but it may omit something published in the weeks before it. `hindsight resolve` prints
the blind window so you can quote it. If you need the other direction, move the cutoff
forward and check what the later release contains.

**RefSeq, not GenBank or nr/nt.** This reconstructs the curated RefSeq viral set,
which is what Kraken2's standard viral library and comparable reference sets are built
from. It is not the whole of GenBank. The original Wu et al. analysis searched much
larger databases, so an arm benchmarked here against `viral.faa` is a faithful
automation of their *method*, not a replay of their *search space*. GenBank has no
equivalent dated catalog, so a time-correct nr is not reconstructible this way. Say so
in the methods.

**Withdrawn records.** A record present in release 97 that NCBI has since suppressed
cannot be fetched today, no matter how the request is phrased. hindsight records every
such accession in `withdrawn.tsv` and in the provenance rather than silently shipping a
smaller database. For release 97 this set is currently empty, but it will not stay that
way forever, so check it.

**Version substitution is rejected, not tolerated.** If E-utilities ever returns a
record at a version other than the one requested, hindsight refuses it, counts it in
`version_mismatches`, and reports it on stderr. A record quietly upgraded to its
current revision is exactly the leak this tool exists to prevent, and it would be
invisible in the output FASTA.

**Monthly taxonomy snapshots.** Taxonomy dumps are archived monthly, so the snapshot
sits up to a month before the release date; the lag is recorded in the provenance
(3 days for release 97). The `taxdmp` archive starts 2014-08-01 and `new_taxdump`
starts 2018-12-01, so cutoffs earlier than that have no contemporaneous taxonomy and
hindsight says so rather than reaching for a later tree.

**NCBI's own products disagree slightly.** The release 97 catalog cites taxid 10658
(*Enterobacteria phage PRD1*), which the 2019-11-01 taxonomy already lists as deleted.
One taxon out of 9,305. hindsight checks every catalog taxid against the dated tree and
splits the misses into merged, deleted and unexplained, so a discrepancy appears in the
provenance instead of surfacing later as an unexplained Kraken2 build warning.

