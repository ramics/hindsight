# What hindsight does not give you

Moved out of the README to keep the front door short. The content is unchanged.


Stated plainly, because a benchmark whose database has unexamined properties is not a
benchmark.

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

## Notes from the archive

Things that are easy to get wrong, found the hard way and now encoded in tests:

- **Release numbering is not contiguous.** RefSeq goes 99 (2 March 2020) → 200
  (4 May 2020). Releases 100–199 do not exist. Any code that interpolates a date
  between release numbers, or walks the range, is wrong.
- **Archived file timestamps are not release dates.** Release 97's catalog carries an
  mtime of 9 January 2020, because that is when it was moved into `archive/` — two
  months after the release it describes. The date comes from the release notes header,
  which is why the release table is built from those and vendored into the repo as
  [`release_dates.tsv`](../src/hindsight/data/release_dates.tsv), reviewable in a diff.
- **The release notes contain typos.** Release 232's header says "Septmeber"; release
  5's says "May 2 , 2004"; release 39's says "January 23,2010". The parser handles
  those but refuses anything looser than an unambiguous three-letter month prefix — a
  fuzzy match would read "Smarch" as March and mis-date a release by a season.
- **E-utilities signals a rate-limit rejection with HTTP 200** and a JSON body. Code
  that checks only the status code will write `{"error":"API rate limit exceeded"}`
  into a FASTA file as a sequence.
