# Notes from the NCBI archive

Properties of NCBI's archived RefSeq releases, taxonomy dumps and E-utilities that this
tool had to be built around. Each was found the hard way and is now encoded in a test.

- **Release numbering is not contiguous.** RefSeq goes 99 (2 March 2020) → 200
  (4 May 2020). Releases 100–199 do not exist. Any code that interpolates a date between
  release numbers, or walks the range, is wrong.
- **Archived file timestamps are not release dates.** Release 97's catalog carries an
  mtime of 9 January 2020, because that is when it was moved into `archive/` — two months
  after the release it describes. The date comes from the release notes header, which is
  why the release table is built from those and vendored into the repo as
  [`release_dates.tsv`](../src/hindsight/data/release_dates.tsv), reviewable in a diff.
- **The release notes contain typos.** Release 232's header says "Septmeber"; release 5's
  says "May 2 , 2004"; release 39's says "January 23,2010". The parser handles those but
  refuses anything looser than an unambiguous three-letter month prefix — a fuzzy match
  would read "Smarch" as March and mis-date a release by a season.
- **E-utilities signals a rate-limit rejection with HTTP 200** and a JSON body. Code that
  checks only the status code will write `{"error":"API rate limit exceeded"}` into a
  FASTA file as a sequence.
- **NCBI's own products disagree slightly.** The release 97 catalog cites taxid 10658
  (*Enterobacteria phage PRD1*), which the 2019-11-01 taxonomy already lists as deleted.
  One taxon out of 9,305, reported in the provenance rather than surfacing later as an
  unexplained Kraken2 build warning.
- **A mapping quality means "unique among the sequences you gave the aligner"**, not "this
  read came from here". The assembly pipeline's read-counting step built its index from 23
  flagged contigs out of an assembly of 7,438, and one negative's read count came out 498
  times too high. It now runs in two passes — subset index as a candidate filter, then
  re-adjudication against the whole assembly — and records the ratio between them for
  every sample.
