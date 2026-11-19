# hindsight

Build a viral reference database as it stood on a date you choose — not an
approximation of it, but the actual record set, at the exact accession versions that
were public at the time.

```console
$ hindsight build --date 2019-12-01 --out db/november-2019
cutoff 2019-12-01 -> release 97 (2019-11-04); first excluded release: release 98 (2020-01-06)
manifest: 12179 nucleotide + 373129 protein records across 9305 taxa
nucleotide: wrote 12179 records to viral.fna (0 withdrawn, 0 version mismatches)
protein: wrote 373129 records to viral.faa (0 withdrawn, 0 version mismatches)
taxonomy: taxdmp snapshot 2019-11-01
```

hindsight exists for one job: making a detection benchmark honest. If you want to ask
whether a pipeline would have flagged a novel pathogen before anyone knew it existed,
every reference the pipeline consults has to predate the discovery. Today's RefSeq
cannot answer that question, and neither can today's RefSeq with the obvious records
deleted.

## The leakage argument

The question this tool serves is: **would an automated, pathogen-agnostic pipeline
have flagged SARS-CoV-2 in the index sample, using only what was knowable before it
existed?**

Answering it requires a database with a specific property: it must contain everything
a real analyst could have consulted in December 2019, and nothing they could not. The
easy approaches all fail:

- **Use today's database and delete SARS-CoV-2.** The genome is gone but its relatives
  are not. RaTG13 (deposited 27 January 2020), the pangolin coronaviruses, and
  thousands of SARS-CoV-2-era sequences remain, and they are the records that make
  the index sample easy to classify. You have deleted the answer and left the
  answer key.
- **Use today's database and filter by deposition date.** Better, but the records
  themselves have moved. Accessions get revised: sequence corrected, annotation
  rebuilt, taxonomy reassigned. A record first deposited in 2007 may now be at
  version 3, and version 3 may postdate your cutoff. Taxon names and the shape of the
  taxonomic tree have also changed, and a classifier's output is a taxon.
- **Rebuild from a snapshot someone happened to keep.** Not reproducible, and not
  available for an arbitrary date.

hindsight takes a different route, which works because NCBI never threw anything away.

**RefSeq release catalogs are archived back to release 1 (June 2003).** Each release
has a dated set of files, and `RefSeq-release<N>.catalog.gz` lists every accession in
that release — crucially, at the **version** it held in that release. The viral set is
a filter on one column:

```
694009  Severe acute respiratory syndrome-related coronavirus  NC_004718.3  complete|viral  na  29751
```

So hindsight:

1. Picks the last release issued **strictly before** your cutoff date.
2. Reads that release's catalog and extracts the viral accession set, versions included.
3. Fetches each record **at that exact version** from E-utilities, which honours the
   version suffix and returns the revision as it was.
4. Pairs it with the NCBI taxonomy dump archived nearest before the same release date,
   so names, ranks and tree structure are contemporaneous too.
5. Writes a `provenance.json` recording the cutoff, the release, the release NCBI
   issued next (the one deliberately excluded), the taxonomy snapshot, every count,
   every withdrawn record, and a SHA-256 of every output.

For a December 2019 cutoff that gives release 97, November 2019:

| | |
| --- | --- |
| Viral nucleotide records | 12,179 |
| Viral protein records | 373,129 |
| Distinct taxa | 9,305 |
| Coronavirus genomes | 45 |
| SARS-CoV | present, `NC_004718.3` |
| Bat SARS-related and other bat CoV genomes | 14, including `NC_014470.1` (BM48-31/BGR/2008) and `NC_025217.1` (Bat Hp-betacoronavirus/Zhejiang2013) |
| SARS-CoV-2 (`NC_045512`, `MN908947`) | absent |
| RaTG13 (`MN996532`) | absent |

That is the shape the benchmark needs. The nearest relatives a December 2019 analyst
could have reached — 79–88% identity to the index sample, which is what Wu et al.
actually had — are in. Everything learned from the pandemic is out. Nothing was
deleted to make it so; those records simply had not been written yet.

These are not claims in a README. `tests/test_leakage.py` asserts them, and the
network-marked test in that file streams the full 385,308-row release 97 viral set
from NCBI and checks the whole thing.

## Install

```bash
uv sync --locked --all-extras     # or: pip install -e ".[dev]"
```

Set `NCBI_API_KEY` if you have one ([free from your NCBI
account](https://www.ncbi.nlm.nih.gov/account/settings/)). It raises the E-utilities
request ceiling from 3/s to 10/s and takes a full build from hours to under one.

## Use

```bash
# Which release does a date select, and what does that exclude?
hindsight resolve 2019-12-01

# The accession set only — fast, no sequence fetched. Good for sanity-checking a cutoff.
hindsight build --date 2019-12-01 --out db/nov2019 --manifest-only

# The nucleotide database, for a Kraken2 arm.
hindsight build --date 2019-12-01 --out db/nov2019 --molecules nucleotide

# Everything: nucleotide, protein, contemporaneous taxonomy.
hindsight build --date 2019-12-01 --out db/nov2019

# Re-check a built database against its own provenance.
hindsight verify db/nov2019

# List the archived releases, or the dated taxonomy snapshots.
hindsight releases --limit 10
hindsight taxdumps
```

Interrupt a build and re-run the same command: it resumes from the last completed
shard. Downloads resume mid-file with HTTP Range, and fetched records are cached by
shard with a digest, so nothing is retrieved twice.

### Outputs

```
db/nov2019/
  manifest.tsv        every accession in the release's viral set, with version,
                      molecule, taxid, organism, RefSeq status and length
  viral.fna           nucleotide records, at their release-97 versions
  viral.faa           protein records, for a translated-search arm
  seqid2taxid.map     accession.version -> taxid, in Kraken2's format
  taxonomy/           nodes.dmp, names.dmp, merged.dmp, delnodes.dmp, as of 2019-11-01
  withdrawn.tsv       accessions the release listed that NCBI will no longer serve
  provenance.json     the full audit record, with SHA-256 of every file above
```

Builds are deterministic: the manifest is sorted, records are emitted in manifest
order regardless of the order NCBI answers in, and FASTA is rewrapped to a fixed
width. Two builds of the same cutoff produce byte-identical files.

## The samples

`hindsight build` is one half of the benchmark. The other half is the reads, and
`hindsight samples` prepares them:

```console
$ scripts/day2.sh
```

That fetches the index sample (`SRR10971381`, the Wu et al. BALF metatranscriptome)
and 20 time-matched pre-COVID Wuhan BALF negatives from `PRJNA699976`, runs read QC
and host filtering over both, labels which fragments in the index sample are really
SARS-CoV-2, and freezes a dilution series — the sample's own viral reads thinned to
50,000 / 10,000 / 1,000 / 100 / 10 fragments and held in the sample's own background.
Every stage is resumable, every seed is recorded, and the counts, digests, tool
versions, timings and machine land in `work/day2-provenance.json`.

The ground-truth label is the one place the pipeline is allowed to know something from
2020: fragments are labelled by aligning to the modern SARS-CoV-2 reference. The
detection arms never see that label or that reference — they see only the November
2019 database above. Knowing the answer is the benchmark's job; it is what the arms
are not allowed to do.

Two properties of the negative set have to travel with every number derived from it:
the deposit holds **non-human reads only**, so a false-alarm rate per raw sequenced
read is not available for those samples, and ENA reports **no collection date** for
any of the 440 runs, so the 2016-17 window is the paper's at cohort level and is not
verifiable per run. [`docs/day2-samples.md`](docs/day2-samples.md) has the full
account, including why the levels are nested, why the negatives are chosen by a seeded
hash rather than by eye, and what the whole thing costs in wall-clock and disk.

## The arms

With the database and the reads in place, an arm is a detection pipeline run against
them and scored on what it found. The first one is Kraken2:

```console
$ scripts/day3.sh
index sample: 25,928 of 62,117 ground-truth viral fragments in Coronaviridae (41.7%)
  viral050000: 20,892 of 50,000 detected
  viral010000:  4,232 of 10,000 detected
  viral001000:    420 of  1,000 detected
  viral000100:     39 of    100 detected
  viral000010:      4 of     10 detected
```

**A standard read classifier, at its own default settings, against a database that had
never seen SARS-CoV-2 or any of its relatives discovered after 4 November 2019, placed
25,928 of the index sample's viral reads in Coronaviridae.** It could not name the
virus — the virus had no name — and it did not need to.

Down the dilution series the detection *rate* is flat at about 41%, because Kraken2
classifies each fragment independently, so what falls is the absolute count: four
coronavirus reads survive at 0.77 viral reads per million. Whether four reads is a
detection depends on the false-alarm rate, which the 20 pre-COVID negatives measure.

One result from that arm is the argument for this whole tool, and it was found rather
than predicted. 4,922 of those reads were assigned to *Bat coronavirus
BM48-31/BGR/2008*, the nearest relative the database holds — and the **November 2019
taxonomy files that genome under "unclassified coronaviruses", not under
Betacoronavirus**, where today's tree puts it. Detection at family level is 25,928
reads; at genus level it is 20,418. A benchmark that dates its sequences but keeps a
modern taxonomy cannot see that difference, because the modern tree has since moved
BM48-31 into Sarbecovirus. A classifier's output is a taxon, and taxa move.

[`docs/day3-arm-a.md`](docs/day3-arm-a.md) has the arm in full: the per-clade curve,
why the negatives run at native depth on a machine that cannot hold them, what the
index sample's own background turns out to contain, and the one check that the
composition shortcut behind the series is exact.

The second arm is not a classifier at all. It is what the original team actually did —
assemble de novo, then match the contigs by protein homology — run against the same
dated database:

```console
$ scripts/day4.sh
index sample: 99.68% of the SARS-CoV-2 genome in one 29,809 bp contig at 100.0% identity
  viral050000: 49,077 of 50,000 detected
  viral010000:  9,824 of 10,000 detected
  viral001000:    936 of  1,000 detected
  viral000100:     41 of    100 detected
  viral000010:      1 of     10 detected
```

**The assembly recovered the genome, not merely a signal, from a database whose nearest
relative matched it at 79.9% amino acid identity.** Its best match was the SARS-CoV
replicase; the lowest common ancestor of everything it tied with is Coronaviridae, which
is the correct answer and not the confident wrong one a single-best-hit rule would give.

And the two arms cross. Above about 1,000 viral fragments assembly is more than twice as
sensitive as the classifier, because a contig recruits reads that would never classify
alone; below 100 it loses, because the assembler runs out of coverage while Kraken2's
per-read rate stays flat. Two methods with opposite failure modes is the most useful
thing a two-arm benchmark can find, and it is the reason Day 5 cannot draw one threshold
across both.

One result in that arm is a caution rather than a finding. **The spike-free background
alarms too**: every assembly yields three short contigs the taxonomy places in
Coronaviridae, whose best matches are Beluga whale coronavirus SW1 and a Rousettus bat
coronavirus at 28–33% identity. Real remote homology, so "found a contig in the family"
is not a detection rule. Bitscore separates it cleanly — the floor tops out at 58 and the
weakest genuine contig scores 270 — so the arm records a signal-strength profile and
leaves the threshold to the negatives, rather than picking one that fits the answer.

On twenty pre-COVID negatives at native depth — 228,407,636 non-host read pairs, the same
cohort and count Arm A used — **the arm places no contig in Sarbecovirus at all**. It
alarms at family level on three of the twenty, and two of those three are correct answers:
real HKU1 and NL63 infections that the November 2019 database contains and names. The
third is the background floor appearing in a real library, one 951 bp contig at a bitscore
of 46.6, and it is gone at any threshold above 100. The NL63 patient is invisible to
Arm A's headline, because NL63 is an alphacoronavirus and Day 3 reported at genus level —
a second, independent argument for scoring at family.

One correction in that arm is worth more than most of its results. The step that converts
contigs into reads built its aligner index from the flagged contigs alone, and **mapping
quality means "unique among the sequences you gave the aligner", not "this read came from
here"**. Given 23 contigs out of an assembly of 7,438, reads belonging to the other 7,415
had nowhere else to go: one negative's endemic coronavirus read count came out 498 times
too high. It now runs in two passes, using the subset index only as a candidate filter and
re-adjudicating against the whole assembly, and records the ratio between the two for
every sample. Nine of the ten were already right; the index sample's 98.1% recovery and
the dilution series' false-positive floor both stand as first measured.

[`docs/day4-arm-b.md`](docs/day4-arm-b.md) has that arm in full: the crossover, the
negatives, why every dilution level costs its own assembly, what `--very-sensitive` turns
out to be worth, the read-mapping error and what generalises from it, and the measurement
behind every runtime decision.

### What the HKU1 negative actually had

One negative's loose end from Arm B turned out to be about the reference set rather than
about the patient. Four of its 23 Coronaviridae contigs had good HKU1 protein hits and no
nucleotide alignment to the dated HKU1 genome, which reads as a second, diverged
betacoronavirus — and is not. RefSeq curates **one genome per species**, and its single
HKU1 genome is genotype A; against a genotype B genome published in 2006, all 23 contigs
align at 99.4% identity. The patient had one ordinary endemic infection whose genotype the
database does not represent, and the "divergence" was in exactly the two genes that define
the genotype.

Two things fell out of chasing it. Arm A's seven non-endemic betacoronavirus reads across
the whole negative cohort are that same infection, misassigned to rodent and bovine
relatives because no genotype B is in the database — so the cohort's non-endemic
false-alarm count at genus level is zero, not seven. And the single largest viral call
anywhere in the cohort, 11,227 reads of a horse lentivirus in a human lung, is a cloned
1.7 kb fragment of its *pol* gene: one window, edges one base wide, no coverage of *gag*,
*env* or either LTR, and the same window in all twenty samples. A read count cannot see
that; where the reads sit on the reference can.

[`docs/hku1-coinfection.md`](docs/hku1-coinfection.md) has it in full, and both checks are
commands: `hindsight arm-b diagnose-strain` for one-strain-or-two, `hindsight arm-a
diagnose-taxon` for organism-or-reagent.

## What this does not give you, and why it matters

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
  [`release_dates.tsv`](src/hindsight/data/release_dates.tsv), reviewable in a diff.
- **The release notes contain typos.** Release 232's header says "Septmeber"; release
  5's says "May 2 , 2004"; release 39's says "January 23,2010". The parser handles
  those but refuses anything looser than an unambiguous three-letter month prefix — a
  fuzzy match would read "Smarch" as March and mis-date a release by a season.
- **E-utilities signals a rate-limit rejection with HTTP 200** and a JSON body. Code
  that checks only the status code will write `{"error":"API rate limit exceeded"}`
  into a FASTA file as a sequence.

## Development

```bash
uv run pytest             # 294 offline tests
uv run pytest -m network  # 3 more against NCBI; downloads ~1.2 GB
uv run ruff check . && uv run ruff format --check .
```

Network tests are deselected by default and excluded from CI: they are rate-limited
per IP and move more than a gigabyte. Run them when changing anything about fetching,
and before trusting a build.

`hindsight refresh-releases` regenerates the vendored release table from NCBI's release
notes. Only needed when a new release lands.
