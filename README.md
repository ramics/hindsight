# hindsight

Build a viral reference database as it stood on a date you choose, at the accession versions
that were public at the time.

```console
$ hindsight build --date 2019-12-01 --out db/november-2019
cutoff 2019-12-01 -> release 97 (2019-11-04); first excluded release: release 98 (2020-01-06)
manifest: 12179 nucleotide + 373129 protein records across 9305 taxa
nucleotide: wrote 12179 records to viral.fna (0 withdrawn, 0 version mismatches)
protein: wrote 373129 records to viral.faa (0 withdrawn, 0 version mismatches)
taxonomy: taxdmp snapshot 2019-11-01
```

hindsight exists to make one kind of benchmark honest: **would an automated,
pathogen-agnostic pipeline have flagged SARS-CoV-2 in the index sample, using only what was
knowable before it existed?** Every reference the pipeline consults has to predate the
discovery, and today's RefSeq cannot provide that, with or without the obvious records deleted.

## Why the easy approaches leak

- **Delete SARS-CoV-2 from today's database.** The genome goes; its relatives stay. RaTG13
  (deposited 27 January 2020), the pangolin coronaviruses and thousands of SARS-CoV-2-era
  sequences are the records that make the index sample easy. You have deleted the answer and
  left the answer key.
- **Filter today's database by deposition date.** The records themselves have moved. Accessions
  get revised — sequence corrected, annotation rebuilt, taxonomy reassigned — so a record
  deposited in 2007 may sit at version 3, and version 3 may postdate your cutoff. The shape of
  the taxonomic tree has changed too, and a classifier's output is a taxon.
- **Use a snapshot someone kept.** Available for one date, reproducible by nobody.

hindsight takes the route NCBI left open: **release catalogs are archived back to release 1
(June 2003)**, and each one lists every accession *at the version it held in that release*. The
viral set is a filter on one column:

```
694009  Severe acute respiratory syndrome-related coronavirus  NC_004718.3  complete|viral  na  29751
```

So hindsight picks the last release issued strictly before your cutoff, extracts the viral
accession set with versions, fetches each record at that exact version from E-utilities, pairs
it with the taxonomy dump archived nearest before the release date, and writes a
`provenance.json` recording the cutoff, the release, the release deliberately excluded, the
taxonomy snapshot, every count, every withdrawn record and a SHA-256 of every output.

December 2019 selects release 97:

| | |
| --- | --- |
| Viral nucleotide / protein records | 12,179 / 373,129 |
| Distinct taxa | 9,305 |
| Coronavirus genomes | 45 |
| SARS-CoV | present, `NC_004718.3` |
| Bat SARS-related and other bat CoV genomes | 14, incl. `NC_014470.1` (BM48-31/BGR/2008), `NC_025217.1` (Bat Hp-betacoronavirus/Zhejiang2013) |
| SARS-CoV-2 (`NC_045512`, `MN908947`), RaTG13 (`MN996532`) | absent |

The nearest relatives a December 2019 analyst could reach — 79–88% identity to the index
sample, which is what Wu et al. had — are in. Everything learned from the pandemic is out, and
nothing was deleted to make it so; those records had yet to be written. `tests/test_leakage.py`
asserts all of it, and its network-marked test streams the full 385,308-row release 97 viral
set from NCBI to check the whole thing.

## Install

```bash
uv sync --locked --all-extras     # or: pip install -e ".[dev]"
```

Set `NCBI_API_KEY` if you have one ([free from your NCBI
account](https://www.ncbi.nlm.nih.gov/account/settings/)). It lifts the E-utilities ceiling
from 3/s to 10/s and takes a full build from hours to under one.

## Use

```bash
hindsight resolve 2019-12-01                    # which release a date selects, and what it excludes
hindsight build --date 2019-12-01 --out db/nov2019 --manifest-only      # accessions only, no fetch
hindsight build --date 2019-12-01 --out db/nov2019 --molecules nucleotide
hindsight build --date 2019-12-01 --out db/nov2019                      # + protein + taxonomy
hindsight verify db/nov2019                     # re-check a build against its own provenance
hindsight releases --limit 10                   # archived releases
hindsight taxdumps                              # dated taxonomy snapshots
```

Interrupt a build and re-run the same command: it resumes from the last completed shard,
downloads resume mid-file with HTTP Range, and fetched records are cached by shard with a
digest. Builds are deterministic — sorted manifest, manifest-order emission, fixed FASTA wrap —
so two builds of one cutoff are byte-identical.

```
db/nov2019/
  manifest.tsv        every accession in the release's viral set, with version, molecule,
                      taxid, organism, RefSeq status and length
  viral.fna           nucleotide records, at their release-97 versions
  viral.faa           protein records, for a translated-search arm
  seqid2taxid.map     accession.version -> taxid, in Kraken2's format
  taxonomy/           nodes.dmp, names.dmp, merged.dmp, delnodes.dmp, as of 2019-11-01
  withdrawn.tsv       accessions the release listed that NCBI will no longer serve
  provenance.json     the audit record, with SHA-256 of every file above
```

## The samples and the arms

`hindsight samples` prepares the other half of the benchmark: the index sample (`SRR10971381`,
the Wu et al. BALF metatranscriptome), 20 time-matched pre-COVID Wuhan BALF negatives from
`PRJNA699976`, read QC and host filtering over both, ground-truth labels for which fragments
are really SARS-CoV-2, and a frozen dilution series — the sample's own viral reads thinned to
50,000 / 10,000 / 1,000 / 100 / 10 fragments in its own background. Every stage resumes and
every seed is recorded. The labels are the one thing allowed to know something from 2020; the
detection arms see only the November 2019 database.
[`docs/day2-samples.md`](docs/day2-samples.md) has the full account, including two properties
of the negative set that travel with every number derived from it.

An arm is a detection pipeline run against the database and the reads. The first is Kraken2:

```console
$ scripts/day3.sh
index sample: 25,928 of 62,117 ground-truth viral fragments in Coronaviridae (41.7%)
  viral050000: 20,892 of 50,000 detected
  viral010000:  4,232 of 10,000 detected
  viral001000:    420 of  1,000 detected
  viral000100:     39 of    100 detected
  viral000010:      4 of     10 detected
```

**A standard read classifier, at default settings, against a database that had never seen
SARS-CoV-2 or any relative discovered after 4 November 2019, placed 25,928 of the index
sample's viral reads in Coronaviridae.** It could not name the virus; the virus had no name.

Detection rate stays flat at about 41% down the series, since Kraken2 classifies each fragment
independently, so the absolute count is what falls: four coronavirus reads survive at 0.77
viral reads per million. Whether four reads counts as a detection depends on the false-alarm
rate, which the negatives measure.

One result from that arm is the argument for this whole tool, and it was found rather than
predicted. 4,922 of those reads were assigned to *Bat coronavirus BM48-31/BGR/2008*, the
nearest relative the database holds — and the **November 2019 taxonomy files that genome under
"unclassified coronaviruses"**, where today's tree has it in Sarbecovirus. Family-level
detection is 25,928 reads; genus-level is 20,418. A benchmark that dates its sequences and
keeps a modern taxonomy cannot see that gap.
[`docs/day3-arm-a.md`](docs/day3-arm-a.md) has the arm in full.

## Before you quote a number

[`docs/limitations.md`](docs/limitations.md) states what this does not give you: release
granularity rather than day granularity, RefSeq rather than GenBank or nr/nt, how withdrawn
records and version substitution are handled, the lag on monthly taxonomy snapshots, and where
NCBI's own products disagree. A benchmark whose database has unexamined properties is not yet a
benchmark, so read it before writing a methods section.

[`docs/archive-notes.md`](docs/archive-notes.md) collects what is easy to get wrong about the
RefSeq archive — non-contiguous release numbering, misleading file timestamps, typos in the
release notes, and E-utilities signalling a rate-limit rejection with HTTP 200 — each now
encoded in a test.

## Development

```bash
uv run pytest             # 243 offline tests
uv run pytest -m network  # 3 more against NCBI; downloads ~1.2 GB
uv run ruff check . && uv run ruff format --check .
```

Network tests are deselected by default and excluded from CI: they are rate-limited per IP and
move more than a gigabyte. Run them when changing anything about fetching, and before trusting
a build. `hindsight refresh-releases` regenerates the vendored release table from NCBI's
release notes, needed only when a new release lands.
