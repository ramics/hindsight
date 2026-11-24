# hindsight

hindsight builds a viral reference database as it stood on a date you choose, with every record
at the accession version that was public on that date.

```console
$ hindsight build --date 2019-12-01 --out db/november-2019
cutoff 2019-12-01 -> release 97 (2019-11-04); first excluded release: release 98 (2020-01-06)
manifest: 12179 nucleotide + 373129 protein records across 9305 taxa
nucleotide: wrote 12179 records to viral.fna (0 withdrawn, 0 version mismatches)
protein: wrote 373129 records to viral.faa (0 withdrawn, 0 version mismatches)
taxonomy: taxdmp snapshot 2019-11-01
```

It is built for retrospective detection benchmarks. The motivating question is whether an
automated, pathogen-agnostic pipeline would have flagged SARS-CoV-2 in the index sample using
only what was known before the virus was described. Answering it requires every reference the
pipeline consults to predate the discovery. Current RefSeq has changed too much since 2019 to
serve as that reference, even with the SARS-CoV-2 records removed.

## Why the common shortcuts leak

- **Removing SARS-CoV-2 from a current database.** The genome goes and its close relatives
  stay. RaTG13 (deposited 27 January 2020), the pangolin coronaviruses and thousands of
  pandemic-era sequences remain, and those are the records that make the index sample easy to
  classify.
- **Filtering a current database by deposition date.** Records are revised after deposition:
  sequences are corrected, annotation is rebuilt, taxonomy is reassigned. A record deposited in
  2007 may now sit at version 3, and version 3 may postdate the cutoff. The taxonomy tree has
  changed as well, and a classifier reports its result as a taxon.
- **Using a snapshot someone kept.** A kept snapshot covers one date, and only its owner can
  reproduce it.

## How a build works

NCBI archives every RefSeq release catalog back to release 1 (June 2003). Each catalog lists
every accession at the version it held in that release, and the viral set is a filter on one
column:

```
694009  Severe acute respiratory syndrome-related coronavirus  NC_004718.3  complete|viral  na  29751
```

Given a cutoff date, hindsight:

1. selects the last release issued strictly before the cutoff;
2. extracts the viral accession set from that release's catalog, with versions;
3. fetches each record at that exact version from E-utilities;
4. pairs the records with the taxonomy dump archived nearest before the release date;
5. writes a `provenance.json` recording the cutoff, the selected release, the first excluded
   release, the taxonomy snapshot, every count, every withdrawn record and a SHA-256 of every
   output.

A cutoff of 1 December 2019 selects release 97:

| | |
| --- | --- |
| Viral nucleotide / protein records | 12,179 / 373,129 |
| Distinct taxa | 9,305 |
| Coronavirus genomes | 45 |
| SARS-CoV | present, `NC_004718.3` |
| Bat SARS-related and other bat CoV genomes | 14, incl. `NC_014470.1` (BM48-31/BGR/2008), `NC_025217.1` (Bat Hp-betacoronavirus/Zhejiang2013) |
| SARS-CoV-2 (`NC_045512`, `MN908947`), RaTG13 (`MN996532`) | absent |

The database holds the nearest relatives available to an analyst in December 2019, at 79–88%
identity to the index sample; these are the relatives Wu et al. had. Pandemic-era records are
absent because release 97 was issued before they existed. `tests/test_leakage.py` asserts each
of these properties, and its network-marked test streams the full 385,308-row release 97 viral
set from NCBI to check them against the complete set.

## Install

```bash
uv sync --locked --all-extras     # or: pip install -e ".[dev]"
```

Set `NCBI_API_KEY` if you have one ([free from your NCBI
account](https://www.ncbi.nlm.nih.gov/account/settings/)). It raises the E-utilities rate limit
from 3 to 10 requests per second, which brings a full build from hours to under one hour.

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

An interrupted build resumes when you re-run the same command. It restarts from the last
completed shard, downloads resume mid-file with HTTP Range, and fetched records are cached by
shard with a digest. Builds are deterministic (sorted manifest, manifest-order emission, fixed
FASTA wrap), so two builds of one cutoff are byte-identical.

```
db/nov2019/
  manifest.tsv        every accession in the release's viral set, with version, molecule,
                      taxid, organism, RefSeq status and length
  viral.fna           nucleotide records, at their release-97 versions
  viral.faa           protein records, for a translated-search arm
  seqid2taxid.map     accession.version -> taxid, in Kraken2's format
  taxonomy/           nodes.dmp, names.dmp, merged.dmp, delnodes.dmp, as of 2019-11-01
  withdrawn.tsv       accessions the release listed that NCBI has since stopped serving
  provenance.json     the audit record, with SHA-256 of every file above
```

## Samples

`hindsight samples` prepares the read side of the benchmark:

- the index sample, `SRR10971381` (the Wu et al. BALF metatranscriptome);
- 20 time-matched pre-COVID Wuhan BALF negatives from `PRJNA699976`;
- read QC and host filtering over both;
- ground-truth labels marking which fragments are SARS-CoV-2;
- a frozen dilution series: the sample's own viral reads thinned to 50,000 / 10,000 / 1,000 /
  100 / 10 fragments in its own background.

Every stage resumes and every seed is recorded. The labels are the only component that uses
knowledge from 2020; the detection arms see only the November 2019 database.
[`docs/day2-samples.md`](docs/day2-samples.md) has the full account, including two properties
of the negative set that affect every number derived from it.

## Detection arms

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

Kraken2 at default settings, run against a database limited to release 97 (4 November 2019),
placed 25,928 of the index sample's viral reads in Coronaviridae. That database predates
SARS-CoV-2 and every relative discovered since, so the assignments land on the family and on
the older relatives it contains.

Kraken2 classifies each fragment independently, so the detection rate holds at about 41% across
the dilution series while the absolute count falls. At 0.77 viral reads per million, four
coronavirus reads remain. Whether four reads count as a detection depends on the false-alarm
rate, which the negatives measure.

This arm also produced an unplanned finding that makes the case for dating the taxonomy along
with the sequences. 4,922 of the reads were assigned to *Bat coronavirus BM48-31/BGR/2008*, the
nearest relative in the database. The November 2019 taxonomy files that genome under
"unclassified coronaviruses"; the current tree places it in Sarbecovirus. Family-level
detection is therefore 25,928 reads and genus-level detection is 20,418. A benchmark that dates
its sequences and keeps a current taxonomy hides that gap.
[`docs/day3-arm-a.md`](docs/day3-arm-a.md) describes the arm in full.

## Limitations

Read [`docs/limitations.md`](docs/limitations.md) before quoting a number or writing a methods
section. It covers:

- date resolution, which is one RefSeq release (coarser than a day);
- coverage, which is RefSeq only (GenBank and nr/nt are out of scope);
- how withdrawn records and version substitution are handled;
- the lag on monthly taxonomy snapshots;
- where NCBI's own products disagree.

[`docs/archive-notes.md`](docs/archive-notes.md) collects the pitfalls of the RefSeq archive:
non-contiguous release numbering, misleading file timestamps, typos in the release notes, and
E-utilities signalling a rate-limit rejection with HTTP 200. Each one is encoded in a test.

## Development

```bash
uv run pytest             # 243 offline tests
uv run pytest -m network  # 3 more against NCBI; downloads ~1.2 GB
uv run ruff check . && uv run ruff format --check .
```

Network tests are deselected by default and excluded from CI because they are rate-limited per
IP and move more than a gigabyte. Run them when changing anything about fetching, and before
trusting a build. `hindsight refresh-releases` regenerates the vendored release table from
NCBI's release notes; run it when a new release lands.
