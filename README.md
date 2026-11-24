# hindsight

Build the viral portion of RefSeq as it stood on a date you choose, with every accession at
the version it held on that date.

```console
$ hindsight build --date 2019-12-01 --out db/november-2019
cutoff 2019-12-01 -> release 97 (2019-11-04); first excluded release: release 98 (2020-01-06)
manifest: 12179 nucleotide + 373129 protein records across 9305 taxa
nucleotide: wrote 12179 records to viral.fna (0 withdrawn, 0 version mismatches)
protein: wrote 373129 records to viral.faa (0 withdrawn, 0 version mismatches)
taxonomy: taxdmp snapshot 2019-11-01
```

## The question

Would a pathogen-agnostic pipeline have flagged SARS-CoV-2 in the index sample, working only
from what was knowable before the virus was described? Answering that requires every reference
the pipeline consults to predate the discovery. Current RefSeq cannot supply such a reference
set, and the three obvious ways to approximate one each fail differently.

**Deleting SARS-CoV-2 from a current database** removes the genome and keeps its relatives.
RaTG13 was deposited on 27 January 2020; the pangolin coronaviruses and thousands of
pandemic-era submissions followed. Those are the records that make the index sample easy, and
they remain after the deletion.

**Filtering a current database by deposition date** keeps records whose contents have since
moved. Accessions are revised — sequence corrected, annotation rebuilt, taxonomy reassigned —
so a record deposited in 2007 may now sit at version 3, and version 3 may postdate the cutoff.
The taxonomic tree has been rearranged over the same period, which matters because a
classifier reports a taxon.

**Reusing a snapshot someone archived** covers whichever date they happened to save, and only
they can regenerate it.

hindsight takes the route NCBI leaves open. Release catalogs are archived back to release 1
(June 2003), and each one lists every accession at the version it held in that release. The
viral set is a filter on one column:

```
694009  Severe acute respiratory syndrome-related coronavirus  NC_004718.3  complete|viral  na  29751
```

So hindsight selects the last release issued strictly before the cutoff, extracts the viral
accession set with its versions, fetches each record at that exact version from E-utilities,
pairs it with the taxonomy dump archived nearest before the release date, and writes a
`provenance.json` recording the cutoff, the selected release, the first release excluded, the
taxonomy snapshot, every count, every withdrawn accession, and a SHA-256 of every output.

A 1 December 2019 cutoff selects release 97:

| | |
| --- | --- |
| Viral nucleotide / protein records | 12,179 / 373,129 |
| Distinct taxa | 9,305 |
| Coronavirus genomes | 45 |
| SARS-CoV | present, `NC_004718.3` |
| Bat SARS-related and other bat CoV genomes | 14, incl. `NC_014470.1` (BM48-31/BGR/2008), `NC_025217.1` (Bat Hp-betacoronavirus/Zhejiang2013) |
| SARS-CoV-2 (`NC_045512`, `MN908947`), RaTG13 (`MN996532`) | absent |

The closest relatives available to a December 2019 analyst are present — SARS-CoV and the
2008–2013 bat sarbecoviruses, at 79–88% nucleotide identity to the index sample, which is the
neighbourhood Wu et al. were working in. The pandemic-era records are absent because they had
yet to be written, so nothing had to be deleted to arrive at this set. `tests/test_leakage.py`
asserts each line of the table, and its network-marked test streams the full 385,308-row
release 97 viral set from NCBI to check the whole thing.

## Install

```bash
uv sync --locked --all-extras     # or: pip install -e ".[dev]"
```

Set `NCBI_API_KEY` if you have one ([free from your NCBI
account](https://www.ncbi.nlm.nih.gov/account/settings/)). It raises the E-utilities request
rate from 2/s to 8/s, which takes a full build from hours to under one.

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
downloads continue mid-file over HTTP Range, and fetched records are cached per shard with a
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

## The samples

`hindsight samples` prepares the other half of the benchmark: the index sample
(`SRR10971381`, the Wu et al. BALF metatranscriptome), 20 time-matched pre-COVID Wuhan BALF
negatives from `PRJNA699976`, read QC and host filtering across both, ground-truth labels for
which fragments are genuinely SARS-CoV-2, and a frozen dilution series — the sample's own viral
reads thinned to 50,000 / 10,000 / 1,000 / 100 / 10 fragments in its own background. Every
stage resumes and every seed is recorded.

The labels are the one component permitted to know something from 2020; they come from aligning
to `NC_045512.2`. The detection arms see the November 2019 database and the reads, and the
labels reach them only as a score at the end. [`docs/day2-samples.md`](docs/day2-samples.md)
gives the full account, including two properties of the negative set that travel with every
number derived from it.

## Arm A: Kraken2 against November 2019

An arm is a detection pipeline run against the dated database and the reads. The first uses
Kraken2 at its own default thresholds.

```console
$ scripts/day3.sh
index sample: 25,928 of 62,117 ground-truth viral fragments in Coronaviridae (41.7%)
  viral050000: 20,892 of 50,000 detected
  viral010000:  4,232 of 10,000 detected
  viral001000:    420 of  1,000 detected
  viral000100:     39 of    100 detected
  viral000010:      4 of     10 detected
```

A standard read classifier, at default settings, against a database that had never seen
SARS-CoV-2 or any relative described after 4 November 2019, placed 25,928 of the index sample's
labelled viral fragments in Coronaviridae. It could not name the virus, which had no name, and
it did not need to: 25,928 coronavirus reads in a clinical metatranscriptome is a signal a
triage process acts on.

Where in the tree those reads land depends on the clade you read the report at, and the spread
is wide enough to matter:

| clade the read was placed in | fragments | of 62,117 |
| --- | --- | --- |
| Coronaviridae (family) | 25,928 | 41.7% |
| Orthocoronavirinae (subfamily) | 20,657 | 33.3% |
| Betacoronavirus (genus) | 20,418 | 32.9% |
| Sarbecovirus (subgenus) | 20,407 | 32.9% |
| SARS-related coronavirus (species) | 20,407 | 32.9% |

The 58% of labelled fragments that went unclassified are the expected outcome against a
12,179-genome viral-only database: a 151 bp read from a virus 12–21% divergent from anything in
the reference set frequently shares no exact 35-mer with it, and exact k-mers are what Kraken2
matches.

### The detection curve

The dilution series holds the sample's own viral reads, thinned, in its own background.
Per-million figures use post-QC pairs, the headline denominator.

| level | viral fragments | per million | Coronaviridae | Betacoronavirus | sensitivity |
| --- | --- | --- | --- | --- | --- |
| viral050000 | 50,000 | 3,853 | 20,892 | 16,412 | 41.8% |
| viral010000 | 10,000 | 771 | 4,232 | 3,315 | 42.3% |
| viral001000 | 1,000 | 77 | 420 | 326 | 42.0% |
| viral000100 | 100 | 7.7 | 39 | 30 | 39.0% |
| viral000010 | 10 | 0.77 | 4 | 3 | 40.0% |

The rate holds at roughly 41% across four orders of magnitude, as it should: Kraken2 assigns
each fragment from that fragment alone, so thinning the viral reads scales the signal and
leaves the per-read probability alone. The absolute count is what falls, and the count is what
decides detectability. Four coronavirus reads at 0.77 viral fragments per million counts as a
detection only against a false-alarm rate, which is what the negatives measure.

### The 2019 taxonomy moves the answer

The gap between the family row and the genus row — 25,928 against 20,418 — has a single cause.
4,922 fragments were assigned to *Bat coronavirus BM48-31/BGR/2008* (`NC_014470.1`, taxid
864596), the nearest relative the database holds, and the November 2019 taxonomy files that
genome under "unclassified coronaviruses" within Coronaviridae. Today's tree places it in
Sarbecovirus.

So a 2019 analyst filtering the report to Betacoronavirus would have seen 20,418 reads, and one
reading it at Coronaviridae would have seen 25,928. The difference comes entirely from where
the contemporaneous tree had put a genome whose relatives had yet to be reorganised. A
benchmark that dates its sequences and keeps a modern taxonomy scores those 4,922 reads as
genus-level hits and never sees the gap.

This is the argument for dating the taxonomy alongside the sequences, and it turned up in the
results rather than in the design: a classifier's output is a taxon, and taxa move.
[`docs/day3-arm-a.md`](docs/day3-arm-a.md) has the arm in full.

### The false-alarm rate

20 pre-COVID BALF samples from the same hospital, at native depth: 228,407,636 non-host pairs,
no cap.

| | reads | samples | per million non-host |
| --- | --- | --- | --- |
| Sarbecovirus, or below | 0 | 0 of 20 | 0 |
| Betacoronavirus, excluding endemic human coronaviruses | 7 | 1 of 20 | 0.031 |
| Coronaviridae, excluding endemic human coronaviruses | 35 | 6 of 20 | 0.153 |
| Betacoronavirus, all | 2,420 | 4 of 20 | 10.6 |
| Coronaviridae, all | 2,495 | 8 of 20 | 10.9 |

Zero reads in 228 million were placed in Sarbecovirus, which is where 20,407 of the index
sample's fragments go. The 2,420 betacoronavirus reads are mostly two patients' real
infections: 2,399 of them are *Human coronavirus HKU1* in `SRR13677922`, and nine more are
*OC43* in `SRR13677878`. Both are endemic cold viruses, both are named species in the November
2019 tree, and a read placed on one of them is the arm working correctly.

Seven reads across the cohort are betacoronaviral without being endemic human coronavirus, and
all seven sit in the same sample as the HKU1 infection, spread across HKU24, rabbit coronavirus
HKU14 and rat coronavirus Parker — embecoviruses, HKU1's neighbours in the reference set. They
are that infection's reads landing one taxon over. Widening to the whole family gives 35 reads
in six samples, including avian, bat and camel coronaviruses in human lung fluid, which are
false alarms on any reading.

Which rate to quote therefore follows from a decision about clade rather than from a
measurement, so the record carries all of them. `docs/day3-arm-a.md` sets out the case for
drawing the threshold at Coronaviridae and reporting 0.153 per million, with the sarbecovirus
zero as the figure a SARS-CoV-2-like novelty is actually measured against.

### Cost

4 cores, 15.7 GiB, no GPU. Measured on the run that produced
[`results/day3-provenance.json`](results/day3-provenance.json):

| | |
| --- | --- |
| Kraken2 database build | 19 min 33 s (12,178 sequences, 318 Mbp, 415 MB hash) |
| Classification, everything | 24.5 min for 250 M read pairs |
| Negatives: fetch, QC, host filter | 163 min, median 7.7 min a sample |
| Whole arm, 5 invocations | 191 min |
| Peak disk, largest single negative | 6.6 GiB |

Classification is the cheap part. The 17 GiB of host-filtered FASTQ that would not fit on this
machine existed to be read once by a classifier that took 23 minutes to read all of it, so the
arm classifies each negative as it comes off the host filter and then deletes the reads,
keeping the report and the classified lines. Peak disk is one negative rather than twenty,
which is how the cohort runs at native depth here at all.

## Before quoting a number

[`docs/limitations.md`](docs/limitations.md) states what this gives you and what it does not:
release granularity rather than day granularity (release 97 to 98 is a 63-day window, and
`hindsight resolve` prints it), RefSeq rather than GenBank or nr/nt, the handling of withdrawn
records and version substitution, the lag on monthly taxonomy snapshots, and the one place
NCBI's own products disagree with each other. Read it before writing a methods section.

[`docs/archive-notes.md`](docs/archive-notes.md) collects what is easy to get wrong about the
RefSeq archive — non-contiguous release numbering, misleading file timestamps, typos in the
release notes, and E-utilities signalling a rate-limit rejection with HTTP 200 — each one now
encoded in a test.

## Development

```bash
uv run pytest             # 243 offline tests
uv run pytest -m network  # 3 more against NCBI; downloads ~1.2 GB
uv run ruff check . && uv run ruff format --check .
```

Network tests are deselected by default and left out of CI: they are rate-limited per IP and
move more than a gigabyte. Run them when changing anything about fetching, and before trusting
a build. `hindsight refresh-releases` regenerates the vendored release table from NCBI's
release notes, which is needed when a new release lands.
