# Day 2: the samples

The database side of this benchmark is built by `hindsight build`. This is the other
half: the reads the detection arms will be run against, and the dilution series that
turns "did it find it" into "down to what abundance".

```console
$ scripts/day2.sh
```

That fetches the index sample and the negatives, QCs them, filters host reads, labels
which fragments are really viral, and writes the frozen dilution series — resumable,
so an interrupted run continues where it stopped. Everything it decided lands in
`work/day2-provenance.json`.

## What it produces

| | |
| --- | --- |
| `work/samples/<run>/raw/` | FASTQ as fetched, size- and MD5-verified |
| `work/samples/<run>/qc/` | after fastp, with `fastp.json` |
| `work/samples/<run>/nonhost/` | fragments where neither mate aligned to GRCh38 |
| `work/samples/<run>/truth/viral_fragment_ids.txt` | index sample only: the ground-truth viral fragments |
| `work/series/<run>_viral<NNNNNN>_{1,2}.fastq.gz` | the dilution series |
| `work/day2-provenance.json` | tool versions, parameters, seeds, counts, digests, timings, machine |

Days 3 and 4 consume `work/series/` and the negatives' `nonhost/`, and nothing else.

## The samples

**The index sample** is `SRR10971381`: the bronchoalveolar lavage fluid
metatranscriptome from Wu et al. 2020 that the first SARS-CoV-2 genome was assembled
from. 28.3 M read pairs, 8.0 Gbases, 2 × 151 bp, Illumina.

**The negatives** come from `PRJNA699976` (Shi et al. 2022, *PLoS Pathogens*
18(2):e1010259) — 408 patient BALF metatranscriptomes collected May 2016 to December
2017 at The Central Hospital of Wuhan, the same hospital and sampling protocol as the
index sample, plus 28 healthy-volunteer and 4 blank controls. The benchmark uses 20 of
them. Two properties of that deposit have to be stated wherever its false-alarm rate
is:

* **The reads are host-depleted.** The authors deposited non-human reads only, by a
  method the submission does not record. So a false-alarm rate per *raw sequenced*
  read is not available for these samples, only one per non-host read — which is the
  denominator the detection arms work in anyway, but it is not the same number, and
  the negatives' composition is not directly comparable with the index sample's.
  The pipeline still runs its own host filter over them and records how much it
  removes, which is a cheap check on how aggressive the submitters' depletion was.
* **There are no per-run collection dates.** ENA reports `collection_date` as
  `missing` for all 440 runs, and the supplementary tables do not map patient to run.
  The 2016-17 window is the paper's, at cohort level. It cannot be verified per run,
  and the runs cannot be ordered or stratified by date.

### Why these 20

Of the 408 patient runs, 175 were deposited as a `_1`/`_2` pair and 233 as a single
merged FASTQ. The eligible pool is the 175, so that one read layout runs through the
whole pipeline; the restriction is a stated criterion rather than an accident of which
runs happened to work.

From those, 20 are chosen by ranking the accessions by `sha256("<seed>:<accession>")`
and taking the lowest. Picking 20 "representative" samples by eye is how a negative
set ends up selected for being easy. The seed is `20191201`, the chosen list is
pinned in `config/day2.toml`, and a test asserts that the pinned list is exactly what
the seed selects — so the cohort cannot drift without the test failing.

```console
$ hindsight samples select-negatives --manifest config/negatives-PRJNA699976.tsv
```

### The depth cap, and why the negatives have one

At native depth the 20 negatives are 18 GiB of FASTQ and 120 million read pairs of
alignment. `config/day2.toml` caps each at 3 million pairs after QC, which also
equalises what each negative contributes: uncapped, a 24-million-read sample would
outweigh a 1.5-million-read one eightfold in a false-alarm rate that is supposed to be
per sample.

The cap is applied by a seeded streaming hash — a fragment is kept when the leading
bytes of `sha256("<seed>:<id>")` fall below the target fraction — so it needs one pass
and no memory, both mates independently keep the same fragments, and a deeper cap at
the same seed is a superset of a shallower one. The count that results is binomial
around the target rather than exact, so the actual number kept is what gets recorded.

The honest cost: the cap lowers sensitivity to genuinely rare background hits, which
is exactly what a false-alarm rate is about. Day 5 should say so, and `max_pairs = 0`
removes the cap on a machine with the disk and the hours for it.

Once a negative's host-filtered pair exists, its raw, QC'd and capped FASTQ are
deleted (`prune_intermediates`). Their digests are recorded first, and the `.sra`
container plus the pinned parameters regenerate them. The index sample is never
pruned, because the dilution series is cut from its QC'd reads.

## The dilution series

The series holds the index sample's **own** SARS-CoV-2 reads, thinned, in the index
sample's **own** background. Nothing is spiked in from elsewhere: the only variable
along the series is viral abundance.

Levels are 50,000 / 10,000 / 1,000 / 100 / 10 viral fragments. The top level sits
below the sample's own viral count, so every rung is a subsample rather than an
extrapolation.

Two design choices are worth defending:

**The selection is a SHA-256 ranking, not a random sample.** `random.sample` is
reproducible for a given CPython, which is not the same thing as reproducible. Here
each fragment identifier is scored by `sha256(f"{seed}:{id}")` and the lowest scores
win, so any implementation in any language reproduces the set from the seed and the
identifier list alone.

**The levels are nested.** They are prefixes of one ranking, so the 10-fragment level
is a subset of the 100-fragment level and so on up. A curve built on nested levels
measures sensitivity falling away as reads are removed; built on independent draws it
also measures the luck of each draw, and at 10 fragments that noise would be the whole
signal. `plan_series(independent_seeds=True)` exists for anyone who wants to quantify
that variance deliberately.

Each level's record carries a `selection_sha256` — the digest of its chosen
identifiers — so two builds can be compared without diffing gigabytes of FASTQ.

## Ground truth, and why it is not leakage

Fragments are labelled viral by aligning to `NC_045512.2`, the modern SARS-CoV-2
reference. That is a 2020 record, used deliberately, and only to decide what the right
answer is. The detection arms see only the November 2019 database `hindsight build`
produces; they never see this reference, this alignment, or the resulting label.
Knowing the answer is the benchmark's job. Knowing it is what the arms are not
allowed to do.

The label is computed **before** host filtering, and the overlap is reported as
`viral_lost_to_host`: how many known-viral fragments the host filter would have
discarded. Host filtering is the step where a careless pipeline quietly loses the
thing it is looking for, and this benchmark would otherwise charge that loss to the
classifier.

## Fetching: two sources, chosen by disk

ENA's ready-made FASTQ is the convenient source — no toolkit, published MD5 per file.
It is not always a fast one. From the machine this was first built on,
`ftp.sra.ebi.ac.uk` delivered about 1 MB/s with the transfer dropped every few
minutes; SRA's Open Data mirror (`sra-pub-run-odp.s3.amazonaws.com`, anonymous)
delivered 24 MB/s. The 2.7 GiB container for the index sample arrived in 2 minutes
there against an estimated 40 for the same data as FASTQ.

The mirror serves `.sra` containers, so the fast path costs a `fasterq-dump`
conversion, and the conversion wants scratch space — about 14 times the container's
size, measured, between the uncompressed pair and `fasterq-dump`'s temporary files.
`--stdout` does not avoid it: piping the interleaved stream into `pigz` looks like it
should stream, but `fasterq-dump` buffers through its temporary directory first, which
filled 20 GB in ten seconds and then died with "storage exhausted" having emitted
nothing.

So the source is chosen per run: the mirror where the conversion fits, ENA where it
does not. Which one served a run is recorded, because the two differ in read *names*
— ENA keeps the submitter's, `fasterq-dump` emits `<accession>.<spot>`. Fragment
identifiers are only ever compared within one sample, so the benchmark is indifferent,
but a reader diffing two builds is not.

One trap worth recording: request the mirror virtual-host style
(`https://sra-pub-run-odp.s3.amazonaws.com/sra/<acc>/<acc>`). The path-style endpoint
`s3.us-east-1.amazonaws.com/sra-pub-run-odp/...` invites any ambient AWS credentials
into a request to a public bucket, which is then refused with `AccessDenied` — a
failure that reads as "the data is gone" rather than "your request was signed by
mistake".

## What it costs

Measured on the machine that produced the committed provenance record — see
`machine` in `work/day2-provenance.json` for the exact CPU, memory and instance type,
which the pipeline records from the host rather than from configuration.

| step | cost |
| --- | --- |
| index sample, fetch | 2.4 GB, ENA at ~1 MB/s on this link |
| negatives, fetch | 20 containers, ~3 GB, minutes from the mirror |
| conversion | ~25 s per negative; needs 14× the container in free disk |
| fastp | a few minutes per sample |
| host filter | the dominant cost; see the provenance record's per-stage timings |
| peak disk | about 20 GB with the human index in place |

The human bowtie2 index is 4.1 GB and is fetched prebuilt by
`scripts/fetch-host-index.sh`. Building it locally instead costs a couple of hours of
CPU for a byte-identical result.

## Pinned versions

`config/day2.toml` pins fastp 0.23.4, bowtie2 2.5.2 and samtools 1.19.2, and the
pipeline probes each tool for its actual version at run time. A mismatch is recorded
in the provenance record as `tool_pin_drift` rather than enforced: trimming and
alignment both change between versions, and a detection limit that moves because the
trimmer changed is a result about the trimmer, which the Day 7 reproducibility pass
needs to be able to see.
