# Day 3: Arm A, Kraken2 against November 2019

The benchmark's first detection arm. Exact k-mer classification, at Kraken2's own
default thresholds, against a database that contains everything a December 2019
analyst could have consulted and nothing they could not.

```console
$ scripts/day3.sh
```

Everything it decided lands in [`results/day3-provenance.json`](../results/day3-provenance.json),
and `hindsight arm-a report` prints the table below from it.

## The headline

**The arm placed 25,928 of the index sample's 62,117 SARS-CoV-2 fragments in
Coronaviridae — 41.7% — using a database that had never seen SARS-CoV-2, RaTG13, or
anything else learned after 4 November 2019.**

That is not a marginal call. 25,928 coronavirus reads in a clinical
metatranscriptome is an alarm no triage process would miss, and the nearest relatives
the database holds are SARS-CoV and the 2008–2013 bat sarbecoviruses at 79–88%
nucleotide identity. The reads were recognisable as coronaviral from what was already
in the reference set.

Two numbers rather than one, because the clade that counts as a detection is not a
free choice here — see [the taxonomy effect](#the-2019-taxonomy-moves-the-answer):

| clade the read was placed in | fragments | of 62,117 |
| --- | --- | --- |
| Coronaviridae (family) | 25,928 | 41.7% |
| Orthocoronavirinae (subfamily) | 20,657 | 33.3% |
| Betacoronavirus (genus) | 20,418 | 32.9% |
| Sarbecovirus (subgenus) | 20,407 | 32.9% |
| SARS-related coronavirus (species) | 20,407 | 32.9% |

The remaining 58% of the labelled viral fragments were not classified at all. Against
a 12,179-genome viral-only database that is the expected outcome rather than a failure:
a 151 bp read from a virus 12–21% divergent from anything in the reference set often
shares no 35-mer with it, and Kraken2 matches exact k-mers. 41.7% is the fraction that
did.

## The detection curve

The dilution series holds the sample's own viral reads, thinned, in its own background
(see [day2-samples.md](day2-samples.md)). Per-million figures are against post-QC
pairs, the headline denominator Day 2 fixed.

| level | viral fragments | per million | Coronaviridae | Betacoronavirus | sensitivity |
| --- | --- | --- | --- | --- | --- |
| viral050000 | 50,000 | 3,853 | 20,892 | 16,412 | 41.8% |
| viral010000 | 10,000 | 771 | 4,232 | 3,315 | 42.3% |
| viral001000 | 1,000 | 77 | 420 | 326 | 42.0% |
| viral000100 | 100 | 7.7 | 39 | 30 | 39.0% |
| viral000010 | 10 | 0.77 | 4 | 3 | 40.0% |

The detection *rate* is flat at about 41% across four orders of magnitude, which is
what it should be: Kraken2 classifies each fragment independently, so thinning the
viral reads scales the signal and changes nothing else. What falls is the absolute
count, and that is what decides detectability — four coronavirus reads at 0.77 per
million is a detection only if four reads beat the false-alarm rate. Drawing that line
is Day 5's job; this day's contribution is the two columns it needs.

## The 2019 taxonomy moves the answer

The gap between the family column and the genus column — 25,928 against 20,418 — is
not noise. 4,922 fragments were assigned to *Bat coronavirus BM48-31/BGR/2008*
(`NC_014470.1`, taxid 864596), the nearest relative in the database, and **the November
2019 taxonomy files that genome under "unclassified coronaviruses" in Coronaviridae,
not under Betacoronavirus**, where today's tree puts it.

So a 2019 analyst running this arm and filtering the report to Betacoronavirus would
have seen 20,418 reads. One looking at Coronaviridae would have seen 25,928. The
difference is entirely an artefact of where the contemporaneous tree had put a genome
whose relatives had not yet been reorganised — and it is invisible to any benchmark
that uses a modern taxonomy with a dated sequence set, because the modern tree has
since moved BM48-31 into Sarbecovirus and would have scored those 4,922 reads as
genus-level hits.

This is the argument for `hindsight`'s whole approach, found rather than assumed:
dating the sequences is not enough, because **a classifier's output is a taxon**, and
taxa move. Day 5 should draw its threshold at Coronaviridae, and say why.

## The false-alarm rate, at native depth

<!-- negatives -->

## What the negatives cost, and why they are not capped

Day 2 capped each negative at 6 million post-QC pairs because 17 GiB of host-filtered
FASTQ did not fit beside the human index on a machine with ~21 GiB of usable disk.
That cap degrades precisely the quantity a false-alarm rate is about: rare background
hits are the thing being counted, and a cap removes reads at random.

Day 3 removes it, and pays for it by coupling the two days rather than by renting a
larger machine. Each negative is fetched, QC'd, host-filtered, classified, and then its
reads are **deleted** — what persists is the Kraken2 report and the classified lines,
which are megabytes. Peak disk is one negative's worth instead of twenty.

```console
$ hindsight arm-a run --stage negatives            # all twenty, in the pinned order
$ hindsight arm-a run --only SRR13677553           # or one at a time
```

The run is resumable at sample granularity, which it has to be: a sample already
classified is skipped without re-fetching the gigabytes it came from.

## What is kept, and what is thrown away

Per classified sample:

| | |
| --- | --- |
| `<name>.kreport` | Kraken2's standard report, against the dated taxonomy |
| `<name>.classified.tsv.gz` | every **classified** read, with its k-mer hit map |
| `<name>.kraken2.log` | the tool's own summary |
| `.<name>.done` | the marker that makes a re-run cheap |

Unclassified reads are dropped on the way out. Against a viral-only database 99.7% of a
clinical metatranscriptome is unclassified, so the full per-read output is gigabytes of
the letter `U`, and the report already carries the count.

Keeping the classified reads *with their hit maps* is deliberate, and it is what lets
Day 5 re-draw a threshold without re-running anything: raising `--confidence` can only
unclassify a read, never classify a new one, so this run's output is a superset of what
any stricter confidence would have produced. Lowering `--minimum-hit-groups` below its
default of 2 is **not** derivable that way and would need a fresh run.

## The composition shortcut, and the check on it

Five dilution levels share one background, so classifying five materialised levels
means classifying the same 10.86 M background pairs five times. Instead the background
is classified once, each level's viral part once, and the level's result is the sum.

That is exact because Kraken2 assigns each fragment from the fragment alone — no
abundance estimation, no normalisation, no dependence on the rest of the file. But
"reads are independent" is the kind of assumption that is true until a tool has a
batching bug, and every dilution number here rests on it, so the run checks it rather
than asserting it: `viral000010` is materialised, classified whole, and compared with
its composed version taxon by taxon.

```
pairs composed       10,860,748
pairs direct         10,860,748
classified composed       4,124
classified direct         4,124
taxa disagreeing              0
```

Identical. The check costs one extra full classification — about a minute — and
`verify_composition = false` turns it off for anyone who does not want to pay it.

## The 78 sarbecovirus reads in the background, which are not false positives

The series' background is the non-host set with every ground-truth viral fragment
removed, so nothing in it should look like a sarbecovirus. 78 fragments do: Kraken2
assigns them to SARS-related coronavirus.

Two explanations fit and they have opposite consequences — either the classifier is
wrong about those reads, which would mean the arm has a false-positive rate on its own
sample, or the *label* is wrong about them, which would mean the sample holds more viral
reads than the benchmark credits it with. So the question the label was asked is asked
again with the threshold removed:

```console
$ hindsight arm-a diagnose-background
71 of 78 fragments align to NC_045512.2, every mate below the label's MAPQ >= 20
cut-off.
```

71 of the 78 align concordantly to the modern SARS-CoV-2 genome, all of them at more
than one position, and the highest mapping quality any mate achieves is 18 — below the
label's cut-off of 20. These are SARS-CoV-2 reads from repetitive stretches of the
genome that the ground-truth label declined to count, and Kraken2 found them because a
k-mer match needs no unique position.

Three consequences, all of which belong in the write-up:

* The arm has **no measured false-positive rate on the index sample's own background**.
  The 20 independent negatives are the only honest estimate, which is what they are for.
* The sample's viral content is marginally larger than 62,117 fragments. The label is
  conservative, and in the direction that understates the arm rather than flattering it.
* `min_truth_mapq = 20` is a defensible choice that is nevertheless doing visible work.
  Day 7 could report the curve against mapping quality instead of one threshold.

## Cost

<!-- cost -->

## For Day 4, Day 5 and Day 7

* **Day 4** (MEGAHIT + DIAMOND) reads the same inputs: `work/series/parts/` and each
  negative's `nonhost/` pair. Note that Day 3 *deletes* the negatives' reads as it
  goes, so Arm B either re-stages them the same way or the two arms share one pass.
  Sharing one pass is the cheaper design and the one this day's structure was built
  for: `classify_negatives` is a loop with one classifier in it.
* **Day 5** should draw the threshold at **Coronaviridae**, for the taxonomy reason
  above, and should take the per-clade columns rather than the headline.
* **Day 7** items this day adds: `bowtie2 --reorder` for byte-identical host-filter
  output (see below), and the mapping-quality curve for the truth label.

## Reproducing Day 2 from the pinned config

This day began by re-staging Day 2 in a container that had none of its outputs — a
reproducibility test nobody designed, and the useful kind. Every count came back
identical:

| | Day 2 | re-staged |
| --- | --- | --- |
| raw spots | 28,282,964 | 28,282,964 |
| redacted | 13,934,971 | 13,934,971 |
| post-QC pairs | 12,975,578 | 12,975,578 |
| non-host pairs | 10,921,957 | 10,921,957 |
| viral fragments | 62,117 | 62,117 |
| viral lost to the host filter | 898 | 898 |
| background fragments | 10,860,738 | 10,860,738 |

The QC output and the 62,117-fragment viral label are byte-identical. The
host-filtered FASTQ and the series parts cut from them are **not**: same reads,
different order, because bowtie2 emits alignments in whatever order its worker threads
finish in unless `--reorder` is given. A SHA-256 cannot tell that apart from a genuinely
different read set, so `reads.unordered_id_digest` adds a digest that depends on the
fragments and not their order — a sum of per-identifier hashes, which is commutative
and needs no sort.

Adding `--reorder` to the host filter would make the whole pipeline byte-reproducible
at a small throughput cost. That is a Day 7 decision, not a Day 3 one, and nothing in
this day's results depends on read order: Kraken2 classifies each fragment
independently, which is the same property the composition shortcut rests on.
