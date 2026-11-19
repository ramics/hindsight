# hindsight

Rebuilds the viral reference set as it stood on any past date, and benchmarks
retrospective pathogen detection against it without leakage.

## Abstract

Would an automated, pathogen-agnostic pipeline have flagged SARS-CoV-2 in the index
sample using only what was knowable before it existed, and down to what abundance? We
rebuilt the viral reference set as RefSeq release 97 (4 November 2019) held it — at the
accession versions public then, with the taxonomy of 2019-11-01. SARS-CoV-2, RaTG13 and
everything else learned since are absent; nothing was deleted to make them so, because
those records had not been written yet. Two pipelines were scored against that database:
Kraken2 at its own defaults, and MEGAHIT assembly followed by DIAMOND blastx. Both were
run on the Wu et al. bronchoalveolar lavage metatranscriptome (`SRR10971381`), on a
621-fold dilution series cut from it, and on twenty time-matched pre-COVID BALF
libraries from which both alarm thresholds were drawn. Both pipelines alarm on the index
sample and at every usable dose of the series, neither raises a false alarm in the
twenty negatives, and the assembly pipeline recovers 99.68% of the genome in a single
29,809 bp contig at 100% nucleotide identity to a reference nobody had yet deposited.
The limit of detection is approximately 8 labelled viral reads per million post-QC pairs
for both pipelines. That figure is set by the benchmark's own floor (§3), not by either
method.

## 1. Results

| | Kraken2 | MEGAHIT + DIAMOND |
| --- | --- | --- |
| index sample | **25,928 of 62,117 viral fragments** in Coronaviridae (41.7%) | **99.68% of the genome in one 29,809 bp contig at 100% identity** |
| detection statistic | reads placed in Coronaviridae | best bitscore of a Coronaviridae contig |
| alarm rule, from the negatives alone | ≥ 12 reads | ≥ 1 contig at bitscore ≥ 250 |
| false alarms, 20 pre-COVID BALF samples | **0 of 20** | **0 of 20** |
| limit of detection | **7.7 labelled viral reads per million** | **≤ 0.77 per million** |

The nearest protein in the dated database sits at 79.9% amino-acid identity to this
virus and the nearest genomes at 79–88% nucleotide identity, yet the assembly pipeline
reconstructed the genome rather than merely registering a signal.

![Two panels. A: each pipeline's evidence as a multiple of its own alarm threshold
across the dilution series, on log axes; both curves lie above the alarm line at every
dose except Kraken2 at the lowest, which is inside the shaded floor region. B: the share
of labelled viral fragments each pipeline recovers; MEGAHIT + DIAMOND falls from 98% to
10% while Kraken2 stays flat near 41%, and the two curves cross between the 10- and
100-fragment rungs.](docs/figures/detection-vs-dilution.svg)

**Figure 1. Detection and recovery across the dilution series.** The index sample's
labelled SARS-CoV-2 fragments, thinned to five doses in that sample's own background and
scored against a November 2019 reference frame. Because the series was diluted into the
index patient's own reads, which already carried ~90 coronavirus fragments of their own,
no rung can hold less virus than that residue: shading marks the doses below 6.9 labelled
reads per million, where most of what is present was never diluted and the rung measures
the benchmark rather than either pipeline.
**(A) Detection.** Each pipeline's evidence — Coronaviridae reads for Kraken2, best
Coronaviridae contig bitscore for MEGAHIT + DIAMOND — as a multiple of that pipeline's own
alarm threshold (§2.3). Both clear the alarm line at every usable dose. The single point
below it is Kraken2 at the lowest rung, inside the floor, where four labelled reads do not
reach its threshold of twelve while MEGAHIT + DIAMOND still alarms on one 462 bp contig at
bitscore 270.
**(B) Recovery.** The share of labelled fragments each pipeline finds. MEGAHIT + DIAMOND
recovers 98% while coverage suffices for assembly and 10% at the lowest dose; Kraken2
recovers a flat ~41% throughout, because it classifies each read independently, so the two
curves cross between the 10- and 100-fragment rungs. That crossing is in recovery and not
in detection: it does not change which pipeline alarms, and at the lowest dose the
pipeline that recovers less is the one that alarms.

The figure is drawn from [`results/day5-provenance.json`](results/day5-provenance.json)
by [`scripts/day6_figure.py`](scripts/day6_figure.py) and introduces no number of its
own; re-score with a different flag and re-run it, and the figure follows.

### 1.1 Detection and recovery come apart

Above ~1,000 fragments the assembly recovers more than twice what the classifier does,
because a contig recruits reads that would never classify alone; below ~100 it collapses
for want of coverage, while Kraken2's per-read rate stays flat. That is the whole of the
recovery crossing. Detection asks a different question — whether one piece of evidence
stands above the cohort's noise — and a 462 bp contig matching a coronavirus replicase
at bitscore 270 does, while four reads, in a cohort whose worst negative reaches eleven,
does not. The method that recovers less can be the method that detects, which is why the
two pipelines carry separate statistics rather than one shared threshold, and why a
single curve of "detection" against dose would be the wrong plot.

### 1.2 Why Kraken2 recovers 41% and not 95%

Kraken2 matches exact k-mers: at its defaults each read is classified from 35-mers
looked up through their 31 bp minimizers. The nearest relatives in the November 2019
database are 79–88% identical at nucleotide level — a mismatch every five to eight bases
on average — so a match survives only where a conserved stretch runs clean for a whole
minimizer. Leaving 58% of the labelled fragments unclassified is the expected cost of
that design against a novel organism, not a defect in the database or its date. It is
also the gap MEGAHIT + DIAMOND closes from the other side: assemble first, then match
protein, where homology remains visible at 30% amino-acid identity.

The exact-match length is a build-time setting, deliberately left alone.
`config/day3.toml` exposes `kmer_length` and `minimizer_length` under `[database]` and
leaves both unset, which is Kraken2's own 35 and 31, with the setting recorded either
way. Turning it is an experiment rather than a free improvement: a shorter k raises
sensitivity and background together, and because the threshold is set from the
negatives, recovering more viral reads does not by itself lower the limit of detection.
Testing it fairly costs a rebuild (~20 min), a reclassification (~25 min) and a
re-measurement of all twenty negatives, whose reads were deleted after classification —
about three hours per value of k. The one knob this benchmark did measure argues for
restraint: DIAMOND's `--very-sensitive` was pre-registered as necessary, matched stock
mode's best bitscore on every genuine signal, and differed only in false positives.

## 2. Methods

Five stages, each with its own provenance record and write-up.

| | what it does | write-up | record |
| --- | --- | --- | --- |
| **Database** | RefSeq release 97 at its own versions, with the 2019-11-01 taxonomy | §2.1 | `db/nov2019/provenance.json` |
| **Reads** | index sample + 20 negatives, QC, host filter, ground-truth label, frozen dilution series | [`docs/day2-samples.md`](docs/day2-samples.md) | [`results/day2-provenance.json`](results/day2-provenance.json) |
| **Kraken2** | a read classifier at its own default settings | [`docs/day3-arm-a.md`](docs/day3-arm-a.md) | [`results/day3-provenance.json`](results/day3-provenance.json) |
| **MEGAHIT + DIAMOND** | de novo assembly, then translated search | [`docs/day4-arm-b.md`](docs/day4-arm-b.md) | [`results/day4-provenance.json`](results/day4-provenance.json) |
| **Scoring** | thresholds from the negatives, limits, 128 operating points | [`docs/day5-scoring.md`](docs/day5-scoring.md) | [`results/day5-provenance.json`](results/day5-provenance.json) |

### 2.1 A dated reference database

The benchmark needs a database containing everything a real analyst could have consulted
in December 2019 and nothing they could not. The obvious routes fail.

- **Today's database with SARS-CoV-2 deleted.** RaTG13 (deposited 27 January 2020), the
  pangolin coronaviruses and thousands of SARS-CoV-2-era sequences remain, and they are
  the records that make the index sample easy to classify. This deletes the answer and
  leaves the answer key.
- **Today's database filtered by deposition date.** The records have moved. Accessions get
  revised — sequence corrected, annotation rebuilt, taxonomy reassigned — so a record
  first deposited in 2007 may now be at version 3, and version 3 may postdate the cutoff.
  Taxon names and the shape of the tree have changed too, and a classifier's output is a
  taxon.
- **A snapshot someone happened to keep.** Not reproducible, and not available for an
  arbitrary date.

hindsight takes a fourth route, which works because NCBI never threw anything away.
RefSeq release catalogs are archived back to release 1 (June 2003), and
`RefSeq-release<N>.catalog.gz` lists every accession in that release at the version it
held then. The viral set is a filter on one column:

```
694009  Severe acute respiratory syndrome-related coronavirus  NC_004718.3  complete|viral  na  29751
```

hindsight picks the last release issued strictly before the cutoff, extracts the viral
accession set with versions, fetches each record at that exact version from E-utilities,
pairs it with the taxonomy dump archived nearest before the release date, and writes a
provenance record with every count and a SHA-256 of every output. A 1 December 2019
cutoff selects release 97:

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

`tests/test_leakage.py` asserts these contents, and its network-marked test streams the
full 385,308-row release 97 viral set from NCBI and checks all of it.

Dating the sequences is not sufficient, because a classifier's output is a taxon. 4,922
of the index sample's reads were assigned to *Bat coronavirus BM48-31/BGR/2008*, the
nearest relative the database holds — and the November 2019 taxonomy files that genome
under "unclassified coronaviruses", not under Betacoronavirus where today's tree puts
it. Detection at family level is 25,928 reads; at genus level, 20,418. A benchmark that
dates its sequences but keeps a modern taxonomy cannot see that difference. It was found
rather than predicted.

### 2.2 Reads

`SRR10971381` (the Wu et al. BALF metatranscriptome) and 20 time-matched pre-COVID Wuhan
BALF negatives from `PRJNA699976`. Half the index sample's spots are NCBI's human-read
redaction, and 15.8% of what survives QC still aligns to GRCh38, so the host filter
earns its place: 12,975,578 post-QC pairs, 10,921,957 non-host, and 62,117 fragments
labelled as really SARS-CoV-2 by alignment to the modern reference, or 0.22% of raw
spots. (Wu et al. report 0.19% against the unredacted original; the two are not directly
comparable.) The dilution series thins the labelled fragments to 50,000, 10,000, 1,000,
100 and 10, holding them in the sample's own background with fixed seeds, stored as parts
rather than as five copies of a 21 GB library.

The ground-truth label is the one place the benchmark is allowed to know something from
2020. The pipelines never see it, and never see the modern reference; they see only the
November 2019 database.

### 2.3 Alarm thresholds

Both thresholds are drawn from the negatives and never from the positive, since a
threshold placed after looking at the answer is not a detection limit. The worst
negative gives Kraken2 11 Coronaviridae reads, so the rule is 12. The assembly rule can
only be picked from the bitscore ladder 50/100/250/500/1000 fixed before the series ran,
and the negatives select 250.

Three scoring choices are argued in full in
[`docs/day5-scoring.md`](docs/day5-scoring.md), and all three are reversible by flag.

1. **Score at the family, Coronaviridae.** The 2019 taxonomy puts the nearest relative
   outside the genus; a novel organism's contig cannot resolve below the family anyway
   (the genome contig's own LCA *is* Coronaviridae); and one negative patient's real NL63
   infection is an alphacoronavirus, invisible below the family.
2. **Two patients' genuine endemic coronavirus infections are not false alarms.** The
   November 2019 database holds all four endemic human coronaviruses and names them
   correctly; a pipeline reporting HKU1 in an HKU1 patient was right. Counting a correct
   answer as a false alarm costs Kraken2 a factor of ten to a hundred and leaves the
   assembly search unthresholdable.
3. **Score DIAMOND's stock mode, not the pre-registered `--very-sensitive`.** Both were
   run on everything. Six comparisons on genuine signal gave six identical best bitscores,
   and the only differences anywhere are false positives — at a sixth of the cost.

Everything the scoring could have decided otherwise is in the record: clade × framing ×
budget × sensitivity × attribution, 128 operating points, of which the headline is two.
A reader who disagrees with a threshold moves it with a flag.

## 3. Limitations

**The series cannot see below its own background: it has four usable rungs, not five.**
The dilution was cut from the index patient's own reads, and that background carries ~90
Coronaviridae reads of its own (78 in Sarbecovirus; 71 verified by alignment as
SARS-CoV-2 reads the ground-truth label declined to count). Those were never diluted, so
they sit at full strength in every rung.

| rung | labelled | labelled /M | resident | effective /M | undiluted share |
| --- | --- | --- | --- | --- | --- |
| 50,000 | 50,000 | 3,853 | 90 | 3,860 | 0.2% |
| 10,000 | 10,000 | 771 | 90 | 778 | 0.9% |
| 1,000 | 1,000 | 77.1 | 90 | 84.0 | 8.3% |
| 100 | 100 | 7.71 | 90 | 14.6 | 47% |
| 10 | 10 | 0.77 | 90 | 7.71 | **90%** |

The benchmark's floor is therefore ~7 viral reads per million, and it belongs to the
benchmark rather than to either method. No sensitivity limit below it should be quoted,
and neither method should be called more sensitive than the other at the limit: ≈15
against ≤7.7 effective reads per million is well inside the floor's own width. The fix
needs no new sequencing — cut a second series into a *negative* library's background,
confirmed free of sarbecovirus reads across 228 million pairs — and two more decades of
dose axis become visible.

**The alarm thresholds are upper bounds, because some "false alarms" are not false.**
Seven of the cohort's non-endemic coronavirus reads have since been traced to one
patient's genuine HKU1 infection, of a genotype RefSeq does not carry — it curates one
genome per species — which Kraken2 scattered onto rodent and bovine relatives because for
those k-mers those really are the nearest thing the database holds. That takes the
cohort's non-endemic count from 35 family-level reads to 28, and to zero at genus level.

The two samples that actually set Kraken2's threshold of 12 — one patient each with 11
*Sparrow coronavirus* reads — have had the same pass and do not survive it. Not one of
those 22 reads aligns to the sparrow genome, none has a protein-level hit anywhere in the
dated set, none is human, and the two samples share no 31-mer, so they are two independent
coincidences rather than one contaminant: an 18 bp match reaching a minimizer, counted
eleven times because the mates overlap and Kraken2 scored the same block twice. Six
non-endemic coronavirus reads then remain cohort-wide, one examined and unaccounted for,
five never examined, and zero contigs.

**The thresholds quoted in §1 are the as-scored ones, and they too are upper bounds.**
With those deductions Kraken2's threshold falls from 12 reads to 3 and its limit from 7.71
to 0.77 viral reads per million, a full decade, and MEGAHIT + DIAMOND's false alarms go to
zero. That figure is in the record as a secondary number rather than the headline, because
every deduction rests on a diagnostic the method does not itself run, each was made after
it was known which samples set the threshold, and all of them move the number in the
method's favour. Each deduction being individually correct is what makes a threshold
fitted to its own cohort persuasive, not what makes it safe; the defensible version
declares the diagnostic in advance and applies it to every alarm. Both columns are in
[`results/day5-provenance.json`](results/day5-provenance.json). Note also where that
decade would land: the rung it buys is the bottom one, whose dose is 90% undiluted
background and which sits below the floor above, so even granted in full it would move
Kraken2's limit onto a rung this series cannot resolve rather than separate the two
methods.

The same correction has now removed a threshold-setting read from **both** methods.
MEGAHIT + DIAMOND's bitscore threshold of 250 is selected by exactly one contig in the
whole cohort, whose best hit is a rat coronavirus at 91.2% amino-acid identity and which
is 100.00% identical at nucleotide level to the HKU1 genotype B genome. Attributing a
contig by its best hit rather than by its LCA cannot catch that one, because the best
hit is itself wrong: the rat protein really is the closest in the dated database, and
HKU1's own replicase ties it at the identical bitscore. A translated search against an
incomplete reference set can be confidently wrong about *which* relative a sequence is
while being right that it is a coronavirus, and no rule applied to its own hit table
fixes that.

**The negatives bound coronavirus false alarms, not pathogen-agnostic ones.** All twenty
carry *Equine infectious anemia virus* — 3,284,717 reads in total, up to 499,690 in one
sample, and the largest viral call in most of them. It is not the virus: the reads are
98.8% identical to the genome but confined to a single 1.7 kb window inside *pol*, with no
coverage of *gag*, *env* or either long terminal repeat, depth that stops within one base
at both ends, and the same window in unrelated patients. It is a cloned gene fragment,
presumably a laboratory reagent. The coronavirus limits are untouched, since they are
drawn per clade, but a pathogen-agnostic pipeline run on this cohort reports a notifiable
equine pathogen in 20 of 20 clinical samples, and any alarm keyed on a read count fires
every time. A read count cannot see the difference; where the reads sit on the reference
can. [`docs/hku1-coinfection.md`](docs/hku1-coinfection.md) has this and the HKU1 case
in full, and both checks are commands rather than arguments:
`hindsight arm-b diagnose-strain` for one diverged strain against two organisms,
`hindsight arm-a diagnose-taxon` for an organism against a laboratory sequence.

**Twenty negatives bound specificity coarsely.** 0 alarms in 20 samples bounds the
per-sample false-alarm rate at under 14% with 95% confidence, and no better. The per-read
bound is much tighter and is a different statement: 0 sarbecovirus reads in 228,407,636
non-host pairs, about 0.013 per million. The cohort is also *cleaner* than the positive —
host-depleted before submission, 0.075% human against the index sample's 15.8% — so these
are specificity figures on unusually clean input, and that flatters a translated search
more than exact k-mer matching.

**n = 1 on the positive.** One index sample, one virus, one library preparation. The
dilution series is the generalisation axis, not the sample, and nothing here says how
either method behaves on a different novel pathogen in a different matrix.

**The read-level translated search is a mechanism result, not a limit.** It detects a flat
~93% of the labelled fragments at every dose, and from 100 fragments downward it is below
its own noise floor of 177 background fragments. Nine true reads among 177 false ones is
not a detection by any operational standard, and the read half has no independent
false-alarm rate (it is off for the negatives), so it cannot be thresholded at all.

**RefSeq, not GenBank or nr/nt.** This reconstructs the curated RefSeq viral set, which is
what Kraken2's standard viral library and comparable reference sets are built from. Wu et
al. searched much larger databases, so the assembly pipeline is a faithful automation of
their *method*, not a replay of their *search space*. GenBank has no equivalent dated
catalog, so a time-correct nr is not reconstructible this way.

**Release granularity, not day granularity.** Releases are roughly two months apart (97 is
4 November 2019, 98 is 6 January 2020 — a 63-day gap), so a 1 December 2019 cutoff yields
the database as of 4 November. The error is one-sided and deliberate: hindsight never
includes something published after the cutoff, but it may omit something published in the
weeks before it. `hindsight resolve` prints the blind window.

**No protein-language-model arm.** The leakage argument for a model trained on today's
sequence databases is not resolvable in seven days.

**This is a study of detection sensitivity.** It makes no claim about the virus's origin,
and it is not a finding that surveillance does not work. It is a measured detection limit
under one leakage-controlled condition.

## 4. Reproduction

```bash
uv sync --locked --all-extras     # or: pip install -e ".[dev]"
```

Set `NCBI_API_KEY` if you have one ([free from your NCBI
account](https://www.ncbi.nlm.nih.gov/account/settings/)). It raises the E-utilities
request ceiling from 3/s to 10/s and takes a full build from hours to under one.

```console
$ hindsight build --date 2019-12-01 --out db/november-2019   # the dated database
$ scripts/day2.sh                                            # reads, label, dilution series
$ scripts/day3.sh                                            # Kraken2
$ scripts/day4.sh                                            # MEGAHIT + DIAMOND
$ scripts/day5.sh                                            # scoring, both pipelines
$ python scripts/day6_figure.py                              # the figure above
```

Each stage is resumable and writes its own provenance; each later stage reads the
earlier ones' records rather than their working files.

```console
$ scripts/day3.sh
index sample: 25,928 of 62,117 ground-truth viral fragments in Coronaviridae (41.7%)
  viral050000: 20,892 of 50,000 detected
  viral010000:  4,232 of 10,000 detected
  viral001000:    420 of  1,000 detected
  viral000100:     39 of    100 detected
  viral000010:      4 of     10 detected

$ scripts/day4.sh
index sample: 99.68% of the SARS-CoV-2 genome in one 29,809 bp contig at 100.0% identity
  viral050000: 49,077 of 50,000 detected
  viral010000:  9,824 of 10,000 detected
  viral001000:    936 of  1,000 detected
  viral000100:     41 of    100 detected
  viral000010:      1 of     10 detected

$ scripts/day5.sh
Kraken2 rule: at least 12 reads placed in coronaviridae
  false alarms: 0 of 20 samples; limit of detection 7.71 viral reads per million
  (viral000100)
MEGAHIT + DIAMOND rule: at least one coronaviridae contig at bitscore >= 250
  false alarms: 0 of 20 samples; limit of detection 0.771 viral reads per million
  (viral000010) (the series' lowest rung; the limit is at or below it)
```

Both limits are the rungs' *labelled* doses, which is how the scorer reports them and is
the conservative choice for both methods: it credits neither with the virus the
background already carried. The 0.771 figure is the bottom rung's label and not a
sensitivity claim, because that rung sits below the floor of §3 — which is why the
scorer prints it as an upper bound rather than a point, and why the limit quoted in the
abstract is ~8 reads per million for both.

### 4.1 Cost

Measured on 4 cores, 15.7 GiB and about 30 GB of writable disk: a small cloud container,
nothing rented, no GPU.

| | |
| --- | --- |
| Database build (release 97, nucleotide + protein + taxonomy) | under 1 h with an API key |
| Reads: fetch, QC, host filter, label, series, 21 samples | 147 min |
| Kraken2, including 20 negatives at native depth | 191 min, of which classification is 24.5 min for 250 M read pairs |
| MEGAHIT + DIAMOND, including 20 negatives at native depth | ~14 h, dominated by one 34.6-min assembly per sample |
| Scoring | ~1 s, no tool invoked |

Classification is never the cost; fetching and host filtering is. A disk ceiling is
therefore a staging problem rather than a compute one, and the negatives run at native
depth on a machine that cannot hold them, by classifying each as it comes off the filter
and keeping only the output. Assembly is always the cost in the assembly pipeline, and
it does not compose: Kraken2 scores each read independently, so a dilution level is its
parts' results added up, while an assembly is a statement about a whole read set.

## 5. The tool on its own

`hindsight` builds the viral reference set as it stood on any date since 2003,
independently of this benchmark.

```bash
hindsight resolve 2019-12-01                 # which release a date selects, and what it excludes
hindsight build --date 2019-12-01 --out db/nov2019 --manifest-only   # accessions only, no sequence
hindsight build --date 2019-12-01 --out db/nov2019 --molecules nucleotide
hindsight build --date 2019-12-01 --out db/nov2019                   # nucleotide, protein, taxonomy
hindsight verify db/nov2019                  # re-check a build against its own provenance
hindsight releases --limit 10                # archived releases
hindsight taxdumps                           # dated taxonomy snapshots
```

A build writes `manifest.tsv` (every accession with version, molecule, taxid, organism,
RefSeq status and length), `viral.fna`, `viral.faa`, `seqid2taxid.map` in Kraken2's
format, `taxonomy/` as of 2019-11-01, `withdrawn.tsv` for accessions the release listed
that NCBI will no longer serve, and `provenance.json` with a SHA-256 of each. Builds are
deterministic — sorted manifest, manifest-order output, fixed FASTA width — so two
builds of the same cutoff are byte-identical, and an interrupted build resumes from the
last completed shard.

Two further properties travel with any database it builds. Version substitution is
rejected rather than tolerated, because a record quietly upgraded to its current
revision is the leak this tool exists to prevent. And taxonomy snapshots are monthly, so
the tree sits up to a month before the release date — three days for release 97 — while
cutoffs before 2014-08-01 have no contemporaneous taxonomy at all, which hindsight
reports rather than reaching for a later tree.

The benchmark's own commands have the same shape:

```bash
hindsight samples ...              # fetch, QC, host filter, label, freeze the series
hindsight arm-a ...                # Kraken2; `diagnose-background` explains its own noise
hindsight arm-b ...                # assembly; `diagnose-contigs` checks a contig's identity
hindsight score run|report         # thresholds and limits, from both records
```

[`docs/archive-notes.md`](docs/archive-notes.md) records the NCBI archive's quirks that
this tool had to be built around, each now encoded in a test.

## 6. Development

```bash
uv run pytest             # 347 offline tests (5 skip without diamond and minimap2)
uv run pytest -m network  # 3 more against NCBI; downloads ~1.2 GB
uv run ruff check . && uv run ruff format --check .
```

Network tests are deselected by default and excluded from CI: they are rate-limited per
IP and move more than a gigabyte. Run them when changing anything about fetching, and
before trusting a build. `hindsight refresh-releases` regenerates the vendored release
table from NCBI's release notes, which is only needed when a new release lands.
