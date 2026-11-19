# Day 4: Arm B, assembly plus translated protein search

The benchmark's second detection arm, and the one that is not a classifier. De novo
assembly with MEGAHIT, then `diamond blastx` of the contigs against the protein records
RefSeq release 97 held on 4 November 2019, scored by a lowest common ancestor over the
November 2019 taxonomy.

```console
$ scripts/day4.sh
```

Everything it decided lands in [`results/day4-provenance.json`](../results/day4-provenance.json),
and `hindsight arm-b report` prints the tables below from it.

This arm exists because Arm A answers a different question from the one the discovery
actually answered. Wu et al. did not run a read classifier. They assembled the
host-filtered reads de novo and matched the longest contig against prior coronavirus
sequences by homology. "Would a classifier have called it?" and "would assembly plus a
homology search have called it?" are separate questions with, as it turns out, opposite
failure modes.

## The headline

**The arm reconstructed 99.68% of the SARS-CoV-2 genome in a single 29,809 bp contig at
100.0% nucleotide identity, from a database whose nearest relative matched it at 79.9%
amino acid identity.**

That is worth stating twice because the two numbers describe different things. The
contig is a perfect reconstruction of a genome the reference set had never seen; the
79.9% is how far away the closest thing the reference set *did* hold sits. The assembly
owes nothing to the database — it is built from the reads alone — and the database's only
job is to say what the contig resembles. It said: a coronavirus replicase, best match
`NP_828850.1` from SARS-CoV Tor2, 4,422 amino acids aligned, E = 0, bitscore 7,120.

The LCA over the top 10% bitscore band lands on **Coronaviridae**, not Sarbecovirus. That
is not a weakness of the method; it is the method being honest. A contig from a genuinely
novel sarbecovirus ties SARS-CoV and several bat coronaviruses at near-identical scores,
and the deepest taxon consistent with all of them is the family. A single-best-hit rule
would have reported "SARS-CoV" — precise, confident and wrong.

| | |
| --- | --- |
| non-host read pairs assembled | 10,921,957 |
| contigs | 27,732 (13,695,075 bases, N50 462) |
| assembly wall clock | 34.6 min on 4 cores |
| contigs placed in Coronaviridae by LCA | 4 |
| of those, above bitscore 100 | 1 — the genome |
| ground-truth viral fragments recovered | 60,923 of 62,117 (**98.1%**) |
| fragments recruited that are not in the label | 303 |
| detection, reads per million post-QC pairs | 4,719 |

The arm's own unit is a contig, but the benchmark's unit is reads, so the non-host reads
are mapped back onto the contigs the arm placed in Coronaviridae. A read that aligns to a
contig the arm called coronaviral is this arm's equivalent of a read Kraken2 called
coronaviral, and the two arms become comparable. That step is done in two passes for a
reason worth reading before trusting any number it produces — see
[the section on it](#the-read-mapping-has-to-see-the-whole-assembly).

## The detection curve, and where it crosses Arm A

Each dilution level is assembled whole — there is no composition shortcut here, for
reasons in [its own section](#why-every-level-costs-a-full-assembly) — searched, and
scored by mapping reads back.

| viral fragments | per million post-QC | Arm B detected | Arm B rate | Arm A rate |
| --- | --- | --- | --- | --- |
| 50,000 | 3,853 | 49,077 | **98.2%** | 41.8% |
| 10,000 | 771 | 9,824 | **98.2%** | 42.3% |
| 1,000 | 77 | 936 | **93.6%** | 42.0% |
| 100 | 7.7 | 41 | 41.0% | 39.0% |
| 10 | 0.77 | 1 | 10.0% | 40.0% |

**Above roughly 1,000 viral fragments the assembly arm is more than twice as sensitive as
the classifier. Below 100 it loses to it.** The two curves cross at about 8 viral reads
per million.

The mechanism is visible in the assemblies themselves. Where coverage is sufficient, a
contig recruits reads that would never have classified individually: a 151 bp read sharing
no exact 35-mer with anything in a November 2019 database still assembles with its
neighbours, and once it is part of a 29 kb contig whose replicase matches SARS-CoV, it is
detected. That is the whole advantage of assembly, and it is worth a factor of 2.4.

But it is an advantage with a floor rather than a slope. Arm A's rate is flat at ~41% all
the way down, because Kraken2 scores each fragment independently of every other; the
absolute count falls and the rate does not. Arm B's rate holds near 98% and then collapses,
because below some coverage the assembler cannot build the contig that does the recruiting.
Watch the longest coronaviral contig fall apart:

| viral fragments | coronaviral contigs | their bases | longest | best bitscore |
| --- | --- | --- | --- | --- |
| 50,000 | 4 | 32,115 | 29,868 | 7,120 |
| 10,000 | 5 | 32,437 | 29,873 | 7,120 |
| 1,000 | 24 | 28,263 | 3,900 | 2,065 |
| 100 | 17 | 8,182 | 975 | 319 |
| 10 | 4 | 2,709 | 975 | 270 |
| background (no spike) | 3 | 2,247 | 975 | 57.8 |

At 10,000 fragments the genome is one contig. At 1,000 it has shattered into 21 pieces
that still total 28 kb — the sequence is all there, in fragments. At 100 only 6 kb of it
survives assembly, in pieces under 600 bp. At 10 there is a single 462 bp contig.

**Note that the contig count is not monotone in dose**, and is therefore the wrong thing to
report. It rises from 4 to 24 as the signal *falls*, because fragmenting one genome into
twenty pieces increases the count. Bases and best bitscore both behave; the count does not.

## The read half: it is the assembler that fails, not the search

The table above raises a question the contig half cannot answer. When the arm stops firing
at 10 fragments, did the **protein search** run out of homology, or did the **assembler**
run out of coverage? Those have opposite implications — the first says translated search
cannot see this far, the second says it can and needs more depth — so the arm runs the same
`diamond blastx` a second way: on the individual reads, with no assembly step at all.

This half is affordable for exactly Arm A's reason. A translated search scores each read
from that read alone, so a dilution level is the shared background's hits plus that level's
own viral reads' hits. One expensive pass over 10.86M background pairs, five trivial ones.

The answer is unambiguous:

| viral fragments | per million | read-level detected | read-level rate | contig rate |
| --- | --- | --- | --- | --- |
| 50,000 | 3,853 | 46,729 | 93.5% | 98.2% |
| 10,000 | 771 | 9,344 | 93.4% | 98.2% |
| 1,000 | 77 | 942 | 94.2% | 93.6% |
| 100 | 7.7 | 92 | 92.0% | 41.0% |
| 10 | 0.77 | 9 | **90.0%** | 10.0% |

**The per-read rate is flat at about 93% all the way down to ten fragments.** Nine of the
ten reads at 0.77 per million are individually recognisable as coronaviral by protein
homology against a November 2019 database. The homology was there the whole way; what
collapsed was the assembler's ability to build a contig out of too few reads. Assembly is
not what makes this virus findable — it is what makes it *nameable*, by producing the
29 kb object an analyst can look at.

That also quantifies what assembly adds where it works: at 50,000 fragments the contig
half recovers 98.2% against the read half's 93.5%, because a contig recruits reads whose
own best protein hit falls below threshold.

### Why a flat rate is not a lower detection limit

It would be easy to stop there and declare the read half the most sensitive arm in the
benchmark. It is not, and the reason is the same one that undercut the contig half: the
background.

The spike-free background's 10.86M fragments contain **177 that this search places in
Coronaviridae** and 22 it places in Sarbecovirus. So a flat 93% detection rate converts
into a detection *limit* only where the signal clears that floor:

| viral fragments | detected | background floor | ratio |
| --- | --- | --- | --- |
| 1,000 | 942 | 177 | 5.3 : 1 |
| 100 | 92 | 177 | 0.5 : 1 |
| 10 | 9 | 177 | 0.05 : 1 |

At a hundred fragments the read half is still detecting 92% of them and is nonetheless
**below its own noise**. Nine true reads among 177 false ones is not a detection by any
operational standard.

Set against Arm A, this is the day's most important cross-arm observation. Kraken2 found
**zero** sarbecovirus reads across 228 million non-host pairs from twenty pre-COVID
negatives, and 35 non-endemic coronaviridae reads — 0.153 per million. Translated search
at 30% identity produces four orders of magnitude more background than exact 35-mer
matching does. **Arm B is markedly more sensitive per read and markedly less specific, and
specificity is what sets a detection limit.**

Two honest caveats on that comparison. The floor above is measured on one sample's own
background, not on the twenty negatives, and a single library is a single observation —
which is what the negatives are for. And Day 3 found that most of the index sample's
apparent sarbecovirus "false positives" were real viral reads its ground-truth label had
missed, so some of these 177 may be too. Both point the same way: the negatives settle it,
not this table.

The read half is off for the negatives (`[negatives] read_level = false`), because
searching 228 million read pairs at `--very-sensitive` is about 49 hours. That is a
measurement, and it is the one real gap in the day: a read-level false-alarm rate needs a
machine with more cores, and it is a Day 7 item.

## The arm alarms on a clean background, and what to do about it

The bottom row of that table is the uncomfortable one. **The spike-free background — the
index sample with every ground-truth viral fragment removed — still yields three contigs
that the LCA places in Coronaviridae.** They are the same three at every level, because
they come from the background, which every level shares:

| contig length | best hit in the dated database | identity | bitscore |
| --- | --- | --- | --- |
| 975 bp | `YP_001039970.1`, Rousettus bat coronavirus HKU9 | 28.4% | 43.9 |
| 862 bp | `YP_001876447.1`, Beluga whale coronavirus SW1 | 32.6% | 57.8 |
| 410 bp | `YP_001876447.1`, Beluga whale coronavirus SW1 | 29.5% | 40.8 |

These are real remote homology, not a bug: a translated search against a database of
distant relatives will always find something at 30% identity over a hundred residues, and
the arm is right to report it. But it means **"the arm found a contig in Coronaviridae" is
not a detection rule**, because that rule fires on a sample containing no SARS-CoV-2 at all.
The bare contig count cannot tell a 29 kb genome from this floor: 3 contigs against 4.

Bitscore separates them, and not narrowly. The floor tops out at **57.8**. The genome
scores **7,120**, and even the single genuine contig at the 10-fragment limit scores **270**
— nearly five times the floor. Above a bitscore of 100 the floor vanishes completely:

| | contigs ≥ bitscore 100 | their bases |
| --- | --- | --- |
| background | 0 | 0 |
| 10 fragments | 1 | 462 |
| 100 fragments | 12 | 5,935 |
| 1,000 fragments | 21 | 26,016 |
| 10,000 fragments | 1 | 29,873 |
| 50,000 fragments | 1 | 29,868 |

**The arm applies no such threshold.** `alarm` in the record is the unfiltered rule, exactly
as the arm behaves, and alongside it every sample carries an `evidence` profile: contigs,
bases, longest contig and best bitscore, each also reported above a fixed ladder of
bitscore cuts (50, 100, 250, 500, 1000). Round numbers, fixed before the series ran, not
fitted to this curve.

Choosing the cut is Day 5's job and it needs the negatives to do it, because a floor
measured on one sample's background is one observation. Twenty pre-COVID samples at native
depth measure what that floor looks like across real libraries, and the threshold belongs
where it separates them from the 270 above. A threshold picked after seeing the positive's
answer is not a detection limit.

## What `--very-sensitive` is actually worth

`config/day4.toml` argues at length that this arm must run at `--very-sensitive`, because
DIAMOND's default mode is tuned for high-identity matches between close homologues and
SARS-CoV-2's proteins sit 10–25% diverged from anything a 2019 database held. Running the
*homology* arm in a mode that cannot see homology at that distance would measure the mode
rather than the method.

The honest way to hold that position is to measure the alternative rather than assert it,
so every assembly is also searched at DIAMOND's stock sensitivity. On the index sample:

| | coronaviral contigs | longest | best bitscore | search time |
| --- | --- | --- | --- | --- |
| `--very-sensitive` | 4 | 29,809 | 7,120 | 0.8 min |
| stock | 1 | 29,809 | 7,120 | ~0.15 min |

**Stock mode found the genome, identically, and found none of the three false positives.**

One sample could not settle it, though — with 62,117 viral fragments either mode finds a
29 kb contig easily. The question is whether the extra sensitivity buys anything at the
*limit*, where the arm is choosing between 462 bases and nothing. So the comparison runs at
every dilution level and on the background too:

| | `--very-sensitive` contigs / bases | stock contigs / bases | best bitscore, either mode |
| --- | --- | --- | --- |
| background (no spike) | 3 / 2,247 | **0 / 0** | 57.8 → none |
| 10 fragments | 4 / 2,709 | **1 / 462** | 270 |
| 100 fragments | 17 / 8,182 | **14 / 5,935** | 319 |
| 1,000 fragments | 24 / 28,263 | **21 / 26,016** | 2,065 |
| 10,000 fragments | 5 / 32,437 | **2 / 30,190** | 7,120 |
| 50,000 fragments | 4 / 32,115 | **1 / 29,868** | 7,120 |

**At every level, stock mode found the same genuine signal at the same bitscore, and none of
the background floor.** The difference is exactly three contigs everywhere: the three
background artefacts, which only `--very-sensitive` sees. At the 10-fragment limit stock
mode returns one contig of 462 bases — the real one — against very-sensitive's four, three
of which are noise.

So `config/day4.toml`'s argument, which this file quoted above, is wrong on this cohort, and
the measurement says so plainly: `--very-sensitive` cost six times as much per search and
contributed nothing but false positives, while stock mode has a **clean zero on a
spike-free background** and still fires at ten fragments.

The configuration has not been changed, because a setting edited after seeing the result is
no longer a pre-registered choice, and because one cohort's background is thin evidence for
a specificity claim. Instead both modes now run on every negative as well — the comparison
is nearly free beside a 35-minute assembly — so the record carries a false-alarm rate for
each, measured on twenty real libraries, and Day 5 chooses on that. On present evidence it
should choose stock.

**Day 5 did choose stock**, from this record and without re-running anything, and swapped
the two values in `config/day4.toml` while leaving the pre-registration above it intact.
See [day5-scoring.md](day5-scoring.md).

## A limitation in the negatives that Arm B is more exposed to than Arm A

Worth stating before the false-alarm rate is quoted, because it bounds what that rate
means. **The twenty pre-COVID negatives are almost free of human reads, and the index
sample is not.**

| | human reads removed by the host filter |
| --- | --- |
| index sample (SRR10971381) | **15.83%** |
| negatives, median of 20 | **0.075%** |
| negatives, worst of 20 | 0.43% |

Two orders of magnitude apart. The negative cohort appears to have been host-depleted or
virome-enriched before submission, so its reads reach the arms already clean, while the
positive sample arrives 16% human even after NCBI's own redaction has stripped half the
spots.

That asymmetry costs Arm B more than Arm A. Kraken2 matches exact 35-mers against a
viral-only database, and a human read has to share one to be miscalled. A translated
protein search at 30% amino-acid identity over a hundred residues is a much looser
instrument, and human transcripts contain plenty of domains distantly related to viral
proteins. **So the false-alarm rate these negatives measure is a rate on unusually clean
input, and it is likely optimistic for this arm in a way it is not for Arm A** — the two
arms' specificity figures are measured on the same cohort but are not equally stressed by
it.

This is a property of the cohort chosen on Day 2 rather than anything Day 4 did, and it is
recorded here because the write-ups should not compare the two arms' specificity without
it. The fix is a Day 7 item: either negatives with native human content, or the index
sample's own host-aligned reads run through both arms as a dedicated host-derived
false-positive test. The second is cheap: it needs one more bowtie2 pass over the index
sample's QC reads, which are on disk, to recover its 2,053,621 host-aligned pairs.

## What twenty pre-COVID negatives say

Twenty of the 408 matched 2016–17 Wuhan BALF libraries, at native depth, **228,407,636
non-host read pairs** — the same cohort and the same count Arm A used, so the two arms'
specificity is measured on identical input.

**No contig in Sarbecovirus, in any of the twenty.** That is the clade the index sample's
genome lands in once the LCA is allowed below family, and it is empty across a quarter of
a billion read pairs.

At family level the arm alarms on three of the twenty, and the interesting part is what
the three are:

| sample | non-host pairs | coronaviral contigs | longest | best bitscore | what it is |
| --- | --- | --- | --- | --- | --- |
| SRR13677922 | 14,751,395 | 23 | 2,652 | 1,188 | *Human coronavirus HKU1*, a real infection |
| SRR13677700 | 10,530,048 | 3 | 664 | 270 | *Human coronavirus NL63*, a real infection |
| SRR13677843 | 9,407,026 | 1 | 951 | 46.6 | Beluga whale coronavirus at 28.5% identity — the floor |

**Two of the three alarms are correct answers.** The 2019 database contains HKU1 and NL63
and names them, at 97.8–100% amino-acid identity in the first case and 97.8–100% in the
second. An operational pipeline that reported them would have been right; they are false
alarms only against the question "is there a novel sarbecovirus here".

The third is the background floor from the previous two sections, appearing in a real
library exactly as predicted: one short contig, remote homology to a whale coronavirus, a
bitscore of 46.6. **At a bitscore cut of 100 it disappears and the arm alarms on two of
twenty, both of them genuinely infected.** And as the previous section predicted, stock
DIAMOND mode does not produce it at all — 26 contigs in two samples against
`--very-sensitive`'s 27 in three. That is the second independent confirmation, now on real
negatives rather than one sample's background, that the tuned setting contributes only
false positives.

Excluding the four endemic human coronaviruses by name, as Day 3 did, leaves **8 contigs**
across the same three samples rather than 27. But six of those eight are in the HKU1
patient and are HKU1's own sequence whose LCA resolved one rank up, to genus or subgenus,
plus one whose best hit is Rat coronavirus Parker at 91.2% — embecovirus neighbours of the
virus actually present. This is the same effect Day 3 reported as "spillover" at read
level, and it means the endemic-exclusion figure is not a clean false-alarm count either.
Day 5 gets all three framings; the record carries every clade per sample.

### One cross-arm difference worth keeping

The NL63 patient is invisible to Arm A's headline. NL63 is an *alpha*coronavirus, and
Day 3 reported at Betacoronavirus level, so Arm A's negatives summary names the HKU1 and
OC43 infections and not this one. Arm B finds it because it is scored at family level.

That is an argument for the family-level threshold independent of the one Day 3 made from
the dated taxonomy, and it points the same way: **a benchmark that scores below family
systematically under-reports what its own arms found.**

## The read mapping has to see the whole assembly

This section exists because the first version of this arm got it wrong, the error was
large, and the shape of it is worth more than the correction.

The read mapping builds a bowtie2 index from the contigs the arm placed in the clade and
aligns the non-host reads to it. The flaw is in the word "from": **mapping quality means
"unique among the sequences you gave the aligner", not "this read came from here".** Hand
the aligner 23 contigs out of an assembly of 7,438 and every read belonging to the other
7,415 has nowhere else to go. With local alignment it lands on whatever thirty-odd bases
it happens to share with a clade contig, and is then scored confidently, because nothing
in the index competes for it.

Measured on the HKU1 patient:

| index given to the aligner | reads mapped onto the 23 coronaviral contigs |
| --- | --- |
| the 23 clade contigs | 1,435,005 |
| the whole 7,438-contig assembly | 2,885 |

A factor of 498, and the larger figure would have made an endemic coronavirus a tenth of
that library. Three independent lines agree on the smaller one: Arm A counted 2,399 HKU1
reads in the same sample, the assembly's own k-mer coverage for those contigs is 2–102,
and the contigs are close enough to the dated HKU1 reference that a 35-mer classifier
should and does find them:

```console
$ hindsight arm-b diagnose-contigs --sample SRR13677922 --dated-reference NC_006577.2
k141_681 is 2,509 bp and covers 8.27% of NC_006577.2 at 92.93% nucleotide identity. No
single contig is the genome, but 19 of them tile 52.31% of NC_006577.2 at 95.05%
weighted identity over 15,744 bases, so the organism was present and the assembler
fragmented it.
```

Individual contigs run 86.2% to 100.0% identity. **The reference is the dated one**, so
unlike the index sample's diagnostic this involves no post-2019 sequence and carries no
leakage caveat: it asks only how far the virus in the patient sat from what the database
actually held. Four of the 23 contigs do not align to it at nucleotide level at all
despite HKU1 protein hits at 76–88% amino-acid identity, which is either a diverged
region of the same strain or a second betacoronavirus, and is a question for Day 7
rather than something this day settles.

**The fix is two passes, and the subset index survives as a candidate filter rather than
as the measurement.** Adding references to an index can take reads away from a contig but
never add them, so the reads the subset index reports are a safe superset of the right
answer. They are then re-adjudicated against the whole assembly. That costs one cheap pass
over the library plus one expensive pass over a few hundred thousand candidates, which on
this sample is 11.5 minutes against the 2 hours 15 the exhaustive version took, and it
reproduces it to 5.4% — 3,041 against 2,885, high because the second pass adjudicates
candidate mates singly rather than as pairs, so it bounds the true count from above.

Every sample now records both counts and the ratio between them, because the ratio is the
thing that tells you whether the shortcut was safe:

| sample | candidates | adjudicated | ratio |
| --- | --- | --- | --- |
| index sample (SRR10971381) | 61,233 | 61,226 | 1.0 |
| each of the six series levels | — | — | 1.0 |
| SRR13677700 | 23 | 23 | 1.0 |
| SRR13677843 | 45 | 45 | 1.0 |
| SRR13677922 | 1,435,005 | 3,041 | **471.9** |

**Nine of the ten were already right and one was wrong by 472-fold**, which localises the
failure instead of spreading it. The index sample's clade contigs are one near-complete
29.8 kb genome with nothing like it in the library, so a subset index costs nothing there;
the series inherit that property; and the HKU1 patient's are 23 short fragments of an
endemic virus in a library full of its relatives, which is the worst case. The index
sample's 98.1% recovery and the series' false-positive floor both stand as first measured
— the 169 spurious reads at the lowest dilution are real, not an artefact, which is the
opposite of what I expected before measuring it.

The general lesson is not about bowtie2. **A confidence score is relative to the
hypotheses you offered**, and subsetting a reference set to the part you are interested in
silently removes the alternatives that score is defined against. Arm A is not exposed to
it, because Kraken2 scores a read against the whole database or not at all.

## Why every level costs a full assembly

Arm A composes a dilution level from its parts: one expensive run over the shared
background, five trivial runs over the viral reads, added together. That is exact, because
Kraken2 assigns each fragment from the fragment alone.

Assembly is the opposite kind of operation. **A contig is a statement about the set of
reads**, so removing 49,990 viral reads changes which contigs exist, how long they are and
what coverage they carry. There is nothing to compose, and the table above is the proof:
the same virus appears as one 29,868 bp contig at 50,000 fragments and as twenty-one
fragments at 1,000. Five levels therefore cost five full assemblies at ~34 minutes each,
plus the background. That is most of this arm's runtime and it is not optimisable away —
it is a property of the method, not of the implementation.

The *per-read* search does compose, for exactly Arm A's reason, which is the only thing
that makes the read half of this arm affordable at all.

## What is kept, and what is thrown away

Disk is the binding constraint in the container (~21 GB usable), and this arm's working
set would not otherwise fit. The same prune-as-you-go discipline Day 3 used:

- each dilution level is materialised from its frozen parts, assembled, searched, and the
  materialised FASTQ deleted before the next level starts, so peak disk is one level
- MEGAHIT's intermediate directory is deleted on success and only `*.contigs.fa` survives
- the read half's query FASTA (~3 GB) is deleted after its search
- both bowtie2 indexes the read mapping builds, over the flagged contigs and over the
  whole assembly, are rebuildable and are not kept, nor is the candidate FASTA between
  the two passes
- each negative's reads are deleted once it has been searched

What is kept is the arm's output: contigs, the gzipped tabular hits, and the record.

## Cost

Measured on 4 cores and 15 GiB, which is what the container has.

| step | time |
| --- | --- |
| one 10.9M-pair MEGAHIT assembly | 34.6 min |
| contig `blastx`, 27,732 contigs, `--very-sensitive` | 0.8 min |
| mapping reads back, both passes, index sample | 11.1 min |
| the same exhaustively, whole library on whole assembly | 2 h 15 min |
| read-level `blastx`, `--very-sensitive` | 78 s per 200,000 reads |

The search is never the cost; the assembly is. The one number that shaped the day's plan is
the last one: a full read-level pass over the 10.86M-pair background is 2.35 hours, and the
same pass over all twenty negatives would be about **49 hours**. Hence `[negatives]
read_level = false` — a measurement, not a guess. The negatives' contig half, at native
depth, is about 14 hours.

## For Day 5 and Day 7

**For Day 5.** The crossover is the finding: the two arms fail in opposite directions, so a
single detection threshold across both is the wrong shape. Both arms record every clade at
every level, so the line can be redrawn without re-running anything. Score at
**Coronaviridae**, as Day 3 recommended — Arm B makes it a requirement rather than a
refinement, because the genome contig's LCA *is* the family and a genus-level headline would
report zero. And the threshold for Arm B is a bitscore, not a contig count, with the
negatives measuring where it goes.

Second, the sensitivity question is now a decision rather than an open problem. Stock mode
matched `--very-sensitive` on genuine signal at every level and produced none of the
background floor, at a sixth of the cost. The negatives measure both. If they agree, Arm B
should be scored in stock mode and the tuned setting dropped.

Third, and least comfortable: the read half's flat 93% is not a sensitivity result to quote
on its own. It is a *mechanism* result — the homology is there at every level — and the
specificity figure beside it is what turns a rate into a limit. The same discipline applies
to the contig half's 98%.

**For Day 7.**

- A read-level false-alarm rate, on a machine with the cores for it. 49 hours at
  `--very-sensitive` on 4 cores is the only measurement this day wanted and could not
  afford, and it is the number that would make the read half's 93% quotable.
- A host-derived false-positive test, because the negatives cannot provide one: they are
  host-depleted and the positive is 16% human. The pipeline keeps only the non-host side,
  so this needs one more bowtie2 pass over the index sample's QC reads — which are on
  disk — to recover its 2,053,621 host-aligned pairs, and then both arms run on them.
  About an hour, and it is the measurement that would say whether a translated search
  miscalls human transcripts as viral.
- ~~Settle what the HKU1 patient's four non-aligning contigs are.~~ **Settled**, and the
  answer was about the reference set rather than the patient: the virus is HKU1
  **genotype B**, and RefSeq's single HKU1 genome is genotype A. All 23 contigs align to a
  2006 genotype B genome at 99.4% identity, the "divergence" is confined to the two genes
  that define the genotype, and one of the four contigs was never diverged at all — it is a
  324 bp inverted-repeat chimera that `minimap2 -x asm20` cannot place. See
  [`hku1-coinfection.md`](hku1-coinfection.md); it also revises this cohort's non-endemic
  betacoronavirus false-alarm count to zero, and identifies the largest viral call in the
  negatives as a cloned gene fragment present in all twenty.
- A combined Day 3 + Day 4 runner would classify and assemble each negative in one visit and
  save a whole fetch-and-host-filter pass (163 min). Day 3 had already deleted its negatives'
  reads and ran in a different container, so Day 4 re-staged them.
- Record order-independent read-ID digests. Day 2's counts reproduce exactly in a fresh
  container but its per-file sha256 digests do not, because bowtie2 runs without `--reorder`.
- Compare bwa-mem against bowtie2 for the host filter, carried over from Day 3.
