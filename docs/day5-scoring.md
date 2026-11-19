# Day 5: scoring both arms, and drawing the line

The benchmark promised one number per method: a limit of detection for a genuinely novel
betacoronavirus in a clinical metatranscriptome, in viral reads per million, at a stated
false-alarm rate. This is the day that produces it, and the day that makes the three
choices Days 3 and 4 deliberately left open.

```console
$ scripts/day5.sh
```

It takes about a second, runs no tool, and reads nothing but
[`results/day3-provenance.json`](../results/day3-provenance.json) and
[`results/day4-provenance.json`](../results/day4-provenance.json). That is the day's
quietest contribution: both arms recorded every clade at every dilution level and every
negative's per-taxon detail, so a reader who disagrees with a threshold can move it
without a machine.

```console
$ hindsight score run --clade betacoronavirus      # disagree with the clade
$ hindsight score run --framing all --budget 1     # or with the false-alarm count
```

Everything this day decided lands in
[`results/day5-provenance.json`](../results/day5-provenance.json), including the
operating points nobody chose: clade × framing × budget for Arm A, and that times both
sensitivity settings and both attribution rules for Arm B. 128 rows, of which the
headline is two.

## The headline

**Both arms would have flagged SARS-CoV-2 in the index sample, and both keep flagging it
down to the bottom of the dilution series, with no false alarm in twenty pre-COVID
samples.** The limit of detection is about **8 viral reads per million post-QC pairs** —
and it is a property of the benchmark's own floor rather than of either method, which is
the first thing this day found and the reason the two arms' headline figures differ.

| | Arm A (Kraken2) | Arm B (MEGAHIT + DIAMOND) |
| --- | --- | --- |
| statistic | reads placed in Coronaviridae | best bitscore of a Coronaviridae contig |
| alarm rule, from the negatives alone | **≥ 12 reads** | **≥ 1 contig at bitscore ≥ 250** |
| false alarms | **0 of 20 samples** | **0 of 20 samples** |
| lowest rung that alarms | 100 labelled fragments | 10 labelled fragments |
| limit, labelled dose | **7.7 reads per million** | **≤ 0.77 per million** (series floor) |
| limit, effective dose | ≈ 15 per million | **≤ 7.7 per million** |

Two columns of limit because the dilution series' background is not virus-free, and
below about a hundred labelled fragments that stops being a detail — see
[the floor](#the-series-cannot-see-below-its-own-background). The honest single number is
the effective one: **neither arm's limit is resolved below ~8 viral reads per million,
because that is as low as this series goes.** The index sample itself carried 4,787 labelled
viral fragments per million, so the benchmark demonstrates detection across a 621-fold
dilution and cannot say where either arm actually fails.

The rules above were set from the twenty negatives and never from the positive. A
threshold placed after looking at the positive's answer is not a detection limit, and
every number in this file follows that order: specificity first, then whatever
sensitivity follows from it.

## What counts as an alarm

Each arm gets its own statistic, because a single threshold across both is the wrong
shape — Day 4's phrasing and it is right. Arm A's output is reads; Arm B's output is a
contig, and its evidence is that contig's bitscore rather than the number of contigs,
which is not monotone in dose.

What *is* shared is the question: did this sample alarm? That turns both arms into the
same binary, at thresholds taken from the cohort of negatives:

| | Arm A, reads in Coronaviridae | Arm B, best bitscore |
| --- | --- | --- |
| worst negative | 11 (`SRR13677975`, sparrow coronavirus HKU17) | 226 (`SRR13677922`, rat coronavirus Parker) |
| second worst | 11 (`SRR13677700`) | — |
| threshold, 0 alarms permitted | 12 reads | 250, the next rung up the ladder |
| threshold, ≤1 alarm permitted | 12 reads (the two worst tie) | the arm's own unfiltered rule |

Arm B's threshold can only be chosen from the ladder `50/100/250/500/1000` that
`config/day4.toml` fixed before the series ran. That constraint matters: a bitscore
fitted to this cohort's negatives would be a free parameter dressed as a threshold, and
the ladder is what keeps it honest. The rung the negatives select, 250, happens to leave
the positive's lowest rung alarming at 270 — close, and recorded as close.

## Decision 1: the clade is Coronaviridae

Three arguments, from three different places, and each sufficient on its own.

**The 2019 taxonomy puts the nearest relative outside the genus.** 4,922 of the index
sample's fragments were assigned to *Bat coronavirus BM48-31/BGR/2008*, and the November
2019 taxonomy files that genome under "unclassified coronaviruses" in Coronaviridae, not
in Betacoronavirus where today's tree puts it. Family and genus differ by 5,510 reads,
and no benchmark that pairs dated sequences with a modern taxonomy can see that at all
(Day 3).

**A novel organism's contig cannot resolve below the family.** The 29,809 bp contig that
is 99.68% of the SARS-CoV-2 genome ties SARS-CoV and several bat sarbecoviruses at
near-identical bitscores, so its lowest common ancestor *is* Coronaviridae. Scoring Arm B
at the genus reports zero on a sample where it reconstructed the genome (Day 4).

**A real coronavirus infection is invisible below the family.** One negative patient had
*HCoV-NL63*, an alphacoronavirus. Arm A's genus-level reporting cannot see it; Arm B's
family-level reporting can.

What the choice is worth, at this day's framing and budget:

| clade | index sample's fragments | Arm A threshold | Arm A limit | Arm B threshold | Arm B limit |
| --- | --- | --- | --- | --- | --- |
| Coronaviridae | 25,928 (41.7%) | 12 reads | 7.7 /M | bits 250 | ≤0.77 /M |
| Orthocoronavirinae | 20,657 (33.3%) | 12 reads | 7.7 /M | bits 250 | **never alarms** |
| Betacoronavirus | 20,418 (32.9%) | 8 reads | 7.7 /M | bits 250 | **never alarms** |
| Sarbecovirus | 20,407 (32.9%) | 1 read | 0.77 /M | unfiltered | 7.7 /M |
| Nidovirales | 25,928 | 64 reads | 77 /M | bits 250 | ≤0.77 /M |
| Riboviria | 25,928 | 14,812,360 reads | never alarms | bits 250 | ≤0.77 /M |

The two ends of that table are the argument in one view. Widen past the family and one
negative's genuine *Enterovirus D68* infection — 6.4 million reads — sets a threshold
nothing can clear. Narrow past the family and Arm B stops alarming altogether, on the
sample where it rebuilt the genome, because the LCA of a novel organism's contig does not
go there.

Sarbecovirus is the seductive row: zero reads and zero contigs in 228,407,636 non-host
pairs, in both arms, so the threshold falls to one read and the limit improves tenfold.
It is not the headline, for the reason the clade decision exists: a 2019 analyst filtering
a report would have filtered it at a clade they could name, and the subgenus a virus
nobody had seen would land in is not knowable in advance. The row is in the record, and
anyone who wants to quote it has to argue for it.

## Decision 2: the false-alarm count excludes the endemic human coronaviruses

The cohort's coronavirus reads are mostly two patients' colds. `SRR13677922` has a
2,399-read *HCoV-HKU1* infection and `SRR13677878` nine reads of *OC43*; Arm B assembles
23 HKU1 contigs and three NL63 contigs at 97.8–100% amino-acid identity. The November
2019 database holds all four endemic human coronaviruses and names them correctly. An
operational pipeline that reported them would have been **right**.

| framing | Arm A threshold | cohort reads | per million | Arm A limit | Arm B threshold | Arm B limit |
| --- | --- | --- | --- | --- | --- | --- |
| **excluding endemic** | 12 | 35 | 0.153 | **7.7 /M** | bits 250 | **≤0.77 /M** |
| all coronavirus reads | 2,409 | 2,495 | 10.9 | 771 /M | no rung works | — |
| all, ≤1 alarm permitted | 58 | 2,495 | 10.9 | 77 /M | bits 500 | 77 /M |

Counting a correct answer as a false alarm costs Arm A a factor of ten to a hundred, and
costs Arm B the ability to be thresholded at all: `SRR13677922`'s HKU1 contig scores
1,188, above every rung of the ladder, so no pre-registered bitscore can exclude it. That
is not a measurement of either method. It is a measurement of how often a BALF patient
has a cold.

The objection Day 4 raised to this framing is real and is handled by a different knob.
Excluding the endemic species by name leaves an infection's *spillover*: reads and contigs
that are genuinely HKU1's but landed one taxon over. At read level the spillover is 8 of
the HKU1 patient's reads and the threshold sits at 12 — the exclusion does not have to be
clean, because the threshold absorbs it, and the two samples that set the threshold at 11
reads are somebody else's problem entirely (see
[what the spillover turned out to be](#the-spillover-is-smaller-than-this-day-scored-it)).
At contig level the fix is to attribute a contig to its **best hit** rather than to its
LCA:

| `SRR13677922`, 23 Coronaviridae contigs | by LCA | by best hit |
| --- | --- | --- |
| counted as the patient's HKU1 | 17 | 22 |
| counted as a false alarm | 6 | **1** (rat coronavirus Parker, 91.2%, bits 226) |
| arm's best "false alarm" bitscore | 494 | 226 |
| threshold the negatives then select | bits 500 | **bits 250** |
| resulting limit of detection | 77 /M | **≤0.77 /M** |

Six of those contigs are HKU1's own sequence whose LCA resolved up to Betacoronavirus,
and calling them false alarms costs the benchmark two decades of sensitivity for a
bookkeeping artefact. The best hit names the organism at 98–100% identity; that is the
question that was meant. Both attributions are scored and both are in the record.

What remains after all of it, per sample, is this — the whole measured false-alarm
surface of the benchmark:

| sample | Arm A: all | endemic | counted | Arm B: contigs | counted | best bits |
| --- | --- | --- | --- | --- | --- | --- |
| SRR13677922 | 2,408 | 2,400 | 8 | 23 | 1 | 226 |
| SRR13677700 | 57 | 46 | 11 | 3 | 0 | — |
| SRR13677975 | 11 | 0 | **11** | 0 | 0 | — |
| SRR13677589 | 3 | 1 | 2 | 0 | 0 | — |
| SRR13677740 | 3 | 1 | 2 | 0 | 0 | — |
| SRR13677744 | 1 | 1 | 1 | 0 | 0 | — |
| SRR13677878 | 11 | 11 | 0 | 0 | 0 | — |
| SRR13677977 | 1 | 1 | 0 | 0 | 0 | — |
| the other twelve | 0 | 0 | 0 | 0 | 0 | — |

35 reads in 228 million, and one contig. `SRR13677843`'s single whale-coronavirus contig
at 28.5% identity — the one Day 4 predicted from the index sample's own background — is
absent from this table because Day 5 scores DIAMOND's stock mode, which does not produce
it. That is the third decision.

### The spillover is smaller than this day scored it

Added after scoring, from the coinfection work in
[`docs/hku1-coinfection.md`](hku1-coinfection.md) (`hindsight arm-b diagnose-strain`),
which examined `SRR13677922` for an unrelated reason and settled where its stray reads
came from. **Seven of that sample's eight counted reads are its own HKU1 infection**: they
are 99.3–100% identical to a 2006 HKU1 **genotype B** genome (`DQ415911`) against
91.9–98.3% to the genotype A reference, which is the only HKU1 the dated database holds.
Kraken2 reached for rodent and bovine embecoviruses because no genotype B was in the
database to reach for — five reads to *Betacoronavirus HKU24*, one to *Betacoronavirus 1*,
one to *Rat coronavirus Parker*.

So the counted figures above are **upper bounds rather than rates**, and the correction
moves two numbers while leaving the headline alone:

| | as scored | after the strain check |
| --- | --- | --- |
| cohort's unexplained family-level reads | 35 | **28** |
| cohort's unexplained genus-level reads | 7 | **0** |
| Arm A's family-level threshold | 12 reads | 12 reads (set by two *other* samples) |
| Arm A's family-level limit | 7.7 /M | 7.7 /M |
| Arm A's genus-level threshold | 8 reads | would fall to 1 read |

The headline survives because the binding negatives were never this sample: the threshold
of 12 comes from `SRR13677700` and `SRR13677975`, which carry 11 reads of *Sparrow
coronavirus HKU17* each — the oddity Day 3 flagged as worth somebody's attention, two
unrelated samples at exactly the same count. Those two have since been chased down as
well, which is [its own section](#what-a-diagnosed-clean-cohort-would-give-and-why-it-is-not-the-headline).
**The arm's specificity is better than this day credited it with, and the quantity that is
actually unexplained is 28 reads in 228 million, not 35.**

**Arm B's one counted contig is the same thing, confirmed.** `k141_852` is **100.00%
identical at nucleotide level to `DQ415911` over its 767 bases**, against 97.42% to the
genotype A reference. So **Arm B's measured false-alarm count across the cohort is zero
contigs**, and its threshold could fall from bitscore 250 to the arm's unfiltered rule.
That does not move Arm B's limit, which already sits at the series floor; it moves how
strong the specificity claim is allowed to be.

**And it exposes the limit of this day's attribution rule.** Decision 2 attributes a
contig to its best hit, which cannot catch this contig, because *the best hit is itself
wrong*: *Rat coronavirus Parker*'s ORF1a genuinely is the closest protein in the dated
database at 91.2%, and HKU1's own ORF1ab comes second at 98.1% over 106 residues with the
**identical bitscore of 226**. Attribution by best hit picks the rat because a coin landed
that way. Best-hit attribution is still the better of the two rules available inside the
arm's own output — it is what takes the HKU1 patient's counted contigs from six to one —
but the thing that settles the last one is a nucleotide check against a genotype the
database does not hold, which is outside the arm entirely. A protein search against an
incomplete reference set can be confidently wrong about *which* relative a sequence is,
while being right that it is a coronavirus, and no rule over its own hit table fixes that.

The mechanism is general: a divergent strain of an organism the database holds only one
genotype of will scatter reads onto whatever relatives *are* there. One read in this
sample — a single *Coronavirus AcCoV-JC34* call — aligns to neither HKU1 genotype and
remains unexplained, which is why it still contributes 1 to the cohort's 28 rather than 0.

### What a diagnosed-clean cohort would give, and why it is not the headline

The sparrow pair turned out to be artefacts too, and not even the same artefact. In each
sample the eleven pairs **tile a single locus** — one coincidence counted eleven times —
and nothing in either sample's reads aligns to the sparrow genome (`NC_016992.1`), to
GRCh38, or to any protein in the dated set at an e-value of 10. The two samples share no
31-mer, so these are two independent coincidences rather than one contaminant. Kraken2's
own output shows the mechanism:

```
SRR13677700.59768  138|138  0:19 1159906:5 0:80 |:| 0:80 1159906:5 0:19
```

Five assigned k-mers in one block per mate, and the mates are reverse complements of each
other, because the insert is shorter than twice the read length. **`--minimum-hit-groups 2`
cannot protect that pair**: the two "independent" hit groups are the same block, read
twice. (From the coinfection thread; `results/hku1-coinfection.json`,
`sparrow_coronavirus_threshold_reads`.)

So every coronavirus alarm in the twenty negatives has now been accounted for, and none of
them is the thing a false-alarm rate is supposed to count — two genuine endemic
infections, one patient's own HKU1 at a genotype the database lacks, and two single-locus
artefacts. `config/day5.toml` records those diagnoses, and the record carries exactly one
number from them, `arm_a.adjudicated`:

| | as scored | adjudicated |
| --- | --- | --- |
| reads deducted | — | 29 of 35 |
| unexplained reads remaining | 35 | **6** |
| worst negative | 11 (`SRR13677975`) | 2 (`SRR13677589`) |
| Arm A threshold | 12 reads | **3 reads** |
| Arm A limit | 7.7 /M | **0.77 /M**, the series' lowest rung |
| Arm B contigs in the cohort | 1 | **0**, threshold falls to the unfiltered rule |

**That is not the headline, and the reason is a rule rather than a preference.** Every
deduction rests on a diagnostic the arm does not run; each was performed after it was
known which samples set the threshold; and all of them move the method's own number in
the method's favour. That is the shape of a threshold fitted to its cohort, and each
deduction being individually correct does not change the shape — it is what makes the
fitting persuasive. The coinfection thread reached the same conclusion from the other
side and deliberately did not apply the exclusion, which is the right instinct.

What would make this number quotable is a **diagnostic declared in advance and applied to
every alarm**, including the ones that are inconvenient to lose: align each alarm's reads
to its assigned taxon's genome, check the hit-group geometry, and exclude on a stated
rule rather than on a verdict. That is a Day 7 design and it is now the second item on
its list. Deleting the `adjudications` section from `config/day5.toml` changes no headline
number, which is the test of whether this has been kept at arm's length.

One thing the adjudicated column does establish, and it does not depend on any threshold:
**after a cohort-wide hunt, the benchmark's twenty negatives contain six coronavirus reads
nobody can explain, and zero contigs.** Both arms' specificity on this cohort is limited
by what twenty samples can show, not by anything either method did wrong.

### One read that does not reconcile, and why it is here

`SRR13677744`'s endemic-excluded count is 1 by enumerating its coronavirus taxa and 0 by
subtracting its endemic reads from its family count, because the clade counts are
cumulative over a clade while the taxon list is direct assignments. The scorer takes the
larger of the two, which can only make the arm look worse, and records both. One read in
228 million moves nothing; an unexplained discrepancy in a published table moves a
reader's confidence in all of it.

## Decision 3: drop `--very-sensitive`, and score DIAMOND's stock mode

`config/day4.toml` argued at length that this arm must run at `--very-sensitive`, because
SARS-CoV-2's proteins sit 10–25% diverged from anything a 2019 database held and running
the homology arm in a mode that cannot see homology at that distance would measure the
mode rather than the method. Day 4 measured the alternative anyway, on everything, and
the argument did not survive it:

| | `--very-sensitive` | stock | same best bitscore |
| --- | --- | --- | --- |
| index sample | 4 contigs, bits 7,120 | 1 contig, bits 7,120 | yes |
| 50,000 fragments | 4, bits 7,120 | 1, bits 7,120 | yes |
| 10,000 fragments | 5, bits 7,120 | 2, bits 7,120 | yes |
| 1,000 fragments | 24, bits 2,065 | 21, bits 2,065 | yes |
| 100 fragments | 17, bits 319 | 14, bits 319 | yes |
| 10 fragments | 4, bits 270 | **1, bits 270** | yes |
| spike-free background | 3, bits 57.8 | **0, nothing** | — |
| negatives alarming | 3 of 20 | **2 of 20** | — |

**Six comparisons on genuine signal, six identical best bitscores, and the only
differences anywhere are false positives.** At the limit — ten fragments — stock mode
returns one contig of 462 bases, the real one, where `--very-sensitive` returns four, of
which three are the background's remote homology to whale and bat coronaviruses at 28–33%
identity. It costs about a sixth as much per search.

Day 4 deliberately did not edit the config after measuring this, on the grounds that a
setting changed after seeing the result is no longer pre-registered. That instinct is
right and it is also why the decision belongs here: Day 5 is scoring, and choosing which
*already recorded* search to score is a scoring choice. Nothing was re-run. What changed:

* `config/day5.toml` sets `arm_b_sensitivity = "comparison"`, so the headline is stock.
* `config/day4.toml` has its two sensitivity values swapped, with the original argument
  left in place above them as the record of what was pre-registered. Both searches still
  run on every sample, so a future run measures the comparison in the other direction and
  the swap is reversible by swapping back.
* Day 4's write-up carries a pointer here.

A caveat the write-up should keep: this is a specificity argument made on twenty
host-depleted libraries, and the mode that found three artefacts in one spike-free
background may find more in dirtier input. The decision is recorded with its evidence
rather than asserted, which is what makes it arguable later.

One consequence for Day 7: the read-level search's 49-hour estimate for the negatives was
measured at `--very-sensitive`. At stock sensitivity it is about eight hours, which turns
"needs a bigger machine" into "fits in an overnight run" — the single measurement the
benchmark wanted and could not afford.

## The crossover that reverses when you ask about alarms

Day 4's headline cross-arm finding was that the curves cross: assembly recovers more than
twice as many reads as the classifier above ~1,000 fragments and collapses below 100,
while the classifier's per-read rate is flat at ~41% throughout. That is correct, and it
is a statement about **recovery**. It is not the benchmark's question.

| level | per million | Arm A reads | alarm | Arm B reads recovered | Arm B best bits | alarm |
| --- | --- | --- | --- | --- | --- | --- |
| 50,000 | 3,853 | 20,892 | yes | 49,077 | 7,120 | yes |
| 10,000 | 771 | 4,232 | yes | 9,824 | 7,120 | yes |
| 1,000 | 77 | 420 | yes | 936 | 2,065 | yes |
| 100 | 7.7 | 39 | yes | 41 | 319 | yes |
| 10 | 0.77 | **4** | **no** | **1** | **270** | **yes** |

At the bottom rung Arm B recovers one read of ten and Arm A recovers four, and Arm B
alarms while Arm A does not. Nothing is contradictory: an alarm needs one piece of
evidence above the cohort's noise, and a 462 bp contig whose protein matches a
coronavirus replicase at bitscore 270 is that, while four reads in a cohort whose
negatives reach eleven is not. **Recovery and detection are different questions, and the
arm that recovers less can be the arm that alarms.**

This is why the two arms get separate statistics rather than a shared one, and it is the
sentence Day 6 should build the cross-arm section around.

## The series cannot see below its own background

The load-bearing finding of this day, and it bounds every number above.

Day 3 reported 78 background fragments that Kraken2 placed in SARS-related coronavirus,
and showed that 71 of them align concordantly to `NC_045512.2` with every mate below the
label's MAPQ 20 cut-off: **real SARS-CoV-2 reads that the ground-truth label declined to
count.** 90 land in Coronaviridae. The dilution series was cut from that background, and
those reads were never diluted, so they sit at full strength in every rung:

| rung | labelled | labelled /M | resident | effective dose | effective /M | undiluted share |
| --- | --- | --- | --- | --- | --- | --- |
| viral050000 | 50,000 | 3,853 | 90 | 50,090 | 3,860 | 0.2% |
| viral010000 | 10,000 | 771 | 90 | 10,090 | 778 | 0.9% |
| viral001000 | 1,000 | 77.1 | 90 | 1,090 | 84.0 | 8.3% |
| viral000100 | 100 | 7.71 | 90 | 190 | 14.6 | 47% |
| viral000010 | 10 | 0.77 | 90 | 100 | 7.71 | **90%** |

The bottom rung's effective dose equals the next rung's labelled dose. In other words
**the series has four usable rungs, not five**, and its floor is about 7 viral reads per
million — the dose its own background already carries.

Three consequences, and they are all in the headline table above:

1. **Scoring uses labelled detections**, which is conservative for both arms: it credits
   neither with the resident virus. Arm A's operational count — what a 2019 analyst's
   report would actually have shown — is 94 Coronaviridae reads at the bottom rung, which
   alarms. Both numbers are recorded; the limit is quoted from the conservative one.
2. **Arm B's limit is unresolved.** It alarms at every rung, so the measurement is an
   upper bound at the floor, not a point.
3. **The two arms' limits differ by less than the floor's width.** 15 against ≤7.7
   effective reads per million is well inside a factor of two, and this series cannot
   separate them. Any claim that one arm is more sensitive than the other at the limit is
   not supported by this benchmark, whatever the recovery curves do.

The fix is cheap and is a Day 7 item: cut a second dilution series into a **negative**
library's background rather than the index sample's own. The negatives are confirmed free
of sarbecovirus reads across 228 million pairs, so a series spiked into one of them has a
dose axis that goes down as far as the spike does, and two more decades become visible.
It needs no new sequencing and no new fetch — one negative's reads and the index sample's
labelled viral fragments, both already defined by Day 2's frozen seeds.

## Arm B's read half is a mechanism result, not a limit

The read-level search is the most sensitive thing in the benchmark per read, and quoting
it as a sensitivity result would be the one serious mistake available on this day.

| level | detected | rate | background floor | signal : floor | clears it |
| --- | --- | --- | --- | --- | --- |
| 50,000 | 46,729 | 93.5% | 177 | 264 : 1 | yes |
| 10,000 | 9,344 | 93.4% | 177 | 52.8 : 1 | yes |
| 1,000 | 942 | 94.2% | 177 | 5.3 : 1 | yes |
| 100 | 92 | 92.0% | 177 | 0.52 : 1 | **no** |
| 10 | 9 | 90.0% | 177 | 0.05 : 1 | **no** |

A flat 93% per-read detection rate, all the way down, **and below its own noise from a
hundred fragments downward.** Nine true reads among 177 false ones is not a detection by
any operational standard, and the flat rate is a statement about homology — the protein
similarity was recognisable at every dose — rather than about detectability. It also has
no false-alarm rate from an independent cohort, because the read half is off for the
negatives, so it cannot be given a threshold at all. It is scored nowhere in this day's
headline and belongs in Day 6's mechanism section.

## What twenty negatives can and cannot establish

The specificity figures above are zero, and zero out of twenty is a weaker statement than
it looks.

* **Per sample**: 0 alarms in 20 samples bounds the per-sample false-alarm rate at
  **under 14%** with 95% confidence, and no better. The benchmark's promise of "≤1 false
  alarm per N samples" is met at N = 20 by measurement; anything stronger needs more of
  the 408 libraries PRJNA699976 holds.
* **Per read**: 0 sarbecovirus reads in 228,407,636 non-host pairs bounds the per-read
  rate at about **0.013 per million**, which is a genuinely tight number and is the one
  worth quoting for the read-level statistic.
* **The cohort is cleaner than the positive.** The negatives were host-depleted before
  submission — 0.075% human against the index sample's 15.8% — so these specificity
  figures are measured on unusually clean input, and they flatter Arm B more than Arm A,
  because a translated search at 30% identity is a looser instrument than exact 35-mer
  matching. Day 4 set this out and it stands; the host-derived false-positive test that
  would settle it is a Day 7 item costing about an hour.
* **n = 1 on the positive.** The dilution series is the generalisation axis, not the
  sample, and the write-up has to say so.

## Cost

| | |
| --- | --- |
| compute | none |
| wall clock | ~1 s |
| inputs read | 2 JSON files, 573 KB |
| operating points computed | 128 |
| tools invoked | 0 |

A day that costs a second is worth stating as a result rather than as a footnote. The
reason this day is free is that Days 3 and 4 wrote down every clade, every level and every
negative's per-taxon detail instead of only their own headline, and the practical test of
whether a benchmark's thresholds are honest is whether moving one costs anybody a machine.
Here it costs a flag.

## For Day 6 and Day 7

**For Day 6**, the write-up:

* The headline is "both arms would have caught it, down to about 8 viral reads per
  million, with no false alarm in twenty pre-COVID samples" — and the limit is the
  series' floor, not either method's. Say that in the abstract rather than in the
  limitations.
* The cross-arm section should be built on recovery-versus-detection, not on the
  crossover. The crossover is real and it does not change who alarms.
* The three decisions here are the methods section's spine: the clade (and the
  BM48-31 taxonomy effect that forces it), the false-alarm framing, and the sensitivity
  setting that was pre-registered and then dropped on its own evidence.
* Quote the per-read specificity bound (0.013 per million) and the per-sample one (<14%)
  as the two different things they are.
* The unexplained family-level figure is **28 reads in 228 million, not 35**, and after
  the sparrow pair was diagnosed as well, **six** — but quote the as-scored threshold and
  limit, not the adjudicated ones, and say why. The diagnosed-clean column (threshold 3
  reads, limit 0.77 /M) belongs in the paper as a sentence about what a pre-declared
  diagnostic pass would be worth, not as the result.

**For Day 7**, in rough order of what the benchmark gains:

1. **A series in a negative's background.** Two more decades of dose axis, no new
   sequencing, and it is the only way to find out where either arm actually fails.
2. **A pre-declared adjudication pass over every alarm.** Each of the cohort's alarms was
   diagnosed one at a time, after the fact, and the result would improve Arm A's limit by
   a decade — which is exactly why it cannot be quoted. Declare the diagnostic first
   (align each alarm's reads to its assigned taxon's genome; check the hit-group
   geometry; exclude on a stated rule), apply it to all twenty samples, and the
   adjudicated column becomes a result instead of a footnote. The `adjudications` table
   in `config/day5.toml` is the scaffolding.
3. **The read-level false-alarm rate**, now about eight hours rather than 49 at stock
   sensitivity — the measurement that would let the read half be scored at all.
4. **The host-derived false-positive test**, one bowtie2 pass over the index sample's QC
   reads to recover its 2,053,621 host-aligned pairs, run through both arms.
5. A mapping-quality curve for the ground-truth label instead of one MAPQ 20 threshold.
   The resident background above *is* that threshold's cost, measured.
6. The items Days 3 and 4 carried: `bowtie2 --reorder`, bwa-mem against bowtie2 for the
   host filter, a combined Day 3 + Day 4 negatives runner, and the HKU1 patient's four
   non-aligning contigs.
