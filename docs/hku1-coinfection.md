# What else did the HKU1 patient have?

Day 4 left a loose end. `SRR13677922`, one of the twenty pre-COVID negative controls, has
a genuine near-genome-wide *Human coronavirus HKU1* infection — and four of its 23
Coronaviridae contigs had good HKU1 protein hits at 76 to 88% amino-acid identity and no
nucleotide alignment to the dated HKU1 reference `NC_006577.2` at all. That is either a
very diverged region of the same virus or a second betacoronavirus in the same patient,
and one reference cannot tell those apart. The question was carried into the Day 7 list
and then pulled forward, in the broader form: what, if anything, else was in this
patient's library?

The answer is **nothing else, and the virus is one strain**: human coronavirus HKU1,
genotype B. Along the way two other things fell out, and both matter more than the
original question. Arm A's non-endemic betacoronavirus false alarms in this sample are not
false alarms, which changes a number Day 5 was going to use. And the single largest viral
call anywhere in the negative-control cohort is not a virus at all.

This day touches no post-cutoff sequence anywhere. Where it reaches outside the dated
database it does so for sequence published in 2004 and 2006, which the database lacks
because of RefSeq's curation policy rather than because of its date — and every reference
used is checked against the cutoff at NCBI, with its publication date recorded in the
output.

---

## The sample, reproduced

The Day 4 container was gone, so the sample was re-staged and re-assembled from scratch in
a fresh one. Every figure reproduced exactly:

| | Day 4 | this run |
| --- | --- | --- |
| raw pairs | 14,756,419 | 14,756,419 |
| post-QC pairs | 14,752,099 | 14,752,099 |
| host pairs | 704 | 704 |
| non-host pairs | 14,751,395 | 14,751,395 |
| contigs | 7,438 | 7,438 |
| assembly bases | 4,539,007 | 4,539,007 |
| longest contig | 14,415 | 14,415 |
| N50 | 653 | 653 |
| DIAMOND hit rows | 10,107 | 10,107 |

**One thing does not reproduce: MEGAHIT's contig names.** The assembly is identical and the
names are shifted — Day 4's `k141_681` is this run's `k141_628`, the same 2,509 bp
sequence. Contig identifiers are therefore not citable across runs, and the Day 7
reproducibility pass should say so rather than treat a name as a reference. This document
uses this run's names and gives lengths alongside them.

---

## One virus: HKU1 genotype B

### Why one reference could not answer it

RefSeq curates one genome per species. Its single HKU1 genome, `NC_006577.2`, is
**genotype A**. HKU1 has three genotypes; A and B differ most in the two surface
glycoprotein genes, hemagglutinin-esterase and spike, and genotype C is a recombinant.
All three were published in 2006, in the same paper that described the genotypes.

So the four unexplained contigs had a mundane possible explanation that the dated database
was structurally incapable of offering: they could be the parts of an ordinary HKU1 strain
that the one curated HKU1 genome happens not to represent. Testing that needs a panel of
references rather than one, which is what `hindsight arm-b diagnose-strain` now does.

### The panel

```console
$ hindsight arm-b diagnose-strain --sample SRR13677922 \
    --reference NC_006577.2 --reference AY597011.1 --reference DQ415914.1 \
    --reference DQ415911.1 --reference DQ415902.1 \
    --reference DQ415913.1 --reference DQ415897.1
```

| reference | genotype | published | contigs | aligned bases | identity | genome covered |
| --- | --- | --- | --- | --- | --- | --- |
| `NC_006577.2` (RefSeq) | A | 2004-12-28 | 22/23 | 16,713 | 94.72% | 55.42% |
| `AY597011.1` | A | 2004-12-21 | 22/23 | 16,726 | 94.73% | 55.43% |
| `DQ415914.1` (N18) | A | 2006-06-30 | 22/23 | 16,737 | 94.56% | 55.57% |
| **`DQ415911.1` (N15)** | **B** | **2006-06-30** | **23/23** | **17,685** | **99.36%** | **58.76%** |
| `DQ415902.1` (N25) | B | 2006-06-30 | 23/23 | 17,714 | 99.10% | 58.96% |
| `DQ415913.1` (N17) | C | 2006-06-30 | 23/23 | 17,741 | 97.18% | 59.17% |
| `DQ415897.1` (N20) | C | 2006-06-30 | 23/23 | 17,681 | 97.47% | 58.85% |

Every contig is accounted for by one genotype B genome at 99.36% nucleotide identity, which
is ordinary strain-level variation between a 2016-17 Wuhan sample and a 2004-05 Hong Kong
isolate. The split is by genotype and not by a single odd reference: all three genotype A
genomes give the same 22 contigs at the same ~94.7%, and both genotype B genomes give 23 at
~99.2%.

### What the four contigs were

Three of them are the genotype-defining genes, and the fourth is not biology at all.

| contig | bp | k-mer coverage | protein hit | vs A | vs B |
| --- | --- | --- | --- | --- | --- |
| `k141_1504` | 1,292 | 7.0 | spike, aa 195-625, 77.5% | 86.18% | 98.99% |
| `k141_1670` | 835 | 38.6 | spike, aa 847-1026, 88.9% | 88.32% | 99.49% |
| `k141_61` | 535 | 10.0 | spike, aa 1205-1356, 76.3% | 85.29% | 98.46% |
| `k141_273` | 490 | 3.0 | spike, aa 1110-1210, 85.1% | 89.93% | 98.89% |
| `k141_1480` | 392 | 9.0 | HE, aa 91-217, 80.8% | *no alignment* | 99.49% |
| `k141_1616` | 324 | 3.0 | nsp3, aa 56-130, 88.0% | 99.48% | 100.00% |

The spike contigs tile the protein in order — residues 195-625, 847-1026, 1110-1210,
1205-1356 — with a six-residue overlap where two of them meet and no region covered twice.
In genome coordinates they fall exactly into the gaps left by the genotype A nucleotide
tiling. Hemagglutinin-esterase is the only part of the genome that will not align to
genotype A at all at any preset.

`k141_1616` is an artefact. Its two halves align to the *same* region of either genotype in
*opposite* orientations — query 0-191 forward onto 2,790-2,981, query 187-324 reverse onto
2,873-3,010 — so it is an inverted-repeat join, from the assembler or from a template
switch during library construction. It is HKU1 ORF1ab at 99 to 100% nucleotide identity and
was never diverged at all; its 88% protein identity over 75 residues is the same artefact
seen in translation, because only one half is in frame.

**Which brings out a second lesson about the Day 4 diagnostic.** `minimap2 -x asm20` could
not place `k141_1616`, `k141_61` or `k141_273` against genotype A, and only one of those
three is genuinely divergent sequence. `asm20` needs a longer seed chain than a 300 to 500
bp contig with a split alignment can offer, so on short contigs **"does not align" is a
preset limit and not a finding**. `diagnose-strain` therefore retries every unplaced contig
with `-x sr` and records which preset placed it; with the fallback, 22 of 23 contigs align
to genotype A, and the discriminator between the genotypes becomes identity — 94.7% against
99.4% — rather than alignability.

### Two checks that would have caught a second organism

Neither fires.

**Coverage classes.** MEGAHIT's own k-mer coverage for the 23 contigs runs from 2 to 102,
and the six contigs above sit at 3, 3, 7, 9, 10 and 38.6 — inside that range, not beside
it. A second virus in the same library would have its own abundance, and two genomes at
coincidentally interleaved depths is a stretch.

**Duplicated coverage.** A fragmented assembly of one genome tiles its reference: contigs
abut, and where they overlap they agree, because they are the same sequence. Two organisms
produce two contigs for the same gene whose identities to any one reference differ. Against
genotype B, no region of the genome is covered twice by contigs that disagree by more than
three points. `diagnose-strain` reports this test explicitly, because it is the one that
actually separates the hypotheses.

And the simplest argument of all: there is no second replicase or nucleocapsid contig
anywhere in the 7,438. A betacoronavirus present deeply enough to assemble 3 kb of spike
would have left those first — they are longer and more conserved.

---

## Arm A's embecovirus false alarms in this sample are not false alarms

Day 3 reported that excluding the endemic human coronaviruses left **7 betacoronavirus
reads** across the whole 20-sample cohort, 0.031 per million non-host pairs, every one of
them in this patient and described as embecovirus spillover. Day 5 counts this sample at
**8** reads once an alphacoronavirus read is included at family level. Every one of them
was checked:

| reads | Kraken2's call | identity to genotype B `DQ415911.1` | to genotype A | verdict |
| --- | --- | --- | --- | --- |
| 5 | *Betacoronavirus HKU24* | 99.3–100% | 91.9–98.3% | the patient's HKU1 |
| 1 | *Rat coronavirus Parker* | 100% | 93.2% | the patient's HKU1 |
| 1 | *Rabbit coronavirus HKU14* | 100% | 93.1–94.5% | the patient's HKU1 |
| 1 | *Coronavirus AcCoV-JC34* | no alignment | no alignment | **unexplained** |

So **7 of the 8 are this patient's own HKU1**, and all 7 are the non-endemic
*betacoronavirus* reads — the genus-level figure. They sit in three places on the genome,
around positions 15,210–15,400, 17,530–18,030 and 21,440–21,620: ORF1ab, and the ORF1ab/HE
boundary. Kraken2 reached for rodent, bovine and rabbit embecoviruses because for those
35-mers a rodent embecovirus genuinely is the closest thing the database holds. The eighth
read is an alphacoronavirus call that aligns to neither HKU1 genotype and stays on the
books as one unexplained read in 14.75 million.

**And Arm B's single counted false-alarm contig in the whole cohort is the same thing.**
Day 5's threshold of bitscore 250 is selected by one contig: `k141_852`, 767 bp, best hit
*Rat coronavirus Parker* ORF1a at 91.2% amino-acid identity, bitscore 226. It is **100.00%
identical at nucleotide level to `DQ415911.1`** and 97.42% to genotype A. Day 5's fix of
attributing a contig by best hit rather than by LCA cannot catch this one, because the best
hit is itself wrong — a rat coronavirus really is the closest protein in the database, and
DIAMOND's second hit is HKU1's own ORF1ab at 98.1% with the identical bitscore of 226. With
it accounted for, **Arm B's measured false-alarm count across the negatives is zero
contigs**, and nothing in the cohort selects a bitscore threshold above the arm's lowest
pre-registered rung. The limit of detection does not move, because it already sits at the
dilution series' floor; the specificity claim gets stronger.

**Consequences for Day 5.** The cohort's non-endemic betacoronavirus count at genus level
is **zero**, not 7, and the family-level count is 28 reads rather than 35; both are upper
bounds rather than rates. Arm A's threshold of 12 reads **does not move**, because it was
set by `SRR13677700` and `SRR13677975` at 11 reads each rather than by this patient — but
those two samples are 11 reads of *Sparrow coronavirus HKU17* each, in two unrelated
patients, which is the same signature. If they are the same mechanism the threshold drops
and Arm A's limit of detection improves. **That is the one remaining thread worth pulling,
and corrections of this kind can only ever lower a threshold**, so the published figures
are conservative in the right direction.

---

## The reads that set Arm A's threshold are not coronavirus reads

Arm A's alarm rule is one read above the worst negative, so the whole threshold rests on
the two samples that tie for worst: `SRR13677700` and `SRR13677975`, at 11 reads each of
*Sparrow coronavirus HKU17*, an avian deltacoronavirus. Two unrelated patients, the same
taxon, the same count — Day 3 flagged it as suspicious at the time. It is worse than
suspicious.

```console
$ hindsight arm-a diagnose-taxon --sample SRR13677700 --taxid 1159906 \
    --reference NC_016992.1
11 fragments were assigned to Sparrow coronavirus HKU17, and **not one of them aligns
to NC_016992.1**. The call rests on shared k-mers rather than on the organism's
sequence being present: whatever these reads are, they are not this.
```

Both samples give that same result. Five independent checks agree:

| check | `SRR13677700` | `SRR13677975` |
| --- | --- | --- |
| reads aligning to `NC_016992.1` (bowtie2 `--local`) | 0 of 22 | 0 of 22 |
| longest exact match to that genome | **18 bp** | **none at 18 bp** |
| DIAMOND `blastx` hits in the whole dated viral protein set, at e-value 10 | 0 | 0 |
| reads aligning to GRCh38 (`--very-sensitive-local`) | 0 of 22 | 0 of 22 |
| 31-mers shared with the *other* sample's reads | 0 | 0 |

So these are not sparrow coronavirus, not any virus the November 2019 database holds, not
human, and not even the same sequence in the two samples. In each sample the 11 pairs
overlap each other and tile one locus — 11 of 11 pairs and 9 of 11 have mates overlapping
by 80 bp or more — so each sample contributes one locus counted eleven times, and the two
loci are unrelated to each other. They are two independent coincidences.

**The mechanism, from Kraken2's own output.** Every one of the 22 pairs looks like this:

```
SRR13677700.59768   138|138   0:19 1159906:5 0:80 |:| 0:80 1159906:5 0:19
```

Exactly **five** k-mers assigned, in one block, in each mate — and the mates are reverse
complements of each other, so that is the *same* five k-mers seen twice. Kraken2's
`--minimum-hit-groups 2` is meant to stop a single stretch of coincidental similarity from
classifying a read, and here it does not, because a read pair whose insert is shorter than
twice the read length presents one hit group twice and clears the threshold. The
underlying similarity is an 18 bp exact match in one sample and less than that in the
other, reaching a 31-base minimizer through Kraken2's spaced seeds.

**What this does to the benchmark.** The arm's unmodified output gives a threshold of 12
reads and a limit of detection of 7.71 viral reads per million. With these 22 reads
recognised as artefacts, the worst remaining negative carries **2** reads, and the bottom
rung of the dilution series produces 3 — so the threshold falls by most of a decade and
the arm detects a rung it previously missed.

**That is not this document's call to make, for a reason worth stating.** The exclusion
rests on a diagnostic the arm does not itself run, applied after seeing which samples set
the threshold. Scoring an arm on its own output is what makes the benchmark honest, and a
post-hoc exclusion that improves a method's headline is exactly the move a reader should
distrust. The defensible form is to report both: the arm as scored, and the arm after a
**pre-declared** diagnostic pass that is applied to every alarm rather than to the
inconvenient ones. The numbers are here, the command is in the repository, and the
decision belongs to Day 5 and to Imo.

---

## The largest viral call in the cohort is a cloned gene fragment

Arm A's single largest viral call in this sample is **11,227 reads of *Equine infectious
anemia virus***, a horse lentivirus, in a human respiratory library — 58% of everything the
classifier called viral here, and more than four times the count of the HKU1 infection the
patient actually had. The natural readings are host-derived retroviral sequence or bare
homology to reverse transcriptase, and both would make it a false positive to count against
the arm.

Both are wrong, and a read count cannot see why. **Where the reads sit on the reference
can**, which is what `hindsight arm-a diagnose-taxon` now reports.

```console
$ hindsight arm-a diagnose-taxon --sample SRR13677922 --taxid 11665 \
    --reference NC_001450.1
```

| | `SRR13677922` | `SRR13677559` |
| --- | --- | --- |
| fragments assigned | 11,227 | 499,690 |
| reads aligning to `NC_001450.1` | 20,857 of 22,454 | 994,175 |
| identity | 98.81% | 98.69% |
| genome covered | 20.19% | 20.28% |
| covered windows | 1 | 1 |
| window | 2,293-3,980 | 2,286-3,980 |
| median depth in window | 1,181 | 57,802 |
| bases for coverage to stop, each end | 1 and 1 | 1 and 1 |

The reads are 98.8% identical to EIAV, so they are not homologues of it — they are its
sequence. But they occupy a single 1,688 bp window and nothing else: zero coverage of *gag*
(525-1,985), of *env* (5,372-7,951) or of either long terminal repeat, and the depth goes
from nothing to over a thousandfold within one base at both ends. An organism replicating
in a patient is not covered that way. This is a cloned fragment — a vector, a plasmid or a
spiked control. The window sits inside *pol* (1,766-5,182) and spans its
reverse-transcriptase region, which makes carry-over from a recombinant enzyme preparation
or a lentiviral vector the obvious suspects, though nothing here identifies the product.

**It is in all twenty negatives.** Re-reading Day 3's record, every one of the 20 samples
has an EIAV call: 3,284,717 reads in total, from 291 reads in the shallowest sample to
499,690 in `SRR13677559`. And in the two samples examined here the window is the *same
window*, to within seven bases at one end and exactly at the other, in two unrelated
patients. That is a shared reagent, not two infections.

**Why this belongs in the write-up.** It does not move the coronavirus detection limits:
the benchmark's thresholds are per clade, and EIAV is a lentivirus. What it does do is
falsify a weaker claim the paper might otherwise have made in passing — that the negatives
establish a false-alarm rate for a *pathogen-agnostic* pipeline. A pathogen-agnostic
pipeline run on this cohort reports a notifiable equine pathogen in 20 of 20 clinical
samples, at up to half a million reads, and any alarm rule keyed on read count fires every
time. The honest form of the claim is that the negatives establish a false-alarm rate for
coronavirus detection, and that a pathogen-agnostic alarm on this cohort is dominated by
laboratory sequence that a shape-of-coverage check removes and a read count does not.

---

## What else was in the library

19,239 of 14,751,395 non-host pairs were classified at all — 0.13%. Against a viral-only
database the rest is host and bacterial RNA, and invisible here.

| reads | call | reading |
| --- | --- | --- |
| 11,227 | *Equine infectious anemia virus* | cloned *pol* fragment, above |
| 2,399 | *Human coronavirus HKU1* | the infection |
| 7 | rodent and bovine embecoviruses | the same infection, misassigned |
| 4,575 | Caudovirales (724 *Escherichia* phage HK639, 1,945 *Hpunavirus*, …) | bacteriophage of the respiratory and oral flora |
| 401 | *Choristoneura fumiferana granulovirus* | an insect virus; environmental or spurious |
| 411 | giant algal and amoebal viruses (*Cafeteria roenbergensis* virus, *Pandoravirus*, …) | spurious; these taxa recur across the cohort |
| 340 | unclassified RNA viruses ShiM-2016, *Wenling narna-like virus 5*, *Hubei picorna-like virus 4* | invertebrate virome reference sequences; spurious |
| 21 | *Dyoupsilonpapillomavirus 1* | a bird papillomavirus; spurious |
| 9 | *Enterovirus D68* | noise-level |
| 2 | *Rhinovirus C* | noise-level |

Nothing else of clinical interest is present at any depth. No influenza, no RSV or
metapneumovirus, no adenovirus, no herpesvirus above single digits. The nine *Enterovirus
D68* reads are worth one line of context: a different negative in this cohort came back 43%
classified and almost entirely *Enterovirus D68*, so that taxon is genuinely circulating in
these samples — nine reads in this one is not an infection.

**The scope limit, stated plainly.** The benchmark's dated database is viral-only curated
RefSeq. This analysis therefore says what viruses were in the patient and says nothing
whatever about bacteria or fungi. For a bronchoalveolar lavage from a patient with
pneumonia that is a real gap: a bacterial co-pathogen is entirely plausible and nothing
here excludes one. Closing it would need a dated bacterial database, which is a different
and much larger build than the viral one — `hindsight build` fetches 12,179 nucleotide
records for the viral set — and is out of scope for the week.

---

## Reproducing this

Both commands are diagnostics. They run after an arm has finished, nothing they compute
feeds back into any arm's score, and both refuse a reference published on or after the
cutoff unless told otherwise in as many words.

```sh
# the strain panel: one diverged strain, or more than one organism?
hindsight arm-b diagnose-strain --sample SRR13677922 \
  --reference NC_006577.2 --reference DQ415911.1 --reference DQ415913.1

# is a classifier's call an organism or a laboratory sequence?
hindsight arm-a diagnose-taxon --sample SRR13677922 --taxid 11665 \
  --reference NC_001450.1
```

Both need the sample's reads on disk. `prune_after_search` and `prune_after_classify` delete
a negative's reads once its arm has run, so diagnosing a negative means re-staging it —
about 10 minutes for this sample, then 22 minutes for the assembly if the contig half is
needed too. The runs behind this document used copies of `config/day3.toml` and
`config/day4.toml` with pruning turned off; the committed configs were deliberately left
alone, since a pre-registered setting changed after seeing a result is no longer
pre-registered.

Results are in `results/hku1-coinfection.json`.

## Carried forward

- **Day 5 should use zero, not 7**, for the cohort's non-endemic betacoronavirus false
  alarms at genus level, and should quote it as a bound.
- **The other 27 family-level non-endemic reads** in other samples deserve the same
  treatment before Day 6 quotes that figure.
- **The EIAV contaminant belongs in the limitations**, as the reason the negatives bound
  coronavirus false alarms rather than pathogen-agnostic ones.
- **RefSeq's one-genome-per-species policy is a limitation in its own right**: an ordinary
  genotype of an endemic virus read as a divergent betacoronavirus of unknown type. The
  arms' false-alarm behaviour partly measures that curation policy, not what was knowable
  in November 2019 — Wu et al. searched nr/nt and did not have this blind spot.
- **Contig names are not stable across MEGAHIT runs** even when the assembly is identical.
