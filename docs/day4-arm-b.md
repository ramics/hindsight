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
| ground-truth viral fragments recovered | 60,926 of 62,117 (**98.1%**) |
| fragments recruited that are not in the label | 307 |
| detection, reads per million post-QC pairs | 4,719 |

The arm's own unit is a contig, but the benchmark's unit is reads, so the contigs the arm
placed in Coronaviridae are made into a bowtie2 index and the non-host reads are mapped
back onto them. A read that aligns to a contig the arm called coronaviral is this arm's
equivalent of a read Kraken2 called coronaviral, and the two arms become comparable.

## The detection curve, and where it crosses Arm A

Each dilution level is assembled whole — there is no composition shortcut here, for
reasons in [its own section](#why-every-level-costs-a-full-assembly) — searched, and
scored by mapping reads back.

| viral fragments | per million post-QC | Arm B detected | Arm B rate | Arm A rate |
| --- | --- | --- | --- | --- |
| 50,000 | 3,853 | 49,078 | **98.2%** | 41.8% |
| 10,000 | 771 | 9,829 | **98.3%** | 42.3% |
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
On this sample the extra sensitivity bought nothing but noise, at six times the cost — and
the write-up says so against its own configuration.

That is not yet the whole answer, because the index sample cannot settle it: with 62,117
viral fragments either mode finds a 29 kb contig easily. The question is whether the extra
sensitivity buys detections at the *limit*, where the arm is choosing between 462 bases and
nothing. The comparison therefore runs at every dilution level too, and the per-level
numbers are in the record under `sensitivity_comparison`. If stock mode holds up at 10 and
100 fragments, Day 7 should drop `--very-sensitive` and save six times on every search.

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
- the bowtie2 index over the flagged contigs is rebuildable and is not kept
- each negative's reads are deleted once it has been searched

What is kept is the arm's output: contigs, the gzipped tabular hits, and the record.

## Cost

Measured on 4 cores and 15 GiB, which is what the container has.

| step | time |
| --- | --- |
| one 10.9M-pair MEGAHIT assembly | 34.6 min |
| contig `blastx`, 27,732 contigs, `--very-sensitive` | 0.8 min |
| mapping reads back onto the flagged contigs | 2.1 min |
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

**For Day 7.**

- Drop `--very-sensitive` if the per-level comparison confirms the index sample's verdict.
  Six times faster, and on present evidence cleaner.
- A combined Day 3 + Day 4 runner would classify and assemble each negative in one visit and
  save a whole fetch-and-host-filter pass (163 min). Day 3 had already deleted its negatives'
  reads and ran in a different container, so Day 4 re-staged them.
- Record order-independent read-ID digests. Day 2's counts reproduce exactly in a fresh
  container but its per-file sha256 digests do not, because bowtie2 runs without `--reorder`.
- Compare bwa-mem against bowtie2 for the host filter, carried over from Day 3.
