# hindsight

**Could we have caught the last pandemic?**

hindsight rebuilds a viral reference database as it actually stood on a date you choose: the real record set, at the accession versions that were public at the time. Then it runs detection pipelines against the index SARS-CoV-2 sample using only that.

```console
$ hindsight build --date 2019-12-01 --out db/november-2019
cutoff 2019-12-01 -> release 97 (2019-11-04); first excluded release: release 98 (2020-01-06)
manifest: 12179 nucleotide + 373129 protein records across 9305 taxa
taxonomy: taxdmp snapshot 2019-11-01
```

## The result so far

Kraken2, at its own default settings, against a database that had never seen SARS-CoV-2, RaTG13 or any relative described after 4 November 2019:

```console
$ scripts/day3.sh
index sample: 25,928 of 62,117 ground-truth viral fragments in Coronaviridae (41.7%)
  viral050000: 20,892 of 50,000 detected      viral001000:    420 of  1,000 detected
  viral010000:  4,232 of 10,000 detected      viral000100:     39 of    100 detected
                                              viral000010:      4 of     10 detected
```

**A standard read classifier placed 25,928 of the index sample's viral reads in Coronaviridae.** It could not name the virus, and it did not need to. Detection rate holds flat at about 41% down the dilution series, because Kraken2 classifies each fragment independently, so what falls is the absolute count: four coronavirus reads survive at 0.77 viral reads per million.

One finding from that arm is the argument for the whole tool, and it was found rather than predicted. 4,922 reads were assigned to *Bat coronavirus BM48-31/BGR/2008*, the nearest relative the database holds. The **November 2019 taxonomy files that genome under "unclassified coronaviruses"**, where today's tree puts it in Sarbecovirus. So family-level detection is 25,928 reads and genus-level is 20,418, and a benchmark that dates its sequences while keeping a modern taxonomy cannot see the difference. A classifier's output is a taxon, and taxa move.

[`docs/day3-arm-a.md`](docs/day3-arm-a.md) has the arm in full.

## Why a time-correct database is hard

The question is whether an automated, pathogen-agnostic pipeline would have flagged SARS-CoV-2 in the index sample, using only what was knowable before it existed. Answering it needs a database holding everything a real analyst could have consulted in December 2019, and only that. The easy routes all leak:

| Route | What leaks |
|---|---|
| Today's database, SARS-CoV-2 deleted | RaTG13, the pangolin coronaviruses and thousands of pandemic-era sequences stay. You have deleted the answer and left the answer key. |
| Today's database, filtered by deposition date | The records themselves have moved. A 2007 accession may now be at version 3, and version 3 may postdate the cutoff. Taxon names and tree shape have changed too. |
| A snapshot someone kept | Reproducible for one date, if that. |

hindsight takes a route that works because NCBI never threw anything away. **RefSeq release catalogs are archived back to release 1 (June 2003)**, and each lists every accession at the version it held in that release. So hindsight picks the last release issued strictly before your cutoff, extracts the viral accession set with versions, fetches each record at that exact version, pairs it with the taxonomy dump archived nearest before the release date, and writes a `provenance.json` with every count, every withdrawn record and a SHA-256 of every output.

For a December 2019 cutoff, release 97:

| | |
| --- | --- |
| Viral nucleotide / protein records | 12,179 / 373,129 |
| Distinct taxa, coronavirus genomes | 9,305 / 45 |
| SARS-CoV | present, `NC_004718.3` |
| Bat SARS-related and other bat CoV genomes | 14, including BM48-31/BGR/2008 |
| SARS-CoV-2, RaTG13 | absent |

The nearest relatives a December 2019 analyst could reach are in, at 79–88% identity to the index sample, which is what Wu et al. actually had. Everything learned from the pandemic is out, because those records had not been written yet. `tests/test_leakage.py` asserts all of it, and its network-marked test streams the full 385,308-row release 97 viral set from NCBI to check.

## The samples

`hindsight samples` prepares the other half of the benchmark: the index sample (`SRR10971381`, the Wu et al. BALF metatranscriptome) and 20 time-matched pre-COVID Wuhan BALF negatives from `PRJNA699976`, with read QC, host filtering, ground-truth labelling and a frozen dilution series down to 10 viral fragments. Every stage resumes, every seed is recorded.

The label is the one place the pipeline may know something from 2020: fragments are labelled against the modern reference. The detection arms see only the November 2019 database. Knowing the answer is the benchmark's job, and it is what the arms are barred from doing. [`docs/day2-samples.md`](docs/day2-samples.md) has the account, including two properties of the negative set that travel with every number derived from it.

## Install and use

```bash
uv sync --locked --all-extras     # or: pip install -e ".[dev]"
```

Set `NCBI_API_KEY` if you have one: it raises the E-utilities ceiling from 3/s to 10/s and takes a full build from hours to under one.

```bash
hindsight resolve 2019-12-01                 # which release a date selects, and what that excludes
hindsight build --date 2019-12-01 --out db/nov2019 --manifest-only   # accession set only, no fetch
hindsight build --date 2019-12-01 --out db/nov2019                   # nucleotide, protein, taxonomy
hindsight verify db/nov2019                  # re-check a build against its own provenance
hindsight releases --limit 10                # archived releases, or: hindsight taxdumps
```

A build writes `manifest.tsv`, `viral.fna`, `viral.faa`, `seqid2taxid.map` in Kraken2's format, a contemporaneous `taxonomy/`, `withdrawn.tsv` and `provenance.json`. Builds are deterministic: two builds of the same cutoff are byte-identical. Interrupt and re-run, and it resumes from the last completed shard.

## Limits

Four properties bear on how a result from this database may be worded: releases land roughly two months apart, so a cutoff selects the last release before it; this is RefSeq's curated viral set and not the whole of GenBank; a record NCBI has since suppressed is reported in `withdrawn.tsv` instead of silently dropped; and a record returned at the wrong version is refused rather than tolerated.

[`docs/caveats.md`](docs/caveats.md) states all four in full, with the archive notes: release numbering that skips from 99 to 200, typos in the release headers, and E-utilities signalling a rate-limit rejection with HTTP 200.

## Development

```bash
uv run pytest             # 243 offline tests
uv run pytest -m network  # 3 more against NCBI; downloads ~1.2 GB
uv run ruff check . && uv run ruff format --check .
```

Network tests are deselected by default and excluded from CI: they are rate-limited per IP and move more than a gigabyte. Run them when changing anything about fetching, and before trusting a build.
