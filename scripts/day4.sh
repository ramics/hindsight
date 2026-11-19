#!/usr/bin/env bash
# One command for Day 4: Arm B, MEGAHIT plus DIAMOND against the November 2019 proteins.
#
#   scripts/day4.sh                        # database, index sample, series, negatives
#   scripts/day4.sh --stage series         # one stage
#   scripts/day4.sh --level viral001000    # one dilution level
#   scripts/day4.sh --only SRR13677553     # one negative
#   CONFIG=config/mine.toml scripts/day4.sh
#
# Resumable, and it has to be: the series is six full assemblies and the negatives are
# twenty more, fetched, assembled and deleted one at a time, so an interrupted run
# continues from the sample or level it was on rather than starting again.
#
# Expects megahit, diamond and bowtie2 on PATH plus the tools config/day2.toml pins,
# since Day 4 re-runs the Day 2 pipeline per negative to get it to native depth. See
# environment/ for a container that has them all.
set -euo pipefail

CONFIG="${CONFIG:-config/day4.toml}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$here"

python -c 'import hindsight' 2>/dev/null || {
  echo "installing hindsight into the current environment" >&2
  pip install -e . >/dev/null
}

# sra-tools is the one pinned tool that is not in the Ubuntu archive, so a container not
# built from environment/Dockerfile tends to have everything except this. Checking here
# rather than letting the fetcher discover it matters because the failure is silent and
# slow: the fetcher downloads and MD5-verifies the .sra container, finds no
# fasterq-dump, falls back to ENA's FASTQ mirror, and then stalls in a retry loop as ENA
# resets a 700 MB transfer — which looks like a network fault and is not one.
if ! command -v fasterq-dump >/dev/null; then
  if [[ -x /opt/sratoolkit/bin/fasterq-dump ]]; then
    export PATH="/opt/sratoolkit/bin:$PATH"
  else
    echo "fasterq-dump is not on PATH. Without it the sample fetcher silently falls" >&2
    echo "back to ENA, which is unreliable for files this size. Install NCBI's build:" >&2
    echo "  scripts/fetch-sra-tools.sh && export PATH=/opt/sratoolkit/bin:\$PATH" >&2
    exit 1
  fi
fi

read -r source_db samples_config < <(python - "$CONFIG" <<'PY'
import sys, tomllib
with open(sys.argv[1], "rb") as fh:
    config = tomllib.load(fh)
print(config["database"]["source_db"], config["arm"]["samples_config"])
PY
)

# Arm B needs `viral.faa`, which Arm A's build does not produce: scripts/day3.sh builds
# with `--molecules nucleotide`. A directory holding only the nucleotide set is the most
# likely way into this script, so the protein set is checked for by name rather than by
# the directory existing.
if [[ ! -f "$source_db/viral.faa" || ! -f "$source_db/seqid2taxid.map" ]]; then
  echo "no dated protein set at $source_db/viral.faa; building it first" >&2
  cutoff="$(python - "$CONFIG" <<'PY'
import sys, tomllib
with open(sys.argv[1], "rb") as fh:
    print(tomllib.load(fh)["database"]["cutoff"])
PY
)"
  hindsight build --date "$cutoff" --out "$source_db" --molecules protein
fi

host_index="$(python - "$samples_config" <<'PY'
import sys, tomllib
with open(sys.argv[1], "rb") as fh:
    print(tomllib.load(fh)["run"]["host_index"])
PY
)"
if [[ ! -f "${host_index}.1.bt2" ]]; then
  echo "host index ${host_index}.1.bt2 is missing; fetching GRCh38_noalt_as" >&2
  scripts/fetch-host-index.sh "$(dirname "$host_index")"
fi

# The index sample and the frozen dilution series come from Day 2 and Arm B reads them
# as they are. The negatives are a different matter — Day 4 re-runs the Day 2 pipeline
# per negative, inside the arm, to reach native depth one sample at a time.
index_run="$(python - "$samples_config" <<'PY'
import sys, tomllib
with open(sys.argv[1], "rb") as fh:
    print(tomllib.load(fh)["samples"]["index_run"])
PY
)"
if [[ ! -f "work/series/parts/background_1.fastq.gz" ]]; then
  echo "no frozen dilution series in work/series; running the Day 2 pipeline first" >&2
  hindsight samples run --config "$samples_config" --only "$index_run"
fi

exec hindsight arm-b run --config "$CONFIG" "$@"
