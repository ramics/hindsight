#!/usr/bin/env bash
# One command for Day 3: Arm A, Kraken2 against the November 2019 database.
#
#   scripts/day3.sh                        # database, index sample, series, negatives
#   scripts/day3.sh --stage series         # one stage
#   scripts/day3.sh --only SRR13677553     # one negative
#   CONFIG=config/mine.toml scripts/day3.sh
#
# Resumable, and it has to be: the negatives are fetched, filtered, classified and
# deleted one at a time, so an interrupted cohort continues from the sample it was on
# rather than starting the twenty again.
#
# Expects kraken2 on PATH plus the tools config/day2.toml pins, since Day 3 re-runs
# the Day 2 pipeline per negative to get it to native depth. See environment/ for a
# container that has them all.
set -euo pipefail

CONFIG="${CONFIG:-config/day3.toml}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$here"

python -c 'import hindsight' 2>/dev/null || {
  echo "installing hindsight into the current environment" >&2
  pip install -e . >/dev/null
}

read -r source_db samples_config < <(python - "$CONFIG" <<'PY'
import sys, tomllib
with open(sys.argv[1], "rb") as fh:
    config = tomllib.load(fh)
print(config["database"]["source_db"], config["arm"]["samples_config"])
PY
)

if [[ ! -f "$source_db/viral.fna" || ! -f "$source_db/seqid2taxid.map" ]]; then
  echo "no dated database at $source_db; building it first" >&2
  cutoff="$(python - "$CONFIG" <<'PY'
import sys, tomllib
with open(sys.argv[1], "rb") as fh:
    print(tomllib.load(fh)["database"]["cutoff"])
PY
)"
  hindsight build --date "$cutoff" --out "$source_db" --molecules nucleotide
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

exec hindsight arm-a run --config "$CONFIG" "$@"
