#!/usr/bin/env bash
# One command for Day 5: score both arms and write the benchmark's headline numbers.
#
#   scripts/day5.sh                             # the recorded decisions
#   scripts/day5.sh --clade betacoronavirus     # disagree with one of them
#   scripts/day5.sh --framing all --budget 1    # or with another
#
# Unlike every other day in this benchmark, this one needs no tools, no database, no
# reads and no network: it reads results/day3-provenance.json and
# results/day4-provenance.json and does arithmetic. It takes about a second, and that is
# a deliberate property of how Days 3 and 4 recorded their output rather than a sign the
# day is trivial — a scoring step a reader can re-run with a different threshold is the
# difference between a published limit of detection and a published opinion.
set -euo pipefail

CONFIG="${CONFIG:-config/day5.toml}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$here"

python -c 'import hindsight' 2>/dev/null || {
  echo "installing hindsight into the current environment" >&2
  pip install -e . >/dev/null
}

for record in results/day3-provenance.json results/day4-provenance.json; do
  if [[ ! -f "$record" ]]; then
    echo "$record is missing; run scripts/day3.sh and scripts/day4.sh first" >&2
    exit 1
  fi
done

exec hindsight score run --config "$CONFIG" "$@"
