#!/usr/bin/env bash
# One command for Day 2: the index sample, the negatives, and the frozen dilution
# series. Safe to interrupt and re-run — every stage is resumable, and a stage that
# finished is skipped rather than redone.
#
#   scripts/day2.sh                      # everything, using config/day2.toml
#   scripts/day2.sh --only SRR10971381   # just the index sample
#   CONFIG=config/mine.toml scripts/day2.sh
#
# Expects the tools the config pins (fastp, bowtie2, samtools, curl) on PATH. See
# environment/ for a container that has them, and docs/day2-samples.md for what it
# costs in wall-clock and disk.
set -euo pipefail

CONFIG="${CONFIG:-config/day2.toml}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$here"

python -c 'import hindsight' 2>/dev/null || {
  echo "installing hindsight into the current environment" >&2
  pip install -e . >/dev/null
}

host_index="$(python - "$CONFIG" <<'PY'
import sys, tomllib
with open(sys.argv[1], "rb") as fh:
    print(tomllib.load(fh)["run"]["host_index"])
PY
)"

if [[ ! -f "${host_index}.1.bt2" ]]; then
  echo "host index ${host_index}.1.bt2 is missing; fetching GRCh38_noalt_as" >&2
  scripts/fetch-host-index.sh "$(dirname "$host_index")"
fi

exec hindsight samples run --config "$CONFIG" "$@"
