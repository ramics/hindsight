#!/usr/bin/env bash
# Fetch the prebuilt bowtie2 index for the GRCh38 no-alt analysis set (~4.1 GB).
#
# Prebuilt rather than built here on purpose: bowtie2-build on the human genome costs
# a couple of hours of CPU for a byte-identical result, and the authors' index is the
# one most published host-filtering steps actually used. The .md5 file NCBI's mirror
# ships beside it is checked after the download.
#
#   scripts/fetch-host-index.sh work/ref/bt2
set -euo pipefail

dest="${1:-work/ref/bt2}"
base="https://genome-idx.s3.amazonaws.com/bt"
name="GRCh38_noalt_as"
mkdir -p "$dest"

files=(1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2 md5 build.txt)
for suffix in "${files[@]}"; do
  f="${name}.${suffix}"
  want="$(curl -sS -I --max-time 60 "$base/$f" | awk '/[Cc]ontent-[Ll]ength/{gsub("\r","");print $2}' | tail -1)"
  for _ in $(seq 1 80); do
    have="$(stat -c%s "$dest/$f" 2>/dev/null || stat -f%z "$dest/$f" 2>/dev/null || echo 0)"
    [[ -n "$want" && "$have" == "$want" ]] && break
    # Resume rather than restart: a 1 GB transfer through a proxy or a hotel wifi gets
    # cut, and a truncated index is not an obvious error until bowtie2 reads it.
    curl -sS -C - -o "$dest/$f" --max-time 900 "$base/$f" || true
  done
  echo "$f $(stat -c%s "$dest/$f" 2>/dev/null || stat -f%z "$dest/$f") (expected ${want:-unknown})"
done

if command -v md5sum >/dev/null && [[ -f "$dest/$name.md5" ]]; then
  # The published .md5 lists the files by the absolute path they had on the machine
  # that built them (/Users/langmead/...), so `md5sum -c` matches nothing and says
  # "no file was verified" — which reads like a pass if you are not looking. Strip
  # the directory, and fail loudly if nothing was checked.
  (
    cd "$dest"
    sed 's#  .*/#  #' "$name.md5" > "$name.md5.local"
    md5sum -c --ignore-missing "$name.md5.local" | tee "$name.md5.checked"
    grep -c ': OK$' "$name.md5.checked" > /dev/null || {
      echo "md5 check verified nothing in $dest" >&2
      exit 1
    }
  )
fi
