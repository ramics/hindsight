#!/usr/bin/env bash
# Install NCBI's sra-tools build into /opt/sratoolkit, at the version
# environment/Dockerfile pins.
#
#   scripts/fetch-sra-tools.sh
#   export PATH=/opt/sratoolkit/bin:$PATH
#
# This exists because sra-tools is the one tool this pipeline needs that is not in the
# Ubuntu archive, so a machine that is not built from environment/Dockerfile — a cloud
# session, a colleague's laptop — tends to have every other pinned tool and not this one.
#
# Its absence is worth a script rather than a line in a README because of how it fails.
# The sample fetcher downloads the .sra container, verifies its MD5, discovers it cannot
# run fasterq-dump, and falls back to ENA's FASTQ mirror. ENA then resets the connection
# partway through a 700 MB transfer, and the retry loop stalls at the same byte offset
# with `curl: (18) transfer closed`, `(35) Recv failure` and `(52) Empty reply`. That
# reads like a broken network or a proxy problem, and it is neither.
set -euo pipefail

VERSION="${SRATOOLKIT_VERSION:-3.2.1}"
PREFIX="${PREFIX:-/opt/sratoolkit}"
URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/sratoolkit.${VERSION}-ubuntu64.tar.gz"

if [[ -x "$PREFIX/bin/fasterq-dump" ]]; then
  echo "sra-tools already at $PREFIX ($("$PREFIX/bin/fasterq-dump" --version | tr -s ' '))" >&2
  exit 0
fi

tarball="$(mktemp -t sratoolkit.XXXXXX.tar.gz)"
trap 'rm -f "$tarball"' EXIT
echo "fetching sra-tools $VERSION" >&2
curl -sSL --retry 4 --retry-delay 5 "$URL" -o "$tarball"

mkdir -p "$PREFIX"
tar -xzf "$tarball" -C "$PREFIX" --strip-components=1

# Without a configuration the tools prompt interactively on first use, which in a
# non-interactive run looks like a hang rather than a question.
PATH="$PREFIX/bin:$PATH" vdb-config --restore-defaults >/dev/null 2>&1 || true

echo "installed $("$PREFIX/bin/fasterq-dump" --version | tr -s ' ')" >&2
echo "add it to PATH:  export PATH=$PREFIX/bin:\$PATH" >&2
