# The benchmark's environment, pinned.
#
# Day 7 is a reproducibility pass, and "pinned versions" has to mean a thing a reader
# can run, not a list in a config file. Ubuntu 24.04 is used because its archive
# happens to carry exactly the versions config/day2.toml pins — fastp 0.23.4,
# bowtie2 2.5.2, samtools 1.19.2 — so the pins are satisfiable from a distribution
# rather than from a bespoke build.
#
#   docker build -t hindsight environment/
#   docker run --rm -v "$PWD:/work" -w /work hindsight scripts/day2.sh
#
# The container carries no data. The human index (4.1 GB) and the samples are fetched
# into the mounted working directory, so they survive a rebuild and can be shared
# between runs.
#
# Day 7 should replace this tag with the digest it actually built from
# (`docker inspect --format '{{index .RepoDigests 0}}' ubuntu:24.04`); a tag is a
# moving target and this file does not claim otherwise.
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONUNBUFFERED=1 \
    PATH=/opt/sratoolkit/bin:$PATH

# The tools whose version changes what the benchmark measures are pinned exactly to
# the versions config/day2.toml names. An apt archive that moves one of these should
# force a deliberate edit here, visible in a diff.
#
# Python, curl and ca-certificates are deliberately *not* pinned: they take security
# updates, they do not touch a read, and pinning them only guarantees the build breaks
# later for no scientific gain.
RUN apt-get update && apt-get install -y --no-install-recommends \
      fastp=0.23.4+dfsg-1 \
      bowtie2=2.5.2-1 \
      samtools=1.19.2-1build2 \
      kraken2=2.1.3-1 \
      megahit=1.2.9-5 \
      diamond-aligner=2.1.9-1 \
      pigz=2.8-1 \
      python3 \
      python3-pip \
      python3-venv \
      curl \
      ca-certificates \
      gawk \
    && rm -rf /var/lib/apt/lists/* \
    && fastp --version && bowtie2 --version | head -1 && samtools --version | head -1 \
    && kraken2 --version | head -1 \
    && megahit --version && diamond --version

# Kraken2 keeps build_db, k2mask and friends in a libexec directory that is not on
# PATH; kraken2-build puts its own directory there before calling them, and Day 3's
# arm looks for them in the same places. Nothing to install, but worth knowing when a
# masking step cannot find k2mask on a different distribution.

# sra-tools is not in the Ubuntu archive; NCBI's own build is used, by version.
ARG SRATOOLKIT_VERSION=3.2.1
RUN curl -fsSL \
      "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRATOOLKIT_VERSION}/sratoolkit.${SRATOOLKIT_VERSION}-ubuntu64.tar.gz" \
      -o /tmp/sratoolkit.tar.gz \
    && mkdir -p /opt/sratoolkit \
    && tar -xzf /tmp/sratoolkit.tar.gz -C /opt/sratoolkit --strip-components=1 \
    && rm /tmp/sratoolkit.tar.gz \
    && fastq-dump --version && fasterq-dump --version

WORKDIR /src
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --break-system-packages --no-cache-dir .

WORKDIR /work
ENTRYPOINT ["/bin/bash", "-lc"]
CMD ["hindsight samples plan"]

# The days this container runs:
#
#   docker run --rm -v "$PWD:/work" -w /work hindsight scripts/day2.sh   # the reads
#   docker run --rm -v "$PWD:/work" -w /work hindsight scripts/day3.sh   # Arm A
#   docker run --rm -v "$PWD:/work" -w /work hindsight scripts/day4.sh   # Arm B
