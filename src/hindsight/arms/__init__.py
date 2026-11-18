"""The detection arms: the pipelines the benchmark actually scores.

Each arm sees the reads Day 2 froze and a database built by ``hindsight build`` for a
cutoff date, and nothing else. In particular no arm sees the ground-truth viral
labels or the modern SARS-CoV-2 reference used to compute them: knowing the answer is
the benchmark's job, not the arm's.

* ``kraken`` — Arm A, exact k-mer classification against the time-correct nucleotide
  database (Day 3).
"""
