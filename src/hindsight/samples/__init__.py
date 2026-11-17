"""Sample acquisition, read QC, host filtering and the frozen dilution series.

Day 1 of the benchmark built the database side: a viral reference set as it stood
before December 2019. This package builds the sample side, and nothing else — it does
no classification and consults no reference database other than the modern SARS-CoV-2
genome, which it uses *only* to label ground truth. The detection arms never see that
label.

The products are:

* the index sample, SRR10971381 (Wu et al. BALF), QC'd and host-filtered;
* a handful of time-matched 2016-17 Wuhan BALF negatives, same treatment;
* a dilution series: the index sample's own SARS-CoV-2 reads thinned to fixed counts
  and held in the index sample's own background, frozen with a recorded seed.
"""

from __future__ import annotations

__all__ = ["dilute", "ena", "pipeline", "reads", "registry", "tools"]
