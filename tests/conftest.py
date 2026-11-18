"""Fixtures shared across the test modules.

The fake dumper lives here rather than in one test module because two modules need
it, and importing a fixture from another test module needs ``tests`` to be a package
— which it is not, and making it one to share one fixture is the wrong trade.
"""

from __future__ import annotations

import gzip
import subprocess

import pytest

from hindsight.samples import sra


@pytest.fixture
def fake_dumper(monkeypatch, tmp_path):
    """Replace fastq-dump and pigz with a dumper that serves a fixed pool of spots."""
    state = {"total": 7, "calls": []}
    real_run = subprocess.run

    def fake_run(argv, **kwargs):
        if argv[0] == "fastq-dump":
            start = int(argv[argv.index("-N") + 1])
            stop = int(argv[argv.index("-X") + 1])
            out = argv[argv.index("-O") + 1]
            sra_path = argv[-1]
            stem = sra_path.rsplit("/", 1)[-1].rsplit(".", 1)[0]
            spots = [s for s in range(start, stop + 1) if s <= state["total"]]
            state["calls"].append((start, stop, len(spots)))
            for mate in (1, 2):
                target = f"{out}/{stem}_{mate}.fastq"
                with open(target, "w") as handle:
                    for spot in spots:
                        handle.write(f"@spot{spot}/{mate}\nACGT\n+\nIIII\n")
            return subprocess.CompletedProcess(argv, 0)
        if argv[0] == "pigz":
            source = argv[-1]
            with open(source, "rb") as handle:
                data = handle.read()
            kwargs["stdout"].write(gzip.compress(data))
            return subprocess.CompletedProcess(argv, 0)
        return real_run(argv, **kwargs)

    monkeypatch.setattr(sra.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(sra.subprocess, "run", fake_run)
    return state
