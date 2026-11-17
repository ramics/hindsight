"""The sample cohort: how the negatives are chosen, and what the config promises."""

from __future__ import annotations

import textwrap

import pytest

from hindsight.samples import registry

MANIFEST_HEADER = (
    "run_accession\tsample_accession\tread_count\tfastq_ftp\tfastq_bytes\t"
    "collection_date\tbiological_replicate\thealth_state\n"
)


def manifest_row(accession, *, replicate="Patient_1", paired=True, orphan=False):
    names = []
    if orphan:
        names.append(f"{accession}.fastq.gz")
    if paired:
        names += [f"{accession}_1.fastq.gz", f"{accession}_2.fastq.gz"]
    if not paired and not orphan:
        names.append(f"{accession}.fastq.gz")
    urls = ";".join(f"ftp.sra.ebi.ac.uk/vol1/fastq/x/{accession}/{n}" for n in names)
    sizes = ";".join("100" for _ in names)
    state = "Disease" if replicate.startswith("Patient") else "Health"
    return f"{accession}\tSAMN1\t1000\t{urls}\t{sizes}\tmissing\t{replicate}\t{state}\n"


def write_manifest(tmp_path, rows):
    path = tmp_path / "negatives.tsv"
    path.write_text(MANIFEST_HEADER + "".join(rows))
    return path


def test_patients_are_separated_from_volunteers_and_blanks(tmp_path):
    path = write_manifest(
        tmp_path,
        [
            manifest_row("SRR1", replicate="Patient_1"),
            manifest_row("SRR2", replicate="Health_control_1"),
            manifest_row("SRR3", replicate="Blank_1"),
            manifest_row("SRR4", replicate="Patient_2"),
        ],
    )
    rows = registry.read_negatives_manifest(path)
    assert len(rows) == 4
    assert [r["run_accession"] for r in registry.patient_rows(rows)] == ["SRR1", "SRR4"]


def test_only_runs_deposited_as_pairs_are_eligible(tmp_path):
    path = write_manifest(
        tmp_path,
        [
            manifest_row("SRR1"),
            manifest_row("SRR2", paired=False),
            manifest_row("SRR3", orphan=True),
        ],
    )
    rows = registry.read_negatives_manifest(path)
    eligible = {r["run_accession"] for r in registry.paired_rows(rows)}
    assert eligible == {"SRR1", "SRR3"}


def test_selection_is_seeded_and_reproducible(tmp_path):
    path = write_manifest(tmp_path, [manifest_row(f"SRR{i:04d}") for i in range(100)])
    rows = registry.read_negatives_manifest(path)
    first = registry.choose_negatives(rows, count=20, seed=5)
    again = registry.choose_negatives(rows, count=20, seed=5)
    other = registry.choose_negatives(rows, count=20, seed=6)
    assert first == again
    assert first != other
    assert len(set(first)) == 20
    assert first == sorted(first), "the chosen list is returned sorted for a readable diff"


def test_selection_refuses_to_pad_a_short_pool(tmp_path):
    path = write_manifest(tmp_path, [manifest_row(f"SRR{i}") for i in range(5)])
    rows = registry.read_negatives_manifest(path)
    with pytest.raises(ValueError, match="only 5 eligible"):
        registry.choose_negatives(rows, count=20, seed=1)


def _config_text(extra: str = "") -> str:
    return textwrap.dedent(
        f"""
        [samples]
        index_run = "SRR10971381"
        truth_reference = "NC_045512.2"

        [negatives]
        project = "PRJNA699976"
        count = 2
        seed = 3
        {extra}

        [dilution]
        levels = [100, 10]
        seed = 11

        [run]
        threads = 2
        host_index = "work/ref/bt2/GRCh38_noalt_as"
        work_dir = "work"

        [tools]
        fastp = "0.23.4"
        """
    )


def test_config_round_trips_the_things_the_methods_section_needs(tmp_path):
    path = tmp_path / "day2.toml"
    path.write_text(_config_text('pinned_runs = ["SRRA", "SRRB"]'))
    config = registry.load_config(path)
    assert config.index_run == "SRR10971381"
    assert config.levels == (100, 10)
    assert config.dilution_seed == 11
    assert config.negatives_seed == 3
    assert config.tool_pins == {"fastp": "0.23.4"}
    specs = registry.samples_from_config(config)
    assert [s.accession for s in specs] == ["SRR10971381", "SRRA", "SRRB"]
    assert specs[0].is_index and not specs[1].is_index
    assert "pinned" in specs[1].notes


def test_pinned_runs_win_over_the_manifest(tmp_path):
    manifest = write_manifest(tmp_path, [manifest_row(f"SRR{i}") for i in range(30)])
    path = tmp_path / "day2.toml"
    path.write_text(
        _config_text(f'pinned_runs = ["SRRX", "SRRY"]\nmanifest = "{manifest.name}"')
    )
    config = registry.load_config(path)
    assert [s.accession for s in registry.samples_from_config(config)[1:]] == ["SRRX", "SRRY"]


def test_without_pins_the_cohort_comes_from_the_manifest_and_seed(tmp_path):
    manifest = write_manifest(tmp_path, [manifest_row(f"SRR{i:03d}") for i in range(30)])
    path = tmp_path / "day2.toml"
    path.write_text(_config_text(f'manifest = "{manifest.name}"'))
    config = registry.load_config(path)
    specs = registry.samples_from_config(config)
    assert len(specs) == 3
    assert "seed 3" in specs[1].notes


def test_a_cohort_with_neither_pins_nor_manifest_is_an_error(tmp_path):
    path = tmp_path / "day2.toml"
    path.write_text(_config_text())
    config = registry.load_config(path)
    with pytest.raises(ValueError, match="no negatives"):
        registry.samples_from_config(config)


def test_the_committed_config_matches_the_committed_manifest():
    """The 20 pinned runs must be exactly what the recorded seed selects.

    This is the test that keeps the methods section true: if someone edits the pinned
    list or the seed without the other, the cohort is no longer the one the paper says
    was chosen blind.
    """
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    config = registry.load_config(root / "config" / "day2.toml")
    rows = registry.read_negatives_manifest(config.negatives_manifest)
    derived = registry.choose_negatives(
        rows, count=config.negatives_count, seed=config.negatives_seed
    )
    assert sorted(config.negatives_pinned) == derived
    assert len(registry.patient_rows(rows)) == 408
