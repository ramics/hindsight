"""Pipeline mechanics: stage resumability, source choice, and the series record.

The external tools are not run here. What is tested is the glue that decides what to
do, because that is where a mistake is invisible: a stage wrongly believed finished, a
source silently switched, a level written with the wrong background.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from hindsight.samples import dilute, ena, pipeline, reads, registry


def test_stage_marker_round_trip(tmp_path):
    stage = pipeline.Stage(tmp_path, "fetch")
    assert not stage.done
    stage.finish({"bytes": 12})
    assert stage.done
    assert pipeline.Stage(tmp_path, "fetch").load() == {"bytes": 12}
    assert not pipeline.Stage(tmp_path, "qc").done


def test_stage_markers_are_per_stage_and_per_sample(tmp_path):
    pipeline.Stage(tmp_path / "a", "host").finish()
    assert pipeline.Stage(tmp_path / "a", "host").done
    assert not pipeline.Stage(tmp_path / "b", "host").done


def _run_file(name, size=10, md5="x"):
    return ena.RunFile(url=f"https://example/{name}", size=size, md5=md5)


def test_fetch_prefers_the_mirror_when_the_conversion_fits(tmp_path, monkeypatch):
    calls = {}
    monkeypatch.setattr(
        pipeline.sra,
        "locate",
        lambda acc: pipeline.sra.RunObject(acc, "https://x", 1000, "md5", "sra-pub-run-odp"),
    )
    monkeypatch.setattr(pipeline.sra, "conversion_fits", lambda size, free: True)
    monkeypatch.setattr(
        pipeline.sra, "download", lambda obj, dest, log=None: _touch(dest / "x.sra")
    )

    def fake_whole(container, out_dir, *, accession, threads, tmp_dir, log_path=None):
        calls["whole"] = True
        return _touch(out_dir / f"{accession}_1.fastq.gz"), _touch(
            out_dir / f"{accession}_2.fastq.gz"
        )

    monkeypatch.setattr(pipeline.sra, "fasterq_dump_paired", fake_whole)
    monkeypatch.setattr(
        pipeline.sra,
        "fastq_dump_chunked",
        lambda *a, **k: pytest.fail("should not chunk when the whole run fits"),
    )
    r1, r2, source = pipeline._fetch_pair(
        "SRRX",
        (_run_file("SRRX_1.fastq.gz"), _run_file("SRRX_2.fastq.gz")),
        tmp_path / "raw",
        source="auto",
        threads=2,
        tmp_dir=tmp_path / "tmp",
        already_done=False,
        log=lambda msg: None,
    )
    assert source == "sra-odp"
    assert calls == {"whole": True}
    assert r1.name.endswith("_1.fastq.gz") and r2.name.endswith("_2.fastq.gz")


def test_fetch_chunks_when_the_whole_conversion_would_not_fit(tmp_path, monkeypatch):
    monkeypatch.setattr(
        pipeline.sra,
        "locate",
        lambda acc: pipeline.sra.RunObject(acc, "https://x", 10**10, "md5", "sra-pub-run-odp"),
    )
    monkeypatch.setattr(pipeline.sra, "conversion_fits", lambda size, free: False)
    monkeypatch.setattr(
        pipeline.sra, "download", lambda obj, dest, log=None: _touch(dest / "x.sra")
    )
    monkeypatch.setattr(
        pipeline.sra,
        "fasterq_dump_paired",
        lambda *a, **k: pytest.fail("should not convert whole when it does not fit"),
    )
    monkeypatch.setattr(
        pipeline.sra,
        "fastq_dump_chunked",
        lambda container, out_dir, *, accession, threads, log=None: (
            _touch(out_dir / f"{accession}_1.fastq.gz"),
            _touch(out_dir / f"{accession}_2.fastq.gz"),
        ),
    )
    _, _, source = pipeline._fetch_pair(
        "SRRX",
        (_run_file("SRRX_1.fastq.gz"), _run_file("SRRX_2.fastq.gz")),
        tmp_path / "raw",
        source="auto",
        threads=2,
        tmp_dir=tmp_path / "tmp",
        already_done=False,
        log=lambda msg: None,
    )
    assert source == "sra-odp-chunked"


def test_fetch_falls_back_to_ena_when_the_mirror_fails(tmp_path, monkeypatch):
    def explode(accession):
        raise RuntimeError("SDL unavailable")

    monkeypatch.setattr(pipeline.sra, "locate", explode)
    monkeypatch.setattr(
        pipeline.ena,
        "download",
        lambda file, dest_dir, verify=True, log=None: _touch(dest_dir / file.url.rsplit("/")[-1]),
    )
    messages = []
    _, _, source = pipeline._fetch_pair(
        "SRRX",
        (_run_file("SRRX_1.fastq.gz"), _run_file("SRRX_2.fastq.gz")),
        tmp_path / "raw",
        source="auto",
        threads=2,
        tmp_dir=tmp_path / "tmp",
        already_done=False,
        log=messages.append,
    )
    assert source == "ena"
    assert any("falling back to ENA" in m for m in messages)


def test_a_finished_fetch_is_not_repeated(tmp_path, monkeypatch):
    raw = tmp_path / "raw"
    raw.mkdir()
    _touch(raw / "SRRX_1.fastq.gz")
    _touch(raw / "SRRX_2.fastq.gz")
    monkeypatch.setattr(
        pipeline.sra, "locate", lambda acc: pytest.fail("no network on a resumed fetch")
    )
    r1, r2, source = pipeline._fetch_pair(
        "SRRX",
        (_run_file("SRRX_1.fastq.gz"), _run_file("SRRX_2.fastq.gz")),
        raw,
        source="auto",
        threads=2,
        tmp_dir=tmp_path / "tmp",
        already_done=True,
        log=lambda msg: None,
    )
    assert (r1.exists(), r2.exists(), source) == (True, True, "auto")


def _touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")
    return path


def _fastq(path, names):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for name in names:
            handle.write(f"@{name}\nACGT\n+\nIIII\n")
    return path


def _series_config(tmp_path, levels=(5, 2)):
    return registry.Config(
        index_run="SRRX",
        negatives_project="PRJNA1",
        negatives_count=0,
        negatives_seed=1,
        negatives_manifest=None,
        negatives_pinned=(),
        negatives_max_pairs=0,
        truth_reference="NC_045512.2",
        levels=tuple(levels),
        dilution_seed=99,
        threads=1,
        fetch_source="auto",
        prune_intermediates=False,
        host_index=tmp_path / "idx",
        work_dir=tmp_path / "work",
    )


def test_build_series_partitions_background_and_viral(tmp_path):
    config = _series_config(tmp_path)
    sample = config.work_dir / "samples" / "SRRX"
    all_ids = [f"frag{i:03d}" for i in range(40)]
    viral = set(all_ids[:12])
    for mate in (1, 2):
        _fastq(sample / "qc" / f"SRRX_qc_{mate}.fastq.gz", all_ids)
        _fastq(sample / "nonhost" / f"SRRX_nonhost_{mate}.fastq.gz", all_ids[:35])

    record = pipeline.build_series(
        config, index_accession="SRRX", viral_ids=viral, log=lambda msg: None
    )

    # 35 non-host fragments, 12 of which are viral -> 23 background.
    assert record["background_fragments"] == 23
    assert record["viral_fragments_available"] == 12
    assert record["nested"] is True
    assert [lvl["viral_fragments"] for lvl in record["levels"]] == [5, 2]

    # The background is stored once and shared; only the viral part is per level.
    parts = config.work_dir / "series" / "parts"
    assert reads.read_ids(parts / "background_1.fastq.gz") == set(all_ids[:35]) - viral
    for level in record["levels"]:
        names = [Path(p["path"]).name for p in level["parts"]["read_1"]]
        assert names == ["background_1.fastq.gz", f"{level['name']}_viral_1.fastq.gz"]
        assert all(p["sha256"] for p in level["parts"]["read_1"])
        got = reads.read_ids(parts / f"{level['name']}_viral_1.fastq.gz")
        assert got <= viral and len(got) == level["viral_fragments"]
        assert level["viral_fragments_per_million"] == pytest.approx(
            level["viral_fragments"] / level["total_fragments"] * 1e6
        )

    # A level is its parts concatenated, and that is a real FASTQ pair.
    out = config.work_dir / "materialised"
    r1, r2 = pipeline.materialise_level(record, record["levels"][0]["name"], out)
    ids = reads.read_ids(r1)
    assert len(ids) == record["levels"][0]["total_fragments"]
    assert len(ids & viral) == record["levels"][0]["viral_fragments"]
    assert reads.read_ids(r2) == ids


def test_the_headline_denominator_names_a_denominator_that_was_recorded(tmp_path):
    """The headline is post-QC reads, and the label has to point at a real column.

    A label naming a key that is not in per_million would make the report render a
    number under the wrong heading, which is worse than no label at all.
    """
    config = _series_config(tmp_path, levels=(3,))
    sample = config.work_dir / "samples" / "SRRX"
    all_ids = [f"frag{i:03d}" for i in range(20)]
    viral = set(all_ids[:6])
    for mate in (1, 2):
        _fastq(sample / "qc" / f"SRRX_qc_{mate}.fastq.gz", all_ids)
        _fastq(sample / "nonhost" / f"SRRX_nonhost_{mate}.fastq.gz", all_ids)

    record = pipeline.build_series(
        config,
        index_accession="SRRX",
        viral_ids=viral,
        denominators={"raw_spots": 100, "non_redacted_spots": 50, "post_qc_pairs": 20},
        log=lambda msg: None,
    )

    assert record["headline_denominator"] == "post_qc_pairs"
    for level in record["levels"]:
        assert record["headline_denominator"] in level["per_million"]
        assert level["per_million"]["post_qc_pairs"] == pytest.approx(
            level["viral_fragments"] / 20 * 1e6
        )


def test_materialise_refuses_a_part_that_changed_since_freezing(tmp_path):
    """The recorded digests are there to be checked, not just published."""
    config = _series_config(tmp_path, levels=(3,))
    sample = config.work_dir / "samples" / "SRRX"
    all_ids = [f"frag{i:03d}" for i in range(20)]
    viral = set(all_ids[:6])
    for mate in (1, 2):
        _fastq(sample / "qc" / f"SRRX_qc_{mate}.fastq.gz", all_ids)
        _fastq(sample / "nonhost" / f"SRRX_nonhost_{mate}.fastq.gz", all_ids)
    record = pipeline.build_series(
        config, index_accession="SRRX", viral_ids=viral, log=lambda msg: None
    )
    _fastq(config.work_dir / "series" / "parts" / "background_1.fastq.gz", all_ids[:4])
    with pytest.raises(RuntimeError, match="changed since the series was frozen"):
        pipeline.materialise_level(
            record, record["levels"][0]["name"], config.work_dir / "materialised"
        )


def test_build_series_levels_are_nested_on_disk(tmp_path):
    config = _series_config(tmp_path, levels=(6, 3, 1))
    sample = config.work_dir / "samples" / "SRRX"
    all_ids = [f"frag{i:03d}" for i in range(30)]
    viral = set(all_ids[:10])
    for mate in (1, 2):
        _fastq(sample / "qc" / f"SRRX_qc_{mate}.fastq.gz", all_ids)
        _fastq(sample / "nonhost" / f"SRRX_nonhost_{mate}.fastq.gz", all_ids)
    pipeline.build_series(config, index_accession="SRRX", viral_ids=viral, log=lambda msg: None)
    parts = config.work_dir / "series" / "parts"
    sets = [reads.read_ids(parts / f"{dilute.level_name(n)}_viral_1.fastq.gz") for n in (6, 3, 1)]
    assert sets[2] <= sets[1] <= sets[0]
    assert all(s <= viral for s in sets)


def test_build_series_is_resumable_and_does_not_double_write(tmp_path):
    config = _series_config(tmp_path, levels=(4,))
    sample = config.work_dir / "samples" / "SRRX"
    all_ids = [f"frag{i:03d}" for i in range(20)]
    viral = set(all_ids[:8])
    for mate in (1, 2):
        _fastq(sample / "qc" / f"SRRX_qc_{mate}.fastq.gz", all_ids)
        _fastq(sample / "nonhost" / f"SRRX_nonhost_{mate}.fastq.gz", all_ids)
    first = pipeline.build_series(
        config, index_accession="SRRX", viral_ids=viral, log=lambda msg: None
    )
    again = pipeline.build_series(
        config, index_accession="SRRX", viral_ids=viral, log=lambda msg: None
    )
    assert first["levels"][0]["parts"] == again["levels"][0]["parts"]
    # 20 fragments, 8 viral -> 12 background, plus a 4-fragment viral part.
    parts = config.work_dir / "series" / "parts"
    assert reads.count_reads(parts / "background_1.fastq.gz") == 12
    assert reads.count_reads(parts / f"{dilute.level_name(4)}_viral_1.fastq.gz") == 4


def test_build_series_refuses_when_a_viral_fragment_is_missing_from_the_qc_fastq(tmp_path):
    """A level that cannot be filled must stop the build, not ship thin.

    This is the failure mode where the truth label and the FASTQ have drifted apart —
    one regenerated without the other — and the resulting series would claim 1,000
    viral fragments while holding 900.
    """
    config = _series_config(tmp_path, levels=(5,))
    sample = config.work_dir / "samples" / "SRRX"
    present = [f"frag{i:03d}" for i in range(10)]
    for mate in (1, 2):
        _fastq(sample / "qc" / f"SRRX_qc_{mate}.fastq.gz", present)
        _fastq(sample / "nonhost" / f"SRRX_nonhost_{mate}.fastq.gz", present)
    viral = set(present[:3]) | {"ghost1", "ghost2", "ghost3"}
    with pytest.raises(RuntimeError, match="wanted 5 fragments"):
        pipeline.build_series(config, index_accession="SRRX", viral_ids=viral, log=lambda msg: None)


def test_machine_record_describes_the_host(tmp_path):
    record = pipeline.machine_record(tmp_path)
    assert record["cpu_count"] and record["cpu_count"] > 0
    assert "platform" in record
    assert "work_dir_disk" in record


def test_fetch_reference_refuses_a_rate_limit_json_body(tmp_path, monkeypatch):
    class FakeResponse:
        text = '{"error":"API rate limit exceeded"}'
        ok = True

        def raise_for_status(self):
            return None

    class FakeSession:
        def get(self, url, params=None, timeout=None):
            return FakeResponse()

    with pytest.raises(RuntimeError, match="did not return FASTA"):
        pipeline.fetch_reference("NC_045512.2", tmp_path / "ref.fna", session=FakeSession())
    assert not (tmp_path / "ref.fna").exists()


def test_fetch_reference_caches(tmp_path):
    dest = tmp_path / "ref.fna"
    dest.write_text(">NC_045512.2\nACGT\n")

    class Boom:
        def get(self, *a, **k):
            raise AssertionError("should not refetch an existing reference")

    assert pipeline.fetch_reference("NC_045512.2", dest, session=Boom()) == dest


def test_sample_record_is_json_serialisable():
    result = pipeline.SampleResult(accession="SRRX", role="index")
    result.timings.append(pipeline.Timing("fetch", 1.5, "sra-odp"))
    json.dumps(pipeline._sample_record(result))


def test_provenance_merges_rather_than_replaces(tmp_path):
    """A partial run must not erase what an earlier one recorded.

    This is a real loss, not a hypothetical: a three-negative run overwrote the record
    describing the index sample and the frozen dilution series, and the numbers had to
    be regenerated.
    """
    path = tmp_path / "day2-provenance.json"
    first = {
        "generated_at": "2026-11-17T17:00:00",
        "command": ["hindsight", "samples", "run", "--only", "SRR_INDEX"],
        "samples": [{"accession": "SRR_INDEX", "role": "index", "viral_pairs": 62117}],
        "dilution_series": {"levels": [{"name": "viral000010"}]},
    }
    path.write_text(json.dumps(pipeline.merge_provenance(path, first)))

    second = {
        "generated_at": "2026-11-17T18:00:00",
        "command": ["hindsight", "samples", "run", "--skip", "SRR_INDEX"],
        "samples": [{"accession": "SRR_NEG1", "role": "negative"}],
        "dilution_series": None,
    }
    merged = pipeline.merge_provenance(path, second)

    assert [s["accession"] for s in merged["samples"]] == ["SRR_INDEX", "SRR_NEG1"]
    assert merged["samples"][0]["viral_pairs"] == 62117
    assert merged["dilution_series"]["levels"] == [{"name": "viral000010"}]
    assert [r["samples"] for r in merged["runs"]] == [["SRR_INDEX"], ["SRR_NEG1"]]
    assert [r["wrote_series"] for r in merged["runs"]] == [True, False]


def test_a_rerun_of_a_sample_replaces_its_own_record(tmp_path):
    path = tmp_path / "day2-provenance.json"
    old = {"samples": [{"accession": "SRR1", "role": "negative", "qc_pairs": 1}]}
    path.write_text(json.dumps(pipeline.merge_provenance(path, old)))
    merged = pipeline.merge_provenance(
        path, {"samples": [{"accession": "SRR1", "role": "negative", "qc_pairs": 2}]}
    )
    assert len(merged["samples"]) == 1
    assert merged["samples"][0]["qc_pairs"] == 2


def _pruned_sample(tmp_path, accession="SRRN", *, capped=True):
    """A sample as it looks after a finished, pruned run: markers and nonhost only."""
    config = _series_config(tmp_path)
    sample = config.work_dir / "samples" / accession
    for mate in (1, 2):
        _fastq(sample / "nonhost" / f"{accession}_nonhost_{mate}.fastq.gz", ["a", "b", "c"])
    (sample / "qc").mkdir(parents=True, exist_ok=True)
    (sample / "qc" / "fastp.json").write_text(
        json.dumps(
            {
                "summary": {
                    "before_filtering": {"total_reads": 2000},
                    "after_filtering": {"total_reads": 1800},
                }
            }
        )
    )
    pipeline.Stage(sample, "redaction").finish({"redacted_pairs": 40})
    if capped:
        pipeline.Stage(sample, "cap").finish({"capped_pairs": 700, "fraction": 0.7})
    pipeline.Stage(sample, "counts").finish({"nonhost_pairs": 650})
    pipeline.Stage(sample, "pruned").finish(
        {
            "freed_bytes": 1234,
            "source": "sra-odp",
            "files": {"raw_1": {"path": "gone", "sha256": "abc", "bytes": 10}},
        }
    )
    return config, sample


class _FakeRun:
    study = "PRJNA1"
    sample = "SAMN1"
    library_strategy = "RNA-Seq"
    library_layout = "PAIRED"
    instrument_model = "Illumina"
    read_count = 1000
    base_count = 100000
    collection_date = "missing"
    layout = "paired"
    files = ()


def test_a_pruned_sample_is_rebuilt_from_its_markers(tmp_path):
    """Re-running a finished sample must not re-fetch what pruning deleted on purpose.

    The fetch stage looks at files rather than markers, so a pruned negative used to
    re-download and re-convert a 1.5 GB container to rebuild files nothing reads and
    the next prune deletes again.
    """
    config, _ = _pruned_sample(tmp_path)
    spec = registry.SampleSpec(accession="SRRN", role="negative", notes="")
    result = pipeline._rebuild_pruned_record(spec, config, run=_FakeRun(), log=lambda m: None)

    assert result is not None
    assert (result.raw_pairs, result.qc_pairs) == (1000, 900)
    assert result.redacted_pairs == 40 and result.non_redacted_pairs == 960
    assert result.capped_pairs == 700 and result.cap_fraction == 0.7
    assert result.nonhost_pairs == 650
    assert result.host_pairs == 50, "host count is against the capped input, not the QC one"
    assert result.pruned_bytes == 1234
    assert result.source == "sra-odp"
    # The deleted files keep their recorded digests; the surviving ones are measured.
    assert result.files["raw_1"]["sha256"] == "abc"
    assert result.files["nonhost_1"]["sha256"] != "abc"


def test_rebuilding_is_refused_when_the_evidence_is_incomplete(tmp_path):
    config, sample = _pruned_sample(tmp_path, accession="SRRM")
    spec = registry.SampleSpec(accession="SRRM", role="negative", notes="")
    (sample / "qc" / "fastp.json").unlink()
    assert pipeline._rebuild_pruned_record(spec, config, run=_FakeRun(), log=lambda m: None) is None

    config2, sample2 = _pruned_sample(tmp_path / "other", accession="SRRO")
    spec2 = registry.SampleSpec(accession="SRRO", role="negative", notes="")
    (sample2 / "nonhost" / "SRRO_nonhost_2.fastq.gz").unlink()
    assert (
        pipeline._rebuild_pruned_record(spec2, config2, run=_FakeRun(), log=lambda m: None) is None
    )
