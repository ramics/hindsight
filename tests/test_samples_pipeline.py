"""Pipeline mechanics: stage resumability, source choice, and the series record.

The external tools are not run here. What is tested is the glue that decides what to
do, because that is where a mistake is invisible: a stage wrongly believed finished, a
source silently switched, a level written with the wrong background.
"""

from __future__ import annotations

import gzip
import json

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

    for level in record["levels"]:
        path = config.work_dir / "series" / f"SRRX_{level['name']}_1.fastq.gz"
        got = reads.read_ids(path)
        assert len(got) == level["total_fragments"]
        assert len(got & viral) == level["viral_fragments"]
        assert level["viral_fragments_per_million"] == pytest.approx(
            level["viral_fragments"] / level["total_fragments"] * 1e6
        )
        assert level["files"]["read_1"]["sha256"]


def test_build_series_levels_are_nested_on_disk(tmp_path):
    config = _series_config(tmp_path, levels=(6, 3, 1))
    sample = config.work_dir / "samples" / "SRRX"
    all_ids = [f"frag{i:03d}" for i in range(30)]
    viral = set(all_ids[:10])
    for mate in (1, 2):
        _fastq(sample / "qc" / f"SRRX_qc_{mate}.fastq.gz", all_ids)
        _fastq(sample / "nonhost" / f"SRRX_nonhost_{mate}.fastq.gz", all_ids)
    pipeline.build_series(config, index_accession="SRRX", viral_ids=viral, log=lambda msg: None)
    series = config.work_dir / "series"
    sets = [
        reads.read_ids(series / f"SRRX_{dilute.level_name(n)}_1.fastq.gz") & viral
        for n in (6, 3, 1)
    ]
    assert sets[2] <= sets[1] <= sets[0]


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
    assert first["levels"][0]["files"] == again["levels"][0]["files"]
    path = config.work_dir / "series" / f"SRRX_{dilute.level_name(4)}_1.fastq.gz"
    assert reads.count_reads(path) == 16


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
