"""The hindsight command line."""

from __future__ import annotations

import datetime as dt
import json
import os
import sys
from pathlib import Path

import click

from . import __version__
from . import provenance as provenance_mod
from . import releases as releases_mod
from . import taxonomy as taxonomy_mod
from .accessions import Molecule
from .build import build as run_build
from .ncbi import NcbiClient

DEFAULT_CACHE = Path(os.environ.get("HINDSIGHT_CACHE", Path.home() / ".cache" / "hindsight"))

MOLECULE_CHOICES = {
    "nucleotide": (Molecule.NUCLEOTIDE,),
    "protein": (Molecule.PROTEIN,),
    "both": (Molecule.NUCLEOTIDE, Molecule.PROTEIN),
}


def _client(ctx: click.Context) -> NcbiClient:
    return NcbiClient(
        api_key=ctx.obj.get("api_key"),
        rate_per_second=ctx.obj.get("rate"),
    )


def _releases(ctx: click.Context) -> list[releases_mod.Release]:
    return releases_mod.load_releases(ctx.obj.get("release_table"))


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="hindsight")
@click.option(
    "--cache",
    type=click.Path(file_okay=False, path_type=Path),
    default=DEFAULT_CACHE,
    show_default=True,
    help="Directory for downloaded catalogs, taxdumps and fetched record shards.",
)
@click.option(
    "--release-table",
    type=click.Path(dir_okay=False, exists=True, path_type=Path),
    default=None,
    help="Override the vendored release-date table.",
)
@click.option(
    "--api-key",
    envvar="NCBI_API_KEY",
    default=None,
    help="NCBI E-utilities API key; raises the request ceiling. [env: NCBI_API_KEY]",
)
@click.option(
    "--rate",
    type=float,
    default=None,
    help="Override E-utilities requests per second (default: 2 without a key, 8 with).",
)
@click.option("--quiet", is_flag=True, help="Suppress progress commentary on stderr.")
@click.pass_context
def main(
    ctx: click.Context,
    cache: Path,
    release_table: Path | None,
    api_key: str | None,
    rate: float | None,
    quiet: bool,
) -> None:
    """Build a time-correct viral reference database for a given cutoff date.

    hindsight reconstructs the viral portion of RefSeq as it stood before a date, from
    NCBI's archived release catalogs, fetching every record at the exact version that
    catalog listed. Nothing published after the cutoff can enter the database.
    """
    ctx.ensure_object(dict)
    ctx.obj.update(
        cache=cache, release_table=release_table, api_key=api_key, rate=rate, quiet=quiet
    )


@main.command("releases")
@click.option("--limit", type=int, default=0, help="Show only the most recent N releases.")
@click.option("--json", "as_json", is_flag=True, help="Emit JSON instead of a table.")
@click.pass_context
def releases_cmd(ctx: click.Context, limit: int, as_json: bool) -> None:
    """List the archived RefSeq releases and their release dates."""
    table = sorted(_releases(ctx), key=lambda r: r.number)
    if limit:
        table = table[-limit:]
    if as_json:
        click.echo(
            json.dumps([{"release": r.number, "date": r.date.isoformat()} for r in table], indent=2)
        )
        return
    click.echo(f"{'release':>8}  date")
    for release in table:
        click.echo(f"{release.number:>8}  {release.date.isoformat()}")


@main.command("resolve")
@click.argument("cutoff", type=click.DateTime(formats=["%Y-%m-%d"]))
@click.pass_context
def resolve_cmd(ctx: click.Context, cutoff: dt.datetime) -> None:
    """Show which release a CUTOFF date (YYYY-MM-DD) selects, and what it excludes."""
    table = _releases(ctx)
    cutoff_date = cutoff.date()
    release = releases_mod.resolve(cutoff_date, table)
    following = releases_mod.next_release_after(release, table)
    click.echo(f"cutoff:            {cutoff_date.isoformat()}")
    click.echo(f"selected release:  {release.number} ({release.date.isoformat()})")
    click.echo(f"catalog:           {release.catalog_url()}")
    if following:
        click.echo(f"first excluded:    release {following.number} ({following.date.isoformat()})")
        click.echo(
            f"blind window:      {(following.date - release.date).days} days between the "
            "selected release and the next one"
        )
    else:
        click.echo("first excluded:    none; this is the most recent archived release")


@main.command("refresh-releases")
@click.option(
    "--out",
    type=click.Path(dir_okay=False, path_type=Path),
    default=Path("src/hindsight/data/release_dates.tsv"),
    show_default=True,
    help="Where to write the regenerated table.",
)
@click.option("--progress/--no-progress", default=True)
@click.pass_context
def refresh_cmd(ctx: click.Context, out: Path, progress: bool) -> None:
    """Rebuild the vendored release-date table from NCBI's release notes.

    Only needed when NCBI issues a new release. The table is checked into the repo so
    that resolving a date needs no network and the dates are reviewable in a diff.
    """
    client = _client(ctx)
    table = releases_mod.refresh_releases(client, progress=progress)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(releases_mod.format_table(table))
    click.echo(f"wrote {len(table)} releases to {out}", err=True)


@main.command("taxdumps")
@click.option(
    "--flavour",
    type=click.Choice(taxonomy_mod.FLAVOURS),
    default="taxdmp",
    show_default=True,
)
@click.option("--limit", type=int, default=12, help="Show only the most recent N snapshots.")
@click.pass_context
def taxdumps_cmd(ctx: click.Context, flavour: str, limit: int) -> None:
    """List the dated taxonomy snapshots NCBI has archived."""
    snapshots = taxonomy_mod.list_snapshots(_client(ctx), flavour)
    click.echo(f"{len(snapshots)} {flavour} snapshots, {snapshots[0].date} .. {snapshots[-1].date}")
    for snapshot in snapshots[-limit:] if limit else snapshots:
        click.echo(f"  {snapshot.date.isoformat()}  {snapshot.url}")


@main.command("build")
@click.option(
    "--date",
    "cutoff",
    type=click.DateTime(formats=["%Y-%m-%d"]),
    required=True,
    help="Cutoff date. Nothing published on or after this date may enter the database.",
)
@click.option(
    "--out",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="Output directory for the database.",
)
@click.option(
    "--molecules",
    type=click.Choice(sorted(MOLECULE_CHOICES)),
    default="both",
    show_default=True,
    help="Which record sets to fetch. 'protein' is what the DIAMOND arm needs.",
)
@click.option(
    "--taxonomy/--no-taxonomy",
    default=True,
    show_default=True,
    help="Include the contemporaneous NCBI taxonomy dump and a seqid2taxid map.",
)
@click.option(
    "--taxdump-flavour",
    type=click.Choice(taxonomy_mod.FLAVOURS),
    default="taxdmp",
    show_default=True,
)
@click.option(
    "--manifest-only",
    is_flag=True,
    help="Resolve the accession set and stop, without fetching any sequence.",
)
@click.option(
    "--stream-catalog",
    is_flag=True,
    help="Decompress the catalog off the wire instead of caching it (saves ~1.3 GB of "
    "disk, loses resumability).",
)
@click.option("--batch-size", type=int, default=None, help="Accessions per efetch request.")
@click.option(
    "--limit",
    type=int,
    default=None,
    help="Fetch only the first N accessions per molecule. For smoke tests, not science.",
)
@click.option("--progress/--no-progress", default=True)
@click.pass_context
def build_cmd(
    ctx: click.Context,
    cutoff: dt.datetime,
    out: Path,
    molecules: str,
    taxonomy: bool,
    taxdump_flavour: str,
    manifest_only: bool,
    stream_catalog: bool,
    batch_size: int | None,
    limit: int | None,
    progress: bool,
) -> None:
    """Build the database for a cutoff date.

    The build is resumable: interrupt it and re-run the same command, and it picks up
    from the last completed shard.
    """
    result = run_build(
        _client(ctx),
        cutoff.date(),
        out,
        ctx.obj["cache"],
        molecules=MOLECULE_CHOICES[molecules],
        releases=_releases(ctx),
        with_taxonomy=taxonomy,
        taxdump_flavour=taxdump_flavour,
        batch_size=batch_size,
        limit=limit,
        stream_catalog=stream_catalog,
        manifest_only=manifest_only,
        progress=progress and not ctx.obj["quiet"],
        quiet=ctx.obj["quiet"],
        command=["hindsight", *sys.argv[1:]],
    )
    click.echo(str(result.out_dir))


@main.command("verify")
@click.argument("out_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
def verify_cmd(out_dir: Path) -> None:
    """Re-check a built database against the digests in its provenance.json."""
    problems = provenance_mod.verify(out_dir)
    if problems:
        for problem in problems:
            click.echo(problem, err=True)
        raise SystemExit(1)
    record = provenance_mod.read(out_dir)
    control = record["time_control"]
    click.echo(
        f"ok: {len(record.get('outputs', {}))} outputs match; "
        f"release {control['release']} ({control['release_date']}), "
        f"cutoff {control['cutoff_date']}"
    )


@main.group("arm-a")
def arm_a_grp() -> None:
    """Arm A: Kraken2 against the time-correct nucleotide database (Day 3).

    The arm sees the reads Day 2 froze and a database `hindsight build` produced for a
    cutoff date. It never sees the ground-truth viral labels, and it never sees a
    post-2019 reference of any kind.
    """


def _arm_config(path: Path):
    from .arms.day3 import load_config

    return load_config(path)


_ARM_CONFIG_OPTION = click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("config/day3.toml"),
    show_default=True,
)


@arm_a_grp.command("build-db")
@_ARM_CONFIG_OPTION
def arm_a_build_db_cmd(config_path: Path) -> None:
    """Assemble the Kraken2 database from a hindsight build directory."""
    from .arms.day3 import prepare_database

    arm = _arm_config(config_path)
    record, clades = prepare_database(arm, log=lambda msg: click.echo(msg, err=True))
    inspect = record.get("inspect", {})
    click.echo(f"{record['db']}: {inspect.get('taxa_in_db', '?')} taxa in the database")
    for name, info in record.get("clades", {}).items():
        click.echo(f"  {name:<22}{info['taxa_in_clade']:>7} taxa  (taxid {info['root_taxid']})")
    for absent in record.get("clades_absent_from_dated_taxonomy", []):
        click.echo(f"  absent from the 2019 taxonomy: {absent}")
    _ = clades


@arm_a_grp.command("run")
@_ARM_CONFIG_OPTION
@click.option(
    "--stage",
    "stages",
    multiple=True,
    type=click.Choice(["database", "index", "series", "negatives"]),
    help="Run only these stages (repeatable). Default: all of them, in order.",
)
@click.option(
    "--only",
    "only_negatives",
    multiple=True,
    help="Classify only these negative runs (repeatable). Implies --stage negatives.",
)
def arm_a_run_cmd(config_path: Path, stages: tuple[str, ...], only_negatives: tuple[str, ...]):
    """Run Arm A: the index sample, the dilution series, and the negatives.

    Resumable. A sample already classified is skipped without re-fetching it, which
    matters because the negatives are fetched and deleted one at a time.
    """
    from .arms.day3 import run_day3

    arm = _arm_config(config_path)
    chosen = list(stages) or (["negatives"] if only_negatives else None)
    record = run_day3(
        arm,
        stages=chosen or ("database", "index", "series", "negatives"),
        only_negatives=list(only_negatives) or None,
        log=lambda msg: click.echo(msg, err=True),
    )
    index = record.get("index_sample") or {}
    if index.get("truth"):
        truth = index["truth"]
        click.echo(
            f"index sample: {truth['detected']:,} of "
            f"{truth.get('truth_fragments_in_nonhost') or truth['truth_fragments']:,} "
            f"ground-truth viral fragments in {truth['clade']}"
        )
    for level in (record.get("series") or {}).get("levels", []):
        click.echo(
            f"  {level['name']}: {level['detected_viral_reads']:,} of "
            f"{level['viral_fragments']:,} detected, "
            f"{level['false_positive_reads']:,} background false positives"
        )
    alarms = record.get("false_alarms") or {}
    if alarms:
        click.echo(
            f"negatives: {alarms['reads_total']:,} {alarms['clade']} reads across "
            f"{alarms['samples']} samples ({alarms['samples_with_at_least_one_hit']} with "
            f"at least one), {alarms['non_host_pairs_total']:,} non-host pairs"
        )


@arm_a_grp.command("diagnose-background")
@_ARM_CONFIG_OPTION
@click.option(
    "--taxid",
    type=int,
    default=694009,
    show_default=True,
    help="The taxon whose background hits to explain. Default: SARS-related coronavirus.",
)
def arm_a_diagnose_cmd(config_path: Path, taxid: int) -> None:
    """Explain the clade hits in the index sample's own background.

    The background has every ground-truth viral fragment removed, so a sarbecovirus
    read in it is either the classifier's mistake or the label's. This re-aligns those
    reads to the truth reference with no mapping-quality threshold and says which.
    """
    from .arms.day3 import diagnose_background_hits, prepare_database

    arm = _arm_config(config_path)
    samples_config = _samples_config(arm.samples_config)
    _, clades = prepare_database(arm, log=lambda msg: click.echo(msg, err=True))
    record = diagnose_background_hits(
        arm, samples_config, clades, taxid=taxid, log=lambda msg: click.echo(msg, err=True)
    )
    click.echo(json.dumps(record, indent=2))


@arm_a_grp.command("report")
@click.option(
    "--provenance",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("results/day3-provenance.json"),
    show_default=True,
)
def arm_a_report_cmd(provenance: Path) -> None:
    """Summarise Arm A: the detection curve, and the false alarms it has to beat."""
    record = json.loads(provenance.read_text())
    database = record.get("database", {})
    release = (database.get("source_release") or {}).get("time_control") or {}
    params = record.get("parameters", {})
    machine = record.get("machine", {})
    runs = record.get("runs") or []
    total_seconds = sum(r.get("wall_clock_seconds", 0) for r in runs) or record.get(
        "wall_clock_seconds", 0
    )
    excluded = release.get("excluded_next_release") or {}
    click.echo(
        f"kraken2 vs RefSeq release {release.get('release', '?')} "
        f"({release.get('release_date', '?')}), cutoff {params.get('cutoff')}, "
        f"first excluded release {excluded.get('release', '?')} "
        f"({excluded.get('release_date', '?')})"
    )
    click.echo(
        f"confidence {params.get('confidence')}, "
        f"minimum-hit-groups {params.get('minimum_hit_groups')}, "
        f"masked library: {params.get('mask_low_complexity')}, "
        f"headline clade: {params.get('headline_clade')}"
    )
    click.echo(
        f"ran on {machine.get('cpu_count')} cores, {machine.get('memory_gib')} GiB, "
        f"{total_seconds / 60:.0f} min wall clock over {len(runs)} invocations; "
        f"negatives at {params.get('negatives_depth')} depth"
    )
    check = (record.get("series") or {}).get("composition_check")
    if check:
        click.echo(
            "composition check: "
            + (
                f"{check['level']} composed == classified whole"
                if check.get("identical")
                else f"MISMATCH on {check['taxa_disagreeing']} taxa"
            )
        )

    index = record.get("index_sample") or {}
    if index:
        truth = index.get("truth", {})
        click.echo("")
        click.echo(
            f"index sample ({index['name']}): {index['pairs']:,} non-host pairs, "
            f"{truth.get('detected', 0):,} of "
            f"{truth.get('truth_fragments_in_nonhost') or 0:,} ground-truth viral "
            f"fragments in {truth.get('clade')} "
            f"({100 * (truth.get('sensitivity') or 0):.1f}%), "
            f"{truth.get('false_positive', 0):,} false positives"
        )

    levels = (record.get("series") or {}).get("levels", [])
    if levels:
        # Both clade columns, because the dated taxonomy puts the nearest relative in
        # the database outside Betacoronavirus and the difference is 5,000 reads.
        click.echo("")
        click.echo(
            f"{'level':<14}{'viral':>8}{'per million':>13}"
            f"{'Coronaviridae':>15}{'Betacoronavirus':>17}{'sensitivity':>13}"
        )
        for level in levels:
            per_million = (level.get("per_million") or {}).get("post_qc_pairs")
            by_clade = level.get("detected_by_clade") or {}
            family = by_clade.get("coronaviridae")
            genus = by_clade.get("betacoronavirus", level["detected_viral_reads"])
            viral = level["viral_fragments"]
            click.echo(
                f"{level['name']:<14}{viral:>8,}"
                f"{per_million if per_million is None else round(per_million, 2):>13}"
                f"{'-' if family is None else format(family, ','):>15}"
                f"{genus:>17,}"
                f"{100 * ((family or genus) / viral if viral else 0):>12.1f}%"
            )
        diagnosis = record.get("background_diagnosis") or {}
        explained = (
            "which the record explains as reads the truth label missed"
            if diagnosis.get("all_below_the_label_threshold")
            else "which the record does not explain"
        )
        background = levels[0]["false_positive_by_clade"]["coronaviridae"]
        click.echo(
            f"sensitivity is at family level; the series' own background contributes "
            f"{background} Coronaviridae reads at every level, {explained}"
        )

    alarms = record.get("false_alarms") or {}
    if alarms:
        click.echo("")
        click.echo(
            f"false alarms in {alarms['samples']} negatives at "
            f"{params.get('negatives_depth')} depth: {alarms['reads_total']:,} "
            f"{alarms['clade']} reads in {alarms['non_host_pairs_total']:,} non-host pairs"
        )
        click.echo(
            f"  {alarms['samples_with_at_least_one_hit']} of {alarms['samples']} samples "
            f"have at least one; most in one sample: {alarms['reads_max_in_one_sample']:,}; "
            f"median {alarms['reads_median_per_sample']:,}"
        )
        novel = alarms.get("excluding_endemic_human_coronaviruses") or {}
        if novel:
            click.echo(
                f"  excluding the endemic human coronaviruses "
                f"({alarms.get('endemic_human_coronavirus_reads_total', 0):,} reads of "
                f"HKU1, OC43, 229E, NL63): {novel['reads_total']:,} reads in "
                f"{novel['samples_with_at_least_one_hit']} of {alarms['samples']} samples, "
                f"{novel['reads_per_million_nonhost']:.3f} per million"
            )
        for negative in record.get("negatives", []):
            taxa = negative.get("coronavirus_taxa") or []
            if not taxa:
                continue
            click.echo(
                f"  {negative['accession']}: "
                + ", ".join(f"{t['reads']} {t['taxon'] or t['taxid']}" for t in taxa)
            )


@main.group("arm-b")
def arm_b_grp() -> None:
    """Arm B: MEGAHIT plus DIAMOND against the time-correct protein set (Day 4).

    The reference-free arm, and the one that automates what Wu et al. actually did:
    assemble the non-host reads, then match the contigs by protein homology against
    what was public before December 2019. Like Arm A it never sees the ground-truth
    viral labels and never sees a post-2019 reference.
    """


def _arm_b_config(path: Path):
    from .arms.day4 import load_config

    return load_config(path)


_ARM_B_CONFIG_OPTION = click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("config/day4.toml"),
    show_default=True,
)


@arm_b_grp.command("build-db")
@_ARM_B_CONFIG_OPTION
def arm_b_build_db_cmd(config_path: Path) -> None:
    """Assemble the DIAMOND protein database from a hindsight build directory."""
    from .arms.day4 import prepare_database

    arm = _arm_b_config(config_path)
    record, _taxonomy, _clades = prepare_database(arm, log=lambda msg: click.echo(msg, err=True))
    click.echo(
        f"{record['db']}: {record.get('protein_records_in_source', '?'):,} protein records, "
        f"{record.get('bytes', 0) / 2**20:.0f} MiB"
    )
    taxonmap = record.get("taxonmap") or {}
    if "records_with_a_taxid" in taxonmap:
        click.echo(
            f"  taxon map: {taxonmap['records_with_a_taxid']:,} of "
            f"{taxonmap['protein_records']:,} proteins carry a taxid"
        )
    proteins = record.get("protein_counts_by_clade") or {}
    for name, info in record.get("clades", {}).items():
        click.echo(
            f"  {name:<22}{info['taxa_in_clade']:>7} taxa  "
            f"{proteins.get(name, 0):>7} proteins  (taxid {info['root_taxid']})"
        )
    for absent in record.get("clades_absent_from_dated_taxonomy", []):
        click.echo(f"  absent from the 2019 taxonomy: {absent}")


@arm_b_grp.command("run")
@_ARM_B_CONFIG_OPTION
@click.option(
    "--stage",
    "stages",
    multiple=True,
    type=click.Choice(["database", "index", "series", "negatives"]),
    help="Run only these stages (repeatable). Default: all of them, in order.",
)
@click.option(
    "--only",
    "only_negatives",
    multiple=True,
    help="Search only these negative runs (repeatable). Implies --stage negatives.",
)
@click.option(
    "--level",
    "only_levels",
    multiple=True,
    help="Assemble only these dilution levels (repeatable). Implies --stage series.",
)
@click.option(
    "--read-level/--no-read-level",
    "read_level",
    default=None,
    help=(
        "Override the config's read_level setting for this run. --no-read-level does "
        "the contig half only, which is how to get the assembly results before paying "
        "for the read half's one expensive pass; re-run with it on afterwards and every "
        "assembly and search is reused from its marker. The recorded parameters say "
        "what actually ran, not what the config asked for."
    ),
)
def arm_b_run_cmd(
    config_path: Path,
    stages: tuple[str, ...],
    only_negatives: tuple[str, ...],
    only_levels: tuple[str, ...],
    read_level: bool | None,
):
    """Run Arm B: the index sample, the dilution series, and the negatives.

    Resumable at the level of one assembly and one search, which matters more here than
    it did for Arm A: a dilution level is a full MEGAHIT run and there are six of them
    before the negatives start.
    """
    import dataclasses

    from .arms.day4 import run_day4

    arm = _arm_b_config(config_path)
    if read_level is not None:
        arm = dataclasses.replace(arm, read_level=read_level)
    chosen = list(stages)
    if not chosen:
        if only_negatives:
            chosen = ["negatives"]
        elif only_levels:
            chosen = ["series"]
    record = run_day4(
        arm,
        stages=chosen or ("database", "index", "series", "negatives"),
        only_negatives=list(only_negatives) or None,
        only_levels=list(only_levels) or None,
        log=lambda msg: click.echo(msg, err=True),
    )
    _arm_b_echo_summary(record)


def _arm_b_echo_summary(record: dict) -> None:
    headline = (record.get("parameters") or {}).get("headline_clade", "coronaviridae")
    index = record.get("index_sample") or {}
    if index:
        clade = (index.get("search") or {}).get("clades", {}).get(headline, {})
        click.echo(
            f"index sample: {index['assembly']['contig_count']:,} contigs, "
            f"{clade.get('contigs', 0)} in {headline}, "
            f"longest {clade.get('longest_contig', 0):,} bp"
        )
        readback = index.get("reads_on_clade_contigs") or {}
        if readback:
            click.echo(
                f"  reads mapping back: {readback.get('detected', 0):,} of "
                f"{readback.get('truth_fragments', 0):,} ground-truth viral fragments, "
                f"{readback.get('false_positive', 0):,} others"
            )
    for level in (record.get("series") or {}).get("levels", []):
        clade = (level.get("search") or {}).get("clades", {}).get(headline, {})
        line = (
            f"  {level['name']}: {level['viral_fragments']:,} viral fragments -> "
            f"{clade.get('contigs', 0)} {headline} contigs"
        )
        if clade.get("longest_contig"):
            line += f", longest {clade['longest_contig']:,} bp"
        reads_half = level.get("read_level") or {}
        if reads_half:
            line += f"; reads: {reads_half.get('detected', 0):,} detected"
        click.echo(line)
    alarms = record.get("false_alarms") or {}
    if alarms:
        click.echo(
            f"negatives: {alarms['samples_with_at_least_one_clade_contig']} of "
            f"{alarms['samples']} samples have a {alarms['clade']} contig "
            f"({alarms['contigs_total']} contigs, "
            f"{alarms['non_host_pairs_total']:,} non-host pairs)"
        )


@arm_b_grp.command("diagnose-contigs")
@_ARM_B_CONFIG_OPTION
@click.option(
    "--clade",
    default="coronaviridae",
    show_default=True,
    help="Which clade's contigs to describe.",
)
@click.option(
    "--sample",
    "sample",
    default=None,
    help="A negative's run accession. Default: the index sample.",
)
@click.option(
    "--dated-reference",
    "dated_reference",
    default=None,
    help=(
        "Align to this accession.version from the dated manifest instead of the modern "
        "reference. Answers how far the thing in the patient sat from what the database "
        "actually held, and involves no post-2019 sequence."
    ),
)
def arm_b_diagnose_cmd(
    config_path: Path, clade: str, sample: str | None, dated_reference: str | None
) -> None:
    """Say what the arm's clade contigs actually are, by aligning them to the virus.

    The arm scores itself on a November 2019 database and reports "something
    coronavirus-shaped that is not in it". That is the honest output and it does not say
    how good the reconstruction was, which is the number a reader wants, since Wu et
    al.'s result was a genome rather than a classification. This aligns those contigs to
    the modern reference afterwards and reports the identity and the genome coverage.
    Nothing it computes feeds back into the arm's score.
    """
    from .arms.day4 import diagnose_contigs, prepare_database

    arm = _arm_b_config(config_path)
    samples_config = _samples_config(arm.samples_config)
    _, taxonomy, clades = prepare_database(arm, log=lambda msg: click.echo(msg, err=True))
    record = diagnose_contigs(
        arm,
        samples_config,
        taxonomy,
        clades,
        clade=clade,
        accession=sample,
        dated_reference=dated_reference,
        log=lambda msg: click.echo(msg, err=True),
    )
    click.echo(json.dumps(record, indent=2))


@arm_b_grp.command("report")
@click.option(
    "--provenance",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("results/day4-provenance.json"),
    show_default=True,
)
def arm_b_report_cmd(provenance: Path) -> None:
    """Summarise Arm B: the assembly-limited detection curve, and the false alarms."""
    record = json.loads(provenance.read_text())
    database = record.get("database", {})
    release = (database.get("source_release") or {}).get("time_control") or {}
    params = record.get("parameters", {})
    machine = record.get("machine", {})
    runs = record.get("runs") or []
    total_seconds = sum(r.get("wall_clock_seconds", 0) for r in runs) or record.get(
        "wall_clock_seconds", 0
    )
    excluded = release.get("excluded_next_release") or {}
    headline = params.get("headline_clade", "coronaviridae")
    click.echo(
        f"megahit + diamond blastx vs RefSeq release {release.get('release', '?')} "
        f"({release.get('release_date', '?')}) proteins, cutoff {params.get('cutoff')}, "
        f"first excluded release {excluded.get('release', '?')} "
        f"({excluded.get('release_date', '?')})"
    )
    click.echo(
        f"{database.get('protein_records_in_source', 0):,} proteins; "
        f"--{params.get('sensitivity')}, e-value {params.get('evalue')}, "
        f"LCA over the top {params.get('lca_top_percent')}% of bitscore; "
        f"headline clade: {headline}"
    )
    click.echo(
        f"ran on {machine.get('cpu_count')} cores, {machine.get('memory_gib')} GiB, "
        f"{total_seconds / 60:.0f} min wall clock over {len(runs)} invocations; "
        f"negatives at {params.get('negatives_depth')} depth"
    )
    click.echo("")
    _arm_b_echo_summary(record)

    index = record.get("index_sample") or {}
    comparison = index.get("sensitivity_comparison") or {}
    if comparison:
        theirs = (comparison.get("search") or {}).get("clades", {}).get(headline, {})
        mine = (index.get("search") or {}).get("clades", {}).get(headline, {})
        click.echo("")
        click.echo(
            f"sensitivity setting: --{index.get('sensitivity')} found "
            f"{mine.get('contigs', 0)} {headline} contigs, "
            f"--{comparison.get('sensitivity')} found {theirs.get('contigs', 0)}"
        )
    reads_half = index.get("read_level") or {}
    if reads_half.get("truth"):
        truth = reads_half["truth"]
        click.echo(
            f"read-level search on the index sample: {truth['detected']:,} of "
            f"{truth['truth_fragments']:,} ground-truth viral fragments in {headline}, "
            f"{truth['false_positive']:,} false positives"
        )
    for negative in record.get("negatives", []):
        detail = (negative.get("search") or {}).get("contig_detail") or []
        if not negative.get("false_alarm_contigs"):
            continue
        click.echo(
            f"  {negative['accession']}: "
            + ", ".join(f"{e['length']} bp {e['lca_taxon'] or e['lca_taxid']}" for e in detail[:4])
        )


@main.group("samples")
def samples_grp() -> None:
    """Fetch, QC and host-filter the benchmark samples, and freeze the dilution series.

    This is the sample side of the benchmark. It never consults the dated database, and
    the only post-2019 reference it touches is the SARS-CoV-2 genome, used solely to
    label which fragments are viral. The detection arms never see that label.
    """


def _samples_config(path: Path):
    from .samples.registry import load_config

    return load_config(path)


@samples_grp.command("plan")
@click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("config/day2.toml"),
    show_default=True,
)
@click.option("--json", "as_json", is_flag=True, help="Emit JSON instead of a table.")
def samples_plan_cmd(config_path: Path, as_json: bool) -> None:
    """Show which runs will be processed, and the dilution levels, without fetching."""
    from .samples.registry import samples_from_config

    config = _samples_config(config_path)
    specs = samples_from_config(config)
    if as_json:
        click.echo(
            json.dumps(
                {
                    "samples": [
                        {"accession": s.accession, "role": s.role, "notes": s.notes} for s in specs
                    ],
                    "levels": list(config.levels),
                    "dilution_seed": config.dilution_seed,
                    "negatives_seed": config.negatives_seed,
                },
                indent=2,
            )
        )
        return
    index = [s for s in specs if s.is_index]
    negatives = [s for s in specs if not s.is_index]
    click.echo(f"index sample:   {', '.join(s.accession for s in index)}")
    click.echo(f"negatives ({len(negatives)}):  {' '.join(s.accession for s in negatives)}")
    if negatives:
        click.echo(f"chosen by:      {negatives[0].notes}")
    click.echo(f"levels:         {', '.join(str(x) for x in config.levels)} viral fragments")
    click.echo(f"dilution seed:  {config.dilution_seed} (nested, sha256 ranking)")


@samples_grp.command("select-negatives")
@click.option(
    "--manifest",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Annotated ENA filereport for the negatives project.",
)
@click.option("--count", type=int, default=20, show_default=True)
@click.option("--seed", type=int, default=20191201, show_default=True)
def samples_select_cmd(manifest: Path, count: int, seed: int) -> None:
    """Re-derive the negative-control subset from the project manifest and a seed."""
    from .samples.registry import choose_negatives, patient_rows, read_negatives_manifest

    rows = read_negatives_manifest(manifest)
    chosen = choose_negatives(rows, count=count, seed=seed)
    click.echo(
        f"# {len(patient_rows(rows))} patient runs in the manifest; "
        f"{len(chosen)} chosen with seed {seed}",
        err=True,
    )
    for accession in chosen:
        click.echo(accession)


@samples_grp.command("run")
@click.option(
    "--config",
    "config_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("config/day2.toml"),
    show_default=True,
)
@click.option(
    "--only",
    multiple=True,
    help="Process only these run accessions (repeatable). Default: everything in the config.",
)
@click.option(
    "--skip",
    multiple=True,
    help=(
        "Leave these run accessions alone (repeatable). Needed once a sample's "
        "intermediates have been pruned: its stage markers remain, but the files they "
        "describe are gone, so re-running it re-fetches gigabytes to rebuild what "
        "nothing downstream reads any more."
    ),
)
@click.option(
    "--skip-series",
    is_flag=True,
    help="Fetch, QC and filter, but do not write the dilution series.",
)
def samples_run_cmd(
    config_path: Path, only: tuple[str, ...], skip: tuple[str, ...], skip_series: bool
) -> None:
    """Run the whole sample pipeline. Resumable: re-running finishes what was left."""
    from .samples.pipeline import run_day2

    config = _samples_config(config_path)
    record = run_day2(
        config,
        only=list(only) or None,
        skip=list(skip) or None,
        skip_series=skip_series,
        log=lambda msg: click.echo(msg, err=True),
    )
    index = next((s for s in record["samples"] if s["role"] == "index"), None)
    if index:
        click.echo(
            f"{index['accession']}: {index['viral_pairs']} viral fragments of "
            f"{index['raw_pairs']} raw ({index['viral_fraction_raw'] * 100:.3f}%)"
        )
    series = record.get("dilution_series")
    if series:
        for level in series["levels"]:
            click.echo(
                f"  {level['name']}: {level['viral_fragments']} viral, "
                f"{level['viral_fragments_per_million']:.2f} per million"
            )


@samples_grp.command("materialise")
@click.option(
    "--provenance",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("work/day2-provenance.json"),
    show_default=True,
)
@click.option("--level", "level_name", required=True, help="Level name, e.g. viral001000.")
@click.option(
    "--out",
    "out_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=Path("work/series/materialised"),
    show_default=True,
)
@click.option(
    "--no-verify",
    is_flag=True,
    help="Skip re-checking each part against its recorded digest. Not recommended.",
)
def samples_materialise_cmd(
    provenance: Path, level_name: str, out_dir: Path, no_verify: bool
) -> None:
    """Write one dilution level out as an ordinary FASTQ pair.

    The series is stored as parts — the background once, the viral fragments once per
    level — because five materialised levels are 21 GB of files that are 99.9%
    identical. A classifier that reads stdin needs no materialisation at all:

    \b
        cat work/series/parts/background_1.fastq.gz \\
            work/series/parts/viral001000_viral_1.fastq.gz | kraken2 --gzip-compressed ...
    """
    from .samples.pipeline import materialise_level

    record = json.loads(provenance.read_text())
    series = record.get("dilution_series")
    if not series:
        raise click.ClickException(f"{provenance} holds no dilution series")
    r1, r2 = materialise_level(series, level_name, out_dir, verify=not no_verify)
    for path in (r1, r2):
        click.echo(f"{path} ({path.stat().st_size:,} bytes)")


@samples_grp.command("report")
@click.option(
    "--provenance",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("work/day2-provenance.json"),
    show_default=True,
)
def samples_report_cmd(provenance: Path) -> None:
    """Summarise a completed run: counts, denominators, the series, and what it cost.

    The provenance record is the authority, but it is 40 KB of JSON. This is the table
    that goes in a methods section.
    """
    record = json.loads(provenance.read_text())
    machine = record.get("machine", {})
    # wall_clock_seconds is the last invocation only. A cohort this size is built
    # over several invocations (a crash, a resume, a rerun of one sample), so the
    # honest figure is the sum of the stamps, with the count said out loud.
    runs = record.get("runs") or []
    total_seconds = sum(r.get("wall_clock_seconds", 0) for r in runs) or record.get(
        "wall_clock_seconds", 0
    )
    invocations = (
        f" over {len(runs)} invocations" if len(runs) > 1 else " in one invocation" if runs else ""
    )
    click.echo(
        f"ran on {machine.get('instance_type') or machine.get('machine', 'unknown')}, "
        f"{machine.get('cpu_count')} cores, {machine.get('memory_gib')} GiB, "
        f"{total_seconds / 60:.0f} min wall clock{invocations}"
    )
    cap = (record.get("parameters") or {}).get("negatives_max_pairs")
    if cap:
        capped = sum(
            1
            for sample in record.get("samples", [])
            if sample.get("role") == "negative" and sample.get("cap_fraction", 1.0) < 1.0
        )
        click.echo(f"negatives capped at {cap:,} pairs each; the cap bound on {capped} of them")
    drift = record.get("tool_pin_drift") or []
    tools_line = ", ".join(f"{t['name']} {t['version']}" for t in record.get("tools", []))
    click.echo(f"tools: {tools_line}" + (f"  [drift: {', '.join(drift)}]" if drift else ""))

    click.echo("")
    click.echo(
        f"{'run':<14}{'role':<10}{'raw':>12}{'redacted':>12}"
        f"{'post-QC':>12}{'non-host':>12}{'viral':>9}"
    )
    for sample in record.get("samples", []):
        click.echo(
            f"{sample['accession']:<14}{sample['role']:<10}"
            f"{sample['raw_pairs']:>12,}{sample['redacted_pairs']:>12,}"
            f"{sample['qc_pairs']:>12,}{sample['nonhost_pairs']:>12,}"
            f"{sample['viral_pairs'] or '':>9}"
        )

    series = record.get("dilution_series")
    if not series:
        return
    click.echo("")
    click.echo(
        f"dilution series: seed {series['seed']}, "
        f"{'nested' if series['nested'] else 'independent draws'}, "
        f"{series['background_fragments']:,} background fragments"
    )
    denominators = series.get("denominators", {})
    keys = list(denominators)
    click.echo(
        f"{'level':<14}{'viral':>8}{'total':>12}"
        + "".join(f"{k.replace('_', ' '):>20}" for k in keys)
    )
    for level in series["levels"]:
        row = f"{level['name']:<14}{level['viral_fragments']:>8,}{level['total_fragments']:>12,}"
        for key in keys:
            value = level.get("per_million", {}).get(key)
            row += f"{value:>20.3f}" if value is not None else f"{'-':>20}"
        click.echo(row)
    click.echo("")
    headline = series.get("headline_denominator", "post_qc_pairs")
    click.echo(
        "per-million columns are viral fragments per million; the headline denominator "
        f"is {headline.replace('_', ' ')}"
    )


if __name__ == "__main__":  # pragma: no cover
    main()


@main.group("score")
def score_grp() -> None:
    """Day 5: score both arms in one frame and draw the detection threshold.

    Reads the two arms' provenance records and writes the benchmark's headline output —
    a limit of detection per method, in viral reads per million, at a stated false-alarm
    rate measured on the twenty pre-COVID negatives. It runs no tool and touches no
    reads, which is the point: the clade, the false-alarm framing and both arms'
    thresholds can be moved by a reader in seconds, because Days 3 and 4 recorded every
    clade at every dilution level rather than only their own headline.
    """


_SCORE_CONFIG_OPTION = click.option(
    "--config",
    "config_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=Path("config/day5.toml"),
    show_default=True,
)


@score_grp.command("run")
@_SCORE_CONFIG_OPTION
@click.option("--clade", default=None, help="Override the clade that counts as a detection.")
@click.option(
    "--framing",
    type=click.Choice(["all", "excluding_endemic"]),
    default=None,
    help="Whether the endemic human coronaviruses count as false alarms.",
)
@click.option(
    "--budget",
    "budget",
    type=int,
    default=None,
    help="How many of the twenty negatives may alarm. Default 0.",
)
@click.option(
    "--arm-b-sensitivity",
    type=click.Choice(["as_run", "comparison"]),
    default=None,
    help="Score Arm B's --very-sensitive run (as_run) or DIAMOND's stock mode (comparison).",
)
@click.option(
    "--out",
    "out",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Where to write the record. Default: the config's provenance path.",
)
def score_run_cmd(
    config_path: Path,
    clade: str | None,
    framing: str | None,
    budget: int | None,
    arm_b_sensitivity: str | None,
    out: Path | None,
) -> None:
    """Score both arms and write results/day5-provenance.json."""
    import dataclasses as _dc

    from .arms.day5 import load_config as load_score_config
    from .arms.day5 import report as score_report
    from .arms.day5 import score as run_score

    config = load_score_config(config_path)
    overrides: dict = {}
    if clade:
        overrides["clade"] = clade
    if framing:
        overrides["framing"] = framing
    if budget is not None:
        overrides["false_alarm_budget"] = budget
    if arm_b_sensitivity:
        overrides["arm_b_sensitivity"] = arm_b_sensitivity
    if out:
        overrides["provenance"] = out
    if overrides:
        config = _dc.replace(config, **overrides)
    record = run_score(config)
    for line in score_report(record):
        click.echo(line)
    click.echo("")
    click.echo(f"wrote {config.provenance}")


@score_grp.command("report")
@click.option(
    "--provenance",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=Path("results/day5-provenance.json"),
    show_default=True,
)
def score_report_cmd(provenance: Path) -> None:
    """Print the cross-arm tables from a record that has already been written."""
    from .arms.day5 import report as score_report

    record = json.loads(provenance.read_text())
    for line in score_report(record):
        click.echo(line)
