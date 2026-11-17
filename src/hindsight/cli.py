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
        click.echo(
            f"first excluded:    release {following.number} " f"({following.date.isoformat()})"
        )
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
                        {"accession": s.accession, "role": s.role, "notes": s.notes}
                        for s in specs
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
    "--skip-series",
    is_flag=True,
    help="Fetch, QC and filter, but do not write the dilution series.",
)
def samples_run_cmd(config_path: Path, only: tuple[str, ...], skip_series: bool) -> None:
    """Run the whole sample pipeline. Resumable: re-running finishes what was left."""
    from .samples.pipeline import run_day2

    config = _samples_config(config_path)
    record = run_day2(
        config,
        only=list(only) or None,
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


if __name__ == "__main__":  # pragma: no cover
    main()
