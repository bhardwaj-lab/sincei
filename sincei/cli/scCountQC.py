from __future__ import annotations

import logging
import re
import sys
import warnings
from pathlib import Path
from typing import Annotated, cast

import anndata as ad
import pandas as pd
import scanpy as sc
import typer

from sincei import _sincei as internal
from sincei.utils import gini

from ._common_args import (
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    configure_logging,
    log_parameters,
    override,
    preprocess_args,
)
from ._parsers import validate_anndata

DESCRIPTION = (
    "Perform quality control and filter cells and regions from a cell-by-feature "
    "matrix.\n\n"
    "``scCountQC`` calculates multiple quality controls metrics on the input .h5ad "
    "file (output of scCountReads) and (optionally) filters the input file based on "
    "filterCellArgs/filterRegionArgs. The output is either an updated .h5ad object (if "
    "filtering is requested) or the filtering metrics (--outMetrics)."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)

_QC = "QC options"

Bounds = dict[str, tuple[float, float | None]]


def parse_filter_args(value: str | None, option: str) -> Bounds:
    """Parse ``"metric: min, max; metric: min"`` into ``{metric: (min, max)}``.

    A metric given a single value has no upper bound.
    """
    bounds: Bounds = {}
    for item in (value or "").split(";"):
        if not item.strip():
            continue
        metric, separator, numbers = item.partition(":")
        try:
            limits = [float(number) for number in numbers.split(",")]
        except ValueError:
            limits = []
        if not separator or not metric.strip() or len(limits) not in {1, 2}:
            msg = (
                f"cannot parse {item.strip()!r}; expected 'metric: min, max' or "
                "'metric: min'"
            )
            raise typer.BadParameter(msg, param_hint=option)
        bounds[metric.strip()] = (limits[0], limits[1] if len(limits) == 2 else None)
    return bounds


def within_bounds(frame: pd.DataFrame, bounds: Bounds) -> pd.Series:
    """Rows of ``frame`` whose metrics are all inside their bounds."""
    keep = pd.Series(data=True, index=frame.index)
    for metric, (low, high) in bounds.items():
        if metric not in frame.columns:
            logging.warning("Filter metric %r is not available. Skipping.", metric)
            continue
        keep &= frame[metric] >= low
        if high is not None:
            keep &= frame[metric] <= high
    return keep


def blacklisted_regions(var: pd.DataFrame, blacklist: list[str]) -> list[str]:
    """Names of the regions in ``var`` that overlap a region of ``blacklist``."""
    annotation = internal.parse_annotation(blacklist)
    return [
        name
        for name, chrom, start, end in zip(
            var.index, var["chrom"].astype(str), var["start"], var["end"], strict=True
        )
        if annotation.find_overlaps(chrom, int(start), int(end))
    ]


def filter_adata(
    adata: ad.AnnData,
    region_filters: Bounds | None = None,
    cell_filters: Bounds | None = None,
    bad_chroms: list[str] | None = None,
    bad_regions: list[str] | None = None,
    bad_cells: list[str] | None = None,
) -> ad.AnnData:
    """Remove the regions and cells that fail a filter, or are blacklisted."""
    var = cast("pd.DataFrame", adata.var)
    obs = cast("pd.DataFrame", adata.obs)

    regions = within_bounds(var, region_filters or {})
    if bad_chroms:
        regions &= ~var["chrom"].astype(str).isin(bad_chroms)
    if bad_regions:
        regions &= ~adata.var_names.isin(bad_regions)

    cells = within_bounds(obs, cell_filters or {})
    if bad_cells:
        cells &= ~adata.obs_names.isin(bad_cells)

    return adata[cells.to_numpy(), regions.to_numpy()].copy()


def describe_metrics(frame: pd.DataFrame, label: str, unit: str) -> None:
    numeric = frame.select_dtypes(include="number")
    typer.echo(f"\n{label} metrics:\nTotal {unit}: {len(frame)}")
    typer.echo(pd.DataFrame({"min": numeric.min(), "max": numeric.max()}).to_string())


@app.callback(invoke_without_command=True)
def main(
    # Input and output options
    input: Annotated[str, INPUT_OUTPUT_OPTS["h5ad_file"]],
    out_file: Annotated[
        str | None,
        override(
            INPUT_OUTPUT_OPTS["out_file"],
            help=(
                "The filtered .h5ad file. Required when a filter or a blacklist is "
                "given."
            ),
        ),
    ] = None,
    # QC options
    describe: Annotated[
        bool,
        typer.Option(
            "-d",
            "--describe",
            rich_help_panel=_QC,
            help="Print a list of cell and region metrics available for QC/filtering.",
        ),
    ] = False,
    out_metrics: Annotated[
        str | None,
        typer.Option(
            "-om",
            "--outMetrics",
            rich_help_panel=_QC,
            help=(
                "Prefix of the output file with calculated QC metrics. If given, the "
                "cell metrics are printed in <prefix>.cells.tsv and region metrics as "
                "<prefix>.regions.tsv."
            ),
        ),
    ] = None,
    filter_cell_args: Annotated[
        str | None,
        typer.Option(
            "-fc",
            "--filterCellArgs",
            rich_help_panel=_QC,
            help=(
                'List of arguments to filter cells. The format is "arg_name: minvalue, '
                'maxvalue; arg_name: minvalue; ...." where arg_name is a cell QC '
                'metric present in the input h5ad file. Run with "--describe" to view '
                "all available metrics. The two values are used as lower and upper "
                "bounds to filter cells. A single value is a lower bound."
            ),
        ),
    ] = None,
    filter_region_args: Annotated[
        str | None,
        typer.Option(
            "-fr",
            "--filterRegionArgs",
            rich_help_panel=_QC,
            help=(
                'List of arguments to filter regions. The format is "arg_name: '
                'minvalue, maxvalue; arg_name: minvalue; ...." where arg_name is a '
                "region QC metric present in the input h5ad file. Run with "
                '"--describe" to view all available metrics. The two values are used '
                "as lower and upper bounds to filter regions. A single value is a "
                "lower bound."
            ),
        ),
    ] = None,
    region_blacklist: Annotated[
        list[str] | None,
        typer.Option(
            "-rb",
            "--regionBlacklist",
            metavar="BED",
            rich_help_panel=_QC,
            help=(
                "A BED or GTF file containing regions that should be excluded from all "
                "analyses. Regions in the anndata object that overlap with blacklisted "
                "regions will be removed."
            ),
        ),
    ] = None,
    cell_blacklist: Annotated[
        str | None,
        typer.Option(
            "-cb",
            "--cellBlacklist",
            metavar=".txt",
            rich_help_panel=_QC,
            help=(
                "A file with the cells to exclude, one per line. The cells must be "
                "named as in the input object (sample::barcode)."
            ),
        ),
    ] = None,
    chrom_blacklist: Annotated[
        list[str] | None,
        typer.Option(
            "-chb",
            "--chromBlacklist",
            metavar="CHR",
            rich_help_panel=_QC,
            help="A space separated list of chromosomes to exclude, e.g. chrM chrUn.",
        ),
    ] = None,
    # Other options
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    if verbose:
        log_parameters(
            input=input,
            out_file=out_file,
            describe=describe,
            out_metrics=out_metrics,
            filter_cell_args=filter_cell_args,
            filter_region_args=filter_region_args,
            region_blacklist=region_blacklist,
            cell_blacklist=cell_blacklist,
            chrom_blacklist=chrom_blacklist,
        )
    else:
        warnings.filterwarnings("ignore")

    cell_filters = parse_filter_args(filter_cell_args, "--filterCellArgs")
    region_filters = parse_filter_args(filter_region_args, "--filterRegionArgs")
    filtering = bool(
        cell_filters
        or region_filters
        or region_blacklist
        or cell_blacklist
        or chrom_blacklist
    )
    if filtering and out_file is None and not describe:
        msg = "is required to write the filtered file"
        raise typer.BadParameter(msg, param_hint="--outFile")

    adata = validate_anndata(ad.read_h5ad(input), input)

    try:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
    except IndexError as exc:
        sys.stderr.write("Too few regions in the input file to perform QC.\n")
        raise typer.Exit(code=1) from exc
    adata.obs["gini_coefficient"] = [gini(i, adata.X) for i in range(adata.n_obs)]
    obs = cast("pd.DataFrame", adata.obs)
    var = cast("pd.DataFrame", adata.var)

    if out_metrics is not None:
        prefix = re.sub(r"\.(txt|tsv|csv)$", "", out_metrics)
        obs.to_csv(f"{prefix}.cells.tsv", sep="\t", index_label="Cell_ID")
        var.to_csv(f"{prefix}.regions.tsv", sep="\t", index_label="Feature_ID")

    if describe:
        describe_metrics(obs, "Cell", "cells")
        describe_metrics(var, "Feature", "features")
        return 0

    if not filtering or out_file is None:
        return 0

    bad_regions = None
    if region_blacklist:
        bad_regions = blacklisted_regions(var, region_blacklist)
        logging.info("Found %d regions overlapping the blacklist.", len(bad_regions))
    bad_cells = None
    if cell_blacklist:
        bad_cells = [
            line.strip()
            for line in Path(cell_blacklist).read_text().splitlines()
            if line.strip()
        ]

    typer.echo("Applying filters")
    filtered = filter_adata(
        adata,
        region_filters=region_filters,
        cell_filters=cell_filters,
        bad_chroms=chrom_blacklist,
        bad_regions=bad_regions,
        bad_cells=bad_cells,
    )
    typer.echo(f"Cells post-filtering: {filtered.n_obs}")
    typer.echo(f"Features post-filtering: {filtered.n_vars}")
    filtered.write_h5ad(out_file)
    return 0


def cli() -> None:
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
