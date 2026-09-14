from __future__ import annotations

import sys
import warnings
from pathlib import Path
from typing import Annotated, cast

import anndata as ad
import numpy as np
import pandas as pd
import typer
from scipy import io, sparse

from sincei.utils import chromosome_to_numeric

from ._common_args import (
    INPUT_OUTPUT_OPTS,
    OTHER_OPTS,
    ExportFormat,
    configure_logging,
    log_parameters,
    override,
    preprocess_args,
)
from ._parsers import parse_region, validate_anndata

DESCRIPTION = (
    "Export .h5ad objects to other formats.\n\n"
    "``scExportSignal`` exports a sincei-supported .h5ad (AnnData) file to other "
    "formats."
)


app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
    help=DESCRIPTION,
    context_settings={"help_option_names": []},
)


def region_mask(var: pd.DataFrame, region: str) -> np.ndarray:
    """Features on the region's chromosome that lie entirely inside the region."""
    chrom, start, end = parse_region(region)
    mask = var["chrom"].astype(str) == chrom
    if start is not None and end is not None:
        mask &= (var["start"] >= start) & (var["end"] <= end)
    return mask.to_numpy()


def dense_counts(adata: ad.AnnData) -> np.ndarray:
    counts = adata.X
    to_dense = getattr(counts, "toarray", None)
    return np.asarray(counts if to_dense is None else to_dense())


def write_bedgraph_matrix(adata: ad.AnnData, path: str) -> None:
    """Write the BedGraphMatrix format: one row per feature, one column per cell."""
    var = cast("pd.DataFrame", adata.var)
    chrom = var["chrom"].to_numpy()
    start = var["start"].to_numpy()
    end = var["end"].to_numpy()
    row_data = dense_counts(adata).T  # features x cells

    # Convert chromosomes and region bounds to numeric values for sorting
    chrom_order = chromosome_to_numeric(chrom)
    frame = pd.DataFrame({
        "chrom_numeric": pd.Series(chrom).map(chrom_order).to_numpy(),
        "start_numeric": start.astype(int),
        "end_numeric": end.astype(int),
        "chrom": chrom,
        "start": start.astype(int),
        "end": end.astype(int),
    })
    frame = pd.concat([frame, pd.DataFrame(row_data.astype(float))], axis=1)
    frame = frame.sort_values(["chrom_numeric", "start_numeric", "end_numeric"])
    frame = frame.drop(["chrom_numeric", "start_numeric", "end_numeric"], axis=1)
    frame.to_csv(path, sep="\t", header=False, index=False)


def write_matrix_market(adata: ad.AnnData, prefix: str) -> None:
    """Write the counts as MatrixMarket, with the cell and feature names beside."""
    Path(f"{prefix}.rownames.txt").write_text("\n".join(adata.obs_names) + "\n")
    Path(f"{prefix}.colnames.txt").write_text("\n".join(adata.var_names) + "\n")
    counts = sparse.csr_matrix(dense_counts(adata))
    io.mmwrite(f"{prefix}.counts.mtx", counts, field="integer")


@app.callback(invoke_without_command=True)
def main(
    # Input / Output options
    input: Annotated[str, INPUT_OUTPUT_OPTS["h5ad_file"]],
    out_prefix: Annotated[
        str,
        override(
            INPUT_OUTPUT_OPTS["out_prefix"],
            help=(
                "Prefix for the output file names: <prefix>.bm with "
                "[bold yellow]bm[/bold yellow], or <prefix>.counts.mtx, "
                "<prefix>.rownames.txt and <prefix>.colnames.txt with "
                "[bold yellow]mtx[/bold yellow]."
            ),
        ),
    ],
    out_file_format: Annotated[
        ExportFormat,
        typer.Option(
            "--outFileFormat",
            metavar="FORMAT",
            rich_help_panel="Export options",
            help=(
                "Output file format.\n\n"
                "[bold yellow]bm[/bold yellow]: BedGraphMatrix format, useful for "
                "single-cell visualization with pyGenomeTracks; stores data densely, "
                "so prefer exporting a region rather than the whole dataset.\n\n"
                "[bold yellow]mtx[/bold yellow]: MatrixMarket sparse format "
                "(<prefix>.counts.mtx plus <prefix>.rownames.txt and "
                "<prefix>.colnames.txt)."
            ),
        ),
    ],
    region: Annotated[str | None, INPUT_OUTPUT_OPTS["region"]] = None,
    # Other options
    verbose: Annotated[bool, OTHER_OPTS["verbose"]] = False,
    help: Annotated[bool, OTHER_OPTS["help"]] = False,
) -> int:
    if verbose:
        log_parameters(
            input=input,
            out_prefix=out_prefix,
            out_file_format=out_file_format,
            region=region,
        )
    else:
        warnings.filterwarnings("ignore")

    adata = validate_anndata(ad.read_h5ad(input), input)

    if region is not None:
        mask = region_mask(cast("pd.DataFrame", adata.var), region)
        if not mask.any():
            sys.stderr.write(f"No features found in region {region}.\n")
            raise typer.Exit(code=1)
        adata = adata[:, mask]

    if out_file_format is ExportFormat.bm:
        write_bedgraph_matrix(adata, f"{out_prefix}.bm")
    else:
        write_matrix_market(adata, out_prefix)

    return 0


def cli() -> None:
    configure_logging()
    preprocess_args()
    app()


if __name__ == "__main__":
    cli()
