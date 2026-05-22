#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

DESCRIPTION = """
``scClusterCells`` clusters cells based on the dimensionality reduction in the input h5ad file.
The result is an updated h5ad object, and (optionally) a plot file and a tsv file with UMAP coordinates
and corresponding cluster id for each cell.
"""

USAGE = "scClusterCells -i cellCounts.h5ad -o clustered.h5ad -op umap.png"


def parse_arguments(args: list[str] | None = None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(opts=["h5adfile", "outFile"], requiredOpts=["outFile"])
    plot_args = ParserCommon.plotOptions()
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, get_args(), plot_args, other_args],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=DESCRIPTION,
        usage=USAGE,
        add_help=False,
    )

    # If no arguments are provided, show help and exit
    if args is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    return parser


def get_args() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)

    general = parser.add_argument_group("Clustering Options")
    general.add_argument(
        "-op",
        "--outFileUMAP",
        metavar="STR",
        type=str,
        required=False,
        help="The output plot file (for UMAP). If you specify this option, a 4-column .tsv file with the same prefix"
        "is also created with the cell IDs, raw UMAP coordinates (UMAP1 and UMAP2) and Leiden cluster number.",
    )

    general.add_argument(
        "-cr",
        "--clusterResolution",
        metavar="FLOAT",
        default=1.0,
        type=float,
        help="Resolution parameter for Leiden clustering. Values lower than 1.0 would result in less clusters, "
        "while higher values lead to splitting of clusters. In most cases, the optimum value would be between "
        "0.8 and 1.2. (Default: %(default)s)",
    )

    general.add_argument(
        "--dimRed",
        metavar="STR",
        type=str,
        help="Dimensionality reduction modality to perform Leiden clustering on. If not given, the program searches "
        "in the ``.obsm`` field of the input h5ad for the output of ``scReduceDims`` in order of preference: "
        "LSA, LDA, logPCA, glmPCA.",
    )

    return parser


def main(args: list[str] | None = None) -> int:
    args = parse_arguments().parse_args(args)

    for arg in args.__dict__:
        logging.debug(f"{arg}: {args.__dict__[arg]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
