#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys

from . import ParserCommon
from . import _parsers as backend

DESCRIPTION = (
    "Cluster cells from a cell-by-feature matrix using the Leiden "
    "algorithm.\n\n``scClusterCells`` clusters cells based on the "
    "dimensionality reduction in the input h5ad file. The result is an "
    "updated h5ad object, and (optionally) a plot and a .tsv file with UMAP "
    "coordinates and corresponding cluster id for each cell."
)

USAGE = "scClusterCells -i reduced.h5ad -o clustered.h5ad"


_CLUSTERING = "Clustering options"


def clustering_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_CLUSTERING)
    group.add_argument(
        "-ou",
        "--outFileUMAP",
        dest="outFileUMAP",
        default=None,
        help="The output plot file (for UMAP). If specified, a 4-column .tsv file with "
        "the same prefix is also created with the cell IDs, raw UMAP coordinates "
        "(UMAP1 and UMAP2) and Leiden cluster number.",
    )
    group.add_argument(
        "-cr",
        "--clusterResolution",
        dest="clusterResolution",
        type=float,
        default=1.0,
        help="Resolution parameter for Leiden clustering. Values lower than 1.0 result "
        "in fewer clusters, while higher values lead to splitting of clusters. In most "
        "cases the optimum is between 0.8 and 1.2.",
    )
    group.add_argument(
        "--dimRed",
        dest="dimRed",
        metavar="METHOD",
        choices=ParserCommon.DIM_RED_METHODS,
        default=None,
        help="Dimensionality reduction modality to perform Leiden clustering on. If "
        "not given, the program searches the ``.obsm`` field of the input h5ad for the "
        "output of ``scReduceDims`` in order of preference: LSA, LDA, logPCA, glmPCA.",
    )
    return parser


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options(["h5ad_file", "out_file"]),
            clustering_options(),
            ParserCommon.plot_options(),
            ParserCommon.other_options(),
        ],
    )


def main(argv: list[str] | None = None) -> int:
    parser = parse_arguments()
    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)

    backend.log_parameters(
        input=args.input,
        out_file=args.outFile,
        out_file_umap=args.outFileUMAP,
        cluster_resolution=args.clusterResolution,
        dim_red=args.dimRed,
        plot_width=args.plotWidth,
        plot_height=args.plotHeight,
        plot_file_format=args.plotFileFormat,
        number_of_processors=args.numberOfProcessors,
        verbose=args.verbose,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
