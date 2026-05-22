#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys

import sincei._sincei as internal

DESCRIPTION = f"""
sincei is a suite of command-line tools to explore single-cell genomics data.
Version: {internal.version()}

Every tool name begins with the prefix `sc` followed by <tool_name>, such as:

    scBulkCoverage -b reads.bam -g groupinfo.txt -o coverage

[ List of tools ]

    scCallBarcodes          Identify cell barcode whitelist de-novo.
    scFilterBarcodes        Filter cell barcodes from BAM file (for droplet-based single-cell seq)
    scFilterStats           Produce per-cell statistics after filtering reads by user-defined criteria.
    scCountReads            Counts reads for each cell barcode on genomic bins or user-defined features.
    scCountQC               Perform quality control and filter cells and regions from a cell-by-feature matrix.
    scFindVCRs              Call variable chromatin regions (VCRs) from binned chromatin data.
    scScoreFeatures         Calculate gene activity scores from chromatin features/bins.
    scCombineCounts         Concatenate/merge AnnData files from different samples/batches or modalities.
    scReduceDims            Perform dimensionality reduction on a cell-by-feature matrix.
    scClusterCells          Cluster cells from a cell-by-feature matrix using the Leiden algorithm.
    scPlotRegion            Plot pseudo-bulk and per cell coverage for a genomic region.
    scBulkCoverage          Get pseudo-bulk coverage per group using a cell->group mapping (output of scClusterCells).
"""


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=DESCRIPTION,
        add_help=False,
    )

    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this message and exit",
    )

    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"%(prog)s {internal.version()}",
        help="Print the program version and exit.",
    )

    # If no arguments are provided, show help and exit
    if args is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    return parser


def main(args=None):
    args = parse_arguments().parse_args(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
