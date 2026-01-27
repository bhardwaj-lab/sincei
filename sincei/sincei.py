#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
from importlib import metadata


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"""
sincei is a suite of command-line tools developed for a user-friendly analysis of single-cell sequencing data.
Version: {metadata.version("sincei")}

Each tool begins with the prefix sc<tool_name>, such as:

 $ scBulkCoverage -b file1.bam -g groupinfo.txt -o coverage

[ Tools for a typical single-cell analysis workflow ]

    scFilterBarcodes        Identify and filter cell barcodes from BAM file (for droplet-based single-cell seq)
    scFilterStats           Produce per-cell statistics after filtering reads by user-defined criteria.
    scFindVCRs              Call variable chromatin regions (VCRs) from binned chromatin data.
    scScoreFeatures         Calculate gene activity scores from chromatin features/bins.
    scCountReads            Counts reads for each barcode on genomic bins or user-defined features.
    scCountQC               Perform quality control and filter the output of scCountReads.
    scCombineCounts         Concatenate/merge the counts from different samples/batches or modalities.
    scClusterCells          Perform dimensionality reduction and clustering on the output of scCountReads.
    scBulkCoverage          Get pseudo-bulk coverage per group using a user-supplied cell->group mapping (output of scClusterCells).
""",
    )

    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    return args


def main(args=None):
    if args is None and len(sys.argv) == 1:
        args = ["--help"]
    process_args(args)
