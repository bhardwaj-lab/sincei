#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import sys
import os

## own functions
# scriptdir=os.path.abspath(os.path.join(__file__, "../../sincei"))
# sys.path.append(scriptdir)
from sincei._version import __version__


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
sincei is a suite of command-line tools developed for a user-friendly analysis of single-cell sequencing data.
Version: {}

Each tool begins with the prefix sc<tool_name>, such as:

 $ scBulkCoverage -b file1.bam -g groupinfo.txt -o coverage

[ Tools for a typical single-cell analysis workflow ] (WIP: work in progress/not available yet)

    scFilterBarcodes        Identify and filter cell barcodes from BAM file (for droplet-based single-cell seq)
    scFilterStats           Produce per-cell statistics after filtering reads by user-defined criteria.
    scCountReads            Counts reads for each barcode on genomic bins or user-defined features.
    scCountQC               Perform quality control and filter the output of scCountReads.
    scCombineCounts         [WIP] Concatenate/merge the counts from different samples/batches or modalities
    scClusterCells          Perform dimensionality reduction and clustering on the output of scCountReads.
    scBulkCoverage          Get pseudo-bulk coverage per group using a user-supplied cell->group mapping (output of scClusterCells).
    scBAMops                Modify a BAM file to group cells (using cell barcodes), or filter/shift mapped reads.
    scFindMarkers           [WIP] Find marker genes per group, given the output of scCountReads and a user-defined group.
    scFeaturePlot           [WIP] Plot the counts for a given feature on a UMAP or on a (IGV-style) genomic-track.

""".format(
            __version__
        ),
    )

    return parser


def process_args(args=None):
    args = parse_arguments().parse_args(args)

    return args


def main(args=None):
    if args is None and len(sys.argv) == 1:
        args = ["--help"]
    process_args(args)
