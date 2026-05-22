#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

import argparse
import sys
import logging

from . import ParserCommon

# from sincei import _sincei

DESCRIPTION = """
``scFindVCRs`` calls variable chromatin regions (VCRs) from binned chromatin data. It takes a
.h5ad file containing single-cell genomic signal in bins, and outputs BED files with genome
segmentations for different sensitivities.

First, a bin-to-bin correlation matrix is computed for each chromosome.

Then, the correlation matrix is turned into a score map by convolving a number of square
Gaussian kernels along its main diagonal. Each kernel has a sigma calculated using a maximum
region size to consider. Each kernel produces a 1-D score for each bin, which are stacked
into a matrix where each row corresponds to a kernel scale and each column to a bin.

Finally, the PELT change-point detection algorithm is applied to the score map to identify
regions with distinct correlation patterns. This step depends on a penalty parameter that
controls the number of detected regions.
"""

USAGE = "scFindVCRs -i binned_signal.h5ad -bs 2000 -mr 100000 -nk 20 -pen 0.05 0.1 0.5 -o detected_VCRs.bed"


def parse_arguments(args=None) -> argparse.ArgumentParser:
    io_args = ParserCommon.inputOutputOptions(opts=["h5adfile", "region"], requiredOpts=["h5adfile"])
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        parents=[io_args, get_args(), other_args],
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

    vcr_options = parser.add_argument_group("VCR detection options")

    vcr_options.add_argument(
        "-bs",
        "--binSize",
        metavar="INT",
        type=int,
        help="The size of the bins in the input Anndata object.",
        required=True,
    )

    vcr_options.add_argument(
        "-mr",
        "--maxRegionSize",
        metavar="INT",
        type=int,
        help="""
        The maximum region size to be considered, in base pairs. Larger regions may increase compute time.
        Defaults to 100 times the bin size.
        """,
        default=None,
    )

    vcr_options.add_argument(
        "-nk",
        "--nKernels",
        metavar="INT",
        type=int,
        help="""
        The number of kernels to use for the score map. More kernels generally lead to a better segmentation,
        but increase the computational cost.""",
        default=20,
    )

    vcr_options.add_argument(
        "-pen",
        "--penalties",
        metavar="INT",
        nargs="+",
        type=float,
        help="""
        Penalty value for change-point detection. Higher values result in fewer segments. Multiple values
        can be provided (separated by space). Each penalty value will produce a separate set of regions within
        which can be seperated from the output BED file by filtering on the "score" column.
        """,
        default=[0.05, 0.1, 0.5],
    )

    vcr_options.add_argument(
        "-o",
        "--outFile",
        metavar=".bed",
        type=str,
        help="""
        Name of the output file (BED format) with genome segmentation result. The penalty threshold that defines
        each segment is saved in the "score" column of the BED file, and the BED file can be filtered based on this
        column to obtain non-overlapping segments.
        """,
        required=True,
    )

    return parser


def main(args: list[str] | None = None) -> int:
    args = parse_arguments().parse_args(args)

    for arg in args.__dict__:
        logging.debug(f"{arg}: {args.__dict__[arg]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
