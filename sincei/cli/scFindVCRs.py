#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys

from sincei import _sincei as internal

from . import ParserCommon
from . import _parsers as backend

DESCRIPTION = (
    "Call variable chromatin regions (VCRs) from binned chromatin data.\n\n"
    "``scFindVCRs`` calls variable chromatin regions (VCRs) from binned chromatin "
    "data. It takes a .h5ad file containing single-cell genomic signal in bins, and "
    "outputs BED files with genome segmentations for different sensitivities.\n"
    "First, a bin-to-bin correlation matrix is computed for each chromosome.\n"
    "Then, the correlation matrix is turned into a score map by convolving a number "
    "of square Gaussian kernels along its main diagonal. Each kernel has a sigma "
    "calculated using a maximum region size to consider. Each kernel produces a 1-D "
    "score for each bin, which are stacked into a matrix where each row corresponds "
    "to a kernel scale and each column to a bin.\n"
    "Finally, the PELT change-point detection algorithm is applied to the score map "
    "to identify regions with distinct correlation patterns. This step depends on a "
    "penalty parameter that controls the number of detected regions."
)

USAGE = "scFindVCRs -i counts.h5ad -bs 10000 -o vcrs.bed"

_VCR = "VCR options"
_DEFAULT_PENALTIES = [0.05, 0.1, 0.5]


def vcr_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_VCR)

    group.add_argument(
        "-bs",
        "--binsize",
        dest="binsize",
        metavar="INT",
        type=int,
        required=True,
        help="The size of the bins in the input AnnData object.",
    )
    group.add_argument(
        "-mr",
        "--maxRegionSize",
        dest="maxRegionSize",
        metavar="INT",
        type=int,
        default=None,
        help="The maximum region size to be considered, in base pairs. Larger regions "
        "may increase compute time. Defaults to 100 times the bin size.",
    )
    group.add_argument(
        "-nk",
        "--nKernels",
        dest="nKernels",
        metavar="INT",
        type=int,
        default=20,
        help="The number of kernels to use for the score map. More kernels generally "
        "lead to a better segmentation, but increase the computational cost.",
    )
    group.add_argument(
        "-pen",
        "--penalties",
        dest="penalties",
        metavar="FLOAT",
        type=float,
        action="extend",
        nargs="+",
        default=None,
        help="Penalty value(s) for change-point detection. Higher values result in "
        "fewer segments. Multiple values can be provided (separated by space); each "
        "produces a separate set of regions, distinguishable in the output BED file by "
        'filtering on the "score" column. (Default: 0.05, 0.1, 0.5)',
    )
    return parser


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options(
                ["h5ad_file", "out_file", "region"],
                defaults={
                    "out_file_metavar": ".bed",
                    "out_file_help": (
                        "Name of the output file (BED format) with the genome "
                        "segmentation result. The penalty threshold that defines each "
                        'segment is saved in the "score" column of the BED file, which '
                        "can be filtered to obtain non-overlapping segments."
                    ),
                },
            ),
            vcr_options(),
            ParserCommon.other_options(),
        ],
    )


def main(argv: list[str] | None = None) -> int:
    parser = parse_arguments()
    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)

    penalties = args.penalties or _DEFAULT_PENALTIES

    if args.verbose:
        backend.log_parameters(
            input=args.input,
            region=args.region,
            bin_size=args.binsize,
            max_region_size=args.maxRegionSize,
            n_kernels=args.nKernels,
            penalties=penalties,
            out_file=args.outFile,
            number_of_processors=args.numberOfProcessors,
            verbose=args.verbose,
        )

    result_path = internal.find_vcrs(
        args.input,
        args.binsize,
        args.outFile,
        max_region_size=args.maxRegionSize,
        n_kernels=args.nKernels,
        penalties=penalties,
        region=args.region,
        num_threads=args.numberOfProcessors,
    )
    print(result_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
