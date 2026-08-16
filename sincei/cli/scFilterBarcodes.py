#!/usr/bin/env python
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

from sincei import _sincei as internal

from . import ParserCommon
from . import _parsers as backend

DESCRIPTION = (
    "Filter cell barcodes from a BAM file (for droplet-based single-cell seq).\n\n"
    "``scFilterBarcodes`` identifies barcodes present in a BAM file and produces a "
    "list. You can optionally filter these barcodes by matching them to a whitelist or "
    "based on total counts. This tool expects single experiment BAM files, not merged "
    "files."
)

USAGE = "scFilterBarcodes -b sample.bam -w whitelist.txt -o barcodes.tsv"

_BARCODE = "Barcode options"


def barcode_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_BARCODE)

    group.add_argument(
        "-d",
        "--minHammingDist",
        dest="minHammingDist",
        metavar="INT",
        type=int,
        default=0,
        help="Minimum hamming distance to match the barcode in whitelist. Note that "
        "increasing the hamming distance really slows down the barcode detection "
        "process.",
    )
    group.add_argument(
        "-mc",
        "--minCount",
        dest="minCount",
        metavar="INT",
        type=int,
        default=0,
        help="Minimum number of bins with non-zero counts in order to report a "
        "barcode. Note that this number ranges from 0 to genome size / bin size.",
    )
    group.add_argument(
        "-rp",
        "--rankPlot",
        dest="rankPlot",
        metavar="STR",
        default=None,
        help="The output file name to plot the ranked counts per barcode (similar to "
        'the "knee plot", but counts are the number of non-zero bins in this case).',
    )
    group.add_argument(
        "-mq",
        "--minMappingQuality",
        dest="minMappingQuality",
        metavar="INT",
        type=int,
        default=None,
        help="If set, only reads that have a mapping quality score of at least this "
        "are considered.",
    )
    return parser


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options([
                "bam_file",
                "out_file",
                "whitelist",
                "region",
            ]),
            barcode_options(),
            ParserCommon.bam_options(
                [
                    "cell_tag",
                    "group_tag",
                    "labels",
                    "smart_labels",
                    "blacklist",
                    "chr_to_skip",
                    "bin_size",
                    "distance_between_bins",
                ],
                defaults={"bin_size": 100000, "distance_between_bins": None},
            ),
            ParserCommon.filter_options([
                "duplicate_filter",
                "motif_filter",
                "genome_2bit",
                "gc_content_filter",
                "min_aligned_fraction",
            ]),
            ParserCommon.read_options([
                "min_fragment_length",
                "max_fragment_length",
                "filter_rna_strand",
                "extend_reads",
                "center_reads",
            ]),
            ParserCommon.other_options(),
        ],
    )


def main(argv: list[str] | None = None) -> int:
    parser = parse_arguments()
    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)

    if args.verbose:
        backend.log_parameters(
            bam_file=args.bamfile, whitelist=args.whitelist, out_file=args.outFile
        )

    # Options exposed for parity but not honored by the barcode-detection backend.
    backend.warn_unsupported(
        region=args.region,
        group_tag=args.groupTag,
        labels=args.labels,
        smart_labels=args.smartLabels,
        distance_between_bins=args.distanceBetweenBins,
        duplicate_filter=args.duplicateFilter,
        motif_filter=args.motifFilter,
        genome_2bit=args.genome2bit,
        gc_content_filter=args.GCcontentFilter,
        min_aligned_fraction=args.minAlignedFraction,
        min_fragment_length=backend.optional_length(args.minFragmentLength),
        max_fragment_length=backend.optional_length(args.maxFragmentLength),
        filter_rna_strand=args.filterRNAstrand,
        extend_reads=args.extendReads,
        center_reads=args.centerReads,
    )

    barcode_counts = internal.filter_barcodes(
        args.bamfile,
        whitelist=backend.read_barcodes(args.whitelist) if args.whitelist else None,
        blacklist_file_name=backend.first_blacklist(args.blacklist),
        cell_tag=args.cellTag,
        min_hamming_dist=args.minHammingDist,
        min_mapping_quality=args.minMappingQuality,
        bin_size=args.binSize,
        chr_to_skip=args.chrToSkip or [],
        num_threads=args.numberOfProcessors,
    )

    # A TSV of every detected barcode, sorted by descending count for a
    # knee-plot-friendly ordering: `count` is the number of non-zero bins the
    # barcode was seen in, `selected` is True when count >= minCount, and
    # `count_log10` and `count_rank` are the knee-plot axes.
    #
    # Two barcodes with the same count share the same rank, and the next rank is
    # incremented by the number of barcodes with that count.
    rank_of_count: dict[int, int] = {}
    for position, (_, count) in enumerate(barcode_counts):
        rank_of_count.setdefault(count, position + 1)

    with Path(args.outFile).open("w") as out:
        out.write("barcode\tcount\tselected\tcount_log10\tcount_rank\n")
        for barcode, count in barcode_counts:
            selected = count >= args.minCount
            # The backend only reports barcodes it saw, so count >= 1 and
            # log10 is always defined.
            out.write(
                f"{barcode}\t{count}\t{selected}\t{math.log10(count)}\t"
                f"{rank_of_count[count]}\n"
            )

    if args.rankPlot:
        # Imported lazily: it pulls in matplotlib, which would otherwise be paid
        # for on every invocation of this command.
        from sincei.plotting._barcode_rank import plot_barcode_rank  # noqa: PLC0415

        plot_barcode_rank(
            barcode_counts, args.rankPlot, min_count=args.minCount or None
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
