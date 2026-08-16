#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys

from sincei import _sincei as internal

from . import ParserCommon
from . import _parsers as backend

DESCRIPTION = (
    "Get pseudo-bulk coverage per group using a cell->group mapping (output of "
    "scClusterCells).\n\n"
    "``scBulkCoverage`` takes alignments of reads or fragments as input (BAM files), "
    "along with cell grouping information, such as barcode -> batch, or barcode -> "
    "cluster, as tsv file, and generates a coverage track (bigWig or bedGraph) per "
    "group as output. The coverage is calculated as the number of reads/fragments per "
    "bin, where bins are short consecutive counting windows of a defined size. It is "
    "possible to extended/change the length of the reads to better reflect the actual "
    "fragment length. ``scBulkCoverage`` offers normalization per cluster using "
    "different methods."
)

USAGE = "scBulkCoverage -b sample.bam -gi groupinfo.tsv -o coverage"

_COVERAGE = "Coverage options"


def coverage_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_COVERAGE)

    group.add_argument(
        "-of",
        "--outFileFormat",
        dest="outFileFormat",
        metavar="FORMAT",
        choices=ParserCommon.OUT_FILE_FORMATS,
        default="bigwig",
        help="Output file type.",
    )
    group.add_argument(
        "-n",
        "--normalizeUsing",
        dest="normalizeUsing",
        metavar="METHOD",
        choices=ParserCommon.NORMALIZE_USING,
        default="CPM",
        help="How to normalize the pseudo-bulk counts.\n\n"
        "CPM: normalize each bin to counts per million mapped reads in that group.\n\n"
        "RPKM: reads per kilobase per million mapped reads in that group.\n\n"
        "Frequency: binarize coverage per bin and normalize to total cells per "
        "group.\n\n"
        "Mean: get mean signal per bin across cells in each group.\n\n"
        "None: simply return the sum of coverage per group.",
    )
    group.add_argument(
        "-ig",
        "--ignoreForNormalization",
        dest="ignoreForNormalization",
        metavar="CHR",
        action="extend",
        nargs="+",
        default=None,
        help="Chromosomes to skip while calculating normalization factors.",
    )
    group.add_argument(
        "-nr",
        "--normalizeByReference",
        dest="normalizeByReference",
        default=None,
        help="NOT IMPLEMENTED: Normalize each group of cells by a reference group "
        "(which must be present in the --groupInfo file). Note that the "
        "--normalizeUsing method is applied beforehand.",
    )
    group.add_argument(
        "--scaleFactor",
        dest="scaleFactor",
        type=float,
        default=1.0,
        help="The computed scaling factor (or 1, if not applicable) will be multiplied "
        "by this.",
    )
    group.add_argument(
        "--mnase",
        action="store_true",
        help="Determine nucleosome positions from MNase-seq/CUTnRUN data. Only 3 "
        "nucleotides at the center of each fragment are counted. Only fragment lengths "
        "between 130 - 200 bp are considered to avoid dinucleosomes or other "
        "artifacts. *NOTE*: Requires paired-end data. A bin size of 1 is recommended.",
    )
    group.add_argument(
        "--offset",
        action="extend",
        nargs="+",
        type=int,
        default=None,
        help="Uses this offset inside of each read as the signal. This is useful in "
        "cases like RiboSeq or GROseq. Negative values indicate offsets from the end "
        "of each read. A value of 1 indicates the first base of the alignment; -1 is "
        "the last base. An offset of 0 is not permitted. If two values are specified, "
        "they define a range of positions.",
    )
    return parser


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options([
                "bam_files",
                "group_info",
                "out_prefix",
                "region",
            ]),
            coverage_options(),
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
                defaults={"bin_size": 100, "distance_between_bins": None},
            ),
            ParserCommon.filter_options([
                "duplicate_filter",
                "motif_filter",
                "genome_2bit",
                "gc_content_filter",
                "min_aligned_fraction",
            ]),
            ParserCommon.read_options([
                "min_mapping_quality",
                "sam_flag_include",
                "sam_flag_exclude",
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
            bam_files=args.bamfiles,
            group_info=args.groupInfo,
            out_prefix=args.outFilePrefix,
        )

    # Options exposed for parity but not honored by the coverage backend.
    backend.warn_unsupported(
        normalize_by_reference=args.normalizeByReference,
        group_tag=args.groupTag,
    )

    min_gc, max_gc = backend.parse_gc_content(args.GCcontentFilter)
    bam_labels = backend.resolve_labels(args.bamfiles, args.labels, args.smartLabels)
    # stepSize = binSize + gap; a zero gap yields contiguous bins.
    step_size = args.binSize + (args.distanceBetweenBins or 0)

    output_files = internal.bulk_coverage(
        args.bamfiles,
        bam_labels,
        args.groupInfo,
        args.outFilePrefix,
        bin_size=args.binSize,
        step_size=step_size,
        bc_tag=args.cellTag,
        umi_tag=None,
        region=args.region,
        min_mapq=args.minMappingQuality,
        sam_flag_include=args.samFlagInclude,
        sam_flag_exclude=args.samFlagExclude,
        chr_to_skip=args.chrToSkip or [],
        ignore_for_normalization=args.ignoreForNormalization or [],
        blacklist_path=backend.first_blacklist(args.blacklist),
        extend_reads=args.extendReads,
        center_reads=args.centerReads,
        dup_method=backend.dup_method(args.duplicateFilter),
        genome_2bit=args.genome2bit,
        motif_filter=backend.parse_motif_filter(args.motifFilter),
        min_gc=min_gc,
        max_gc=max_gc,
        min_aligned_fraction=args.minAlignedFraction,
        min_fragment_length=backend.optional_length(args.minFragmentLength),
        max_fragment_length=backend.optional_length(args.maxFragmentLength),
        normalize_using=args.normalizeUsing,
        scale_factor=args.scaleFactor,
        out_format=args.outFileFormat,
        mnase=args.mnase,
        offset=args.offset or None,
        filter_rna_strand=args.filterRNAstrand,
        num_threads=args.numberOfProcessors,
    )
    for path in output_files:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
