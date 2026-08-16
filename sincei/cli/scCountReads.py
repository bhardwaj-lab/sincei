#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys

from sincei import _sincei as internal

from . import ParserCommon
from . import _parsers as backend

DESCRIPTION = (
    "Counts reads for each cell barcode on genomic bins or user-defined features.\n\n"
    "``scCountReads`` computes the read coverages per cell barcode for genomic regions "
    "in the provided BAM file(s). The analysis can be performed for the entire genome "
    "by running the program in ``bins`` mode. If you want to count the read coverage "
    "for specific regions only, use the ``features`` mode instead. The standard output "
    'of ``scCountReads`` is a ".h5ad" file with counts, along with rowName (features) '
    "and colNames (cell barcodes)."
)

USAGE = "scCountReads {bins,features} -b sample.bam -bc barcodes.txt -o counts.h5ad"

_COUNTING = "Counting options"


def counting_options() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(_COUNTING)

    group.add_argument(
        "--valuetag",
        dest="valuetag",
        metavar="XX",
        default=None,
        help='Instead of counting each read/fragment as "1", add the values from a '
        "given BAM tag to the count matrix. For example, this can be used to count the "
        "number of methylated CpG per read.",
    )
    group.add_argument(
        "--genomeChunkSize",
        dest="genomeChunkSize",
        metavar="INT",
        type=int,
        default=None,
        help="Manually specify the size of the genome provided to each processor. The "
        "default (None) determines this from the read density of the BAM file.",
    )
    return parser


def _shared_parents() -> list[argparse.ArgumentParser]:
    """Option groups common to both the `bins` and `features` subcommands."""
    return [
        ParserCommon.bam_options(
            [
                "labels",
                "smart_labels",
                "blacklist",
                "chr_to_skip",
                "bin_size",
                "distance_between_bins",
                "cell_tag",
                "group_tag",
            ],
            defaults={"bin_size": 10000, "distance_between_bins": 0},
        ),
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
        ParserCommon.filter_options([
            "duplicate_filter",
            "motif_filter",
            "genome_2bit",
            "gc_content_filter",
            "min_aligned_fraction",
        ]),
        counting_options(),
    ]


def parse_arguments() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=DESCRIPTION,
        usage=USAGE,
        add_help=False,
    )
    parser.add_argument(
        "-h", "--help", action="help", help="Show this message and exit."
    )
    subparsers = parser.add_subparsers(dest="mode", metavar="{bins,features}")

    bins = subparsers.add_parser(
        "bins",
        parents=[
            ParserCommon.input_output_options([
                "bam_files",
                "barcodes",
                "out_file",
                "region",
                "compression",
            ]),
            *_shared_parents(),
            ParserCommon.other_options(),
        ],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Count reads in fixed-size bins.",
        help="Count reads in fixed-size bins.",
        add_help=False,
        conflict_handler="resolve",
    )
    bins.set_defaults(mode="bins", bed=None)

    features = subparsers.add_parser(
        "features",
        parents=[
            ParserCommon.input_output_options([
                "bam_files",
                "barcodes",
                "out_file",
                "bed",
                "region",
                "compression",
            ]),
            *_shared_parents(),
            ParserCommon.gtf_gff_options([
                "transcript_id",
                "exon_id",
                "transcript_id_tag",
                "metagene",
            ]),
            ParserCommon.other_options(),
        ],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Count reads on provided features.",
        help="Count reads on provided features.",
        add_help=False,
        conflict_handler="resolve",
    )
    features.set_defaults(mode="features")

    return parser


def _count_reads(args: argparse.Namespace) -> int:
    if args.verbose:
        backend.log_parameters(
            mode=args.mode, bam_files=args.bamfiles, out_file=args.outFile
        )

    # Options accepted by the CLI for parity with deepTools but not (yet) implemented
    backend.warn_unsupported(
        group_tag=args.groupTag,
        labels=args.labels,
        smart_labels=args.smartLabels,
    )

    min_gc, max_gc = backend.parse_gc_content(args.GCcontentFilter)

    shared = {
        "barcodes": backend.read_barcodes(args.barcodes),
        "output_path": args.outFile,
        "bc_tag": args.cellTag,
        "umi_tag": None,
        "count_tag": args.valuetag,
        "min_mapq": args.minMappingQuality,
        "sam_flag_include": args.samFlagInclude,
        "sam_flag_exclude": args.samFlagExclude,
        "chr_to_skip": args.chrToSkip or [],
        "region": args.region,
        "blacklist_path": backend.first_blacklist(args.blacklist),
        "extend_reads": args.extendReads,
        "center_reads": args.centerReads,
        "dup_method": backend.dup_method(args.duplicateFilter),
        "filter_rna_strand": args.filterRNAstrand,
        "genome_2bit": args.genome2bit,
        "motif_filter": backend.parse_motif_filter(args.motifFilter),
        "min_gc": min_gc,
        "max_gc": max_gc,
        "min_aligned_fraction": args.minAlignedFraction,
        "min_fragment_length": backend.optional_length(args.minFragmentLength),
        "max_fragment_length": backend.optional_length(args.maxFragmentLength),
        "compression": args.compression,
        "compression_level": args.compressionLevel,
        "num_threads": args.numberOfProcessors,
    }
    if args.genomeChunkSize:
        shared["chunk_size"] = args.genomeChunkSize

    if args.mode == "bins":
        # stepSize = binSize + gap; a zero gap yields contiguous bins.
        step_size = args.binSize + (args.distanceBetweenBins or 0)
        internal.count_bins(
            args.bamfiles, bin_size=args.binSize, step_size=step_size, **shared
        )
    else:
        # GTF/GFF field selection (ignored for BED inputs). `transcriptIDtag`
        # is None by default so the backend picks the per-format name attribute.
        internal.count_features(
            args.bamfiles,
            args.bed,
            feature_type=args.transcriptID,
            exon_type=args.exonID,
            name_attr=args.transcriptIDtag,
            metagene=args.metagene,
            **shared,
        )

    return 0


def main(argv: list[str] | None = None) -> int:
    parser = parse_arguments()
    if argv is None and len(sys.argv) == 1:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)
    if not getattr(args, "mode", None):
        parser.print_help()
        return 0
    return _count_reads(args)


if __name__ == "__main__":
    raise SystemExit(main())
