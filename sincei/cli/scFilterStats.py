#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

from sincei import _sincei as internal

from . import ParserCommon
from . import _parsers as backend

if TYPE_CHECKING:
    import argparse
    from collections.abc import Sequence

# Column order produced by the Rust `filter_stats` backend (BarcodeStat::to_vec).
_STAT_COLUMNS = [
    "Total_sampled",
    "Filtered",
    "Blacklisted",
    "Low_MAPQ",
    "Missing_Flags",
    "Excluded_Flags",
    "Internal_Duplicates",
    "Marked_Duplicates",
    "Singletons",
    "Wrong_strand",
    "Wrong_motif",
    "Unwanted_GC_content",
    "Low_aligned_fraction",
]


def _format_row(row: Sequence[float]) -> list[str]:
    """Format one backend stat row for the TSV.

    Mirrors the original ``sincei``: ``total`` is the absolute number of sampled
    reads, while every other column is reported as a percentage of ``total``
    (``count / total * 100``).  When ``total`` is 0 the percentage is undefined,
    so those cells are left blank (matching the original's NaN -> empty output).
    """
    total = row[0]
    return [
        str(total),
        *("" if total == 0 else str(value / total * 100) for value in row[1:]),
    ]


DESCRIPTION = (
    "Produce per-cell statistics after filtering reads by user-defined criteria.\n\n"
    "``scFilterStats`` estimates the number of reads that would be filtered given a "
    "set of criteria and prints it to the terminal. Furthermore, it tracks the number "
    "of singleton reads. The following metrics will always be tracked regardless of "
    "what you specify (the order output also matches this):\n"
    "* Total sampleed reads (including unmapped)\n"
    "* Mapped reads\n"
    "* Reads in blacklisted regions (--blacklist)\n\n"
    "The following metrics are estimated according to the --binSize and "
    "--distanceBetweenBins parameters:\n"
    "* Estimated mapped reads filtered (the total number of mapped reads filtered for "
    "any reason)\n"
    "* Alignments with a below threshold MAPQ (--minMappingQuality)\n"
    "* Alignments with at least one missing flag (--samFlagInclude)\n"
    "* Alignments with undesirable flags (--samFlagExclude)\n"
    "* Duplicates determined by sincei (--duplicateFilter)\n"
    "* Duplicates marked externally (e.g., by picard)\n"
    "* Singletons (paired-end reads with only one mate aligning)\n"
    "* Wrong strand (due to --filterRNAstrand)\n\n"
    "The sum of these may be more than the total number of reads. Note that alignments "
    "are sampled from bins of size --binSize spaced --distanceBetweenBins apart."
)

USAGE = "scFilterStats -b sample.bam -bc barcodes.txt -o stats.tsv"


def parse_arguments() -> argparse.ArgumentParser:
    return ParserCommon.build_parser(
        DESCRIPTION,
        USAGE,
        [
            ParserCommon.input_output_options(["bam_files", "barcodes", "out_file"]),
            ParserCommon.bam_options(
                [
                    "cell_tag",
                    "bin_size",
                    "distance_between_bins",
                    "labels",
                    "smart_labels",
                    "blacklist",
                    "chr_to_skip",
                ],
                defaults={"bin_size": 100_000, "distance_between_bins": 1_000_000},
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
                "filter_rna_strand",
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
            bam_files=args.bamfiles, barcodes=args.barcodes, out_file=args.outFile
        )

    min_gc, max_gc = backend.parse_gc_content(args.GCcontentFilter)
    barcode_list = backend.read_barcodes(args.barcodes)
    sample_labels = backend.resolve_labels(args.bamfiles, args.labels, args.smartLabels)

    kwargs = {
        "barcodes": barcode_list,
        "bc_tag": args.cellTag,
        "umi_tag": None,
        "bin_size": args.binSize,
        "distance_between_bins": args.distanceBetweenBins or 0,
        "min_mapq": args.minMappingQuality,
        "sam_flag_include": args.samFlagInclude,
        "sam_flag_exclude": args.samFlagExclude,
        "filter_rna_strand": args.filterRNAstrand,
        "chr_to_skip": args.chrToSkip or [],
        "blacklist_path": backend.first_blacklist(args.blacklist),
        "dup_method": backend.dup_method(args.duplicateFilter),
        "genome_2bit": args.genome2bit,
        "motif_filter": backend.parse_motif_filter(args.motifFilter),
        "min_gc": min_gc,
        "max_gc": max_gc,
        "min_aligned_fraction": args.minAlignedFraction,
        "num_threads": args.numberOfProcessors,
    }

    # The backend processes one BAM at a time; aggregate the per-sample rows
    # into a single TSV keyed by `Cell_ID` (`{sample}::{barcode}`, the same
    # identifier the counting tools use for AnnData's obs_names).
    with Path(args.outFile).open("w") as out:
        out.write("\t".join(["Cell_ID", *_STAT_COLUMNS]) + "\n")
        for bam, label in zip(args.bamfiles, sample_labels, strict=True):
            result_barcodes, stat_rows = internal.filter_stats(bam, **kwargs)
            out.writelines(
                "\t".join([f"{label}::{barcode}", *_format_row(row)]) + "\n"
                for barcode, row in zip(result_barcodes, stat_rows, strict=True)
            )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
