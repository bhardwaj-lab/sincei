#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from importlib import import_module

from sincei import _sincei as internal

# The order here is the order the subcommands appear in the help message.
_TOOLS = [
    (
        "scFilterBarcodes",
        "Filter cell barcodes from a BAM file (for droplet-based single-cell seq).",
    ),
    (
        "scFilterStats",
        "Produce per-cell statistics after filtering reads by user-defined criteria.",
    ),
    (
        "scJSD",
        "Compare read coverages on sampled regions using the Jensen-Shannon distance.",
    ),
    (
        "scCountReads",
        "Counts reads for each cell barcode on genomic bins or user-defined features.",
    ),
    (
        "scCountQC",
        (
            "Perform quality control and filter cells and regions from a "
            "cell-by-feature matrix."
        ),
    ),
    (
        "scFindVCRs",
        "Call variable chromatin regions (VCRs) from binned chromatin data.",
    ),
    (
        "scReduceDims",
        (
            "Perform dimensionality reduction and UMAP projection on a "
            "cell-by-feature matrix."
        ),
    ),
    (
        "scClusterCells",
        "Cluster cells from a cell-by-feature matrix using the Leiden algorithm.",
    ),
    (
        "scScoreFeatures",
        "Aggregate a binned chromatin count matrix into per-feature scores.",
    ),
    ("scCombineSamples", "Concatenate/merge AnnData files from different samples."),
    (
        "scCombineMods",
        "Merge AnnData files from different modalities into a MuData object.",
    ),
    ("scPlotRegion", "Plot pseudo-bulk and per cell coverage for a genomic region."),
    (
        "scBulkCoverage",
        (
            "Get pseudo-bulk coverage per group using a cell->group mapping "
            "(output of scClusterCells)."
        ),
    ),
    ("scExportSignal", "Export .h5ad objects to other formats."),
]

DESCRIPTION = f"""
sincei is a suite of command-line tools to explore single-cell genomics data.
Every tool name begins with the prefix `sc` followed by <tool_name>, such as:
    scBulkCoverage -b reads.bam -g groupinfo.txt -o coverage

    sincei {internal.version()}
"""

USAGE = "sincei <tool_name> [options]"


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
    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=f"sincei {internal.version()}",
        help="Print the program version and exit.",
    )
    subparsers = parser.add_subparsers(dest="tool", metavar="<tool_name>")
    for name, summary in _TOOLS:
        subparsers.add_parser(name, help=summary, add_help=False)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Dispatch to a subcommand's own ``main``.

    The subcommand modules are imported lazily so that ``sincei <tool> --help``
    only pays for the tool it is about to run, and ``sincei --help`` pays for
    none of them.  ``ParserCommon`` is imported eagerly and is what the root
    help needs.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = parse_arguments()

    if not argv or argv[0] in ("-h", "--help"):
        parser.print_help()
        return 0
    if argv[0] in ("-V", "--version"):
        parser.parse_args(argv)
        return 0

    known = {name for name, _ in _TOOLS}
    if argv[0] not in known:
        parser.parse_args(argv)  # argparse reports the unknown subcommand
        return 2

    module = import_module(f"sincei.cli.{argv[0]}")
    sys.argv = [argv[0], *argv[1:]]
    return module.main(argv[1:] or None)


if __name__ == "__main__":
    raise SystemExit(main())
