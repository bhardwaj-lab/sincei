#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import anndata as ad
from scipy import sparse, io
import numpy as np
import pandas as pd
from sincei import ParserCommon

# ignore anndata .var modification warning
import warnings
from anndata import ImplicitModificationWarning

warnings.simplefilter(action="ignore", category=ImplicitModificationWarning)


# Helper function to convert chromosome labels to a numeric order
def chromosome_to_numeric(chrom):
    chrom = chrom.upper().strip("CHR")
    try:
        return int(chrom)
    except ValueError:
        # Handle special cases ("X", "Y", "MT") by giving them a higher order value
        if chrom == "X":
            return 23
        elif chrom == "Y":
            return 24
        elif chrom == "MT" or chrom == "M":
            return 25
        return 100  # for any unexpected case, assign a high number


def get_args():
    parser = argparse.ArgumentParser(add_help=False)

    group = parser.add_argument_group("Common Options")

    group.add_argument(
        "--outFileFormat",
        type=str,
        default=None,
        choices=["bm", "mtx"],
        help="Output file format for `scExportSignal`. "
        '"bm" refers to the BedGraphMatrix format, useful for single-cell data visualization '
        "with pyGenomneTracks. It stores data in dense format, so we recommend choosing a "
        "region to export instead of the whole dataset."
        '"mtx" refers to the MatrixMarket sparse-matrix format. The output in this case would '
        "be <prefix>.counts.mtx, along with <prefix>.rownames.txt and <prefix>.colnames.txt",
        metavar="FORMAT",
        required=True,
    )

    group.add_argument(
        "--region",
        "-r",
        help="Region of the genome to export. The format is chr:start:end, for example "
        "``--region chr10`` or ``--region chr10:456700:891000``.",
        type=ParserCommon.genomicRegion,
        metavar="CHR:START:END",
        required=False,
    )

    return parser


def parseArguments(args=None):
    io_args = ParserCommon.inputOutputOptions(
        opts=["h5adfile", "outFilePrefix"], requiredOpts=["input", "outFilePrefix"]
    )
    other_args = ParserCommon.otherOptions()

    parser = argparse.ArgumentParser(
        prog="scExportSignal",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        parents=[io_args, get_args(), other_args],
        description="""
``scExportSignal`` exports sincei-supportred .h5ad (anndata) file to other formats.
""",
        usage="""
scExportSignal -i INPUT.h5ad -f FORMAT -o OUTPUT
""",
        add_help=False,
    )

    # If no arguments are provided, show help and exit
    if args is None and len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    return parser


def main(args=None):
    """Main entry point for scExportSignal."""

    args = parseArguments().parse_args(args)
    adata = ad.read_h5ad(args.input)

    # Subset region
    if args.region is not None:
        chrom, start, end = ParserCommon.parse_region(args.region)
        adata = adata[:, adata.var["chrom"] == chrom]
        adata.var[["start", "end"]] = adata.var[["start", "end"]].apply(pd.to_numeric)
        if start is not None and end is not None:
            adata = adata[:, (adata.var["start"] >= int(start)) & (adata.var["end"] <= int(end))]

    # Export to formats
    if args.outFileFormat == "bm":
        chrom = adata.var["chrom"].to_numpy()  # assuming 'chrom' field
        start = adata.var["start"].to_numpy()  # assuming 'start' field
        end = adata.var["end"].to_numpy()  # assuming 'end' field
        row_data = adata.X.toarray()  # Extract the data matrix
        row_data = np.rot90(row_data, k=3)  # Rotate the matrix 270 degrees

        # Convert chromosomes and region starts to numeric values for sorting
        chrom_numeric = np.array([chromosome_to_numeric(c) for c in chrom])
        start_numeric = np.array([int(s) for s in start])

        # Create a polars DataFrame and sort by chromosome and region start
        df = pd.DataFrame(
            {
                "chrom_numeric": chrom_numeric,
                "start_numeric": start_numeric,
                "chrom": chrom,
                "start": start.astype(int),
                "end": end.astype(int),
            }
        )

        row_data = pd.DataFrame(row_data.astype(float))
        df = pd.concat([df, row_data], axis=1)
        df = df.sort_values(["chrom_numeric", "start_numeric"])
        # Drop numeric columns
        df = df.drop(["chrom_numeric", "start_numeric"], axis=1)
        # Write to .bm file
        df.to_csv(args.outFilePrefix + ".bm", sep="\t", header=False)

    elif args.outFileFormat == "mtx":
        mtxFile = args.outFilePrefix + ".counts.mtx"
        rowNamesFile = args.outFilePrefix + ".rownames.txt"
        colNamesFile = args.outFilePrefix + ".colnames.txt"

        # write row labels
        row_labels = adata.obs_names
        f = open(rowNamesFile, "w")
        f.write("\n".join(row_labels))
        f.write("\n")
        f.close()

        # write column labels
        col_labels = adata.var_names
        f = open(colNamesFile, "w")
        f.write("\n".join(col_labels))
        f.write("\n")
        f.close()

        # write the matrix as .mtx
        sp = sparse.csr_matrix(adata.X)
        io.mmwrite(mtxFile, sp, field="integer")


if __name__ == "__main__":
    main()
