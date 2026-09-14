from __future__ import annotations

from itertools import compress
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    import anndata as ad
    import pandas as pd
    from deeptoolsintervals import GTF

    Overlaps = dict[str, list[tuple[str, str]] | None]


## get overlap of GTF object (deeptoolsintervals) with anndata object (from sincei)
## output: dict (region->gene mapping)
def get_gtf_adata_olaps(adata: ad.AnnData, gtf: GTF) -> Overlaps:
    r"""Get overlaps between AnnData features and GTF regions.

    Parameters
    ----------
    adata : AnnData
        AnnData with regions to overlap.

    gtf : GTF
        GTF object with regions to overlap.

    Returns
    -------
    dict
        Dictionary with overlaps for each feature in adata.

    Examples
    --------
    >>> test = Tester()
    >>> gtf = GTF(test.gtfFile)
    >>> adata = sc.read_10x_mtx(
    ...     test.input_matrix_dir, var_names="gene_symbols", cache=True
    ... )
    >>> olaps = get_gtf_adata_olaps(adata, gtf)
    >>> olaps["Gm37381"]
    [('ENSMUSG00000064372', '+'), ('ENSMUSG00000064372', '-')]
    """
    var = cast("pd.DataFrame", adata.var)
    olaps: Overlaps = dict.fromkeys(var.index)
    for i, key in enumerate(var.index):
        try:
            chrom, start, end = (
                var["chrom"].iloc[i],
                int(var["start"].iloc[i]),
                int(var["end"].iloc[i]),
            )
            ol = gtf.findOverlaps(chrom, start, end, includeStrand=True)
            if ol:
                genelist = [(x[2], x[5]) for x in ol]
                olaps[key] = genelist
        except ValueError:  # noqa: PERF203
            olaps[key] = None
            continue
    return olaps


## Search for bins by gene name, return either the first bin (promoter) or all
## overlapping bins
def get_bins_by_gene(
    dict: Overlaps, gene: str, firstBin: bool = False
) -> str | list[str]:
    r"""
    Returns the bins for a given gene.

    Parameters
    ----------
    dict : dict
        Dictionary of bins and genes.
    gene : str
        Gene name.
    firstBin : bool
        If true, return only the first bin of the gene.

    Returns
    -------
    list
        List of bins.

    Examples
    --------
    >>> dict = {"chr1_1": [("gene1", "+"), ("gene2", "-")], "chr1_2": [("gene1", "+")]}
    >>> get_bins_by_gene(dict, "gene1")
    ['chr1_1', 'chr1_2']
    >>> get_bins_by_gene(dict, "gene1", firstBin=True)
    'chr1_1'
    """
    klist = []
    strand = None
    for k, v in dict.items():
        if v:
            vlist = [x[0] for x in v]  # overlapping genes
            slist = [x[1] for x in v]  # overlapping gene strands
            match = [x.lower() == gene.lower() for x in vlist]
            if any(match):
                klist.append(k)
                # get strand of the gene
                strand = next(compress(slist, match))
            else:
                strand = None

    # if firstBin is asked, sort the bins by start pos and
    # return only the firstBin by strand
    if klist and firstBin:
        spos = [x.split("_")[1] for x in klist]
        first_bin = spos.index(min(spos)) if strand == "+" else spos.index(max(spos))
        return klist[first_bin]
    return klist
