from __future__ import annotations

from itertools import compress
from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    import anndata as ad
    import pandas as pd

    from sincei._sincei import GenomeAnnotation

    Overlaps = dict[str, list[tuple[str, str]] | None]


def get_gtf_adata_olaps(adata: ad.AnnData, annotation: GenomeAnnotation) -> Overlaps:
    r"""Get overlaps between AnnData features and annotation regions.

    Parameters
    ----------
    adata : AnnData
        AnnData with regions to overlap.

    annotation : GenomeAnnotation
        Parsed BED/GTF/GFF file, from ``sincei._sincei.parse_annotation``.

    Returns
    -------
    dict
        For each feature in adata, the (name, strand) of the overlapping annotation
        features, or None if there are none.

    Examples
    --------
    >>> from sincei._sincei import parse_annotation
    >>> annotation = parse_annotation(["Chrna9.gtf"])
    >>> olaps = get_gtf_adata_olaps(adata, annotation)
    >>> olaps["chr5_100000_200000::None"]
    [('ENSMUSG00000029205', '+')]
    """
    features = annotation.features()
    names, strands = features["name"], features["strand"]
    var = cast("pd.DataFrame", adata.var)
    olaps: Overlaps = dict.fromkeys(var.index)
    for i, key in enumerate(var.index):
        try:
            chrom, start, end = (
                str(var["chrom"].iloc[i]),
                int(var["start"].iloc[i]),
                int(var["end"].iloc[i]),
            )
        except ValueError:
            continue
        hits = annotation.find_overlaps(chrom, start, end)
        if hits:
            olaps[key] = [(names[j], strands[j]) for j in hits]
    return olaps


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

    # if firstBin, sort the bins by start pos and return only the firstBin by strand
    if klist and firstBin:
        spos = [x.split("_")[1] for x in klist]
        first_bin = spos.index(min(spos)) if strand == "+" else spos.index(max(spos))
        return klist[first_bin]
    return klist
