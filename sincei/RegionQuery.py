from itertools import compress


## get overlap of GTF object (deeptoolsintervals) with anndata object (from sincei)
## output: dict (region->gene mapping)
def get_gtf_adata_olaps(adata, gtf):
    r"""Get overlaps between adata and gtf

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    gtf : GTF
        GTF object

    Returns
    -------
    dict
        Dictionary with overlaps for each gene in adata.

    Examples
    --------
    >>> test = Tester()
    >>> gtf = GTF(test.gtfFile)
    >>> adata = sc.read_10x_mtx(test.input_matrix_dir, var_names='gene_symbols', cache=True)
    >>> olaps=get_gtf_adata_olaps(adata, gtf)
    >>> olaps['Gm37381']
    [('ENSMUSG00000064372', '+'), ('ENSMUSG00000064372', '-')]
    """
    var = adata.var
    olaps = dict.fromkeys(var.index)
    for i, key in enumerate(var.index):
        try:
            chrom, start, end = (
                var["chrom"][i],
                int(var["start"][i]),
                int(var["end"][i]),
            )
            ol = gtf.findOverlaps(chrom, start, end, includeStrand=True)
            if ol:
                genelist = [(x[2], x[5]) for x in ol]
                olaps[key] = genelist
        except ValueError:
            olaps[key] = None
            continue
    return olaps


## Search for bins by gene name, return either the first bin (promoter) or all overlapping bins
def get_bins_by_gene(dict, gene, firstBin=False):
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
    >>> dict = {'chr1_1': [('gene1', '+'), ('gene2', '-')], 'chr1_2': [('gene1', '+')]}
    >>> get_bins_by_gene(dict, 'gene1')
    ['chr1_1', 'chr1_2']
    >>> get_bins_by_gene(dict, 'gene1', firstBin=True)
    'chr1_1'
    """
    klist = []
    for k, v in dict.items():
        if v:
            vlist = [x[0] for x in v]  # overlapping genes
            slist = [x[1] for x in v]  # overlapping gene strands
            match = [x.lower() == gene.lower() for x in vlist]
            if any(match):
                klist.append(k)
                # get strand of the gene
                strand = list(compress(slist, match))[0]
            else:
                strand = None

    # if firstBin is asked, sort the bins by start pos and
    # return only the firstBin by strand
    if klist and firstBin:
        spos = [x.split("_")[1] for x in klist]
        if strand == "+":
            first_bin = spos.index(min(spos))
        else:
            first_bin = spos.index(max(spos))
        return klist[first_bin]
    else:
        return klist
