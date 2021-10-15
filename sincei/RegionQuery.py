from itertools import compress

## get overlap of GTF object (deeptoolsintervals) with anndata object (from sincei)
## output: dict (region->gene mapping)
def get_gtf_adata_olaps(adata, gtf):
    var=adata.var
    olaps=dict.fromkeys(var.index)
    for i, key in enumerate(var.index):
        try:
            chrom, start, end = var['chrom'][i], int(var['start'][i]), int(var['end'][i])
            ol=gtf.findOverlaps(chrom, start, end,  includeStrand=True)
            if ol:
                genelist=[(x[2], x[5]) for x in ol]
                olaps[key]=genelist
        except ValueError:
            olaps[key]=None
            continue
    return olaps

## Search for bins by gene name, return either the first bin (promoter) or all overlapping bins
def get_bins_by_gene(dict, gene, firstBin=False):
    klist=[]
    for k, v in dict.items():
        if v:
            vlist=[x[0] for x in v]# overlapping genes
            slist=[x[1] for x in v]# overlapping gene strands
            match=[x==gene for x in vlist]
            if any(match):
                klist.append(k)
                # get strand of the gene
                strand=list(compress(slist, match))[0]
            else:
                strand=None

    # if firstBin is asked, sort the bins by start pos and
    # return only the firstBin by strand
    if klist and firstBin:
        spos = [x.split('_')[1] for x in klist]
        if strand == '+':
            first_bin = spos.index(min(spos))
        else:
            first_bin = spos.index(max(spos))
        return klist[first_bin]
    else:
        return klist
