from concurrent.futures import ThreadPoolExecutor
import sys
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
from deeptoolsintervals import GTF
from tqdm import tqdm


def _parse_gtf_genes(gtf_path):
    """
    Parse a GTF/BED file using deeptoolsintervals and extract gene/feature information.

    Returns a DataFrame with gene_name, chrom, start, end, strand, and score.
    For BED files, score corresponds to the 5th column (e.g. penalty value).
    For GTF files, score is typically the file name and can be ignored.
    """
    gtf = GTF(
        gtf_path, exonID="exon", transcriptID="transcript", transcript_id_designator="transcript_id", keepExons=False
    )

    genes = []
    for chrom in gtf.chroms:
        # Get all features on chromosome (avoid overflow for int32)
        for i, gene in enumerate(gtf.findOverlaps(chrom, 0, 2**31 - 1)):
            # gene is a tuple: (start, end, name, source/strand, exons, score)
            gene_start = gene[0]
            gene_end = gene[1]
            gene_name = gene[2] if len(gene) > 2 else f"Feature_{i}"
            gene_strand = gene[3] if len(gene) > 3 else "+"
            gene_score = gene[5] if len(gene) > 5 else None

            genes.append(
                {
                    "gene_name": gene_name,
                    "chrom": chrom,
                    "start": gene_start,
                    "end": gene_end,
                    "strand": gene_strand,
                    "score": gene_score,
                }
            )

    return pd.DataFrame(genes)


def get_indices_overlapping(
    adata,
    chrom,
    start,
    end,
):
    """
    This function takes an AnnData object and a region defined by chromosome, start, and end positions.
    It returns the overlap indices of features overlapping with the region.

    Parameters
    ----------
    adata : AnnData
        The input AnnData object containing the data.
    chrom : str
        The chromosome of the region.
    start : int
        The start position of the region.
    end : int
        The end position of the region.

    Returns
    -------
    overlap_indices : np.ndarray or None
        Array of global feature indices that overlap with the region, or None if no overlaps.
    """
    # Filter to the chromosome
    chrom_mask = adata.var["chrom"] == chrom
    if not chrom_mask.any():
        return None

    chrom_var = adata.var[chrom_mask]

    # Find overlapping features: feature_start < end AND feature_end > start
    overlap_mask = (chrom_var["start"].values < end) & (chrom_var["end"].values > start)

    if not overlap_mask.any():
        return None

    # Get the overlapping feature indices within the chromosome subset
    overlap_indices = np.where(overlap_mask)[0]

    # Get the overlap indices within the whole anndata
    chrom_indices = np.where(chrom_mask.values)[0]
    overlap_indices = chrom_indices[overlap_indices]

    return overlap_indices


def get_decay_weights(
    gene_start,
    gene_end,
    feature_starts,
    feature_ends,
    strand="+",
    decay=0.75,
    gene_body=True,
    excluded_regions=[],
):
    """
    This function computes a vector of weights for calculating the gene activity of a particular
    gene in a given region. The weights are the average exponential decay weight across each
    feature body, assuming uniform count distribution within features.
    Features in ``excluded_regions`` are assigned a weight of 0.

    The weights are computed as the average of: np.exp(-decay * distance / 10000) across each feature.

    Parameters
    ----------
    gene_start : int
        The start position of the gene of interest.
    gene_end : int
        The end position of the gene of interest.
    feature_starts : np.ndarray
        Array of feature start positions.
    feature_ends : np.ndarray
        Array of feature end positions.
    strand : str, optional
        The strand of the gene ('+' or '-'), by default '+'.
    decay : float, optional
        Decay parameter for weighting, by default 1.0. Higher values lead to faster decay.
    gene_body : bool, optional
        Whether the weight of the gene body is considered as 1 like the TSS, by default True.
        If True, the decay starts beyond the gene body.
    excluded_regions : list of tuples, optional
        List of (start, end) tuples defining regions to exclude from contributing to the activity score (weight 0).

    Returns
    -------
    weights : np.ndarray
        Array of average weights for each feature.
    """
    feature_starts = np.asarray(feature_starts, dtype=np.float64)
    feature_ends = np.asarray(feature_ends, dtype=np.float64)

    if decay == 0.0:
        # No decay - all weights are 1
        return np.ones(len(feature_starts), dtype=np.float64)

    # Scale distance decay per kilobase
    lam = decay / 10000.0

    if gene_body:
        # Vectorized computation for gene_body=True
        # Determine overlap with gene body
        overlap_start = np.maximum(feature_starts, gene_start)
        overlap_end = np.minimum(feature_ends, gene_end)

        # Initialize weight sums
        weights = np.zeros(len(feature_starts), dtype=np.float64)

        # Case 1: Features entirely outside gene body
        no_overlap = overlap_start >= overlap_end

        # Upstream features (feature_end <= gene_start)
        upstream = no_overlap & (feature_ends <= gene_start)
        if np.any(upstream):
            weights[upstream] = np.exp(-lam * (gene_start - feature_ends[upstream])) - np.exp(
                -lam * (gene_start - feature_starts[upstream])
            )

        # Downstream features (feature_start >= gene_end)
        downstream = no_overlap & (feature_starts >= gene_end)
        if np.any(downstream):
            weights[downstream] = np.exp(-lam * (feature_starts[downstream] - gene_end)) - np.exp(
                -lam * (feature_ends[downstream] - gene_end)
            )

        # Case 2: Features with some overlap with gene body
        has_overlap = ~no_overlap
        if np.any(has_overlap):
            # Upstream part
            upstream_part = feature_starts[has_overlap] < gene_start
            if np.any(upstream_part):
                idx = np.where(has_overlap)[0][upstream_part]
                weights[idx] += 1.0 - np.exp(-lam * (gene_start - feature_starts[idx]))

            # Inside gene body (weight = 1)
            weights[has_overlap] += overlap_end[has_overlap] - overlap_start[has_overlap]

            # Downstream part
            downstream_part = feature_ends[has_overlap] > gene_end
            if np.any(downstream_part):
                idx = np.where(has_overlap)[0][downstream_part]
                weights[idx] += 1.0 - np.exp(-lam * (feature_ends[idx] - gene_end))
    else:
        # Vectorized computation for gene_body=False (decay from TSS)
        tss = gene_start if strand == "+" else gene_end

        # Initialize weights
        weights = np.zeros(len(feature_starts), dtype=np.float64)

        # Features entirely left of TSS
        left_of_tss = feature_ends <= tss
        if np.any(left_of_tss):
            weights[left_of_tss] = np.exp(-lam * (tss - feature_ends[left_of_tss])) - np.exp(
                -lam * (tss - feature_starts[left_of_tss])
            )

        # Features entirely right of TSS
        right_of_tss = feature_starts >= tss
        if np.any(right_of_tss):
            weights[right_of_tss] = np.exp(-lam * (feature_starts[right_of_tss] - tss)) - np.exp(
                -lam * (feature_ends[right_of_tss] - tss)
            )

        # Features overlapping TSS
        overlap_tss = (feature_starts < tss) & (feature_ends > tss)
        if np.any(overlap_tss):
            weights[overlap_tss] = (
                2.0
                - np.exp(-lam * (tss - feature_starts[overlap_tss]))
                - np.exp(-lam * (feature_ends[overlap_tss] - tss))
            )

    for exclude_start, exclude_end in excluded_regions:
        # Vectorized exclusion of specified regions
        exclude_mask = (feature_starts < exclude_end) & (feature_ends > exclude_start)
        weights[exclude_mask] = 0.0

    return weights


def _compute_gene_activity_single(
    adata,
    gene_row,
    max_region,
    decay,
    gene_body,
    gene_size_factor,
    overlap_policy="partial",
    exclude_in_range=None,
    genes_arrays=None,
):
    """
    Compute gene activity for a single gene.

    Returns (gene_name, activity_vector) or None if no overlapping features.
    """
    chrom = gene_row["chrom"]
    gene_start = int(gene_row["start"])
    gene_end = int(gene_row["end"])
    strand = gene_row.get("strand", "+")
    gene_name = gene_row["gene_name"]

    # Define region to consider (gene body + max_region upstream/downstream)
    max_region = max_region * 1000  # convert from kb to base pairs
    region_start = max(0, gene_start - max_region)
    region_end = gene_end + max_region

    # Get indices of features overlapping with the region
    overlap_indices = get_indices_overlapping(adata, chrom, region_start, region_end)
    if overlap_indices is None:
        return None

    # Get feature coordinates for decay calculation
    feature_starts = adata.var["start"].values[overlap_indices]
    feature_ends = adata.var["end"].values[overlap_indices]

    if overlap_policy not in ["partial", "all", "none"]:
        sys.stderr.write(f"WARNING: Invalid overlap_policy '{overlap_policy}'. Defaulting to 'partial'.")
        overlap_policy = "partial"

    # Apply overlap policy: filter features based on how they overlap the search region
    if overlap_policy == "none":
        # Only keep features fully contained within the region
        fully_contained = (feature_starts >= region_start) & (feature_ends <= region_end)
        if not np.any(fully_contained):
            return None
        overlap_indices = overlap_indices[fully_contained]
        feature_starts = feature_starts[fully_contained]
        feature_ends = feature_ends[fully_contained]

    # Fetch excluded regions for this gene if requested
    excluded_regions = []
    if exclude_in_range in ("TSS", "genes") and genes_arrays is not None:
        # Filter genes using pre-converted arrays for performance
        chrom_mask = genes_arrays["chrom"] == chrom
        name_mask = genes_arrays["gene_name"] != gene_name
        region_mask = (genes_arrays["start"] < region_end) & (genes_arrays["end"] > region_start)
        other_genes_mask = chrom_mask & name_mask & region_mask

        # Extract excluded regions for the gene of interest
        if exclude_in_range == "TSS":
            strand_mask = genes_arrays["strand"][other_genes_mask] == strand
            excluded_regions = np.where(
                strand_mask, genes_arrays["start"][other_genes_mask], genes_arrays["end"][other_genes_mask]
            )
            excluded_regions = list(zip(excluded_regions, excluded_regions))
        elif exclude_in_range == "genes":
            excluded_regions = list(zip(genes_arrays["start"][other_genes_mask], genes_arrays["end"][other_genes_mask]))
        else:
            excluded_regions = []

    # Calculate decay weights (average weight across each feature body)
    weights = get_decay_weights(
        gene_start=gene_start,
        gene_end=gene_end,
        feature_starts=feature_starts,
        feature_ends=feature_ends,
        strand=strand,
        decay=decay,
        gene_body=gene_body,
        excluded_regions=excluded_regions,
    )

    # Apply overlap policy: scale weights for partially overlapping features
    if overlap_policy == "partial":
        clip_start = np.maximum(feature_starts, region_start)
        clip_end = np.minimum(feature_ends, region_end)
        feat_lengths = np.maximum(feature_ends - feature_starts, 1.0)
        overlap_fractions = (clip_end - clip_start) / feat_lengths
        weights = weights * overlap_fractions

    # Apply gene size factor if requested
    if gene_size_factor:
        gene_length = gene_end - gene_start
        size_factor = gene_length
        weights = weights * size_factor

    # Get counts for overlapping features and compute weighted sum
    counts = adata.X[:, overlap_indices]

    if sparse.issparse(counts):
        # Efficient sparse matrix multiplication with weights
        activity = np.asarray(counts.dot(weights)).ravel()
    else:
        activity = counts.dot(weights)

    return (gene_name, activity)


def FeatureScorer(
    adata,
    gtf,
    mode,
    overlap_policy="partial",
    penalty=None,
    decay=0.75,
    max_region=100,
    gene_body=True,
    gene_size_factor=True,
    exclude_in_range=None,
    center_scores=False,
    verbose=False,
    n_threads=1,
):
    """
    This function calculates a cell x gene matrix with gene activity scores.
    First, it parses the input BED/GTF file to get gene/feature annotations, then it identifies
    the relevant genomic region (including upstream/downstream regions if specified),
    retrieves the counts of features overlapping with that region, applies decay weights if specified,
    computes the weighted sum of counts to obtain the gene activity scores for each cell, and
    L1-normalizes the scores row-wise (per cell).

    Parameters
    ----------
    adata : AnnData
        The input AnnData object containing the data.
    gtf : str
        Path to the BED/GTF file with region annotations.
    mode : str
        Scoring mode. Options are 'aggregate' or 'activities'.
        ``aggregate`` calculates the total counts of the genomic features in the input BED/GTF file from the
        input anndata.
        ``activities`` mode calculates the weighted sum of counts based on distance to TSS of the genes
        in the input GTF file. The weights are calculated using an exponential decay function.
    overlap_policy: str, optional
        Policy for handling adata features that only partially overlap regions in the BED/GTF provided.
        Options are:
            - ``partial``: count reads in anndata feature proportionally to the overlap fraction.
              counts_considered = feature_counts * overlap_length / region_length.
            - ``all``: count all reads in the partially overlapping anndata feature.
            - ``none``: exclude reads from partially overlapping anndata features, in other words, only
              count reads in anndata features fully contained within BED/GTF regions.
        Default is 'partial'.
    center_scores : bool, optional
        Whether to scale the scores to unit variance and center them around zero, by default False.
        This destroys the sparsity of the output matrix and can lead to increased memory usage.
        Use with caution for large datasets.
    penalty : float, optional
        Optional parameter to select VCRs of a particular penalty value from a BED file with VCRs
        calculated using multiple penalties.
    decay : float, optional
        Decay parameter for calculating the decay weights, by default 0.75. Higher values lead to
        faster decay. Weights are calculated as ``exp(-decay * distance_in_kb / 10)``. This parameter
        is ignored in ``aggregate`` mode.
    max_region : int, optional
        Maximum region size around the gene (upstream and downstream) to consider (in kilobases),
        by default 100 Kb.
    gene_body : bool, optional
        Whether the weight of the gene body is considered as 1 like the TSS, by default True.
        If True, the weights the decay starts beyond the gene body.
    gene_size_factor : bool, optional
        Whether to divide scores by gene length to account for gene length bias, by default True.
    exclude_in_range : str, optional
        Whether to exclude regions of other genes from contributing to this gene's activity score.
        Options are:
        - None: No exclusion (default)
        - "TSS": Exclude features overlapping the TSS of other genes
        - "genes": Exclude features overlapping the bodies of other genes
        Invalid values default to None.
    center_scores : bool, optional
        Whether to scale the scores to unit variance and center them around zero, by default False.
        This destroys the sparsity of the output matrix and can lead to increased memory usage.
        Use with caution for large datasets.
    verbose : bool, optional
        Print progress messages and warnings. Default is False.
    n_threads : int, optional
        Number of threads to use for parallel processing, by default 1.

    Returns
    -------
    adata_out : AnnData
        AnnData object with cells as obs and genes as var, containing gene activity scores.
    """
    # Parse BED/GTF file to get gene annotations
    sys.stdout.write("Parsing BED/GTF file...\n")
    genes_df = _parse_gtf_genes(gtf)

    if genes_df.empty:
        raise ValueError("No genes/features found in the input file.")

    # Filter VCR BED by penalty value
    if penalty is not None and "gene_name" in genes_df.columns:
        genes_df = genes_df[genes_df["gene_name"].str.contains(f"_pen{penalty}", na=False)]
        if genes_df.empty:
            raise ValueError(
                f"No VCRs found with penalty value {penalty} in the VCR BED file. "
                f"Check the 5th column of the BED file for available penalty values."
            )

    # Ensure adata.var coordinate columns are numeric (may be categorical from h5ad)
    for col in ["start", "end"]:
        if col in adata.var.columns and hasattr(adata.var[col], "cat"):
            adata.var[col] = adata.var[col].astype(int)

    # Keep only chromosomes present in both data and BED/GTF
    common_chroms = set(adata.var["chrom"].unique()) & set(genes_df["chrom"].unique())
    genes_df = genes_df[genes_df["chrom"].isin(common_chroms)]

    if genes_df.empty:
        raise ValueError("No common chromosomes between data and BED/GTF")

    # Remove duplicate gene names (keep first occurrence)
    genes_df = genes_df.drop_duplicates(subset="gene_name", keep="first")

    # Validate exclude_in_range parameter
    if exclude_in_range is not None and exclude_in_range not in ("TSS", "genes"):
        sys.stderr.write(f"WARNING: Invalid exclude_in_range value '{exclude_in_range}'. Defaulting to None.\n")
        exclude_in_range = None

    sys.stdout.write(f"Processing {len(genes_df)} features across {len(common_chroms)} chromosomes\n")

    n_cells = adata.n_obs

    if mode == "aggregate":
        # aggregate mode: simple sum of counts within VCR
        effective_decay = 0.0
        effective_max_region = 0
        gene_body = True
        gene_size_factor = False
        exclude_in_range = None
    elif mode == "activities":
        effective_decay = decay
        effective_max_region = max_region
    else:
        raise ValueError(f"Unknown mode: {mode}. Must be 'aggregate' or 'activities'")

    # Pre-convert genes_df to numpy arrays for faster access
    genes_arrays = None
    if exclude_in_range in ("TSS", "genes"):
        genes_arrays = {
            "chrom": genes_df["chrom"].values,
            "gene_name": genes_df["gene_name"].values,
            "start": genes_df["start"].values.astype(np.int64),
            "end": genes_df["end"].values.astype(np.int64),
            "strand": genes_df["strand"].values,
        }

    # Prepare gene rows for processing
    gene_rows = [row for _, row in genes_df.iterrows()]

    def process_gene(gene_row):
        return _compute_gene_activity_single(
            adata,
            gene_row,
            effective_max_region,
            effective_decay,
            gene_body,
            gene_size_factor,
            overlap_policy=overlap_policy,
            exclude_in_range=exclude_in_range,
            genes_arrays=genes_arrays,
        )

    # Accumulate results using COO format for efficiency
    all_rows = []
    all_cols = []
    all_data = []
    gene_names = []
    gene_col_idx = 0

    sys.stdout.write("Computing features...\n")

    if n_threads > 1:
        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            results = list(
                tqdm(
                    executor.map(process_gene, gene_rows),
                    total=len(gene_rows),
                    desc="Processing features",
                    disable=not verbose,
                )
            )
    else:
        results = [process_gene(g) for g in tqdm(gene_rows, desc="Processing features", disable=not verbose)]

    # Collect results into COO components
    for result in results:
        if result is None:
            continue

        gene_name, activity = result

        # Find non-zero entries
        nonzero_mask = activity != 0
        nonzero_rows = np.where(nonzero_mask)[0]

        if len(nonzero_rows) > 0:
            all_rows.append(nonzero_rows)
            all_cols.append(np.full(len(nonzero_rows), gene_col_idx, dtype=np.int32))
            all_data.append(activity[nonzero_mask].astype(np.float32))

        gene_names.append(gene_name)
        gene_col_idx += 1

    # Build sparse matrix from COO components
    if all_rows:
        all_rows = np.concatenate(all_rows)
        all_cols = np.concatenate(all_cols)
        all_data = np.concatenate(all_data)
        activity_matrix = sparse.csr_matrix(
            (all_data, (all_rows, all_cols)),
            shape=(n_cells, len(gene_names)),
            dtype=np.float32,
        )
    else:
        activity_matrix = sparse.csr_matrix((n_cells, 0), dtype=np.float32)
        sys.stderr.write("WARNING: No gene activities computed - check chromosome naming consistency\n")

    # Create output AnnData
    var_df = pd.DataFrame({"gene_name": gene_names}, index=gene_names)

    # Add gene coordinates to var
    gene_info = genes_df.set_index("gene_name").loc[gene_names, ["chrom", "start", "end", "strand"]]
    var_df = var_df.join(gene_info)

    adata_out = ad.AnnData(
        X=activity_matrix,
        obs=adata.obs.copy(),
        var=var_df,
    )
    if mode == "activities":
        if center_scores:
            from scanpy.pp import scale

            sys.stderr.write(
                "WARNING: Centering the scores destroys the sparsity of the output matrix and can lead "
                "to increased memory usage. Use with caution for large datasets.\n"
            )
            scale(adata_out, zero_center=True)
        else:
            from sklearn.preprocessing import normalize

            normalize(adata_out.X, norm="l1", copy=False)

    sys.stdout.write(f"Created AnnData with {adata_out.n_obs} cells and {adata_out.n_vars} features.\n")

    return adata_out
