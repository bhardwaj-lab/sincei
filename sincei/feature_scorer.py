from concurrent.futures import ThreadPoolExecutor
import sys
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
from deeptoolsintervals import GTF
from tqdm import tqdm


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

    # Get the overlap indices for subsetting
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
):
    """
    This function computes a vector of weights for calculating the gene activity of a particular
    gene in a given region. The weights are the average exponential decay weight across each
    feature body, assuming uniform count distribution within features.

    The weights are computed as the average of: np.exp(-decay * distance / 5000) across each feature.

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

    Returns
    -------
    weights : np.ndarray
        Array of average weights for each feature.
    """
    feature_starts = np.asarray(feature_starts, dtype=np.float64)
    feature_ends = np.asarray(feature_ends, dtype=np.float64)
    feature_lengths = feature_ends - feature_starts

    # Avoid division by zero for zero-length features
    feature_lengths = np.maximum(feature_lengths, 1.0)

    if decay == 0.0:
        # No decay - all weights are 1
        return np.ones(len(feature_starts), dtype=np.float64)

    # Scale distance decay per kilobase
    lam = decay / 1000.0

    if gene_body:
        # Vectorized computation for gene_body=True
        # Determine overlap with gene body
        overlap_start = np.maximum(feature_starts, gene_start)
        overlap_end = np.minimum(feature_ends, gene_end)

        # Initialize weight sums
        weight_sums = np.zeros(len(feature_starts), dtype=np.float64)

        # Case 1: Features entirely outside gene body
        no_overlap = overlap_start >= overlap_end

        # Upstream features (feature_end <= gene_start)
        upstream = no_overlap & (feature_ends <= gene_start)
        if np.any(upstream):
            weight_sums[upstream] = np.exp(-lam * (gene_start - feature_ends[upstream])) - np.exp(
                -lam * (gene_start - feature_starts[upstream])
            )

        # Downstream features (feature_start >= gene_end)
        downstream = no_overlap & (feature_starts >= gene_end)
        if np.any(downstream):
            weight_sums[downstream] = np.exp(-lam * (feature_starts[downstream] - gene_end)) - np.exp(
                -lam * (feature_ends[downstream] - gene_end)
            )

        # Case 2: Features with some overlap with gene body
        has_overlap = ~no_overlap
        if np.any(has_overlap):
            # Upstream part
            upstream_part = feature_starts[has_overlap] < gene_start
            if np.any(upstream_part):
                idx = np.where(has_overlap)[0][upstream_part]
                weight_sums[idx] += 1.0 - np.exp(-lam * (gene_start - feature_starts[idx]))

            # Inside gene body (weight = 1)
            weight_sums[has_overlap] += overlap_end[has_overlap] - overlap_start[has_overlap]

            # Downstream part
            downstream_part = feature_ends[has_overlap] > gene_end
            if np.any(downstream_part):
                idx = np.where(has_overlap)[0][downstream_part]
                weight_sums[idx] += 1.0 - np.exp(-lam * (feature_ends[idx] - gene_end))
    else:
        # Vectorized computation for gene_body=False (decay from TSS)
        tss = gene_start if strand == "+" else gene_end

        # Initialize weight sums
        weight_sums = np.zeros(len(feature_starts), dtype=np.float64)

        # Features entirely left of TSS
        left_of_tss = feature_ends <= tss
        if np.any(left_of_tss):
            weight_sums[left_of_tss] = np.exp(-lam * (tss - feature_ends[left_of_tss])) - np.exp(
                -lam * (tss - feature_starts[left_of_tss])
            )

        # Features entirely right of TSS
        right_of_tss = feature_starts >= tss
        if np.any(right_of_tss):
            weight_sums[right_of_tss] = np.exp(-lam * (feature_starts[right_of_tss] - tss)) - np.exp(
                -lam * (feature_ends[right_of_tss] - tss)
            )

        # Features spanning TSS
        spans_tss = (feature_starts < tss) & (feature_ends > tss)
        if np.any(spans_tss):
            weight_sums[spans_tss] = (
                2.0 - np.exp(-lam * (tss - feature_starts[spans_tss])) - np.exp(-lam * (feature_ends[spans_tss] - tss))
            )

    return weight_sums / feature_lengths


def _compute_weight_sum(
    region_start,
    region_end,
    gene_start,
    gene_end,
    decay,
    gene_body=True,
    strand="+",
):
    """
    Compute the finite sum of decay weights over a specific region.

    This is used to calculate the weight contribution that a sub-region would make,
    so it can be subtracted from the total weight when excluding regions.

    Parameters
    ----------
    region_start : float
        Start position of the region.
    region_end : float
        End position of the region.
    gene_start : int
        Start position of the gene of interest.
    gene_end : int
        End position of the gene of interest.
    decay : float
        Decay parameter.
    gene_body : bool
        Whether gene body has weight 1.
    strand : str
        Strand of the gene ('+' or '-').

    Returns
    -------
    float
        The finite sum of weights over the region.
    """
    if region_start >= region_end:
        return 0.0

    lam = decay / 5000.0

    if lam == 0:
        # No decay - weight is 1 everywhere
        return region_end - region_start

    rs, re = region_start, region_end

    if gene_body:
        # Weight = 1 inside gene body, exponential decay outside
        weight_sum = 0.0

        # Determine overlap with gene body
        overlap_start = max(rs, gene_start)
        overlap_end = min(re, gene_end)

        if overlap_start >= overlap_end:
            # No overlap with gene body - entirely outside
            if re <= gene_start:
                # Entirely upstream of gene
                weight_sum = np.exp(-lam * (gene_start - re)) - np.exp(-lam * (gene_start - rs))
            else:
                # Entirely downstream of gene (rs >= gene_end)
                weight_sum = np.exp(-lam * (rs - gene_end)) - np.exp(-lam * (re - gene_end))
        else:
            # Some overlap with gene body
            # Part upstream of gene body (if any)
            if rs < gene_start:
                weight_sum += 1.0 - np.exp(-lam * (gene_start - rs))

            # Part inside gene body (weight = 1)
            weight_sum += overlap_end - overlap_start

            # Part downstream of gene body (if any)
            if re > gene_end:
                weight_sum += 1.0 - np.exp(-lam * (re - gene_end))
    else:
        # Decay from TSS
        tss = gene_start if strand == "+" else gene_end

        if re <= tss:
            # Region entirely left of TSS
            weight_sum = np.exp(-lam * (tss - re)) - np.exp(-lam * (tss - rs))
        elif rs >= tss:
            # Region entirely right of TSS
            weight_sum = np.exp(-lam * (rs - tss)) - np.exp(-lam * (re - tss))
        else:
            # Region spans TSS
            weight_sum = 2.0 - np.exp(-lam * (tss - rs)) - np.exp(-lam * (re - tss))

    return weight_sum


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
        for i, gene in enumerate(gtf.findOverlaps(chrom, 0, 2**31 - 1)):  # Get all genes on chromosome
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


def _compute_gene_activity_single(
    adata,
    gene_row,
    max_region,
    decay,
    gene_body,
    gene_size_factor,
    avg_gene_length,
    overlap_policy="partial",
    exclude_in_range=None,
    extend_TSS=2000,
    genes_arrays=None,
    chrom_indices_cache=None,
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

    # Define region to search (gene body + max_region upstream/downstream)
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

    # Calculate decay weights (average weight across each feature body)
    weights = get_decay_weights(
        gene_start,
        gene_end,
        feature_starts,
        feature_ends,
        strand=strand,
        decay=decay,
        gene_body=gene_body,
    )

    # Apply overlap policy: scale weights for partially overlapping features
    if overlap_policy == "partial":
        clip_start = np.maximum(feature_starts, region_start)
        clip_end = np.minimum(feature_ends, region_end)
        feat_lengths = np.maximum(feature_ends - feature_starts, 1.0)
        overlap_fractions = (clip_end - clip_start) / feat_lengths
        weights = weights * overlap_fractions

    # Apply exclusion weight reduction if requested
    if exclude_in_range in ("TSS", "genes") and genes_arrays is not None:
        # Filter genes using pre-converted arrays for performance
        chrom_mask = genes_arrays["chrom"] == chrom
        name_mask = genes_arrays["gene_name"] != gene_name
        region_mask = (genes_arrays["start"] < region_end) & (genes_arrays["end"] > region_start)
        other_genes_mask = chrom_mask & name_mask & region_mask

        # Early termination if no other genes in region
        if not np.any(other_genes_mask):
            pass  # No exclusion needed
        else:
            # Extract relevant gene arrays
            other_starts = genes_arrays["start"][other_genes_mask]
            other_ends = genes_arrays["end"][other_genes_mask]
            other_strands = genes_arrays["strand"][other_genes_mask]

            # Compute TSS for all other genes at once
            other_tss = np.where(other_strands == "+", other_starts, other_ends)

            # Compute exclusion regions for all other genes
            if exclude_in_range == "TSS":
                exclude_starts = other_tss - extend_TSS
                exclude_ends = other_tss + extend_TSS
            else:  # exclude_in_range == "genes"
                # Vectorized computation of gene body extensions
                is_plus = other_strands == "+"
                exclude_starts = np.where(is_plus, other_starts - extend_TSS, other_starts)
                exclude_ends = np.where(is_plus, other_ends, other_ends + extend_TSS)

            # Vectorized overlap computation: broadcast features × exclusion regions
            # Shape: (n_features, n_other_genes)
            fs_broadcast = feature_starts[:, np.newaxis]  # (n_features, 1)
            fe_broadcast = feature_ends[:, np.newaxis]

            # Check overlaps for all feature-exclusion pairs at once
            overlaps = (fs_broadcast < exclude_ends) & (fe_broadcast > exclude_starts)

            # Early termination if no overlaps
            if not np.any(overlaps):
                pass  # No exclusion needed
            else:
                # Compute overlap regions (vectorized)
                overlap_starts = np.maximum(fs_broadcast, exclude_starts)
                overlap_ends = np.minimum(fe_broadcast, exclude_ends)

                # Only compute weights for actual overlaps
                overlap_starts_flat = overlap_starts[overlaps]
                overlap_ends_flat = overlap_ends[overlaps]

                # Compute weight sums for all overlaps at once
                overlap_sums = np.array(
                    [
                        _compute_weight_sum(
                            overlap_starts_flat[i],
                            overlap_ends_flat[i],
                            gene_start,
                            gene_end,
                            decay,
                            gene_body=gene_body,
                            strand=strand,
                        )
                        for i in range(len(overlap_starts_flat))
                    ]
                )

                # Pre-compute feature lengths once
                feature_lengths = feature_ends - feature_starts

                # Accumulate reductions per feature
                weight_reductions = np.zeros(len(feature_starts), dtype=np.float64)
                feature_indices, _ = np.where(overlaps)

                # Vectorized accumulation
                np.add.at(weight_reductions, feature_indices, overlap_sums / feature_lengths[feature_indices])

                # Subtract accumulated reductions and clamp to zero minimum
                weights = np.maximum(0.0, weights - weight_reductions)

    # Apply gene size factor if requested
    if gene_size_factor and avg_gene_length > 0:
        gene_length = gene_end - gene_start
        size_factor = gene_length / avg_gene_length
        weights = weights * size_factor

    # Get counts for overlapping features and compute weighted sum
    counts = adata.X[:, overlap_indices]

    if sparse.issparse(counts):
        # Efficient sparse matrix multiplication with weights
        activity = np.asarray(counts.dot(weights)).ravel()
    else:
        activity = counts.dot(weights)

    return (gene_name, activity)


def feature_scorer(
    adata,
    gtf,
    mode,
    overlap_policy="partial",
    penalty=None,
    decay=0.75,
    max_region=200000,
    gene_body=True,
    gene_size_factor=True,
    exclude_in_range=None,
    extend_TSS=2000,
    chrs_to_skip=None,
    verbose=False,
    n_threads=1,
):
    """
    This function calculates a cell x gene matrix with gene activity scores each cell.
    The function first parses the BED/GTF file to get gene/feature annotations, then it identifies
    the relevant genomic region (including upstream/downstream regions if specified),
    retrieves the counts of features overlapping with that region, applies decay weights if specified,
    and computes the weighted sum of counts to obtain the gene activity scores for each cell.

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
    penalty : float, optional
        Optional parameter to select VCRs of a particular penalty value from a BED file with VCRs
        calculated using multiple penalties.
    decay : float, optional
        Decay parameter for calculating the decay weights, by default 0.75. Higher values lead to
        faster decay. Weights are calculated as ``exp(-decay * distance_in_kb)``. This parameter is
        ignored in ``aggregate`` mode.
    max_region : int, optional
        Maximum region size around the gene (upstream and downstream) to consider, by default
        200000 bp.
    gene_body : bool, optional
        Whether the weight of the gene body is considered as 1 like the TSS, by default True.
        If True, the weights the decay starts beyond the gene body.
    gene_size_factor : bool, optional
        Whether to apply a size factor based on gene length, by default True. If True, the weights
        are multiplied by a size factor calculated as (gene_length / average_gene_length) to account
        for gene length bias.
    exclude_in_range : str, optional
        Whether to exclude regions of other genes from contributing to this gene's activity score.
        Options are:
        - None: No exclusion (default)
        - "TSS": Exclude regions around TSS of other genes (TSS ± extend_TSS)
        - "genes": Exclude gene bodies of other genes, extended upstream by extend_TSS
        Invalid values default to None.
    extend_TSS : int, optional
        Number of base pairs to extend around TSS for exclusion regions, by default 2000.
        Used when exclude_in_range is "TSS" or "genes".
    chrs_to_skip : list, optional
        List of chromosomes to skip, by default None.
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
    if penalty is not None and "name" in genes_df.columns:
        genes_df = genes_df[genes_df["name"].str.contains(f"_pen{penalty}", na=False)]
        if genes_df.empty:
            raise ValueError(
                f"No VCRs found with penalty value {penalty} in the VCR BED file. "
                f"Check the 5th column of the BED file for available penalty values."
            )

    # Ensure adata.var coordinate columns are numeric (may be categorical from h5ad)
    for col in ["start", "end"]:
        if col in adata.var.columns and hasattr(adata.var[col], "cat"):
            adata.var[col] = adata.var[col].astype(int)

    # Filter chromosomes if requested
    if chrs_to_skip is not None:
        genes_df = genes_df[~genes_df["chrom"].isin(chrs_to_skip)]
        adata = adata[:, ~adata.var["chrom"].isin(chrs_to_skip)].copy()

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

    # Calculate average gene length for size factor normalization
    avg_gene_length = (genes_df["end"] - genes_df["start"]).mean()

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

    # Pre-convert genes_df to numpy arrays for faster access (opt 4)
    genes_arrays = None
    if exclude_in_range in ("TSS", "genes"):
        genes_arrays = {
            "chrom": genes_df["chrom"].values,
            "gene_name": genes_df["gene_name"].values,
            "start": genes_df["start"].values.astype(np.int64),
            "end": genes_df["end"].values.astype(np.int64),
            "strand": genes_df["strand"].values,
        }

    # Pre-compute chromosome indices cache (opt 7)
    chrom_indices_cache = {}
    for chrom in common_chroms:
        chrom_mask = adata.var["chrom"] == chrom
        chrom_indices_cache[chrom] = np.where(chrom_mask)[0]

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
            avg_gene_length,
            overlap_policy=overlap_policy,
            exclude_in_range=exclude_in_range,
            extend_TSS=extend_TSS,
            genes_arrays=genes_arrays,
            chrom_indices_cache=chrom_indices_cache,
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
                    desc="Processing genes",
                    disable=not verbose,
                )
            )
    else:
        results = [process_gene(g) for g in tqdm(gene_rows, desc="Processing genes", disable=not verbose)]

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

    sys.stdout.write(f"Created AnnData with {adata_out.n_obs} cells and {adata_out.n_vars} genes\n")

    return adata_out
