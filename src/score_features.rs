//! Aggregation of binned single-cell chromatin data into feature scores.
//!
//! For each feature in a BED/GTF/GFF file, sums the counts of the bins that
//! overlap it, giving a per-cell score. Bins that only partly overlap a feature
//! are handled according to the overlap policy.
//!
//! The output is an AnnData h5ad file (cells × features) whose obs mirrors the
//! input and whose var carries `chrom`, `start`, `end`, `strand`, and
//! `gene_name` columns.
use std::path::{Path, PathBuf};

use ahash::{AHashMap, AHashSet};
use anndata::data::DataFrameIndex;
use anndata::{AnnData, AnnDataOp, Backend};
use anndata_hdf5::H5;
use anyhow::{Context, Result, bail};
use nalgebra_sparse::{CscMatrix, CsrMatrix};
use polars::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::annotation::parse_annotation::parse_annotation_files;
use crate::annotation::region_index::{ChromIndex, Feature, Interval};
use crate::counting::count_utils::{df_i64_col, df_str_col, read_x_f64};

/// Parse an annotation file (BED, GTF, or GFF3) into [`Feature`]s, reusing the
/// shared `parse_annotation_files` machinery (and thus the same format
/// detection as `scCountReads`).
///
/// `feature_types` selects the column-3 types kept as a region (e.g. `["gene"]`,
/// or the several transcript biotypes an Ensembl GFF3 uses); `name_attr` is the
/// column-9 attribute used to name each feature (e.g. `"gene_id"`). `None` falls
/// back to the per-format defaults. When `metagene` is set, rows of
/// `exon_types` are merged per transcript (or per gene, with `name_attr`), so
/// each output feature is one transcript. All only affect GTF/GFF inputs and are
/// ignored for BED.
///
/// When `penalty` is `Some`, only features whose name contains `_pen{penalty}`
/// are kept, used for VCR BED files produced by `scFindVCRs` with a range of
/// penalties encoded in the feature name.
fn parse_annotation_features(
    path: &Path,
    penalty: Option<f64>,
    feature_types: Option<&[&str]>,
    exon_types: Option<&[&str]>,
    name_attr: Option<&str>,
    metagene: bool,
) -> Result<Vec<Feature>> {
    let (_index, var) =
        parse_annotation_files([path], feature_types, exon_types, name_attr, metagene)?;
    let penalty_suffix = penalty.map(|p| format!("_pen{}", p));

    Ok(var
        .into_iter()
        .filter(|v| penalty_suffix.as_deref().is_none_or(|s| v.name.contains(s)))
        .collect())
}

/// Compute the score vector (length `n_obs`) for a single feature, or `None`
/// when no overlapping bin carries any signal.
///
/// The score of a cell is the sum of its counts in the bins overlapping the
/// feature. `overlap_policy` decides what to do with a bin that only partly
/// overlaps it:
///
/// * `"all"`: count the bin in full;
/// * `"partial"`: count the fraction of the bin that lies inside the feature;
/// * `"none"`: ignore the bin unless it is wholly inside the feature.
fn compute_feature_score(
    feat: &Feature,
    bin_chrom_idx: &AHashMap<String, ChromIndex>,
    bin_starts: &[i64],
    bin_ends: &[i64],
    x_csc: &CscMatrix<f64>,
    n_obs: usize,
    overlap_policy: &str,
) -> Option<Vec<f32>> {
    let chrom_idx = bin_chrom_idx.get(&feat.chrom)?;

    let feat_start = feat.start;
    let feat_end = feat.end;

    // Bins overlapping the feature, and the weight each contributes with.
    let mut bins: Vec<(usize, f64)> = Vec::new();
    for iv in chrom_idx.find(feat_start, feat_end) {
        let bin_start = bin_starts[iv.var_idx] as f64;
        let bin_end = bin_ends[iv.var_idx] as f64;

        let weight = match overlap_policy {
            "none" => {
                if bin_start < feat_start as f64 || bin_end > feat_end as f64 {
                    continue;
                }
                1.0
            }
            "partial" => {
                let inside = bin_end.min(feat_end as f64) - bin_start.max(feat_start as f64);
                let bin_len = (bin_end - bin_start).max(1.0);
                inside / bin_len
            }
            _ => 1.0,
        };
        bins.push((iv.var_idx, weight));
    }

    if bins.is_empty() {
        return None;
    }

    // Weighted column sum: score[cell] += X[cell, bin] * weight[bin].
    let mut score = vec![0.0f64; n_obs];
    for (col, weight) in bins {
        if weight == 0.0 {
            continue;
        }
        let c = x_csc.col(col);
        for (&row, &val) in c.row_indices().iter().zip(c.values().iter()) {
            score[row] += val * weight;
        }
    }

    // No counts anywhere in the feature: nothing to record.
    if score.iter().all(|&v| v == 0.0) {
        return None;
    }
    Some(score.into_iter().map(|v| v as f32).collect())
}

/// Build a CSR matrix (cells × features) from the per-feature score vectors.
///
/// `scores[j]` is `Some(vec)` for feature `j` when it has signal, or `None`,
/// in which case the feature still occupies an all-zero column.
///
/// Returns `(row_offsets, col_indices, values)` in CSR form.
fn build_score_csr(
    scores: &[Option<Vec<f32>>],
    n_obs: usize,
) -> (Vec<usize>, Vec<usize>, Vec<f32>) {
    // Transpose from feature-major to cell-major.
    let mut coo: Vec<(usize, usize, f32)> = Vec::new();
    for (j, score) in scores.iter().enumerate() {
        let Some(v) = score else { continue };
        for (i, &val) in v.iter().enumerate() {
            if val != 0.0 {
                coo.push((i, j, val));
            }
        }
    }
    coo.sort_unstable_by_key(|&(r, c, _)| (r, c));

    let mut row_offsets = vec![0usize; n_obs + 1];
    for &(r, _, _) in &coo {
        row_offsets[r + 1] += 1;
    }
    for i in 0..n_obs {
        row_offsets[i + 1] += row_offsets[i];
    }

    let col_indices = coo.iter().map(|&(_, c, _)| c).collect();
    let values = coo.iter().map(|&(_, _, v)| v).collect();

    (row_offsets, col_indices, values)
}

// Main entry point
fn run_score_features(
    h5ad_path: &Path,
    features_path: &Path,
    out_path: &Path,
    overlap_policy: &str,
    penalty: Option<f64>,
    feature_types: Option<&[&str]>,
    exon_types: Option<&[&str]>,
    name_attr: Option<&str>,
    metagene: bool,
    num_threads: usize,
) -> Result<()> {
    if overlap_policy != "partial" && overlap_policy != "all" && overlap_policy != "none" {
        bail!(
            "overlap_policy must be 'partial', 'all', or 'none', got {:?}",
            overlap_policy
        );
    }

    // 1. Parse annotation features.
    eprintln!("Parsing annotation file {}…", features_path.display());
    let mut features = parse_annotation_features(
        features_path,
        penalty,
        feature_types,
        exon_types,
        name_attr,
        metagene,
    )?;
    if features.is_empty() {
        bail!("no features found in {}", features_path.display());
    }

    // 2. Read adata.
    eprintln!("Reading AnnData {}…", h5ad_path.display());
    let store =
        H5::open(h5ad_path).with_context(|| format!("failed to open {}", h5ad_path.display()))?;
    let adata = AnnData::<H5>::open(store)?;

    let var = adata.read_var()?;
    let bin_chrom = df_str_col(&var, "chrom")?;
    let bin_start = df_i64_col(&var, "start")?;
    let bin_end = df_i64_col(&var, "end")?;
    let n_obs = adata.n_obs();
    let obs_names = adata.obs_names();
    let obs_df = adata.read_obs()?;
    let x = read_x_f64(&adata)?;
    adata.close()?;
    let x_csc = CscMatrix::from(&x);

    // 3. Filter features to chromosomes present in adata.var.
    let adata_chroms: AHashSet<&str> = bin_chrom.iter().map(|s| s.as_str()).collect();
    features.retain(|f| adata_chroms.contains(f.chrom.as_str()));
    if features.is_empty() {
        bail!("no common chromosomes between adata and annotation file");
    }

    // Remove duplicate feature names (keep first occurrence, matching Python).
    {
        let mut seen: AHashSet<String> = AHashSet::new();
        features.retain(|f| seen.insert(f.name.clone()));
    }
    eprintln!("Processing {} features…", features.len());

    // 4. Build bin interval index (bin col index per chromosome).
    let mut bin_by_chrom: AHashMap<String, Vec<Interval>> = AHashMap::new();
    for (col, chrom) in bin_chrom.iter().enumerate() {
        bin_by_chrom
            .entry(chrom.clone())
            .or_default()
            .push(Interval {
                start: bin_start[col] as usize,
                end: bin_end[col] as usize,
                var_idx: col,
            });
    }
    let bin_chrom_idx: AHashMap<String, ChromIndex> = bin_by_chrom
        .into_iter()
        .map(|(c, ivs)| (c, ChromIndex::build(ivs)))
        .collect();

    // 5. Thread pool.
    let n_threads = if num_threads == 0 {
        rayon::current_num_threads()
    } else {
        num_threads
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .context("failed to build thread pool")?;

    // 6. Score each feature in parallel.
    let scores: Vec<Option<Vec<f32>>> = pool.install(|| {
        features
            .par_iter()
            .map(|feat| {
                compute_feature_score(
                    feat,
                    &bin_chrom_idx,
                    &bin_start,
                    &bin_end,
                    &x_csc,
                    n_obs,
                    overlap_policy,
                )
            })
            .collect()
    });

    let n_features = features.len();
    eprintln!(
        "Scored features: {}/{} have non-zero signal",
        scores.iter().filter(|s| s.is_some()).count(),
        n_features
    );

    // 7. Build and write the output matrix.
    let (row_offsets, col_indices, values) = build_score_csr(&scores, n_obs);
    let mat = CsrMatrix::try_from_csr_data(n_obs, n_features, row_offsets, col_indices, values)
        .map_err(|e| anyhow::anyhow!("failed to build CSR: {:?}", e))?;
    write_output_anndata(out_path, mat, obs_names, obs_df, &features)?;
    Ok(())
}

fn write_output_anndata(
    out_path: &Path,
    mat: CsrMatrix<f32>,
    obs_names: DataFrameIndex,
    obs_df: DataFrame,
    features: &[Feature],
) -> Result<()> {
    let n_genes = features.len();
    let gene_names: Vec<String> = features.iter().map(|f| f.name.clone()).collect();
    let chrom_col: Vec<String> = features.iter().map(|f| f.chrom.clone()).collect();
    let start_col: Vec<i64> = features.iter().map(|f| f.start as i64).collect();
    let end_col: Vec<i64> = features.iter().map(|f| f.end as i64).collect();
    let strand_col: Vec<String> = features.iter().map(|f| f.strand.to_string()).collect();

    let var_df = DataFrame::new(
        n_genes,
        vec![
            Column::new("gene_name".into(), gene_names.clone()),
            Column::new("chrom".into(), chrom_col),
            Column::new("start".into(), start_col),
            Column::new("end".into(), end_col),
            Column::new("strand".into(), strand_col),
        ],
    )?;

    let out = AnnData::<H5>::new(out_path)
        .with_context(|| format!("failed to create {}", out_path.display()))?;
    out.set_x(mat)?;
    out.set_obs_names(obs_names)?;
    out.set_obs(obs_df)?;
    out.set_var_names(gene_names.into_iter().collect())?;
    out.set_var(var_df)?;
    out.close()?;

    eprintln!("Output written to {}", out_path.display());
    Ok(())
}

/// Aggregate binned single-cell chromatin data into per-feature scores.
///
/// Reads `h5ad_path` (an AnnData file produced by `count_bins` or similar),
/// whose `var` carries `chrom`/`start`/`end` bin coordinates. For each feature
/// in `features_path` (BED/GTF/GFF), sums the counts of the overlapping bins to
/// produce a per-cell score, weighting partly-overlapping bins according to
/// `overlap_policy`.
///
/// Returns `out_path` on success.
#[pyfunction(signature = (
    h5ad_path,
    features_path,
    out_path,
    overlap_policy = String::from("partial"),
    penalty = None,
    feature_type = None,
    exon_type = None,
    name_attr = None,
    metagene = false,
    num_threads = 0,
))]
#[allow(clippy::too_many_arguments)]
pub fn score_features(
    h5ad_path: PathBuf,
    features_path: PathBuf,
    out_path: PathBuf,
    overlap_policy: String,
    penalty: Option<f64>,
    feature_type: Option<Vec<String>>,
    exon_type: Option<Vec<String>>,
    name_attr: Option<String>,
    metagene: bool,
    num_threads: usize,
) -> PyResult<String> {
    let feature_types: Option<Vec<&str>> = feature_type
        .as_ref()
        .map(|t| t.iter().map(String::as_str).collect());
    let exon_types: Option<Vec<&str>> = exon_type
        .as_ref()
        .map(|t| t.iter().map(String::as_str).collect());
    run_score_features(
        &h5ad_path,
        &features_path,
        &out_path,
        &overlap_policy,
        penalty,
        feature_types.as_deref(),
        exon_types.as_deref(),
        name_attr.as_deref(),
        metagene,
        num_threads,
    )
    .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?;
    Ok(out_path.display().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn feature(chrom: &str, start: usize, end: usize) -> Feature {
        Feature {
            chrom: chrom.to_string(),
            start,
            end,
            name: format!("{chrom}:{start}-{end}"),
            strand: '*',
        }
    }

    /// Three adjacent 100 bp bins on `chr1`, as columns 0, 1 and 2.
    fn bin_layout() -> (AHashMap<String, ChromIndex>, Vec<i64>, Vec<i64>) {
        let intervals = (0..3)
            .map(|i| Interval {
                start: i * 100,
                end: (i + 1) * 100,
                var_idx: i,
            })
            .collect();
        let mut index = AHashMap::new();
        index.insert("chr1".to_string(), ChromIndex::build(intervals));
        (index, vec![0, 100, 200], vec![100, 200, 300])
    }

    /// Dense (cells x bins) matrix as CSC, which is how the scorer reads it.
    fn csc(rows: usize, cols: usize, entries: &[(usize, usize, f64)]) -> CscMatrix<f64> {
        let coo = nalgebra_sparse::CooMatrix::try_from_triplets(
            rows,
            cols,
            entries.iter().map(|e| e.0).collect(),
            entries.iter().map(|e| e.1).collect(),
            entries.iter().map(|e| e.2).collect(),
        )
        .unwrap();
        CscMatrix::from(&coo)
    }

    // compute_feature_score

    #[test]
    fn a_feature_sums_the_counts_of_every_bin_it_covers() {
        let (index, starts, ends) = bin_layout();
        // Cell 0 has 1 in each of the three bins; cell 1 has 10 in bin 1 only.
        let x = csc(2, 3, &[(0, 0, 1.0), (0, 1, 1.0), (0, 2, 1.0), (1, 1, 10.0)]);

        let score = compute_feature_score(
            &feature("chr1", 0, 300),
            &index,
            &starts,
            &ends,
            &x,
            2,
            "all",
        )
        .expect("the feature covers bins that carry signal");

        assert_eq!(score, vec![3.0, 10.0]);
    }

    #[test]
    fn the_all_policy_counts_a_partly_covered_bin_in_full() {
        let (index, starts, ends) = bin_layout();
        let x = csc(1, 3, &[(0, 0, 4.0), (0, 1, 4.0)]);

        // The feature covers bin 0 fully and half of bin 1.
        let score = compute_feature_score(
            &feature("chr1", 0, 150),
            &index,
            &starts,
            &ends,
            &x,
            1,
            "all",
        )
        .unwrap();

        assert_eq!(score, vec![8.0]);
    }

    #[test]
    fn the_partial_policy_weights_a_bin_by_the_fraction_inside_the_feature() {
        let (index, starts, ends) = bin_layout();
        let x = csc(1, 3, &[(0, 0, 4.0), (0, 1, 4.0)]);

        // Bin 0 lies wholly inside; bin 1 contributes half of its 4.
        let score = compute_feature_score(
            &feature("chr1", 0, 150),
            &index,
            &starts,
            &ends,
            &x,
            1,
            "partial",
        )
        .unwrap();

        assert_eq!(score, vec![6.0]);
    }

    #[test]
    fn the_none_policy_drops_a_bin_that_is_not_wholly_inside_the_feature() {
        let (index, starts, ends) = bin_layout();
        let x = csc(1, 3, &[(0, 0, 4.0), (0, 1, 4.0)]);

        let score = compute_feature_score(
            &feature("chr1", 0, 150),
            &index,
            &starts,
            &ends,
            &x,
            1,
            "none",
        )
        .unwrap();

        // Only bin 0 survives.
        assert_eq!(score, vec![4.0]);
    }

    #[test]
    fn a_feature_on_an_unknown_chromosome_scores_nothing() {
        let (index, starts, ends) = bin_layout();
        let x = csc(1, 3, &[(0, 0, 4.0)]);

        assert!(
            compute_feature_score(
                &feature("chrZ", 0, 300),
                &index,
                &starts,
                &ends,
                &x,
                1,
                "all"
            )
            .is_none()
        );
    }

    #[test]
    fn a_feature_overlapping_no_bin_scores_nothing() {
        let (index, starts, ends) = bin_layout();
        let x = csc(1, 3, &[(0, 0, 4.0)]);

        assert!(
            compute_feature_score(
                &feature("chr1", 5_000, 6_000),
                &index,
                &starts,
                &ends,
                &x,
                1,
                "all"
            )
            .is_none()
        );
    }

    #[test]
    fn a_feature_whose_bins_are_all_empty_scores_nothing() {
        let (index, starts, ends) = bin_layout();
        // Signal exists, but not in the bins this feature covers.
        let x = csc(1, 3, &[(0, 2, 9.0)]);

        assert!(
            compute_feature_score(
                &feature("chr1", 0, 100),
                &index,
                &starts,
                &ends,
                &x,
                1,
                "all"
            )
            .is_none()
        );
    }

    // build_score_csr

    #[test]
    fn the_csr_offsets_count_the_entries_of_each_row() {
        // Two cells, three features; feature 1 has no signal at all.
        let scores = vec![Some(vec![1.0f32, 0.0]), None, Some(vec![2.0f32, 3.0])];

        let (offsets, cols, values) = build_score_csr(&scores, 2);

        // Cell 0 has features 0 and 2; cell 1 has feature 2 only.
        assert_eq!(offsets, vec![0, 2, 3]);
        assert_eq!(cols, vec![0, 2, 2]);
        assert_eq!(values, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn a_feature_without_signal_still_occupies_its_column() {
        let scores = vec![None, Some(vec![5.0f32])];
        let (offsets, cols, values) = build_score_csr(&scores, 1);

        assert_eq!(offsets, vec![0, 1]);
        assert_eq!(cols, vec![1], "the surviving entry keeps column index 1");
        assert_eq!(values, vec![5.0]);
    }

    #[test]
    fn scores_of_zero_are_left_out_of_the_sparse_matrix() {
        let scores = vec![Some(vec![0.0f32, 7.0])];
        let (offsets, cols, values) = build_score_csr(&scores, 2);

        assert_eq!(offsets, vec![0, 0, 1], "cell 0 contributes nothing");
        assert_eq!(cols, vec![0]);
        assert_eq!(values, vec![7.0]);
    }

    #[test]
    fn an_empty_score_set_yields_an_all_zero_matrix() {
        let (offsets, cols, values) = build_score_csr(&[], 3);

        assert_eq!(offsets, vec![0, 0, 0, 0]);
        assert!(cols.is_empty());
        assert!(values.is_empty());
    }

    #[test]
    fn csr_entries_come_out_ordered_by_cell_then_feature() {
        let scores = vec![
            Some(vec![1.0f32, 1.0]),
            Some(vec![1.0f32, 1.0]),
            Some(vec![1.0f32, 1.0]),
        ];
        let (offsets, cols, _) = build_score_csr(&scores, 2);

        assert_eq!(offsets, vec![0, 3, 6]);
        assert_eq!(cols, vec![0, 1, 2, 0, 1, 2]);
    }
}
