//! Gene and region activity scoring from binned single-cell chromatin data.
//!
//! Ports `FeatureScorer` / `scScoreFeatures` from the original Python `sincei`
//! package. Two modes:
//!
//! - **aggregate**: simple weighted sum of bin counts within each BED/GTF
//!   feature, with uniform weight 1 and no extension beyond the feature body.
//! - **activities**: exponential-decay weighted sum extending `max_region` kb
//!   upstream and downstream of each feature.
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

use crate::counting::count_utils::{df_i64_col, df_str_col, read_x_f64};
use crate::counting::parse_annotation::parse_annotation_files;
use crate::counting::region_index::{ChromIndex, Interval, VarMeta};

/// One genomic feature from the annotation file (BED/GTF/GFF).
#[derive(Clone)]
struct Feature {
    name: String,
    chrom: String,
    start: i64,
    end: i64,
    strand: char,
}

impl From<VarMeta> for Feature {
    fn from(v: VarMeta) -> Self {
        Feature {
            name: v.name,
            chrom: v.chrom,
            start: v.start as i64,
            end: v.end as i64,
            strand: v.strand,
        }
    }
}

/// Parse an annotation file (BED / GTF / GFF3) into [`Feature`]s, reusing the
/// shared `parse_annotation_files` machinery (and thus the same format
/// detection as `scCountReads`).
///
/// `feature_type` selects the column-3 type kept as a region (e.g. `"gene"`,
/// `"transcript"`); `name_attr` is the column-9 attribute used to name each
/// feature (e.g. `"gene_id"`). `None` falls back to the per-format defaults.
/// When `metagene` is set, exon rows of `feature_type` are merged per gene
/// (grouped by `name_attr`) so each output feature is one gene. All three only
/// affect GTF/GFF inputs and are ignored for BED.
///
/// When `penalty` is `Some`, only features whose name contains `_pen{penalty}`
/// are kept — for VCR BED files produced by `scFindVCRs` with a range of
/// penalties encoded in the feature name.
fn parse_annotation_features(
    path: &Path,
    penalty: Option<f64>,
    feature_type: Option<&str>,
    name_attr: Option<&str>,
    metagene: bool,
) -> Result<Vec<Feature>> {
    let (_index, var) = parse_annotation_files([path], feature_type, name_attr, metagene)?;
    let penalty_suffix = penalty.map(|p| format!("_pen{}", p));

    Ok(var
        .into_iter()
        .filter(|v| penalty_suffix.as_deref().is_none_or(|s| v.name.contains(s)))
        .map(Feature::from)
        .collect())
}

/// Build a per-chromosome `ChromIndex` from a slice of features, mapping
/// each feature's genomic interval to its index in the slice.
fn build_feature_index(features: &[Feature]) -> AHashMap<String, ChromIndex> {
    let mut by_chrom: AHashMap<String, Vec<Interval>> = AHashMap::new();
    for (i, f) in features.iter().enumerate() {
        by_chrom.entry(f.chrom.clone()).or_default().push(Interval {
            start: f.start as usize,
            end: f.end as usize,
            val: i,
            name: f.name.clone(),
        });
    }
    by_chrom
        .into_iter()
        .map(|(c, ivs)| (c, ChromIndex::build(ivs)))
        .collect()
}

/// Compute the decay weight for a single bin `[fs, fe)` relative to a gene
/// `[gs, ge)`, porting `get_decay_weights` from `FeatureScorer.py` exactly.
///
/// `lam = decay / 10_000` (decay per bp).  When `lam == 0` all weights are 1.
///
/// When `gene_body = true` the gene body is weighted as a constant 1 and decay
/// starts beyond its boundary.  When `gene_body = false` decay starts from the
/// TSS (which is `gs` for '+' and `ge` for '-'/`*`).
fn decay_weight(
    gs: f64,
    ge: f64,
    strand: char,
    fs: f64,
    fe: f64,
    lam: f64,
    gene_body: bool,
) -> f64 {
    if lam == 0.0 {
        return 1.0;
    }

    if gene_body {
        let overlap_s = fs.max(gs);
        let overlap_e = fe.min(ge);

        if overlap_s >= overlap_e {
            // Feature has no overlap with the gene body.
            if fe <= gs {
                // Upstream: integral of exp(-lam*(gs-x)) over [fs, fe]
                (-(lam * (gs - fe))).exp() - (-(lam * (gs - fs))).exp()
            } else {
                // Downstream: integral of exp(-lam*(x-ge)) over [fs, fe]
                (-(lam * (fs - ge))).exp() - (-(lam * (fe - ge))).exp()
            }
        } else {
            let mut w = 0.0;
            // Upstream flap before gene start
            if fs < gs {
                w += 1.0 - (-(lam * (gs - fs))).exp();
            }
            // Inside gene body: raw overlap length
            w += overlap_e - overlap_s;
            // Downstream flap after gene end
            if fe > ge {
                w += 1.0 - (-(lam * (fe - ge))).exp();
            }
            w
        }
    } else {
        // Decay from TSS
        let tss = if strand == '-' { ge } else { gs };
        if fe <= tss {
            (-(lam * (tss - fe))).exp() - (-(lam * (tss - fs))).exp()
        } else if fs >= tss {
            (-(lam * (fs - tss))).exp() - (-(lam * (fe - tss))).exp()
        } else {
            // Feature straddles TSS
            2.0 - (-(lam * (tss - fs))).exp() - (-(lam * (fe - tss))).exp()
        }
    }
}

/// Compute the activity vector (length `n_obs`) for a single feature, or
/// `None` when no overlapping bins exist.
///
/// Returns the raw (un-normalized) float32 activity per cell.
fn compute_activity(
    feat: &Feature,
    bin_chrom_idx: &AHashMap<String, ChromIndex>,
    bin_starts: &[i64],
    bin_ends: &[i64],
    x_csc: &CscMatrix<f64>,
    n_obs: usize,
    lam: f64,
    max_region_bp: i64,
    gene_body: bool,
    normalize_gene_lengths: bool,
    overlap_policy: &str,
    excl_kind: Option<&str>,
    gene_chrom_idx: Option<&AHashMap<String, ChromIndex>>,
    all_genes: &[Feature],
) -> Option<Vec<f32>> {
    let chrom_idx = bin_chrom_idx.get(&feat.chrom)?;

    let region_start = (feat.start - max_region_bp).max(0) as usize;
    let region_end = (feat.end + max_region_bp) as usize;

    // Find bins overlapping the extended search region.
    let mut bin_cols: Vec<usize> = Vec::new();
    let mut bin_fs: Vec<f64> = Vec::new();
    let mut bin_fe: Vec<f64> = Vec::new();

    for iv in chrom_idx.find(region_start, region_end) {
        let fs = bin_starts[iv.val] as f64;
        let fe = bin_ends[iv.val] as f64;

        // "none" policy: only fully-contained bins
        if overlap_policy == "none" && (fs < region_start as f64 || fe > region_end as f64) {
            continue;
        }

        bin_cols.push(iv.val);
        bin_fs.push(fs);
        bin_fe.push(fe);
    }

    if bin_cols.is_empty() {
        return None;
    }

    let gs = feat.start as f64;
    let ge = feat.end as f64;

    // Compute per-bin decay weights.
    let mut weights: Vec<f64> = bin_fs
        .iter()
        .zip(bin_fe.iter())
        .map(|(&fs, &fe)| decay_weight(gs, ge, feat.strand, fs, fe, lam, gene_body))
        .collect();

    // Scale weights for partially-overlapping bins ("partial" policy).
    if overlap_policy == "partial" {
        let rgs = region_start as f64;
        let rge = region_end as f64;
        for (w, (&fs, &fe)) in weights.iter_mut().zip(bin_fs.iter().zip(bin_fe.iter())) {
            let clip_s = fs.max(rgs);
            let clip_e = fe.min(rge);
            let feat_len = (fe - fs).max(1.0);
            *w *= (clip_e - clip_s) / feat_len;
        }
    }

    // Gene-size factor: multiply all weights by gene length.
    if normalize_gene_lengths {
        let gene_len = ge - gs;
        for w in &mut weights {
            *w *= gene_len;
        }
    }

    // Exclude bins overlapping TSS / bodies of other genes.
    if let (Some(kind), Some(Some(cidx))) = (excl_kind, gene_chrom_idx.map(|m| m.get(&feat.chrom)))
    {
        // Build the set of excluded intervals from other genes in the region.
        let mut excluded: Vec<(f64, f64)> = Vec::new();
        for iv in cidx.find(region_start, region_end) {
            if all_genes[iv.val].name == feat.name {
                continue;
            }
            let og = &all_genes[iv.val];
            if kind == "TSS" {
                let tss = if og.strand == '-' {
                    og.end as f64
                } else {
                    og.start as f64
                };
                excluded.push((tss, tss));
            } else {
                excluded.push((og.start as f64, og.end as f64));
            }
        }
        // Zero-out weights of bins overlapping any excluded interval.
        for (w, (&fs, &fe)) in weights.iter_mut().zip(bin_fs.iter().zip(bin_fe.iter())) {
            for &(es, ee) in &excluded {
                if fs < ee && fe > es {
                    *w = 0.0;
                    break;
                }
            }
        }
    }

    // Compute weighted column sum: activity[cell] += X[cell, bin] * weight[bin]
    let mut activity = vec![0.0f64; n_obs];
    for ((&col, &fs), &w) in bin_cols.iter().zip(bin_fs.iter()).zip(weights.iter()) {
        let _ = fs; // used only for overlap_policy above
        if w == 0.0 {
            continue;
        }
        let c = x_csc.col(col);
        for (&row, &val) in c.row_indices().iter().zip(c.values().iter()) {
            activity[row] += val * w;
        }
    }

    // Return None when the entire activity vector is zero (no counts in region).
    if activity.iter().all(|&v| v == 0.0) {
        return None;
    }
    Some(activity.into_iter().map(|v| v as f32).collect())
}

/// L1-normalize the rows of a CSR matrix in-place: each cell's activity
/// vector sums to 1 after normalization.
fn l1_normalize_rows(n_obs: usize, row_offsets: &[usize], values: &mut [f32]) {
    for i in 0..n_obs {
        let start = row_offsets[i];
        let end = row_offsets[i + 1];
        let row_sum: f32 = values[start..end].iter().map(|v| v.abs()).sum();
        if row_sum > 0.0 {
            for v in &mut values[start..end] {
                *v /= row_sum;
            }
        }
    }
}

/// Build a CSR matrix from collected per-gene activity vectors.
///
/// `activities[j]` is `Some(Vec<f32>)` for gene `j` when there is signal, or
/// `None` to skip (the gene is still included in the output column count as
/// an all-zero column).
///
/// Returns `(row_offsets, col_indices, values, n_obs, n_cols)` in CSR form.
fn build_activity_csr(
    activities: &[Option<Vec<f32>>],
    n_obs: usize,
) -> (Vec<usize>, Vec<usize>, Vec<f32>) {
    // Transpose from col-major (gene-first) to row-major (cell-first).
    let n_genes = activities.len();
    // Collect (row, col, val) triples.
    let mut coo: Vec<(usize, usize, f32)> = Vec::new();
    for (j, act) in activities.iter().enumerate() {
        let Some(v) = act else { continue };
        for (i, &val) in v.iter().enumerate() {
            if val != 0.0 {
                coo.push((i, j, val));
            }
        }
    }
    coo.sort_unstable_by_key(|&(r, c, _)| (r, c));

    let mut row_offsets = vec![0usize; n_obs + 1];
    let mut col_indices = Vec::with_capacity(coo.len());
    let mut values = Vec::with_capacity(coo.len());

    for &(r, c, v) in &coo {
        row_offsets[r + 1] = row_offsets[r + 1].max(col_indices.len() + 1);
        col_indices.push(c);
        values.push(v);
    }
    // Fix offsets to be properly cumulative.
    let mut nnz = 0usize;
    for i in 0..n_obs {
        let row_nnz = row_offsets[i + 1];
        row_offsets[i + 1] = nnz + row_nnz;
        nnz += row_nnz;
    }
    // Rebuild using sorted COO directly.
    let mut row_offsets2 = vec![0usize; n_obs + 1];
    for &(r, _, _) in &coo {
        row_offsets2[r + 1] += 1;
    }
    for i in 0..n_obs {
        row_offsets2[i + 1] += row_offsets2[i];
    }
    let col_indices2: Vec<usize> = coo.iter().map(|&(_, c, _)| c).collect();
    let values2: Vec<f32> = coo.iter().map(|&(_, _, v)| v).collect();

    let _ = (n_genes, row_offsets, col_indices, values);
    (row_offsets2, col_indices2, values2)
}

/// Apply column-wise z-score normalization (center + unit variance) to a CSR
/// matrix returned as `(row_offsets, col_indices, values)`.
///
/// This produces a dense-like result (most entries nonzero after centering).
fn center_and_scale(
    n_obs: usize,
    n_genes: usize,
    row_offsets: Vec<usize>,
    col_indices: Vec<usize>,
    values: Vec<f32>,
) -> (Vec<usize>, Vec<usize>, Vec<f32>) {
    // Materialize as dense column-major to compute per-gene mean and std.
    let mut dense = vec![0.0f32; n_obs * n_genes];
    for i in 0..n_obs {
        for idx in row_offsets[i]..row_offsets[i + 1] {
            let j = col_indices[idx];
            dense[i * n_genes + j] = values[idx];
        }
    }

    for j in 0..n_genes {
        let mut sum = 0.0f64;
        for i in 0..n_obs {
            sum += dense[i * n_genes + j] as f64;
        }
        let mean = sum / n_obs as f64;
        let mut var = 0.0f64;
        for i in 0..n_obs {
            let d = dense[i * n_genes + j] as f64 - mean;
            var += d * d;
        }
        let std = (var / n_obs as f64).sqrt();
        for i in 0..n_obs {
            let v = (dense[i * n_genes + j] as f64 - mean) / if std > 0.0 { std } else { 1.0 };
            dense[i * n_genes + j] = v as f32;
        }
    }

    // Convert back to CSR (all nonzero after centering, except exact zeros).
    let mut row_off = vec![0usize; n_obs + 1];
    let mut col_idx = Vec::with_capacity(n_obs * n_genes);
    let mut vals = Vec::with_capacity(n_obs * n_genes);

    for i in 0..n_obs {
        for j in 0..n_genes {
            let v = dense[i * n_genes + j];
            if v != 0.0 {
                col_idx.push(j);
                vals.push(v);
            }
        }
        row_off[i + 1] = col_idx.len();
    }

    (row_off, col_idx, vals)
}

// Main entry point
fn run_score_features(
    h5ad_path: &Path,
    features_path: &Path,
    out_path: &Path,
    mode: &str,
    overlap_policy: &str,
    center_scores: bool,
    penalty: Option<f64>,
    decay: f64,
    max_region: i64,
    gene_body: bool,
    normalize_gene_lengths: bool,
    exclude_in_range: Option<&str>,
    feature_type: Option<&str>,
    name_attr: Option<&str>,
    metagene: bool,
    num_threads: usize,
) -> Result<()> {
    if mode != "aggregate" && mode != "activities" {
        bail!("mode must be 'aggregate' or 'activities', got {:?}", mode);
    }
    if overlap_policy != "partial" && overlap_policy != "all" && overlap_policy != "none" {
        bail!(
            "overlap_policy must be 'partial', 'all', or 'none', got {:?}",
            overlap_policy
        );
    }
    if matches!(exclude_in_range, Some(excl) if excl != "TSS" && excl != "genes") {
        bail!(
            "exclude_in_range must be 'TSS' or 'genes', got {:?}",
            exclude_in_range.unwrap()
        );
    }

    // Effective parameters for aggregate mode (override caller args).
    let (eff_lam, eff_max_region_bp, eff_gene_body, eff_norm_len, eff_excl) = if mode == "aggregate"
    {
        (0.0f64, 0i64, true, false, None)
    } else {
        (
            decay / 10_000.0,
            max_region * 1000,
            gene_body,
            normalize_gene_lengths,
            exclude_in_range,
        )
    };

    // 1. Parse annotation features.
    eprintln!("Parsing annotation file {}…", features_path.display());
    let mut features =
        parse_annotation_features(features_path, penalty, feature_type, name_attr, metagene)?;
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
                val: col,
                name: String::new(),
            });
    }
    let bin_chrom_idx: AHashMap<String, ChromIndex> = bin_by_chrom
        .into_iter()
        .map(|(c, ivs)| (c, ChromIndex::build(ivs)))
        .collect();

    // 5. Build gene interval index (for exclude_in_range).
    let gene_chrom_idx = eff_excl.map(|_| build_feature_index(&features));

    // 6. Thread pool.
    let n_threads = if num_threads == 0 {
        rayon::current_num_threads()
    } else {
        num_threads
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .context("failed to build thread pool")?;

    // 7. Compute per-feature activities in parallel.
    let activities: Vec<Option<Vec<f32>>> = pool.install(|| {
        features
            .par_iter()
            .map(|feat| {
                compute_activity(
                    feat,
                    &bin_chrom_idx,
                    &bin_start,
                    &bin_end,
                    &x_csc,
                    n_obs,
                    eff_lam,
                    eff_max_region_bp,
                    eff_gene_body,
                    eff_norm_len,
                    overlap_policy,
                    eff_excl,
                    gene_chrom_idx.as_ref(),
                    &features,
                )
            })
            .collect()
    });

    let n_genes = features.len();
    eprintln!(
        "Computed activities: {}/{} features have non-zero signal",
        activities.iter().filter(|a| a.is_some()).count(),
        n_genes
    );

    // 8. Build output CSR matrix + apply normalization.
    let (row_offsets, col_indices, mut values) = build_activity_csr(&activities, n_obs);

    if mode == "activities" {
        if center_scores {
            let (ro, ci, va) = center_and_scale(n_obs, n_genes, row_offsets, col_indices, values);
            let mat = CsrMatrix::try_from_csr_data(n_obs, n_genes, ro, ci, va)
                .map_err(|e| anyhow::anyhow!("failed to build CSR: {:?}", e))?;
            write_output_anndata(out_path, mat, obs_names, obs_df, &features)?;
            return Ok(());
        }
        l1_normalize_rows(n_obs, &row_offsets, &mut values);
    }

    let mat = CsrMatrix::try_from_csr_data(n_obs, n_genes, row_offsets, col_indices, values)
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
    let start_col: Vec<i64> = features.iter().map(|f| f.start).collect();
    let end_col: Vec<i64> = features.iter().map(|f| f.end).collect();
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

/// Compute gene/region activity scores from binned single-cell chromatin data.
///
/// Reads `h5ad_path` (an AnnData file produced by `count_bins` or similar),
/// whose `var` carries `chrom`/`start`/`end` bin coordinates.  For each
/// feature in `features_path` (BED/GTF/GFF), computes the weighted sum of
/// bin counts (cells × bins) to produce a per-cell activity score.
///
/// `mode="aggregate"` sums raw counts within each feature with uniform weight.
/// `mode="activities"` uses exponential-decay weights extending `max_region`
/// kb up/downstream of each feature and L1-normalizes the result per cell.
///
/// Returns `out_path` on success.
#[pyfunction(signature = (
    h5ad_path,
    features_path,
    out_path,
    mode,
    overlap_policy = String::from("partial"),
    center_scores = false,
    penalty = None,
    decay = 0.75,
    max_region = 100,
    gene_body = false,
    normalize_gene_lengths = false,
    exclude_in_range = None,
    feature_type = None,
    name_attr = None,
    metagene = false,
    num_threads = 0,
))]
#[allow(clippy::too_many_arguments)]
pub fn score_features(
    h5ad_path: PathBuf,
    features_path: PathBuf,
    out_path: PathBuf,
    mode: String,
    overlap_policy: String,
    center_scores: bool,
    penalty: Option<f64>,
    decay: f64,
    max_region: i64,
    gene_body: bool,
    normalize_gene_lengths: bool,
    exclude_in_range: Option<String>,
    feature_type: Option<String>,
    name_attr: Option<String>,
    metagene: bool,
    num_threads: usize,
) -> PyResult<String> {
    run_score_features(
        &h5ad_path,
        &features_path,
        &out_path,
        &mode,
        &overlap_policy,
        center_scores,
        penalty,
        decay,
        max_region,
        gene_body,
        normalize_gene_lengths,
        exclude_in_range.as_deref(),
        feature_type.as_deref(),
        name_attr.as_deref(),
        metagene,
        num_threads,
    )
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    Ok(out_path.display().to_string())
}
