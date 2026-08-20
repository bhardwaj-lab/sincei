//! Variable chromatin region (VCR) detection.
//!
//! Ports the `VCRfinder` algorithm from the original Python `sincei`
//! package (`VCRfinder.py` / `scFindVCRs.py`). For each chromosome:
//!
//! 1. Compute a banded bin-to-bin Pearson correlation matrix of the binned
//!    signal (only the first `k` diagonals on each side of the main
//!    diagonal are kept).
//! 2. Convolve that banded matrix with a family of Gaussian kernels of
//!    increasing scale to build a per-bin score map (one score per kernel
//!    per bin).
//! 3. Run PELT changepoint detection on the score map, once per penalty
//!    value, to call segment boundaries.
//!
//! The changepoint step is mathematically equivalent to
//! `ruptures.KernelCPD(kernel="linear", min_size=1, jump=1)` (no Rust crate
//! covers this, so it's implemented directly as classic PELT with an L2
//! mean-shift cost).
use ahash::AHashMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use anndata::{AnnData, AnnDataOp, Backend};
use anndata_hdf5::H5;
use anyhow::{Context, Result, bail};
use nalgebra_sparse::CscMatrix;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::counting::count_utils::{df_i64_col, df_str_col, read_x_f64};

/// One row of the output BED file.
struct BedRow {
    chrom: String,
    start: i64,
    end: i64,
    name: String,
    score: f64,
}

/// One genomic bin: its column index into `X` and its genomic extent.
#[derive(Clone, Copy)]
struct Bin {
    col: usize,
    start: i64,
    end: i64,
}

/// Banded correlation
///
/// Compute the first `k` diagonals on each side of the main diagonal of the
/// Pearson correlation matrix of `cols`, in banded storage: `band[k + d]`
/// holds the diagonal at offset `d` (`band[k]` is the all-ones main
/// diagonal).
///
/// `cols[j]` is bin `j`'s column of `X`, as `(row_indices, values)` sorted
/// by row index; `n` is the number of rows (cells) `X` was drawn from.
fn sparse_band_corr(cols: &[(&[usize], &[f64])], n: usize, k: usize) -> Vec<Vec<f64>> {
    let p = cols.len();
    let mut band = vec![vec![0.0f64; p]; 2 * k + 1];
    band[k].fill(1.0);

    let n_f = n as f64;
    let col_means: Vec<f64> = cols
        .iter()
        .map(|(_, vals)| vals.iter().sum::<f64>() / n_f)
        .collect();
    let col_norms: Vec<f64> = cols
        .iter()
        .zip(&col_means)
        .map(|((_, vals), &mean)| {
            let sq_sum: f64 = vals.iter().map(|v| v * v).sum();
            let var = sq_sum - n_f * mean * mean;
            let norm = var.max(0.0).sqrt();
            if norm == 0.0 { f64::INFINITY } else { norm }
        })
        .collect();

    for d in 1..=k {
        let n_pairs = p - d;
        for j in 0..n_pairs {
            let (rows_a, vals_a) = cols[j];
            let (rows_b, vals_b) = cols[j + d];
            let raw_dot = sparse_dot(rows_a, vals_a, rows_b, vals_b);
            let centered = raw_dot - n_f * col_means[j] * col_means[j + d];
            let corr = centered / (col_norms[j] * col_norms[j + d]);
            band[k - d][j + d] = corr;
            band[k + d][j] = corr;
        }
    }
    band
}

/// Dot product of two columns given as sorted `(row_index, value)` pairs.
#[inline]
fn sparse_dot(rows_a: &[usize], vals_a: &[f64], rows_b: &[usize], vals_b: &[f64]) -> f64 {
    let (mut i, mut j) = (0, 0);
    let mut sum = 0.0;
    while i < rows_a.len() && j < rows_b.len() {
        match rows_a[i].cmp(&rows_b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                sum += vals_a[i] * vals_b[j];
                i += 1;
                j += 1;
            }
        }
    }
    sum
}

/// Gaussian kernel + diagonal-convolution scoring
///
/// Square Gaussian kernel from `distance_kernel` in `VCRfinder.py`: despite
/// its `np.eye`/`np.rot90` construction, the value at `(r, c)` reduces to a
/// Gaussian of the Chebyshev distance from the kernel's center, normalized
/// to sum to 1.
fn distance_kernel(sigma: f64, radius: usize) -> Vec<Vec<f64>> {
    let width = 1 + 2 * radius;
    let mut kernel = vec![vec![0.0f64; width]; width];
    let mut total = 0.0;
    for (r, row) in kernel.iter_mut().enumerate() {
        for (c, cell) in row.iter_mut().enumerate() {
            let dr = (r as i64 - radius as i64).unsigned_abs();
            let dc = (c as i64 - radius as i64).unsigned_abs();
            let d = dr.max(dc) as f64;
            let v = (-(d * d) / (2.0 * sigma * sigma)).exp();
            *cell = v;
            total += v;
        }
    }
    for row in &mut kernel {
        for v in row.iter_mut() {
            *v /= total;
        }
    }
    kernel
}

/// The `k`-th diagonal of a square matrix, matching `numpy.diag(a, k)`:
/// `a[i, i + k]` for `k >= 0`, `a[i - k, i]` for `k < 0`.
fn matrix_diagonal(m: &[Vec<f64>], width: usize, k: i64) -> Vec<f64> {
    let mut out = Vec::with_capacity(width);
    if k >= 0 {
        let k = k as usize;
        for i in 0..width - k {
            out.push(m[i][i + k]);
        }
    } else {
        let k = (-k) as usize;
        for i in 0..width - k {
            out.push(m[i + k][i]);
        }
    }
    out
}

/// 1-D convolution of `a` and `b`, truncated to `numpy.convolve(...,
/// mode="same")` semantics: output length `max(len(a), len(b))`, centered
/// within the full convolution.
fn convolve_same(a: &[f64], b: &[f64]) -> Vec<f64> {
    let (m, n) = (a.len(), b.len());
    let out_len = m.max(n);
    let full_len = m + n - 1;
    let mut full = vec![0.0f64; full_len];
    for (i, &av) in a.iter().enumerate() {
        if av == 0.0 {
            continue;
        }
        for (j, &bv) in b.iter().enumerate() {
            full[i + j] += av * bv;
        }
    }
    let start = (full_len - out_len) / 2;
    full[start..start + out_len].to_vec()
}

/// Decompose the 2-D Gaussian-kernel convolution of the banded correlation
/// matrix into 1-D convolutions per diagonal (`_score_sigma` in
/// `VCRfinder.py`), producing one score per bin for a single kernel `sigma`.
fn score_sigma(band: &[Vec<f64>], k: usize, p: usize, sigma: f64) -> Vec<f64> {
    let radius = ((4.0 * sigma) as usize).min(k);
    let kernel = distance_kernel(sigma, radius);
    let k_width = kernel.len();
    let k_radius = (k_width / 2).min(k) as i64;

    let mut score = vec![0.0f64; p];
    for d in -k_radius..=k_radius {
        let mut kernel_diag = matrix_diagonal(&kernel, k_width, d);
        kernel_diag.reverse();
        let band_idx = (k as i64 - d) as usize;
        let conv = convolve_same(&band[band_idx], &kernel_diag);
        for (s, c) in score.iter_mut().zip(&conv) {
            *s += c;
        }
    }
    score
}

/// PELT changepoint detection (L2 / mean-shift cost)
///
/// O(1)-query cost of segment `[a, b)` under the L2 ("mean-shift") cost: the
/// sum of squared deviations from the segment mean, summed over all score
/// dimensions. Mathematically equivalent to the cost
/// `ruptures.KernelCPD(kernel="linear")` minimizes.
struct SegmentCost {
    n_dims: usize,
    prefix_sum: Vec<f64>,
    prefix_sq_sum: Vec<f64>,
}

impl SegmentCost {
    fn new(scores: &[Vec<f64>]) -> Self {
        let n = scores.len();
        let n_dims = scores.first().map_or(0, |row| row.len());
        let mut prefix_sum = vec![0.0f64; (n + 1) * n_dims];
        let mut prefix_sq_sum = vec![0.0f64; (n + 1) * n_dims];
        for (t, row) in scores.iter().enumerate() {
            for (d, &v) in row.iter().enumerate() {
                prefix_sum[(t + 1) * n_dims + d] = prefix_sum[t * n_dims + d] + v;
                prefix_sq_sum[(t + 1) * n_dims + d] = prefix_sq_sum[t * n_dims + d] + v * v;
            }
        }
        Self {
            n_dims,
            prefix_sum,
            prefix_sq_sum,
        }
    }

    #[inline]
    fn cost(&self, a: usize, b: usize) -> f64 {
        let len = (b - a) as f64;
        let mut total = 0.0;
        for d in 0..self.n_dims {
            let s1 = self.prefix_sum[b * self.n_dims + d] - self.prefix_sum[a * self.n_dims + d];
            let s2 =
                self.prefix_sq_sum[b * self.n_dims + d] - self.prefix_sq_sum[a * self.n_dims + d];
            total += s2 - s1 * s1 / len;
        }
        total
    }
}

/// PELT changepoint detection with the L2 cost above, `min_size = 1`,
/// `jump = 1` (every split point considered), the same minimization that
/// `ruptures.KernelCPD(kernel="linear", min_size=1, jump=1).predict(pen=...)`
/// solves, for penalty `pen`.
///
/// Returns the detected breakpoints (segment end positions); the last
/// breakpoint always equals `scores.len()`.
fn pelt_predict(scores: &[Vec<f64>], pen: f64) -> Vec<usize> {
    let n = scores.len();
    let cost = SegmentCost::new(scores);

    let mut f = vec![0.0f64; n + 1];
    let mut last_cp = vec![0usize; n + 1];
    f[0] = -pen;
    let mut candidates = vec![0usize];

    for t in 1..=n {
        let (best_s, best_f) = candidates
            .iter()
            .map(|&s| (s, f[s] + cost.cost(s, t) + pen))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .expect("candidates always contains at least 0");
        f[t] = best_f;
        last_cp[t] = best_s;
        candidates.retain(|&s| f[s] + cost.cost(s, t) <= f[t]);
        candidates.push(t);
    }

    let mut bkps = Vec::new();
    let mut t = n;
    while t > 0 {
        bkps.push(t);
        t = last_cp[t];
    }
    bkps.reverse();
    bkps
}

/// Parse a `chrom[:start[:end]]` region string against the full bin table,
/// matching `VCRfinder`'s region-clamping logic.
fn parse_region(
    region: &str,
    chrom: &[String],
    start: &[i64],
    end: &[i64],
) -> Result<(String, i64, i64)> {
    let parts: Vec<&str> = region.split(':').collect();
    let chrom_name = parts[0].to_string();
    let in_chrom = || chrom.iter().enumerate().filter(|(_, c)| **c == chrom_name);

    match parts.len() {
        3 => Ok((chrom_name, parts[1].parse()?, parts[2].parse()?)),
        2 => {
            let region_start: i64 = parts[1].parse()?;
            let region_end = in_chrom()
                .map(|(i, _)| end[i])
                .max()
                .with_context(|| format!("region chromosome {chrom_name:?} not found"))?;
            Ok((chrom_name, region_start, region_end))
        }
        1 => {
            let bounds: Vec<(i64, i64)> = in_chrom().map(|(i, _)| (start[i], end[i])).collect();
            anyhow::ensure!(
                !bounds.is_empty(),
                "region chromosome {chrom_name:?} not found"
            );
            let region_start = bounds.iter().map(|(s, _)| *s).min().unwrap();
            let region_end = bounds.iter().map(|(_, e)| *e).max().unwrap();
            Ok((chrom_name, region_start, region_end))
        }
        _ => bail!("invalid region {region:?}; expected CHROM[:START[:END]]"),
    }
}

fn process_chromosome(
    chrom: &str,
    mut bins: Vec<Bin>,
    x_csc: &CscMatrix<f64>,
    n_obs: usize,
    bin_size: i64,
    max_region_size: i64,
    n_kernels: usize,
    penalties: &[f64],
) -> Result<Vec<BedRow>> {
    if bins.len() == 1 {
        eprintln!("Skipping chromosome {chrom} with only one bin.");
        let (start, end) = (bins[0].start, bins[0].end);
        return Ok(penalties
            .iter()
            .map(|&pen| BedRow {
                chrom: chrom.to_string(),
                start,
                end,
                name: format!("pen-{pen}_brkpoint-1"),
                score: pen,
            })
            .collect());
    }

    let region_start = bins.iter().map(|b| b.start).min().unwrap();
    bins.sort_unstable_by_key(|b| (b.start, b.end));

    // The last bin's width often differs from `bin_size` (chromosome end not
    // evenly divisible); extend it up to `bin_size` if it's at least half a
    // bin, otherwise drop it. Every remaining bin must then be exactly
    // `bin_size` (or `bin_size - 1`, for inclusive/exclusive off-by-ones).
    {
        let last = bins.last_mut().expect("checked len > 1 above");
        let width = last.end - last.start;
        if width != bin_size {
            if width >= bin_size / 2 {
                last.end = last.start + bin_size;
            } else {
                eprintln!("Feature removed due to difference in binsize in chromosome {chrom}");
                bins.pop();
            }
        }
    }
    for b in &bins {
        let width = b.end - b.start;
        anyhow::ensure!(
            width == bin_size || width == bin_size - 1,
            "variable bin sizes detected in chromosome {chrom}"
        );
    }

    let p = bins.len();
    let reg_to_consider =
        max_region_size.min((p as i64 * bin_size / n_kernels as i64).max(bin_size));
    let k_factor = (reg_to_consider as f64 / bin_size as f64).powf(1.0 / n_kernels as f64);
    let sigmas: Vec<f64> = (1..=n_kernels)
        .map(|i| 0.25 * k_factor.powi(i as i32))
        .collect();
    let max_sigma = sigmas.iter().cloned().fold(0.0f64, f64::max);
    let k = ((4.0 * max_sigma).round() as usize).min(p);

    let csc_cols: Vec<_> = bins.iter().map(|b| x_csc.col(b.col)).collect();
    let cols: Vec<(&[usize], &[f64])> = csc_cols
        .iter()
        .map(|c| (c.row_indices(), c.values()))
        .collect();
    let band = sparse_band_corr(&cols, n_obs, k);

    let score_cols: Vec<Vec<f64>> = sigmas
        .par_iter()
        .map(|&sigma| score_sigma(&band, k, p, sigma))
        .collect();
    let mut scores = vec![vec![0.0f64; n_kernels]; p];
    for (j, col) in score_cols.iter().enumerate() {
        for (t, &v) in col.iter().enumerate() {
            scores[t][j] = v;
        }
    }

    let mut out = Vec::new();
    for &pen in penalties {
        let bkps = pelt_predict(&scores, pen);
        eprintln!(
            "Detected {} VCRs in chromosome {chrom} with penalty {pen}",
            bkps.len()
        );
        let mut prev = 0usize;
        for &bkp in &bkps {
            out.push(BedRow {
                chrom: chrom.to_string(),
                start: region_start + prev as i64 * bin_size,
                end: region_start + bkp as i64 * bin_size,
                name: format!("{chrom}_VCR_{bkp}_pen{pen}"),
                score: pen,
            });
            prev = bkp;
        }
    }
    Ok(out)
}

fn run_find_vcrs(
    h5ad_path: &Path,
    bin_size: i64,
    out_file: &Path,
    max_region_size: i64,
    n_kernels: usize,
    penalties: &[f64],
    region: Option<&str>,
    num_threads: usize,
) -> Result<()> {
    anyhow::ensure!(bin_size > 0, "bin_size must be positive");
    anyhow::ensure!(n_kernels > 0, "n_kernels must be positive");
    anyhow::ensure!(
        !penalties.is_empty(),
        "at least one penalty value is required"
    );

    let store = H5::open(h5ad_path)
        .with_context(|| format!("failed to open AnnData file: {}", h5ad_path.display()))?;
    let adata = AnnData::<H5>::open(store)?;

    let var = adata.read_var()?;
    let chrom = df_str_col(&var, "chrom")?;
    let start = df_i64_col(&var, "start")?;
    let end = df_i64_col(&var, "end")?;
    let n_obs = adata.n_obs();
    let x = read_x_f64(&adata)?;
    adata.close()?;
    let x_csc = CscMatrix::from(&x);

    // Original matrix column indices to keep, after region filtering (or
    // all of them, in their original order, when no region is set).
    let mut cols: Vec<usize> = (0..chrom.len()).collect();
    if let Some(region) = region {
        let (rchrom, rstart, rend) = parse_region(region, &chrom, &start, &end)?;
        cols.retain(|&i| chrom[i] == rchrom && start[i] >= rstart && end[i] <= rend);
        anyhow::ensure!(!cols.is_empty(), "region {region:?} matched no bins");
    }

    // Group the retained bins by chromosome, preserving first-seen order
    // (matching pandas' `Series.unique()`).
    let mut chrom_order: Vec<String> = Vec::new();
    let mut chrom_bins: AHashMap<String, Vec<Bin>> = AHashMap::new();
    for &i in &cols {
        chrom_bins
            .entry(chrom[i].clone())
            .or_insert_with(|| {
                chrom_order.push(chrom[i].clone());
                Vec::new()
            })
            .push(Bin {
                col: i,
                start: start[i],
                end: end[i],
            });
    }

    let n_threads = if num_threads == 0 {
        rayon::current_num_threads()
    } else {
        num_threads
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .context("failed to build thread pool")?;

    let rows: Vec<BedRow> = pool
        .install(|| {
            chrom_order
                .par_iter()
                .map(|c| {
                    eprintln!("Processing chromosome {c}...");
                    let bins = chrom_bins.get(c).expect("just inserted above").clone();
                    process_chromosome(
                        c,
                        bins,
                        &x_csc,
                        n_obs,
                        bin_size,
                        max_region_size,
                        n_kernels,
                        penalties,
                    )
                })
                .collect::<Result<Vec<Vec<BedRow>>>>()
        })?
        .into_iter()
        .flatten()
        .collect();

    let mut out = std::fs::File::create(out_file)
        .with_context(|| format!("failed to create output file: {}", out_file.display()))?;
    for r in &rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t*",
            r.chrom, r.start, r.end, r.name, r.score
        )?;
    }
    Ok(())
}

/// Detect variable chromatin regions (VCRs) from binned single-cell
/// chromatin-accessibility data.
///
/// Reads `h5ad_path` (an AnnData file such as one produced by
/// `count_bins`), whose `var` carries `chrom`/`start`/`end` columns for each
/// bin. For each chromosome, computes a banded bin-to-bin correlation
/// matrix, convolves it with `n_kernels` Gaussian kernels of increasing
/// scale (up to `max_region_size`, default `bin_size * 100`) to build a
/// per-bin score map, then runs PELT changepoint detection on the score map
/// once per `penalty` value. Writes the detected segments to `out_file` as
/// a headerless BED file with columns `chrom, start, end, name, score,
/// strand`, where `score` is the penalty that produced that row, filter on
/// it to recover one non-overlapping segmentation. Returns `out_file`.
#[pyfunction(signature = (
    h5ad_path,
    bin_size,
    out_file,
    max_region_size = None,
    n_kernels = 20,
    penalties = vec![0.05, 0.1, 0.5],
    region = None,
    num_threads = 0,
))]
pub fn find_vcrs(
    h5ad_path: PathBuf,
    bin_size: u32,
    out_file: PathBuf,
    max_region_size: Option<u32>,
    n_kernels: usize,
    penalties: Vec<f64>,
    region: Option<String>,
    num_threads: usize,
) -> PyResult<String> {
    let bin_size = bin_size as i64;
    let max_region_size = max_region_size.map(|v| v as i64).unwrap_or(bin_size * 100);
    run_find_vcrs(
        &h5ad_path,
        bin_size,
        &out_file,
        max_region_size,
        n_kernels,
        &penalties,
        region.as_deref(),
        num_threads,
    )
    .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?;
    Ok(out_file.display().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Absolute tolerance for values that go through a few float operations.
    const EPS: f64 = 1e-9;

    fn assert_close(got: f64, want: f64) {
        assert!((got - want).abs() < EPS, "expected {want}, got {got}");
    }

    // Sparse dot product

    #[test]
    fn sparse_dot_multiplies_only_the_rows_the_two_columns_share() {
        let rows_a = [0usize, 2, 5];
        let vals_a = [1.0, 2.0, 3.0];
        let rows_b = [2usize, 3, 5];
        let vals_b = [10.0, 20.0, 30.0];

        // Rows 2 and 5 are shared: 2*10 + 3*30.
        assert_close(sparse_dot(&rows_a, &vals_a, &rows_b, &vals_b), 110.0);
    }

    #[test]
    fn sparse_dot_of_columns_with_no_shared_rows_is_zero() {
        assert_close(sparse_dot(&[0, 2], &[1.0, 1.0], &[1, 3], &[1.0, 1.0]), 0.0);
    }

    #[test]
    fn sparse_dot_with_an_empty_column_is_zero() {
        assert_close(sparse_dot(&[], &[], &[1, 2], &[1.0, 1.0]), 0.0);
        assert_close(sparse_dot(&[1, 2], &[1.0, 1.0], &[], &[]), 0.0);
    }

    // Banded correlation

    #[test]
    fn the_main_diagonal_of_the_correlation_band_is_all_ones() {
        let rows = [0usize, 1, 2, 3];
        let a = [1.0, 2.0, 3.0, 4.0];
        let cols: Vec<(&[usize], &[f64])> = vec![(&rows, &a), (&rows, &a), (&rows, &a)];

        let band = sparse_band_corr(&cols, 4, 2);

        assert_eq!(band.len(), 2 * 2 + 1, "band holds 2k+1 diagonals");
        assert_eq!(band[2], vec![1.0; 3], "band[k] is the main diagonal");
    }

    #[test]
    fn two_identical_columns_correlate_at_one() {
        let rows = [0usize, 1, 2, 3];
        let a = [1.0, 2.0, 3.0, 4.0];
        let cols: Vec<(&[usize], &[f64])> = vec![(&rows, &a), (&rows, &a)];

        let band = sparse_band_corr(&cols, 4, 1);

        // k = 1, so band[0] is the -1 diagonal and band[2] the +1 diagonal.
        assert_close(band[2][0], 1.0);
        assert_close(band[0][1], 1.0);
    }

    #[test]
    fn two_opposed_columns_correlate_at_minus_one() {
        let rows = [0usize, 1, 2, 3];
        let up = [1.0, 2.0, 3.0, 4.0];
        let down = [4.0, 3.0, 2.0, 1.0];
        let cols: Vec<(&[usize], &[f64])> = vec![(&rows, &up), (&rows, &down)];

        let band = sparse_band_corr(&cols, 4, 1);

        assert_close(band[2][0], -1.0);
    }

    #[test]
    fn a_column_with_no_variance_correlates_with_nothing() {
        // A constant column has zero norm; the code substitutes infinity so the
        // correlation collapses to zero rather than dividing by zero.
        let rows = [0usize, 1, 2, 3];
        let varying = [1.0, 2.0, 3.0, 4.0];
        let constant = [2.0, 2.0, 2.0, 2.0];
        let cols: Vec<(&[usize], &[f64])> = vec![(&rows, &varying), (&rows, &constant)];

        let band = sparse_band_corr(&cols, 4, 1);

        assert_close(band[2][0], 0.0);
        assert!(band[2][0].is_finite(), "must not be NaN or infinite");
    }

    #[test]
    fn the_band_is_symmetric_about_the_main_diagonal() {
        let rows = [0usize, 1, 2, 3];
        let a = [1.0, 5.0, 2.0, 8.0];
        let b = [3.0, 1.0, 9.0, 2.0];
        let c = [7.0, 2.0, 4.0, 1.0];
        let cols: Vec<(&[usize], &[f64])> = vec![(&rows, &a), (&rows, &b), (&rows, &c)];

        let k = 2;
        let band = sparse_band_corr(&cols, 4, k);

        for d in 1..=k {
            for j in 0..cols.len() - d {
                assert_close(band[k - d][j + d], band[k + d][j]);
            }
        }
    }

    // Gaussian kernel

    #[test]
    fn the_distance_kernel_sums_to_one() {
        for radius in [0usize, 1, 3, 5] {
            let kernel = distance_kernel(1.5, radius);
            let total: f64 = kernel.iter().flatten().sum();
            assert_close(total, 1.0);
        }
    }

    #[test]
    fn the_distance_kernel_is_square_and_odd_sized() {
        let kernel = distance_kernel(2.0, 3);
        assert_eq!(kernel.len(), 1 + 2 * 3);
        for row in &kernel {
            assert_eq!(row.len(), kernel.len());
        }
    }

    #[test]
    fn a_zero_radius_kernel_is_a_single_cell_holding_all_the_weight() {
        assert_eq!(distance_kernel(1.0, 0), vec![vec![1.0]]);
    }

    #[test]
    fn the_kernel_peaks_at_its_centre_and_falls_off_symmetrically() {
        let radius = 3;
        let kernel = distance_kernel(2.0, radius);
        let centre = kernel[radius][radius];

        for (r, row) in kernel.iter().enumerate() {
            for (c, &v) in row.iter().enumerate() {
                assert!(v <= centre + EPS, "cell ({r},{c}) exceeded the centre");
                // Chebyshev distance decides the value, so opposite cells match.
                assert_close(v, kernel[2 * radius - r][2 * radius - c]);
            }
        }
    }

    // Matrix diagonals

    #[test]
    fn the_zeroth_diagonal_is_the_main_diagonal() {
        let m = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        assert_eq!(matrix_diagonal(&m, 3, 0), vec![1.0, 5.0, 9.0]);
    }

    #[test]
    fn positive_and_negative_offsets_walk_opposite_sides() {
        let m = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        // numpy.diag(a, 1) -> a[i, i+1]; numpy.diag(a, -1) -> a[i+1, i].
        assert_eq!(matrix_diagonal(&m, 3, 1), vec![2.0, 6.0]);
        assert_eq!(matrix_diagonal(&m, 3, -1), vec![4.0, 8.0]);
        assert_eq!(matrix_diagonal(&m, 3, 2), vec![3.0]);
        assert_eq!(matrix_diagonal(&m, 3, -2), vec![7.0]);
    }

    // Convolution

    #[test]
    fn convolving_with_a_unit_impulse_returns_the_input() {
        let a = [1.0, 2.0, 3.0];
        assert_eq!(convolve_same(&a, &[1.0]), vec![1.0, 2.0, 3.0]);
        // A centred delta of odd width is also the identity.
        assert_eq!(convolve_same(&a, &[0.0, 1.0, 0.0]), vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn convolve_same_matches_numpy_same_mode() {
        // numpy.convolve([1,1,1], [1,1], mode="same") == [1, 2, 2]
        assert_eq!(
            convolve_same(&[1.0, 1.0, 1.0], &[1.0, 1.0]),
            vec![1.0, 2.0, 2.0]
        );
    }

    #[test]
    fn the_convolution_is_as_long_as_its_longer_input() {
        assert_eq!(convolve_same(&[1.0, 2.0, 3.0, 4.0], &[1.0, 1.0]).len(), 4);
        assert_eq!(convolve_same(&[1.0, 1.0], &[1.0, 2.0, 3.0, 4.0]).len(), 4);
    }

    // Segment cost

    #[test]
    fn a_segment_of_identical_values_costs_nothing() {
        let scores = vec![vec![3.0], vec![3.0], vec![3.0], vec![3.0]];
        let cost = SegmentCost::new(&scores);
        assert_close(cost.cost(0, 4), 0.0);
    }

    #[test]
    fn the_segment_cost_is_the_sum_of_squared_deviations_from_the_mean() {
        // Values 1..4: mean 2.5, deviations 1.5/0.5/0.5/1.5, squares sum to 5.
        let scores = vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0]];
        let cost = SegmentCost::new(&scores);
        assert_close(cost.cost(0, 4), 5.0);
    }

    #[test]
    fn segment_costs_add_across_score_dimensions() {
        let scores = vec![
            vec![1.0, 1.0],
            vec![2.0, 2.0],
            vec![3.0, 3.0],
            vec![4.0, 4.0],
        ];
        let cost = SegmentCost::new(&scores);
        assert_close(cost.cost(0, 4), 10.0);
    }

    // PELT changepoints

    #[test]
    fn the_last_breakpoint_is_always_the_end_of_the_series() {
        let scores = vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0]];
        for pen in [0.0, 1.0, 100.0] {
            let bkps = pelt_predict(&scores, pen);
            assert_eq!(*bkps.last().unwrap(), scores.len(), "pen {pen}");
        }
    }

    #[test]
    fn a_flat_series_is_one_segment() {
        let scores = vec![vec![5.0]; 8];
        assert_eq!(pelt_predict(&scores, 1.0), vec![8]);
    }

    #[test]
    fn a_step_change_is_found_where_the_level_shifts() {
        let mut scores = vec![vec![0.0]; 5];
        scores.extend(vec![vec![10.0]; 5]);

        // A small penalty makes the split worth its cost.
        assert_eq!(pelt_predict(&scores, 1.0), vec![5, 10]);
    }

    #[test]
    fn a_large_penalty_suppresses_the_split() {
        let mut scores = vec![vec![0.0]; 5];
        scores.extend(vec![vec![10.0]; 5]);

        // The step costs 250 to leave unsplit, so a penalty above that wins.
        assert_eq!(pelt_predict(&scores, 1_000.0), vec![10]);
    }

    // Region parsing

    fn bin_table() -> (Vec<String>, Vec<i64>, Vec<i64>) {
        (
            vec!["chr1".to_string(), "chr1".to_string(), "chr2".to_string()],
            vec![0, 100, 500],
            vec![100, 200, 600],
        )
    }

    #[test]
    fn a_bare_chromosome_spans_all_of_its_bins() {
        let (chrom, start, end) = bin_table();
        assert_eq!(
            parse_region("chr1", &chrom, &start, &end).unwrap(),
            ("chr1".to_string(), 0, 200)
        );
    }

    #[test]
    fn a_start_only_region_runs_to_the_last_bin_of_the_chromosome() {
        let (chrom, start, end) = bin_table();
        assert_eq!(
            parse_region("chr1:50", &chrom, &start, &end).unwrap(),
            ("chr1".to_string(), 50, 200)
        );
    }

    #[test]
    fn an_explicit_range_is_taken_as_given() {
        let (chrom, start, end) = bin_table();
        assert_eq!(
            parse_region("chr1:10:90", &chrom, &start, &end).unwrap(),
            ("chr1".to_string(), 10, 90)
        );
    }

    #[test]
    fn a_chromosome_with_no_bins_is_rejected() {
        let (chrom, start, end) = bin_table();
        let err = parse_region("chrX", &chrom, &start, &end)
            .unwrap_err()
            .to_string();
        assert!(err.contains("chrX"), "{err}");
    }

    #[test]
    fn a_region_with_too_many_fields_is_rejected() {
        let (chrom, start, end) = bin_table();
        assert!(parse_region("chr1:1:2:3", &chrom, &start, &end).is_err());
    }

    #[test]
    fn a_non_numeric_coordinate_is_rejected() {
        let (chrom, start, end) = bin_table();
        assert!(parse_region("chr1:abc", &chrom, &start, &end).is_err());
    }
}
