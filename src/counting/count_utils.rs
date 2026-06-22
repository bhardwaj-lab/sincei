use std::path::Path;

use ahash::AHashMap;

use anndata::backend::{Compression, WriteConfig, set_default_write_config};
use anndata::data::DynCsrMatrix;
use anndata::{AnnData, AnnDataOp, ArrayElemOp};
use anndata_hdf5::H5;
use anyhow::{Context, Result, bail};
use nalgebra_sparse::CsrMatrix;
use polars::prelude::*;

use super::region_index::VarMeta;

/// Read `adata.X` as `f64`, regardless of its on-disk numeric dtype.
///
/// Bool / string matrices are rejected.  Shared by the downstream commands
/// (`find_vcrs`, `score_features`) that need a dense-friendly float matrix.
pub(crate) fn read_x_f64(adata: &AnnData<H5>) -> Result<CsrMatrix<f64>> {
    let dyn_csr: DynCsrMatrix = adata
        .x()
        .get::<DynCsrMatrix>()?
        .context("AnnData has no X matrix")?;

    fn convert<T: Copy>(m: CsrMatrix<T>, f: impl Fn(T) -> f64) -> CsrMatrix<f64> {
        let (nrows, ncols) = (m.nrows(), m.ncols());
        let (offsets, indices, values) = m.disassemble();
        CsrMatrix::try_from_csr_data(
            nrows,
            ncols,
            offsets,
            indices,
            values.into_iter().map(f).collect(),
        )
        .expect("disassembled CSR data is valid by construction")
    }

    Ok(match dyn_csr {
        DynCsrMatrix::F64(m) => m,
        DynCsrMatrix::F32(m) => convert(m, |v| v as f64),
        DynCsrMatrix::I8(m) => convert(m, |v| v as f64),
        DynCsrMatrix::I16(m) => convert(m, |v| v as f64),
        DynCsrMatrix::I32(m) => convert(m, |v| v as f64),
        DynCsrMatrix::I64(m) => convert(m, |v| v as f64),
        DynCsrMatrix::U8(m) => convert(m, |v| v as f64),
        DynCsrMatrix::U16(m) => convert(m, |v| v as f64),
        DynCsrMatrix::U32(m) => convert(m, |v| v as f64),
        DynCsrMatrix::U64(m) => convert(m, |v| v as f64),
        DynCsrMatrix::Bool(_) | DynCsrMatrix::String(_) => {
            bail!("AnnData X matrix must be numeric (found bool/string data)")
        }
    })
}

/// Read an integer column from a polars `DataFrame` (e.g. an AnnData `var`
/// table), casting to `i64` and erroring on missing column or null values.
pub(crate) fn df_i64_col(df: &DataFrame, name: &str) -> Result<Vec<i64>> {
    let col = df
        .column(name)
        .with_context(|| format!("adata.var is missing required column {name:?}"))?
        .cast(&DataType::Int64)
        .with_context(|| format!("adata.var[{name:?}] must be numeric"))?;
    col.i64()?
        .iter()
        .enumerate()
        .map(|(i, v)| v.ok_or_else(|| anyhow::anyhow!("null value in var[{name:?}] at row {i}")))
        .collect()
}

/// Read a string column from a polars `DataFrame`, erroring on missing column
/// or null values.
pub(crate) fn df_str_col(df: &DataFrame, name: &str) -> Result<Vec<String>> {
    let col = df
        .column(name)
        .with_context(|| format!("adata.var is missing required column {name:?}"))?
        .cast(&DataType::String)
        .with_context(|| format!("adata.var[{name:?}] must be string-typed"))?;
    col.str()?
        .iter()
        .enumerate()
        .map(|(i, v)| {
            v.map(str::to_string)
                .ok_or_else(|| anyhow::anyhow!("null value in var[{name:?}] at row {i}"))
        })
        .collect()
}

/// Build a sparse count matrix in CSR format from a HashMap COO accumulator.
pub(super) fn build_csr(
    accumulator: &AHashMap<(usize, usize), u32>,
    n_rows: usize,
    n_cols: usize,
) -> Result<CsrMatrix<u32>> {
    let mut row_data: Vec<Vec<(usize, u32)>> = vec![vec![]; n_rows];
    for (&(row, col), &val) in accumulator {
        row_data[row].push((col, val));
    }
    for row in &mut row_data {
        row.sort_unstable_by_key(|&(col, _)| col);
    }

    let mut row_offsets = vec![0usize; n_rows + 1];
    let mut col_indices: Vec<usize> = Vec::new();
    let mut values: Vec<u32> = Vec::new();

    for (i, row) in row_data.iter().enumerate() {
        row_offsets[i + 1] = row_offsets[i] + row.len();
        for &(col, val) in row {
            col_indices.push(col);
            values.push(val);
        }
    }

    CsrMatrix::try_from_csr_data(n_rows, n_cols, row_offsets, col_indices, values)
        .map_err(|e| anyhow::anyhow!("failed to build CSR matrix: {:?}", e))
}

/// Write a cell × feature count matrix to an AnnData HDF5 file.
///
/// Cells (`obs`) are the cartesian product of the input BAM samples and the
/// barcode whitelist (`obs_names = "{sample}::{barcode}"`), with `sample` and
/// `barcode` columns.  Features (`var`) carry `chrom`, `start`, `end`, and
/// `name` columns; `var_names` are the feature names.
pub(crate) fn write_counts_anndata(
    output_path: &Path,
    matrix: CsrMatrix<u32>,
    bam_paths: &[(&Path, &str)],
    barcodes: &[String],
    var: &[VarMeta],
    compression: &str,
    compression_level: u8,
) -> Result<()> {
    // Choose the HDF5 dataset compression.  anndata-rs defaults to blosc-zstd,
    // which standard h5py / scanpy cannot read without an external filter
    // plugin, so only support `none` (anndata's modern default) or gzip
    // (built-in deflate, universally readable).
    let compression = match compression {
        "none" => None,
        "gzip" => Some(Compression::Gzip(compression_level)),
        other => anyhow::bail!("unknown compression {:?}; expected 'none' or 'gzip'", other),
    };
    set_default_write_config(WriteConfig {
        compression,
        block_size: None,
    });

    // obs: one row per (sample, barcode), in the same order as the matrix rows.
    let n_cells = bam_paths.len() * barcodes.len();
    let mut obs_index: Vec<String> = Vec::with_capacity(n_cells);
    let mut sample_col: Vec<String> = Vec::with_capacity(n_cells);
    let mut barcode_col: Vec<String> = Vec::with_capacity(n_cells);
    for (_, sample) in bam_paths {
        for bc in barcodes {
            obs_index.push(format!("{}::{}", sample, bc));
            sample_col.push((*sample).to_string());
            barcode_col.push(bc.clone());
        }
    }
    let obs_df = DataFrame::new(
        n_cells,
        vec![
            Column::new("sample".into(), sample_col),
            Column::new("barcode".into(), barcode_col),
        ],
    )?;

    // var: chrom / start / end / name, in feature-index order.
    let var_index: Vec<String> = var.iter().map(|v| v.name.clone()).collect();
    let chrom_col: Vec<String> = var.iter().map(|v| v.chrom.clone()).collect();
    let start_col: Vec<i64> = var.iter().map(|v| v.start as i64).collect();
    let end_col: Vec<i64> = var.iter().map(|v| v.end as i64).collect();
    let name_col: Vec<String> = var.iter().map(|v| v.name.clone()).collect();
    let var_df = DataFrame::new(
        var.len(),
        vec![
            Column::new("chrom".into(), chrom_col),
            Column::new("start".into(), start_col),
            Column::new("end".into(), end_col),
            Column::new("name".into(), name_col),
        ],
    )?;

    let adata = AnnData::<H5>::new(output_path)
        .with_context(|| format!("failed to create AnnData file: {}", output_path.display()))?;
    adata.set_x(matrix)?;
    // Index first (creates the obs/var elements), then the columns.
    adata.set_obs_names(obs_index.into_iter().collect())?;
    adata.set_obs(obs_df)?;
    adata.set_var_names(var_index.into_iter().collect())?;
    adata.set_var(var_df)?;
    adata.close()?;
    Ok(())
}
