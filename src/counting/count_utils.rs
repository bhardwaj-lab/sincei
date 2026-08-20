//! Builds the count matrix, and moves it in and out of AnnData.
//!
//! `build_csr` turns the counting loops' sparse (cell, feature) accumulator
//! into a CSR matrix, and `write_counts_anndata` writes it alongside the obs/var
//! tables naming its cells and regions.
//!
//! The readers (`read_x_f64` and the column helpers) go the other way, for the
//! downstream commands that take a sincei matrix back in.

use std::path::Path;

use ahash::AHashMap;

use anndata::backend::{Compression, WriteConfig, set_default_write_config};
use anndata::data::DynCsrMatrix;
use anndata::{AnnData, AnnDataOp, ArrayElemOp};
use anndata_hdf5::H5;
use anyhow::{Context, Result, bail};
use nalgebra_sparse::CsrMatrix;
use polars::prelude::*;

use crate::annotation::region_index::Feature;

/// Read `adata.X` as `f64`, regardless of its on-disk numeric dtype.
///
/// Bool / string matrices are rejected. Shared by the downstream commands
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
/// `barcode` columns. Features (`var`) carry `chrom`, `start`, `end`, and
/// `name` columns; `var_names` are the feature names.
pub(crate) fn write_counts_anndata(
    output_path: &Path,
    matrix: CsrMatrix<u32>,
    // One label per row block, in row order. Normally the input BAMs' labels;
    // under `--groupTag` the merged BAM's `@RG` IDs.
    samples: &[String],
    barcodes: &[String],
    var: &[Feature],
    compression: &str,
    compression_level: u8,
) -> Result<()> {
    // Choose the HDF5 dataset compression. anndata-rs defaults to blosc-zstd,
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
    let n_cells = samples.len() * barcodes.len();
    let mut obs_index: Vec<String> = Vec::with_capacity(n_cells);
    let mut sample_col: Vec<String> = Vec::with_capacity(n_cells);
    let mut barcode_col: Vec<String> = Vec::with_capacity(n_cells);
    for sample in samples {
        for bc in barcodes {
            obs_index.push(format!("{}::{}", sample, bc));
            sample_col.push(sample.clone());
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

    // var: chrom, start, end, name, in feature-index order.
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

#[cfg(test)]
mod tests {
    use super::*;
    // `H5::open` comes from this trait, which the module itself does not need.
    use anndata::Backend;

    /// Dense view of a CSR matrix, for readable assertions.
    fn dense(m: &CsrMatrix<u32>) -> Vec<Vec<u32>> {
        let mut out = vec![vec![0u32; m.ncols()]; m.nrows()];
        for (r, c, &v) in m.triplet_iter() {
            out[r][c] = v;
        }
        out
    }

    #[test]
    fn build_csr_places_every_entry_at_its_coordinate() {
        let acc: AHashMap<(usize, usize), u32> = [((0, 1), 5), ((1, 0), 3), ((1, 2), 7)]
            .into_iter()
            .collect();

        let m = build_csr(&acc, 2, 3).unwrap();

        assert_eq!((m.nrows(), m.ncols()), (2, 3));
        assert_eq!(m.nnz(), 3);
        assert_eq!(dense(&m), vec![vec![0, 5, 0], vec![3, 0, 7]]);
    }

    #[test]
    fn build_csr_sorts_column_indices_within_each_row() {
        // CSR requires ascending column indices per row; the accumulator is an
        // unordered hash map, so the sort has to happen in build_csr.
        let acc: AHashMap<(usize, usize), u32> = (0..8).map(|c| ((0, 7 - c), 1u32)).collect();

        let m = build_csr(&acc, 1, 8).unwrap();

        let cols: Vec<usize> = m.col_indices().to_vec();
        assert_eq!(cols, (0..8).collect::<Vec<_>>());
    }

    #[test]
    fn build_csr_keeps_empty_rows_and_the_declared_shape() {
        // Cells with no counts must still occupy a row: obs is the full
        // sample × barcode product regardless of coverage.
        let acc: AHashMap<(usize, usize), u32> = [((2, 0), 4)].into_iter().collect();

        let m = build_csr(&acc, 4, 2).unwrap();

        assert_eq!((m.nrows(), m.ncols()), (4, 2));
        assert_eq!(m.nnz(), 1);
        assert_eq!(m.row_offsets(), &[0, 0, 0, 1, 1]);
        assert_eq!(
            dense(&m),
            vec![vec![0, 0], vec![0, 0], vec![4, 0], vec![0, 0]]
        );
    }

    #[test]
    fn build_csr_on_an_empty_accumulator_yields_an_all_zero_matrix() {
        let m = build_csr(&AHashMap::new(), 3, 5).unwrap();

        assert_eq!((m.nrows(), m.ncols()), (3, 5));
        assert_eq!(m.nnz(), 0);
    }

    #[test]
    fn df_i64_col_reads_and_casts_integer_columns() {
        let df = DataFrame::new(
            3,
            vec![
                Column::new("start".into(), [10i32, 20, 30]),
                Column::new("name".into(), ["a", "b", "c"]),
            ],
        )
        .unwrap();

        // Cast from i32 to i64 happens transparently.
        assert_eq!(df_i64_col(&df, "start").unwrap(), vec![10, 20, 30]);
        assert!(df_i64_col(&df, "missing").is_err());
    }

    #[test]
    fn df_i64_col_rejects_null_values() {
        let df = DataFrame::new(2, vec![Column::new("start".into(), [Some(1i64), None])]).unwrap();

        assert!(df_i64_col(&df, "start").is_err());
    }

    #[test]
    fn df_str_col_reads_string_columns_and_rejects_nulls() {
        let df = DataFrame::new(2, vec![Column::new("name".into(), ["gene1", "gene2"])]).unwrap();
        assert_eq!(
            df_str_col(&df, "name").unwrap(),
            vec!["gene1".to_string(), "gene2".to_string()]
        );

        let with_null =
            DataFrame::new(2, vec![Column::new("name".into(), [Some("gene1"), None])]).unwrap();
        assert!(df_str_col(&with_null, "name").is_err());
    }

    // AnnData round trip

    fn feature(chrom: &str, start: usize, end: usize) -> Feature {
        Feature {
            chrom: chrom.to_string(),
            start,
            end,
            name: format!("{chrom}:{start}-{end}"),
            strand: '*',
        }
    }

    fn read_back(path: &Path) -> CsrMatrix<f64> {
        let store = H5::open(path).unwrap();
        let adata = AnnData::<H5>::open(store).unwrap();
        read_x_f64(&adata).unwrap()
    }

    #[test]
    fn a_written_matrix_reads_back_with_the_same_values() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("counts.h5ad");

        let acc: AHashMap<(usize, usize), u32> = [((0, 0), 1), ((0, 2), 5), ((1, 1), 3)]
            .into_iter()
            .collect();
        let matrix = build_csr(&acc, 2, 3).unwrap();

        let var = vec![
            feature("chr1", 0, 100),
            feature("chr1", 100, 200),
            feature("chr2", 0, 100),
        ];
        let barcodes = vec!["AAA".to_string(), "CCC".to_string()];

        write_counts_anndata(
            &path,
            matrix,
            &["s1".to_string()],
            &barcodes,
            &var,
            "none",
            0,
        )
        .unwrap();
        assert!(path.exists(), "the h5ad file was not created");

        // X comes back as f64 whatever it was stored as.
        let x = read_back(&path);
        assert_eq!((x.nrows(), x.ncols()), (2, 3));
        let mut got = vec![vec![0.0f64; 3]; 2];
        for (r, c, &v) in x.triplet_iter() {
            got[r][c] = v;
        }
        assert_eq!(got, vec![vec![1.0, 0.0, 5.0], vec![0.0, 3.0, 0.0]]);
    }

    #[test]
    fn obs_rows_are_one_per_sample_and_barcode_in_matrix_order() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("counts.h5ad");

        // Two samples x two barcodes = four rows, samples varying slowest.
        let matrix = build_csr(&AHashMap::new(), 4, 1).unwrap();
        let barcodes = vec!["AAA".to_string(), "CCC".to_string()];

        write_counts_anndata(
            &path,
            matrix,
            &["s1".to_string(), "s2".to_string()],
            &barcodes,
            &[feature("chr1", 0, 100)],
            "none",
            0,
        )
        .unwrap();

        let store = H5::open(&path).unwrap();
        let adata = AnnData::<H5>::open(store).unwrap();
        assert_eq!(
            adata.obs_names().into_vec(),
            vec!["s1::AAA", "s1::CCC", "s2::AAA", "s2::CCC"]
        );
        assert_eq!(adata.var_names().into_vec(), vec!["chr1:0-100"]);
    }

    #[test]
    fn the_var_table_keeps_the_feature_coordinates() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("counts.h5ad");

        let var = vec![feature("chr1", 0, 100), feature("chr2", 500, 750)];
        let matrix = build_csr(&AHashMap::new(), 1, 2).unwrap();

        write_counts_anndata(
            &path,
            matrix,
            &["s1".to_string()],
            &["AAA".to_string()],
            &var,
            "none",
            0,
        )
        .unwrap();

        let store = H5::open(&path).unwrap();
        let adata = AnnData::<H5>::open(store).unwrap();
        let var_df = adata.read_var().unwrap();

        assert_eq!(df_str_col(&var_df, "chrom").unwrap(), vec!["chr1", "chr2"]);
        assert_eq!(df_i64_col(&var_df, "start").unwrap(), vec![0, 500]);
        assert_eq!(df_i64_col(&var_df, "end").unwrap(), vec![100, 750]);
    }

    #[test]
    fn gzip_is_accepted_and_produces_a_readable_file() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("gzipped.h5ad");

        let acc: AHashMap<(usize, usize), u32> = [((0, 0), 42)].into_iter().collect();
        let matrix = build_csr(&acc, 1, 1).unwrap();

        write_counts_anndata(
            &path,
            matrix,
            &["s1".to_string()],
            &["AAA".to_string()],
            &[feature("chr1", 0, 100)],
            "gzip",
            6,
        )
        .unwrap();

        let x = read_back(&path);
        assert_eq!(x.get_entry(0, 0).map(|e| e.into_value()), Some(42.0));
    }

    #[test]
    fn an_unsupported_compression_is_rejected_by_name() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("nope.h5ad");

        // blosc-zstd is anndata-rs's own default but h5py cannot read it.
        let err = write_counts_anndata(
            &path,
            build_csr(&AHashMap::new(), 1, 1).unwrap(),
            &["s1".to_string()],
            &["AAA".to_string()],
            &[feature("chr1", 0, 100)],
            "blosc",
            0,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("blosc"), "{err}");
        assert!(err.contains("gzip"), "{err}");
    }
}
