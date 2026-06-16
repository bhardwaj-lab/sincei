use std::fs::File;
use std::path::Path;

use ahash::AHashMap;

use anndata::backend::{Compression, WriteConfig, set_default_write_config};
use anndata::{AnnData, AnnDataOp};
use anndata_hdf5::H5;
use anyhow::{Context, Result};
use nalgebra_sparse::CsrMatrix;
use noodles::bam;
use noodles::bgzf;
use polars::prelude::*;

use super::filters::QcFilter;
use super::params::CountingParams;
use super::region_index::VarMeta;
use super::sc_record::{ScRecord, ScRecordOptions};

/// Concrete type of a BAI-indexed BAM reader opened from a file path.
pub(crate) type BamReader = bam::io::IndexedReader<bgzf::io::Reader<File>>;

pub(crate) fn effective_interval(rec: &ScRecord<'_>, params: &CountingParams) -> (usize, usize) {
    let (start, end) = match params.extend_reads {
        None => (rec.alignment_start, rec.alignment_end),
        Some(frag_len) => {
            let max_paired = 4 * frag_len;
            let tlen_abs = rec.template_length.unsigned_abs() as usize;
            let use_tlen = rec.is_proper_pair && tlen_abs > 0 && tlen_abs <= max_paired;

            if use_tlen {
                if rec.is_reverse {
                    let mate = rec.next_alignment_start.unwrap_or(rec.alignment_start);
                    (mate, rec.alignment_end)
                } else {
                    (rec.alignment_start, rec.alignment_start + tlen_abs)
                }
            } else if rec.is_reverse {
                (
                    rec.alignment_end.saturating_sub(frag_len),
                    rec.alignment_end,
                )
            } else {
                (rec.alignment_start, rec.alignment_start + frag_len)
            }
        }
    };

    if params.center_reads && rec.read_length > 0 {
        let center = (start + end) / 2;
        let half = rec.read_length / 2;
        let cs = center.saturating_sub(half);
        (cs, cs + rec.read_length)
    } else {
        (start, end)
    }
}

pub(crate) fn derive_record_opts(qc: Option<&QcFilter>, has_motif: bool) -> ScRecordOptions {
    ScRecordOptions {
        compute_gc: qc.is_some_and(|f| f.needs_gc()),
        compute_aligned_fraction: qc.is_some_and(|f| f.needs_aligned_fraction()),
        store_sequence: has_motif,
    }
}

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
    // plugin, so we never use it: only `none` (anndata's modern default) or
    // gzip (built-in deflate, universally readable).
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
    let obs = DataFrame::new(
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
    adata.set_obs(obs)?;
    adata.set_var_names(var_index.into_iter().collect())?;
    adata.set_var(var_df)?;
    adata.close()?;
    Ok(())
}
