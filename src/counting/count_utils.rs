//! Builds the count matrix, and moves it in and out of AnnData.
//!
//! `count_into_anndata` counts one BAM at a time: `build_csr` turns the
//! counting loops' (cell, feature, count) entries into a CSR block, and
//! `CountsWriter` adds each block to the output file, then writes the obs/var
//! tables naming its cells and regions.
//!
//! The readers (`read_x_f64` and the column helpers) go the other way, reading
//! a written matrix back. They are test-only: nothing in the crate reads a
//! count matrix yet.

use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, MutexGuard, PoisonError};

use ahash::AHashMap;

#[cfg(test)]
use anndata::ArrayElemOp;
use anndata::backend::{
    AttributeOp, Backend, Compression, DatasetOp, GroupOp, StoreOp, WriteConfig,
    set_default_write_config,
};
#[cfg(test)]
use anndata::data::DynCsrMatrix;
use anndata::data::SelectInfoElem;
use anndata::{AnnData, AnnDataOp};
use anndata_hdf5::H5;
#[cfg(test)]
use anyhow::bail;
use anyhow::{Context, Result};
use nalgebra_sparse::CsrMatrix;
use ndarray::{Array1, CowArray};
use polars::prelude::*;
use rayon::ThreadPool;

use crate::annotation::region_index::Feature;
use crate::bam::bam_io::{BamWorker, Chunk, Samples};

/// Read `adata.X` as `f64`, regardless of its on-disk numeric dtype.
///
/// Bool / string matrices are rejected.
///
/// Test-only. The downstream commands that would want a dense-friendly float
/// matrix (`scFindVCRs`, `scScoreFeatures`) have no Rust backend yet, so this
/// exists to read back what the counting tests write.
#[cfg(test)]
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

/// Obs names and dense rows of a written count matrix, in row order.
///
/// Test-only, like [`read_x_f64`].
#[cfg(test)]
pub(crate) fn read_rows(path: &Path) -> Vec<(String, Vec<f64>)> {
    let adata = AnnData::<H5>::open(H5::open(path).unwrap()).unwrap();
    let x = read_x_f64(&adata).unwrap();
    let mut rows = vec![vec![0.0; x.ncols()]; x.nrows()];
    for (row, col, &value) in x.triplet_iter() {
        rows[row][col] = value;
    }
    adata.obs_names().into_vec().into_iter().zip(rows).collect()
}

/// The rows a count without a whitelist must give, from a count whose
/// whitelist lists every barcode: the rows with counts, sorted by name.
#[cfg(test)]
pub(crate) fn observed_part(listed: Vec<(String, Vec<f64>)>) -> Vec<(String, Vec<f64>)> {
    let mut rows: Vec<_> = listed
        .into_iter()
        .filter(|(_, row)| row.iter().any(|&v| v != 0.0))
        .collect();
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    rows
}

/// Read an integer column from a polars `DataFrame` (e.g. an AnnData `var`
/// table), casting to `i64` and erroring on missing column or null values.
#[cfg(test)]
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
#[cfg(test)]
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

/// Counts keyed by `(cell, feature)`.
pub(super) type CellCounts = AHashMap<(usize, usize), u32>;

/// One count matrix entry: `(cell, feature, count)`.
pub(super) type Entry = (u32, u32, u32);

fn to_u32(index: usize) -> Result<u32> {
    u32::try_from(index).context("the count matrix has more than 2^32 rows or columns")
}

/// A chunk's counts as matrix entries, with each cell moved to `row(cell)`.
pub(super) fn to_entries(acc: CellCounts, row: impl Fn(usize) -> usize) -> Result<Vec<Entry>> {
    let mut entries = Vec::with_capacity(acc.len());
    for ((cell, feature), count) in acc {
        entries.push((to_u32(row(cell))?, to_u32(feature)?, count));
    }
    Ok(entries)
}

/// Every sample x barcode pair, samples varying slowest: the rows of a count
/// with a barcode whitelist.
pub(super) fn product_cells(samples: &[String], barcodes: &[String]) -> Vec<(String, String)> {
    samples
        .iter()
        .flat_map(|sample| barcodes.iter().map(move |bc| (sample.clone(), bc.clone())))
        .collect()
}

/// Barcodes numbered in the order they are met, for a count without a
/// whitelist.
///
/// There is one per work chunk and one for the whole run: a chunk numbers its
/// barcodes without a lock, and moves its counts to the run's numbers once, at
/// its end ([`to_run_numbers`]).
#[derive(Default)]
pub(super) struct BarcodeNumbers {
    numbers: AHashMap<Vec<u8>, usize>,
}

impl BarcodeNumbers {
    pub(super) fn number(&mut self, barcode: &[u8]) -> usize {
        if let Some(&number) = self.numbers.get(barcode) {
            return number;
        }
        let number = self.numbers.len();
        self.numbers.insert(barcode.to_vec(), number);
        number
    }

    /// The barcodes, indexed by their number.
    pub(super) fn into_barcodes(self) -> Vec<String> {
        let mut barcodes = vec![String::new(); self.numbers.len()];
        for (barcode, number) in self.numbers {
            barcodes[number] = String::from_utf8_lossy(&barcode).into_owned();
        }
        barcodes
    }
}

/// Move a chunk's counts, keyed `(barcode * n_samples + sample, feature)` by
/// the chunk's barcode numbers, to the run's numbers, as matrix entries.
pub(super) fn to_run_numbers(
    acc: CellCounts,
    n_samples: usize,
    chunk: &BarcodeNumbers,
    run: &Mutex<BarcodeNumbers>,
) -> Result<Vec<Entry>> {
    let mut to_run = vec![0usize; chunk.numbers.len()];
    {
        let mut run = run.lock().unwrap_or_else(PoisonError::into_inner);
        for (barcode, &number) in &chunk.numbers {
            to_run[number] = run.number(barcode);
        }
    }
    to_entries(acc, |cell| {
        to_run[cell / n_samples] * n_samples + cell % n_samples
    })
}

/// Keep only the rows with counts, ordered by sample, then barcode.
///
/// The entries' cells are `barcode * n_samples + sample` by the run's barcode
/// numbers, and `barcodes` holds each number's barcode. Moves the entries to
/// the new row numbers, and returns the (sample, barcode) of each row.
pub(super) fn observed_rows(
    entries: &mut [Vec<Entry>],
    samples: &[String],
    barcodes: &[String],
) -> Vec<(String, String)> {
    let n_samples = samples.len();
    let mut seen = vec![false; barcodes.len() * n_samples];
    for &(key, _, _) in entries.iter().flatten() {
        seen[key as usize] = true;
    }
    let mut rows: Vec<(usize, &str, usize)> = (0..seen.len())
        .filter(|&key| seen[key])
        .map(|key| (key % n_samples, barcodes[key / n_samples].as_str(), key))
        .collect();
    rows.sort_unstable();

    let mut row_of = vec![0u32; seen.len()];
    for (row, &(_, _, key)) in rows.iter().enumerate() {
        row_of[key] = row as u32;
    }
    for entry in entries.iter_mut().flatten() {
        entry.0 = row_of[entry.0 as usize];
    }
    rows.iter()
        .map(|&(sample, barcode, _)| (samples[sample].clone(), barcode.to_string()))
        .collect()
}

/// Add two chunks' `(cell, feature) -> count` maps together.
///
/// The smaller map is drained into the larger: merging costs one hash lookup
/// per entry moved, so moving the shorter side does strictly less work.
pub(super) fn merge_counts(
    a: AHashMap<(usize, usize), u32>,
    b: AHashMap<(usize, usize), u32>,
) -> AHashMap<(usize, usize), u32> {
    let (mut keep, drain) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    for (key, val) in drain {
        *keep.entry(key).or_insert(0) += val;
    }
    keep
}

/// Build a sparse count matrix in CSR format from the chunks' entries.
///
/// Two chunks can each hold an entry for the same cell and feature, as a read
/// can reach past the end of its chunk; such entries are added together. Each
/// chunk's entries are freed as soon as they are placed.
pub(super) fn build_csr(
    entries: Vec<Vec<Entry>>,
    n_rows: usize,
    n_cols: usize,
) -> Result<CsrMatrix<u32>> {
    let mut row_offsets = vec![0usize; n_rows + 1];
    for &(row, _, _) in entries.iter().flatten() {
        row_offsets[row as usize + 1] += 1;
    }
    for i in 0..n_rows {
        row_offsets[i + 1] += row_offsets[i];
    }

    let n_entries = row_offsets[n_rows];
    let mut col_indices = vec![0usize; n_entries];
    let mut values = vec![0u32; n_entries];
    let mut next = row_offsets.clone();
    for chunk in entries {
        for (row, col, count) in chunk {
            let slot = &mut next[row as usize];
            col_indices[*slot] = col as usize;
            values[*slot] = count;
            *slot += 1;
        }
    }
    drop(next);

    let mut row: Vec<(usize, u32)> = Vec::new();
    let mut nnz = 0;
    for r in 0..n_rows {
        let (start, end) = (row_offsets[r], row_offsets[r + 1]);
        row.clear();
        row.extend(
            col_indices[start..end]
                .iter()
                .copied()
                .zip(values[start..end].iter().copied()),
        );
        row.sort_unstable_by_key(|&(col, _)| col);
        row_offsets[r] = nnz;
        for &(col, count) in &row {
            if nnz > row_offsets[r] && col_indices[nnz - 1] == col {
                values[nnz - 1] += count;
            } else {
                col_indices[nnz] = col;
                values[nnz] = count;
                nnz += 1;
            }
        }
    }
    row_offsets[n_rows] = nnz;
    col_indices.truncate(nnz);
    col_indices.shrink_to_fit();
    values.truncate(nnz);
    values.shrink_to_fit();

    CsrMatrix::try_from_csr_data(n_rows, n_cols, row_offsets, col_indices, values)
        .map_err(|e| anyhow::anyhow!("failed to build CSR matrix: {:?}", e))
}

/// Tag the root group as an AnnData object.
///
/// The AnnData spec requires every element to carry `encoding-type` and
/// `encoding-version`, including the root (`"anndata"` / `"0.1.0"`). Note that
/// the encoding version is per-element and fixed; it does not track the anndata
/// package release. anndata-rs writes them for every sub-element but never for
/// the root, so an external anndata reader may fail without it.
///
/// The root cannot be reached through `AnnData`: it keeps its store private, and
/// `Backend::Store` is bound by `StoreOp + GroupOp` but not `AttributeOp`.
/// Reopening and asking for `"/"` yields a group, which does implement it.
///
/// Must run after the `AnnData` handle is closed.
fn tag_anndata_root(path: &Path) -> Result<()> {
    let store = H5::open_rw(path)
        .with_context(|| format!("failed to reopen {} to tag its root", path.display()))?;
    let mut root = store.open_group("/")?;
    root.new_attr("encoding-type", "anndata")?;
    root.new_attr("encoding-version", "0.1.0")?;
    Ok(())
}

/// Rows per HDF5 chunk of the growing `data` and `indices` datasets.
const X_BLOCK: usize = 16384;

/// Writes a cell × feature count matrix to a new AnnData HDF5 file, one block
/// of rows at a time, so only the block being added has to be in memory.
pub(crate) struct CountsWriter {
    path: PathBuf,
    store: <H5 as Backend>::Store,
    x: <H5 as Backend>::Group,
    data: <H5 as Backend>::Dataset,
    indices: <H5 as Backend>::Dataset,
    indptr: Vec<usize>,
    n_cols: usize,
    compression: Option<Compression>,
}

impl CountsWriter {
    /// Create the file, with an empty `X` of `n_cols` columns.
    pub(crate) fn create(
        output_path: &Path,
        n_cols: usize,
        compression: &str,
        compression_level: u8,
    ) -> Result<Self> {
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
            compression: compression.clone(),
            block_size: None,
        });
        anyhow::ensure!(
            i32::try_from(n_cols.saturating_sub(1)).is_ok(),
            "the count matrix has more than 2^31 columns"
        );

        let store = H5::new(output_path)
            .with_context(|| format!("failed to create AnnData file: {}", output_path.display()))?;
        let mut x = store.new_group("X")?;
        x.new_attr("encoding-type", "csr_matrix")?;
        x.new_attr("encoding-version", "0.1.0")?;
        let growing = || WriteConfig {
            compression: compression.clone(),
            block_size: Some(X_BLOCK.into()),
        };
        let data = x.new_empty_dataset::<u32>("data", &0.into(), growing())?;
        let indices = x.new_empty_dataset::<i32>("indices", &0.into(), growing())?;
        Ok(Self {
            path: output_path.to_path_buf(),
            store,
            x,
            data,
            indices,
            indptr: vec![0],
            n_cols,
            compression,
        })
    }

    /// Add `block`'s rows below the rows already written.
    pub(crate) fn append(&mut self, block: CsrMatrix<u32>) -> Result<()> {
        anyhow::ensure!(
            block.ncols() == self.n_cols,
            "a block of {} columns cannot be added to a matrix of {}",
            block.ncols(),
            self.n_cols
        );
        let start = self.indptr[self.indptr.len() - 1];
        let (offsets, indices, values) = block.disassemble();
        self.indptr
            .extend(offsets[1..].iter().map(|&offset| start + offset));
        let end = start + values.len();
        if end == start {
            return Ok(());
        }

        let slice = [SelectInfoElem::from(start..end)];
        self.data.reshape(&end.into())?;
        self.data
            .write_array_slice(CowArray::from(Array1::from_vec(values)), &slice)?;
        let indices: Array1<i32> = indices.into_iter().map(|col| col as i32).collect();
        self.indices.reshape(&end.into())?;
        self.indices
            .write_array_slice(CowArray::from(indices), &slice)?;
        Ok(())
    }

    /// Write the rest of `X` and the obs/var tables, and close the file.
    ///
    /// Cells (`obs`) are one `(sample, barcode)` per matrix row, in row order
    /// (`obs_names = "{sample}::{barcode}"`), with `sample` and `barcode` columns.
    /// The sample is an input BAM's label, or under `--groupTag` a merged BAM's
    /// `@RG` ID. Features (`var`) carry `chrom`, `start`, `end`, and `name`
    /// columns; `var_names` are the feature names.
    pub(crate) fn finish(mut self, cells: &[(String, String)], var: &[Feature]) -> Result<()> {
        let n_rows = self.indptr.len() - 1;
        anyhow::ensure!(
            cells.len() == n_rows,
            "{} cells were given for a matrix of {} rows",
            cells.len(),
            n_rows
        );
        anyhow::ensure!(
            var.len() == self.n_cols,
            "{} features were given for a matrix of {} columns",
            var.len(),
            self.n_cols
        );

        // Use i32 or i64 as indptr type in order to be compatible with scipy
        let config = WriteConfig {
            compression: self.compression.clone(),
            block_size: None,
        };
        if i32::try_from(self.indptr[n_rows]).is_ok() {
            let indptr: Array1<i32> = self.indptr.iter().map(|&o| o as i32).collect();
            self.x
                .new_array_dataset("indptr", CowArray::from(indptr), config)?;
        } else {
            let indptr: Array1<i64> = self.indptr.iter().map(|&o| o as i64).collect();
            self.x
                .new_array_dataset("indptr", CowArray::from(indptr), config)?;
        }
        self.x
            .new_attr("shape", [n_rows as u64, self.n_cols as u64].as_slice())?;
        let Self {
            path,
            store,
            x,
            data,
            indices,
            ..
        } = self;
        drop((x, data, indices));
        store.close()?;

        // obs: one row per (sample, barcode), in the same order as the matrix rows.
        let n_cells = cells.len();
        let mut obs_index: Vec<String> = Vec::with_capacity(n_cells);
        let mut sample_col: Vec<String> = Vec::with_capacity(n_cells);
        let mut barcode_col: Vec<String> = Vec::with_capacity(n_cells);
        for (sample, bc) in cells {
            obs_index.push(format!("{}::{}", sample, bc));
            sample_col.push(sample.clone());
            barcode_col.push(bc.clone());
        }
        let obs_df = DataFrame::new(
            n_cells,
            vec![
                Column::new("sample".into(), sample_col)
                    .cast(&DataType::from_categories(Categories::global()))?,
                Column::new("barcode".into(), barcode_col)
                    .cast(&DataType::from_categories(Categories::global()))?,
            ],
        )?;

        // var: chrom, start, end, name, in feature-index order.
        //
        // `var_names` are `{chrom}_{start}_{end}::{name}`. Bins and unnamed
        // features render the name as the literal `None`.
        let locus = |v: &Feature| format!("{}_{}_{}", v.chrom, v.start, v.end);
        let var_index: Vec<String> = var
            .iter()
            .map(|v| format!("{}::{}", locus(v), v.name.as_deref().unwrap_or("None")))
            .collect();
        let chrom_col: Vec<String> = var.iter().map(|v| v.chrom.clone()).collect();
        let start_col: Vec<i64> = var.iter().map(|v| v.start as i64).collect();
        let end_col: Vec<i64> = var.iter().map(|v| v.end as i64).collect();
        let name_col: Vec<String> = var.iter().map(locus).collect();
        let var_df = DataFrame::new(
            var.len(),
            vec![
                Column::new("chrom".into(), chrom_col),
                Column::new("start".into(), start_col),
                Column::new("end".into(), end_col),
                Column::new("name".into(), name_col),
            ],
        )?;

        let adata = AnnData::<H5>::open(H5::open_rw(&path)?)
            .with_context(|| format!("failed to reopen AnnData file: {}", path.display()))?;
        // Index first (creates the obs/var elements), then the columns.
        adata.set_obs_names(obs_index.into_iter().collect())?;
        adata.set_obs(obs_df)?;
        adata.set_var_names(var_index.into_iter().collect())?;
        adata.set_var(var_df)?;
        adata.close()?;
        tag_anndata_root(&path)?;
        Ok(())
    }
}

/// Write a whole cell × feature count matrix to an AnnData HDF5 file.
///
/// Test-only: the counting commands add one block per BAM through
/// [`count_into_anndata`].
#[cfg(test)]
pub(crate) fn write_counts_anndata(
    output_path: &Path,
    matrix: CsrMatrix<u32>,
    cells: &[(String, String)],
    var: &[Feature],
    compression: &str,
    compression_level: u8,
) -> Result<()> {
    let mut writer = CountsWriter::create(output_path, var.len(), compression, compression_level)?;
    writer.append(matrix)?;
    writer.finish(cells, var)
}

/// Count every BAM, and add each BAM's rows to the output in BAM order.
///
/// The threads take the chunks BAM after BAM, so the next BAM starts as soon
/// as the last chunk of the one before it is taken. The thread that counts the
/// last chunk of a BAM adds that BAM's rows, while the other threads count the
/// next BAM, so usually only two BAMs' counts are in memory.
///
/// The rows are ordered by sample, and a BAM holds one sample, or under
/// `--groupTag` all of them, so the BAMs' rows follow each other in BAM order.
/// `count_chunk` counts one chunk into entries whose cells are numbered as for
/// the whole run, except that without a whitelist each BAM numbers its
/// barcodes in its own `run_barcodes[bam_idx]`. Returns the number of cells
/// written.
pub(super) fn count_into_anndata<'a, F>(
    output_path: &Path,
    compression: &str,
    compression_level: u8,
    var: &[Feature],
    pool: &ThreadPool,
    work: &[Chunk<'a>],
    samples: &Samples,
    barcodes: Option<&[String]>,
    run_barcodes: &[Mutex<BarcodeNumbers>],
    count_chunk: F,
) -> Result<usize>
where
    F: Fn(&mut BamWorker<'a>, &Chunk<'a>) -> Result<Vec<Entry>> + Sync,
{
    let writer = CountsWriter::create(output_path, var.len(), compression, compression_level)?;
    let labels = samples.labels();
    let blocks: Vec<(usize, Range<usize>)> = if samples.by_read_group() {
        vec![(0, 0..labels.len())]
    } else {
        (0..labels.len())
            .map(|bam_idx| (bam_idx, bam_idx..bam_idx + 1))
            .collect()
    };

    let chunks: Vec<(usize, &Chunk<'a>)> = blocks
        .iter()
        .enumerate()
        .flat_map(|(block, &(bam_idx, _))| {
            work.iter()
                .filter(move |c| c.bam_idx == bam_idx)
                .map(move |c| (block, c))
        })
        .collect();
    let mut left = vec![0usize; blocks.len()];
    for &(block, _) in &chunks {
        left[block] += 1;
    }

    let add_block = |(writer, cells): &mut (CountsWriter, Vec<(String, String)>),
                     block: usize,
                     mut entries: Vec<Vec<Entry>>|
     -> Result<()> {
        let (bam_idx, bam_samples) = blocks[block].clone();
        // With a whitelist every sample x barcode is a row; without one, only
        // the (sample, barcode) pairs that have counts.
        let bam_cells = match barcodes {
            Some(barcodes) => {
                let first_row = (bam_samples.start * barcodes.len()) as u32;
                for entry in entries.iter_mut().flatten() {
                    entry.0 -= first_row;
                }
                product_cells(&labels[bam_samples], barcodes)
            }
            None => {
                let found = std::mem::take(&mut *lock(&run_barcodes[bam_idx])).into_barcodes();
                observed_rows(&mut entries, labels, &found)
            }
        };
        writer.append(build_csr(entries, bam_cells.len(), var.len())?)?;
        cells.extend(bam_cells);
        Ok(())
    };

    let next = AtomicUsize::new(0);
    let progress = Mutex::new(Progress {
        counts: vec![Vec::new(); blocks.len()],
        left,
        first: 0,
        adding: false,
        failed: None,
    });
    let output = Mutex::new((writer, Vec::new()));
    // Add every BAM whose chunks are all counted, in BAM order. One thread at
    // a time does this, while the others count on.
    let add_ready = || {
        let mut state = lock(&progress);
        if state.adding {
            return;
        }
        state.adding = true;
        while state.failed.is_none() && state.first < blocks.len() && state.left[state.first] == 0 {
            let block = state.first;
            let entries = std::mem::take(&mut state.counts[block]);
            drop(state);
            let added = add_block(&mut lock(&output), block, entries);
            state = lock(&progress);
            state.first += 1;
            if let Err(err) = added {
                state.failed = Some(err);
                next.store(chunks.len(), Ordering::Relaxed);
            }
        }
        state.adding = false;
    };

    pool.broadcast(|_| {
        let mut worker = BamWorker::default();
        while let Some(&(block, chunk)) = chunks.get(next.fetch_add(1, Ordering::Relaxed)) {
            let counted = count_chunk(&mut worker, chunk);
            let mut state = lock(&progress);
            match counted {
                Ok(entries) => {
                    state.counts[block].push(entries);
                    state.left[block] -= 1;
                }
                Err(err) => {
                    state.failed.get_or_insert(err);
                    next.store(chunks.len(), Ordering::Relaxed);
                }
            }
            drop(state);
            add_ready();
        }
    });
    // A run whose BAMs have no chunks at all is added here.
    add_ready();

    let (writer, cells) = output.into_inner().unwrap_or_else(PoisonError::into_inner);
    let written = match progress
        .into_inner()
        .unwrap_or_else(PoisonError::into_inner)
        .failed
    {
        Some(err) => Err(err),
        None => writer.finish(&cells, var),
    };
    if written.is_err() {
        let _ = std::fs::remove_file(output_path);
    }
    written.map(|()| cells.len())
}

/// Counted chunks waiting for the rest of their BAM, with the number of chunks
/// each BAM still waits for, and the first BAM whose rows are not written yet.
struct Progress {
    counts: Vec<Vec<Vec<Entry>>>,
    left: Vec<usize>,
    first: usize,
    adding: bool,
    failed: Option<anyhow::Error>,
}

fn lock<T>(mutex: &Mutex<T>) -> MutexGuard<'_, T> {
    mutex.lock().unwrap_or_else(PoisonError::into_inner)
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let m = build_csr(vec![vec![(0, 1, 5), (1, 0, 3)], vec![(1, 2, 7)]], 2, 3).unwrap();

        assert_eq!((m.nrows(), m.ncols()), (2, 3));
        assert_eq!(m.nnz(), 3);
        assert_eq!(dense(&m), vec![vec![0, 5, 0], vec![3, 0, 7]]);
    }

    #[test]
    fn build_csr_sorts_column_indices_within_each_row() {
        // CSR requires ascending column indices per row; the entries come from
        // unordered hash maps, so the sort has to happen in build_csr.
        let entries = vec![(0..8).map(|c| (0, 7 - c, 1)).collect()];

        let m = build_csr(entries, 1, 8).unwrap();

        let cols: Vec<usize> = m.col_indices().to_vec();
        assert_eq!(cols, (0..8).collect::<Vec<_>>());
    }

    #[test]
    fn build_csr_keeps_empty_rows_and_the_declared_shape() {
        // Cells with no counts must still occupy a row: obs is the full
        // sample × barcode product regardless of coverage.
        let m = build_csr(vec![vec![(2, 0, 4)]], 4, 2).unwrap();

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
        let m = build_csr(vec![], 3, 5).unwrap();

        assert_eq!((m.nrows(), m.ncols()), (3, 5));
        assert_eq!(m.nnz(), 0);
    }

    #[test]
    fn build_csr_adds_entries_of_one_cell_and_feature_from_two_chunks() {
        let entries = vec![vec![(0, 2, 1), (1, 0, 3)], vec![(0, 2, 4), (0, 1, 2)]];

        let m = build_csr(entries, 2, 3).unwrap();

        assert_eq!(m.nnz(), 3);
        assert_eq!(m.row_offsets(), &[0, 2, 3]);
        assert_eq!(dense(&m), vec![vec![0, 2, 5], vec![3, 0, 0]]);
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
            // Bins and nameless annotation records both arrive unnamed.
            name: None,
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

        let matrix = build_csr(vec![vec![(0, 0, 1), (0, 2, 5), (1, 1, 3)]], 2, 3).unwrap();

        let var = vec![
            feature("chr1", 0, 100),
            feature("chr1", 100, 200),
            feature("chr2", 0, 100),
        ];
        let barcodes = vec!["AAA".to_string(), "CCC".to_string()];

        write_counts_anndata(
            &path,
            matrix,
            &product_cells(&["s1".to_string()], &barcodes),
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
    fn blocks_added_one_after_another_read_back_as_one_matrix() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("blocks.h5ad");
        let var = vec![
            feature("chr1", 0, 100),
            feature("chr1", 100, 200),
            feature("chr2", 0, 100),
        ];
        let barcodes = vec!["AAA".to_string(), "CCC".to_string()];

        let mut writer = CountsWriter::create(&path, 3, "none", 0).unwrap();
        writer
            .append(build_csr(vec![vec![(0, 2, 5), (1, 0, 1)]], 2, 3).unwrap())
            .unwrap();
        writer.append(build_csr(vec![], 0, 3).unwrap()).unwrap();
        writer.append(build_csr(vec![], 1, 3).unwrap()).unwrap();
        writer
            .append(build_csr(vec![vec![(0, 1, 7)]], 1, 3).unwrap())
            .unwrap();
        writer
            .finish(
                &product_cells(&["s1".to_string(), "s2".to_string()], &barcodes),
                &var,
            )
            .unwrap();

        let x = read_back(&path);
        let mut got = vec![vec![0.0f64; 3]; x.nrows()];
        for (r, c, &v) in x.triplet_iter() {
            got[r][c] = v;
        }
        assert_eq!(
            got,
            vec![
                vec![0.0, 0.0, 5.0],
                vec![1.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0],
                vec![0.0, 7.0, 0.0]
            ]
        );
        let adata = AnnData::<H5>::open(H5::open(&path).unwrap()).unwrap();
        assert_eq!(
            adata.obs_names().into_vec(),
            vec!["s1::AAA", "s1::CCC", "s2::AAA", "s2::CCC"]
        );
    }

    #[test]
    fn a_writer_rejects_cells_that_do_not_match_its_rows() {
        let dir = tempfile::TempDir::new().unwrap();
        let mut writer =
            CountsWriter::create(&dir.path().join("short.h5ad"), 1, "none", 0).unwrap();
        writer.append(build_csr(vec![], 2, 1).unwrap()).unwrap();

        let err = writer
            .finish(
                &product_cells(&["s1".to_string()], &["AAA".to_string()]),
                &[feature("chr1", 0, 100)],
            )
            .unwrap_err()
            .to_string();

        assert!(err.contains("1 cells") && err.contains("2 rows"), "{err}");
    }

    #[test]
    fn obs_rows_are_one_per_sample_and_barcode_in_matrix_order() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("counts.h5ad");

        // Two samples x two barcodes = four rows, samples varying slowest.
        let matrix = build_csr(vec![], 4, 1).unwrap();
        let barcodes = vec!["AAA".to_string(), "CCC".to_string()];

        write_counts_anndata(
            &path,
            matrix,
            &product_cells(&["s1".to_string(), "s2".to_string()], &barcodes),
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
        assert_eq!(adata.var_names().into_vec(), vec!["chr1_0_100::None"]);
    }

    #[test]
    fn the_var_table_keeps_the_feature_coordinates() {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("counts.h5ad");

        let var = vec![feature("chr1", 0, 100), feature("chr2", 500, 750)];
        let matrix = build_csr(vec![], 1, 2).unwrap();

        write_counts_anndata(
            &path,
            matrix,
            &product_cells(&["s1".to_string()], &["AAA".to_string()]),
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

        let matrix = build_csr(vec![vec![(0, 0, 42)]], 1, 1).unwrap();

        write_counts_anndata(
            &path,
            matrix,
            &product_cells(&["s1".to_string()], &["AAA".to_string()]),
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
            build_csr(vec![], 1, 1).unwrap(),
            &product_cells(&["s1".to_string()], &["AAA".to_string()]),
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
