use std::path::PathBuf;

use anyhow::Result;
use pyo3::prelude::*;

use crate::bam::filters::{DupMethod, QcFilter, RawRecordFilter, RnaStrand};
use crate::bam::fragment_length::resolve_extend_reads;
use crate::bam::sc_record::AdjustRead;
use crate::counting::params::CountingParams;
use crate::counting::{count_bam_bins, count_bam_features};
use crate::to_py_err;

/// Pair each BAM with the sample name its matrix rows are labelled by.
///
/// An empty `labels` means derive the name from the file stem, which is what
/// both the default and `--smartLabels` ask for. Otherwise the supplied labels
/// are used verbatim, one per BAM.
fn sample_names(bam_paths: &[PathBuf], labels: Vec<String>) -> Result<Vec<(PathBuf, String)>> {
    if labels.is_empty() {
        return Ok(bam_paths
            .iter()
            .map(|p| {
                let stem = p
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("unknown")
                    .to_string();
                (p.clone(), stem)
            })
            .collect());
    }

    // The CLI checks this too; repeated here so the backend cannot be driven
    // into a mislabelled matrix by a direct call.
    anyhow::ensure!(
        labels.len() == bam_paths.len(),
        "got {} labels for {} BAM files; there must be one label per file",
        labels.len(),
        bam_paths.len()
    );
    Ok(bam_paths.iter().cloned().zip(labels).collect())
}

/// Count reads into a cell × genomic-bin matrix and write the result as an
/// AnnData HDF5 file.
///
/// `bam_paths` may contain multiple BAM files; the resulting AnnData will have
/// one row per (sample × barcode) combination. Requires a BAI index alongside
/// each BAM. `num_threads = 0` uses all available cores.
#[pyfunction(signature = (
    bam_paths,
    barcodes,
    output_path,
    bin_size = 10_000,
    step_size = 10_000,
    bc_tag = "CB",
    labels = vec![],
    umi_tag = None,
    count_tag = None,
    group_tag = None,
    min_mapq = None,
    sam_flag_include = None,
    sam_flag_exclude = None,
    chr_to_skip = vec![],
    region = None,
    blacklist_path = None,
    extend_reads = None,
    center_reads = false,
    dup_method = None,
    genome_2bit = None,
    motif_filter = None,
    min_gc = None,
    max_gc = None,
    min_aligned_fraction = None,
    min_fragment_length = None,
    max_fragment_length = None,
    filter_rna_strand = None,
    compression = "none",
    compression_level = 4,
    num_threads = 0,
    chunk_size = 1_000_000,
))]
pub fn count_bins(
    bam_paths: Vec<PathBuf>,
    barcodes: Option<Vec<String>>,
    output_path: PathBuf,
    bin_size: usize,
    step_size: usize,
    bc_tag: &str,
    labels: Vec<String>,
    umi_tag: Option<String>,
    count_tag: Option<String>,
    group_tag: Option<String>,
    min_mapq: Option<u8>,
    sam_flag_include: Option<u16>,
    sam_flag_exclude: Option<u16>,
    chr_to_skip: Vec<String>,
    region: Option<String>,
    blacklist_path: Option<PathBuf>,
    extend_reads: Option<usize>,
    center_reads: bool,
    dup_method: Option<String>,
    genome_2bit: Option<PathBuf>,
    motif_filter: Option<Vec<(String, String)>>,
    min_gc: Option<f32>,
    max_gc: Option<f32>,
    min_aligned_fraction: Option<f32>,
    min_fragment_length: Option<usize>,
    max_fragment_length: Option<usize>,
    filter_rna_strand: Option<String>,
    compression: &str,
    compression_level: u8,
    num_threads: usize,
    chunk_size: usize,
) -> PyResult<()> {
    let params = CountingParams {
        chr_to_skip,
        region,
        blacklist_path,
        feature_type: None,
        exon_type: None,
        name_attr: None,
        metagene: false,
    };

    let qc = QcFilter::from_bounds(
        min_fragment_length,
        max_fragment_length,
        min_gc,
        max_gc,
        min_aligned_fraction,
    );

    let filter_rna_strand = filter_rna_strand
        .as_deref()
        .map(str::parse::<RnaStrand>)
        .transpose()
        .map_err(to_py_err)?;

    let record_filter = RawRecordFilter::from_options(
        min_mapq,
        sam_flag_include,
        sam_flag_exclude,
        filter_rna_strand,
    );

    let dup = dup_method
        .as_deref()
        .map(str::parse::<DupMethod>)
        .transpose()
        .map_err(to_py_err)?;

    let path_sample: Vec<(PathBuf, String)> =
        sample_names(&bam_paths, labels).map_err(to_py_err)?;
    let bam_path_refs: Vec<(&std::path::Path, &str)> = path_sample
        .iter()
        .map(|(p, s)| (p.as_path(), s.as_str()))
        .collect();

    let adjust = AdjustRead {
        extend_reads: resolve_extend_reads(extend_reads, &bam_path_refs).map_err(to_py_err)?,
        center_reads,
        max_paired_fragment_length: max_fragment_length,
    };

    count_bam_bins(
        &bam_path_refs,
        bin_size,
        step_size,
        barcodes.as_deref(),
        bc_tag,
        umi_tag.as_deref(),
        count_tag.as_deref(),
        group_tag.as_deref(),
        &params,
        &adjust,
        record_filter.as_ref(),
        qc.as_ref(),
        dup,
        genome_2bit.as_deref(),
        motif_filter.as_deref(),
        output_path.as_path(),
        compression,
        compression_level,
        num_threads,
        chunk_size,
    )
    .map_err(to_py_err)
}

/// Count reads into a cell × genomic-feature matrix and write the result as an
/// AnnData HDF5 file.
///
/// `bam_paths` may contain multiple BAM files; the resulting AnnData will have
/// one row per (sample × barcode) combination. Requires a BAI index alongside
/// each BAM. `num_threads = 0` uses all available cores.
#[pyfunction(signature = (
    bam_paths,
    annotation_path,
    barcodes,
    output_path,
    bc_tag = "CB",
    labels = vec![],
    umi_tag = None,
    count_tag = None,
    group_tag = None,
    min_mapq = None,
    sam_flag_include = None,
    sam_flag_exclude = None,
    chr_to_skip = vec![],
    region = None,
    blacklist_path = None,
    extend_reads = None,
    center_reads = false,
    feature_type = None,
    exon_type = None,
    name_attr = None,
    metagene = false,
    dup_method = None,
    filter_rna_strand = None,
    genome_2bit = None,
    motif_filter = None,
    min_gc = None,
    max_gc = None,
    min_aligned_fraction = None,
    min_fragment_length = None,
    max_fragment_length = None,
    compression = "none",
    compression_level = 4,
    num_threads = 0,
    chunk_size = 1_000_000,
))]
pub fn count_features(
    bam_paths: Vec<PathBuf>,
    annotation_path: PathBuf,
    barcodes: Option<Vec<String>>,
    output_path: PathBuf,
    bc_tag: &str,
    labels: Vec<String>,
    umi_tag: Option<String>,
    count_tag: Option<String>,
    group_tag: Option<String>,
    min_mapq: Option<u8>,
    sam_flag_include: Option<u16>,
    sam_flag_exclude: Option<u16>,
    chr_to_skip: Vec<String>,
    region: Option<String>,
    blacklist_path: Option<PathBuf>,
    extend_reads: Option<usize>,
    center_reads: bool,
    feature_type: Option<Vec<String>>,
    exon_type: Option<Vec<String>>,
    name_attr: Option<String>,
    metagene: bool,
    dup_method: Option<String>,
    filter_rna_strand: Option<String>,
    genome_2bit: Option<PathBuf>,
    motif_filter: Option<Vec<(String, String)>>,
    min_gc: Option<f32>,
    max_gc: Option<f32>,
    min_aligned_fraction: Option<f32>,
    min_fragment_length: Option<usize>,
    max_fragment_length: Option<usize>,
    compression: &str,
    compression_level: u8,
    num_threads: usize,
    chunk_size: usize,
) -> PyResult<()> {
    let params = CountingParams {
        chr_to_skip,
        region,
        blacklist_path,
        feature_type,
        exon_type,
        name_attr,
        metagene,
    };

    let qc = QcFilter::from_bounds(
        min_fragment_length,
        max_fragment_length,
        min_gc,
        max_gc,
        min_aligned_fraction,
    );

    let filter_rna_strand = filter_rna_strand
        .as_deref()
        .map(str::parse::<RnaStrand>)
        .transpose()
        .map_err(to_py_err)?;

    let record_filter = RawRecordFilter::from_options(
        min_mapq,
        sam_flag_include,
        sam_flag_exclude,
        filter_rna_strand,
    );

    let dup = dup_method
        .as_deref()
        .map(str::parse::<DupMethod>)
        .transpose()
        .map_err(to_py_err)?;

    let path_sample: Vec<(PathBuf, String)> =
        sample_names(&bam_paths, labels).map_err(to_py_err)?;
    let bam_path_refs: Vec<(&std::path::Path, &str)> = path_sample
        .iter()
        .map(|(p, s)| (p.as_path(), s.as_str()))
        .collect();

    let adjust = AdjustRead {
        extend_reads: resolve_extend_reads(extend_reads, &bam_path_refs).map_err(to_py_err)?,
        center_reads,
        max_paired_fragment_length: max_fragment_length,
    };

    count_bam_features(
        &bam_path_refs,
        annotation_path.as_path(),
        barcodes.as_deref(),
        bc_tag,
        umi_tag.as_deref(),
        count_tag.as_deref(),
        group_tag.as_deref(),
        &params,
        &adjust,
        record_filter.as_ref(),
        qc.as_ref(),
        dup,
        genome_2bit.as_deref(),
        motif_filter.as_deref(),
        output_path.as_path(),
        compression,
        compression_level,
        num_threads,
        chunk_size,
    )
    .map_err(to_py_err)
}
