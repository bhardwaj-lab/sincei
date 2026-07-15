use std::path::PathBuf;

use anyhow::Result;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

use crate::counting::filters::{DupMethod, QcFilter, RawRecordFilter};
use crate::counting::params::CountingParams;
use crate::counting::{count_bam_bins, count_bam_features};

fn parse_dup_method(s: &str) -> Result<DupMethod> {
    match s {
        "barcode_start" => Ok(DupMethod::BarcodeStart),
        "barcode_start_end" => Ok(DupMethod::BarcodeStartEnd),
        "barcode_umi_start" => Ok(DupMethod::BarcodeUmiStart),
        "barcode_umi_start_end" => Ok(DupMethod::BarcodeUmiStartEnd),
        _ => anyhow::bail!(
            "unknown dup_method {:?}; expected one of: \
             barcode_start, barcode_start_end, barcode_umi_start, barcode_umi_start_end",
            s
        ),
    }
}

fn build_qc_filter(
    min_fragment_length: Option<usize>,
    max_fragment_length: Option<usize>,
    min_gc: Option<f32>,
    max_gc: Option<f32>,
    min_aligned_fraction: Option<f32>,
) -> Option<QcFilter> {
    if min_fragment_length.is_none()
        && max_fragment_length.is_none()
        && min_gc.is_none()
        && max_gc.is_none()
        && min_aligned_fraction.is_none()
    {
        return None;
    }
    Some(QcFilter {
        min_fragment_length,
        max_fragment_length,
        min_gc,
        max_gc,
        min_aligned_fraction,
    })
}

fn build_record_filter(
    min_mapq: Option<u8>,
    sam_flag_include: Option<u16>,
    sam_flag_exclude: Option<u16>,
    filter_rna_strand: Option<String>,
) -> Option<RawRecordFilter> {
    if min_mapq.is_none()
        && sam_flag_include.is_none()
        && sam_flag_exclude.is_none()
        && filter_rna_strand.is_none()
    {
        return None;
    }
    Some(RawRecordFilter {
        min_mapq,
        sam_flag_include,
        sam_flag_exclude,
        filter_rna_strand,
    })
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
    umi_tag = None,
    count_tag = None,
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
    barcodes: Vec<String>,
    output_path: PathBuf,
    bin_size: usize,
    step_size: usize,
    bc_tag: &str,
    umi_tag: Option<String>,
    count_tag: Option<String>,
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
        extend_reads,
        center_reads,
        feature_type: None,
        exon_type: None,
        name_attr: None,
        metagene: false,
    };

    let qc = build_qc_filter(
        min_fragment_length,
        max_fragment_length,
        min_gc,
        max_gc,
        min_aligned_fraction,
    );

    let record_filter = build_record_filter(
        min_mapq,
        sam_flag_include,
        sam_flag_exclude,
        filter_rna_strand,
    );

    let dup = dup_method
        .as_deref()
        .map(parse_dup_method)
        .transpose()
        .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?;

    // Build (path, sample_name) pairs. sample name is the file stem.
    let path_sample: Vec<(PathBuf, String)> = bam_paths
        .iter()
        .map(|p| {
            let stem = p
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string();
            (p.clone(), stem)
        })
        .collect();
    let bam_path_refs: Vec<(&std::path::Path, &str)> = path_sample
        .iter()
        .map(|(p, s)| (p.as_path(), s.as_str()))
        .collect();

    count_bam_bins(
        &bam_path_refs,
        bin_size,
        step_size,
        &barcodes,
        bc_tag,
        umi_tag.as_deref(),
        count_tag.as_deref(),
        &params,
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
    .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))
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
    umi_tag = None,
    count_tag = None,
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
    barcodes: Vec<String>,
    output_path: PathBuf,
    bc_tag: &str,
    umi_tag: Option<String>,
    count_tag: Option<String>,
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
        extend_reads,
        center_reads,
        feature_type,
        exon_type,
        name_attr,
        metagene,
    };

    let qc = build_qc_filter(
        min_fragment_length,
        max_fragment_length,
        min_gc,
        max_gc,
        min_aligned_fraction,
    );

    let record_filter = build_record_filter(
        min_mapq,
        sam_flag_include,
        sam_flag_exclude,
        filter_rna_strand,
    );

    let dup = dup_method
        .as_deref()
        .map(parse_dup_method)
        .transpose()
        .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?;

    let path_sample: Vec<(PathBuf, String)> = bam_paths
        .iter()
        .map(|p| {
            let stem = p
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string();
            (p.clone(), stem)
        })
        .collect();
    let bam_path_refs: Vec<(&std::path::Path, &str)> = path_sample
        .iter()
        .map(|(p, s)| (p.as_path(), s.as_str()))
        .collect();

    count_bam_features(
        &bam_path_refs,
        annotation_path.as_path(),
        &barcodes,
        bc_tag,
        umi_tag.as_deref(),
        count_tag.as_deref(),
        &params,
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
    .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))
}
