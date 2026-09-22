use std::path::{Path, PathBuf};

use anyhow::Result;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

use crate::bam::filters::{DupMethod, QcFilter, RawRecordFilter, RnaStrand};
use crate::counting::coverage::{NormalizeMethod, OutputFormat, ReadMode, run_bulk_coverage};
use crate::to_py_err;

fn parse_normalize_method(s: &str) -> Result<NormalizeMethod> {
    match s {
        "CPM" => Ok(NormalizeMethod::Cpm),
        "RPKM" => Ok(NormalizeMethod::Rpkm),
        "Frequency" => Ok(NormalizeMethod::Frequency),
        "Mean" => Ok(NormalizeMethod::Mean),
        "None" => Ok(NormalizeMethod::None),
        _ => anyhow::bail!(
            "unknown normalize_using {:?}; expected one of: CPM, RPKM, Frequency, Mean, None",
            s
        ),
    }
}

/// Compute pseudo-bulk coverage tracks, one bigWig (or bedGraph) per cell group.
///
/// `bam_files` and `bam_labels` must be the same length; each label must match
/// the `sample` column in `group_info`.
///
/// `group_info` is the path to a TSV file with a header line, in either of two
/// layouts, told apart by the header's column count: `sample`, `barcode`,
/// `group`; or `sample::barcode`, `UMAP1`, `UMAP2`, `group` as
/// `scClusterCells` writes it. `None` pools every read that carries `bc_tag`
/// into one track, `{output_prefix}.{ext}`, whatever its barcode (a read without
/// the tag is not counted), and rejects Mean and Frequency normalization.
///
/// Returns the list of output file paths created.
#[pyfunction(signature = (
    bam_files,
    bam_labels,
    group_info,
    output_prefix,
    bin_size = 100,
    step_size = 100,
    bc_tag = "CB",
    umi_tag = None,
    group_tag = None,
    region = None,
    min_mapq = None,
    sam_flag_include = None,
    sam_flag_exclude = None,
    chr_to_skip = vec![],
    ignore_for_normalization = vec![],
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
    normalize_using = "CPM",
    scale_factor = 1.0,
    out_format = "bigwig",
    mnase = false,
    offset = None,
    filter_rna_strand = None,
    num_threads = 0,
    chunk_size = 1_000_000,
))]
pub fn bulk_coverage(
    bam_files: Vec<PathBuf>,
    bam_labels: Vec<String>,
    group_info: Option<PathBuf>,
    output_prefix: String,
    bin_size: usize,
    step_size: usize,
    bc_tag: &str,
    umi_tag: Option<String>,
    group_tag: Option<String>,
    region: Option<String>,
    min_mapq: Option<u8>,
    sam_flag_include: Option<u16>,
    sam_flag_exclude: Option<u16>,
    chr_to_skip: Vec<String>,
    ignore_for_normalization: Vec<String>,
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
    normalize_using: &str,
    scale_factor: f64,
    out_format: &str,
    mnase: bool,
    offset: Option<Vec<i32>>,
    filter_rna_strand: Option<String>,
    num_threads: usize,
    chunk_size: usize,
) -> PyResult<Vec<PathBuf>> {
    if bam_files.len() != bam_labels.len() {
        return Err(PyRuntimeError::new_err(
            "bam_files and bam_labels must have the same length",
        ));
    }

    let normalize = parse_normalize_method(normalize_using).map_err(to_py_err)?;

    let format = match out_format {
        "bigwig" | "bw" => OutputFormat::BigWig,
        "bedgraph" | "bg" => OutputFormat::BedGraph,
        _ => {
            return Err(PyRuntimeError::new_err(
                "out_format must be 'bigwig' or 'bedgraph'",
            ));
        }
    };

    let dup = dup_method
        .as_deref()
        .map(str::parse::<DupMethod>)
        .transpose()
        .map_err(to_py_err)?;

    // MNase defaults: fragment length 130–200 bp unless explicitly overridden.
    let min_fragment_length = if mnase && min_fragment_length.is_none() {
        Some(130)
    } else {
        min_fragment_length
    };
    let max_fragment_length = if mnase && max_fragment_length.is_none() {
        Some(200)
    } else {
        max_fragment_length
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

    // Validate and parse offset.
    let parsed_offset: Option<(i32, Option<i32>)> = match offset.as_deref() {
        None => None,
        Some([s]) => {
            if *s == 0 {
                return Err(PyRuntimeError::new_err(
                    "offset value 0 is not allowed (offsets are 1-based)",
                ));
            }
            Some((*s, None))
        }
        Some([s, e]) => {
            if *s == 0 || *e == 0 {
                return Err(PyRuntimeError::new_err(
                    "offset value 0 is not allowed (offsets are 1-based)",
                ));
            }
            if *e > 0 && *e < *s {
                return Err(PyRuntimeError::new_err(
                    "offset end must be >= offset start",
                ));
            }
            Some((*s, Some(*e)))
        }
        Some(_) => return Err(PyRuntimeError::new_err("offset must have 1 or 2 elements")),
    };

    if mnase && parsed_offset.is_some() {
        return Err(PyRuntimeError::new_err(
            "--mnase and --offset are mutually exclusive",
        ));
    }

    let read_mode = if mnase {
        ReadMode::MNase
    } else if let Some((s, e)) = parsed_offset {
        ReadMode::Offset(s, e)
    } else {
        ReadMode::Normal
    };

    let path_label: Vec<(PathBuf, String)> = bam_files.into_iter().zip(bam_labels).collect();
    let bam_path_refs: Vec<(&Path, &str)> = path_label
        .iter()
        .map(|(p, l)| (p.as_path(), l.as_str()))
        .collect();

    run_bulk_coverage(
        &bam_path_refs,
        group_info.as_deref(),
        &output_prefix,
        bin_size,
        step_size,
        bc_tag,
        umi_tag.as_deref(),
        group_tag.as_deref(),
        region.as_deref(),
        &chr_to_skip,
        &ignore_for_normalization,
        blacklist_path.as_deref(),
        extend_reads,
        center_reads,
        dup,
        genome_2bit.as_deref(),
        motif_filter.as_deref(),
        record_filter.as_ref(),
        qc.as_ref(),
        normalize,
        scale_factor,
        format,
        read_mode,
        num_threads,
        chunk_size,
    )
    .map_err(to_py_err)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_documented_normalize_method_parses() {
        assert!(matches!(
            parse_normalize_method("CPM").unwrap(),
            NormalizeMethod::Cpm
        ));
        assert!(matches!(
            parse_normalize_method("RPKM").unwrap(),
            NormalizeMethod::Rpkm
        ));
        assert!(matches!(
            parse_normalize_method("Frequency").unwrap(),
            NormalizeMethod::Frequency
        ));
        assert!(matches!(
            parse_normalize_method("Mean").unwrap(),
            NormalizeMethod::Mean
        ));
        assert!(matches!(
            parse_normalize_method("None").unwrap(),
            NormalizeMethod::None
        ));
    }

    #[test]
    fn an_unknown_normalize_method_names_the_valid_ones() {
        // The match is case-sensitive, so the lowercase spelling is an error.
        let err = parse_normalize_method("cpm").unwrap_err().to_string();
        assert!(err.contains("CPM"), "{err}");
        assert!(err.contains("RPKM"), "{err}");
    }
}
