use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write as IoWrite};
use std::path::{Path, PathBuf};

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result};
use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BigWigWrite, Value};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::counting::bam_io::read_bam_header;
use crate::counting::count_utils::{BamWorker, derive_record_opts, effective_interval};
use crate::counting::filters::{DupMethod, DuplicateFilter, QcFilter};
use crate::counting::params::CountingParams;
use crate::counting::parse_annotation::parse_bed_file;
use crate::counting::region_index::build_bin_index;
use crate::counting::sc_record::{ScRecord, parse_tag};

#[derive(Clone, Copy, Debug)]
pub enum NormalizeMethod {
    Cpm,
    Rpkm,
    Frequency,
    Mean,
    None,
}

#[derive(Clone, Copy, Debug)]
pub enum OutputFormat {
    BigWig,
    BedGraph,
}

/// How to derive the signal interval from each read.
#[derive(Clone, Copy, Debug)]
pub enum ReadMode {
    /// Standard behavior: apply extend_reads / center_reads from CountingParams.
    Normal,
    /// MNase / CUT&RUN: use only the 2–3 central bp of the paired-end fragment.
    /// Reads that are not proper-pair forward reads return `None` and are skipped.
    MNase,
    /// Report only a specific offset position (or range) within the read.
    /// Values are 1-based; negative values count from the read end.
    Offset(i32, Option<i32>),
}

/// Which RNA strand to keep (for strand-specific RNA-seq libraries).
#[derive(Clone, Copy, Debug)]
pub enum RnaStrand {
    Forward,
    Reverse,
}

// Read-mode helpers

/// Returns the 2–3 central bp of a paired-end fragment (MNase / CUT&RUN mode).
/// Only proper-pair, non-reverse reads are considered; returns `None` otherwise.
fn apply_mnase(rec: &ScRecord) -> Option<(usize, usize)> {
    let tlen = rec.template_length.unsigned_abs() as usize;
    if !rec.is_proper_pair || rec.is_reverse || tlen <= 1 {
        return None;
    }
    // Python: fragment_start = read.pos + read.tlen / 2 - 1
    // Integer tlen/2 matches Python 3 float division here because both round down.
    let frag_start = rec.alignment_start + tlen / 2 - 1;
    let frag_end = frag_start + if tlen % 2 == 0 { 2 } else { 3 };
    Some((frag_start, frag_end))
}

/// Returns the genomic interval corresponding to an offset (or offset range)
/// within the read.  Offsets are 1-based; negative values index from the read
/// end.  A single value returns a 1-bp interval; two values return a range.
///
/// Simplified: treats each read as a contiguous block (no CIGAR-block
/// traversal), which is correct for unspliced data.
fn apply_offset(rec: &ScRecord, start: i32, end: Option<i32>) -> Option<(usize, usize)> {
    let len = rec.read_length as i32;
    if len == 0 { return None; }

    // Convert 1-based possibly-negative offset to 0-based index into the
    // conceptual stretch (list of read positions from 5′ to 3′).
    let to_idx = |o: i32| -> i32 {
        if o > 0 { o - 1 } else { len + o }
    };

    let idx_s = to_idx(start);
    let idx_e = match end {
        None => idx_s + 1,
        Some(e) => to_idx(e) + 1,
    };

    if idx_s < 0 || idx_e <= idx_s || idx_s >= len {
        return None;
    }
    let idx_s = idx_s as usize;
    let idx_e = (idx_e.min(len)) as usize;

    // For reverse reads the stretch is reversed (3′→5′ in genomic order),
    // so index 0 corresponds to the genomic end of the alignment.
    Some(if rec.is_reverse {
        let geno_s = rec.alignment_end.saturating_sub(idx_e);
        let geno_e = rec.alignment_end.saturating_sub(idx_s);
        (geno_s, geno_e)
    } else {
        (rec.alignment_start + idx_s, rec.alignment_start + idx_e)
    })
}

/// Returns `true` if the read passes the RNA-strand filter.
///
/// For paired-end reads this follows the dUTP convention (deepTools):
/// "forward" = R2 maps forward OR R1 mate maps forward.
/// For single-end reads the mapping strand is inverted relative to the RNA
/// strand (dUTP), so "forward" = read maps to reverse strand.
fn passes_rna_strand_filter(rec: &ScRecord, strand: RnaStrand) -> bool {
    if rec.is_paired {
        match strand {
            RnaStrand::Forward => {
                (!rec.is_read1 && !rec.is_reverse) || (rec.is_read1 && !rec.mate_is_reverse)
            }
            RnaStrand::Reverse => {
                (!rec.is_read1 && rec.is_reverse) || (rec.is_read1 && rec.mate_is_reverse)
            }
        }
    } else {
        match strand {
            RnaStrand::Forward => rec.is_reverse,
            RnaStrand::Reverse => !rec.is_reverse,
        }
    }
}

fn get_effective_interval(
    rec: &ScRecord,
    params: &CountingParams,
    mode: ReadMode,
) -> Option<(usize, usize)> {
    match mode {
        ReadMode::Normal => Some(effective_interval(rec, params)),
        ReadMode::MNase => apply_mnase(rec),
        ReadMode::Offset(start, end) => apply_offset(rec, start, end),
    }
}

// Group-info parsing

struct ParsedGroups {
    /// Ordered unique group names.
    groups: Vec<String>,
    /// Number of cells per group (indexed by group_idx).
    n_cells_per_group: Vec<usize>,
    /// All cells: (bam_idx, barcode, group_idx).
    cells: Vec<(usize, String, usize)>,
}

fn parse_group_info(
    path: &Path,
    bam_labels: &[&str],
) -> Result<ParsedGroups> {
    let file = File::open(path)
        .with_context(|| format!("failed to open group info file: {}", path.display()))?;
    let reader = BufReader::new(file);

    // label → bam index
    let label_to_bam: AHashMap<&str, usize> = bam_labels
        .iter()
        .enumerate()
        .map(|(i, l)| (*l, i))
        .collect();

    let mut group_index: AHashMap<String, usize> = AHashMap::new();
    let mut groups: Vec<String> = Vec::new();
    let mut cells: Vec<(usize, String, usize)> = Vec::new();

    let mut lines = reader.lines();
    // Skip header line.
    lines.next();

    for (line_no, line) in lines.enumerate() {
        let line = line.with_context(|| format!("failed to read group info line {}", line_no + 2))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let sample = fields.next().context("missing sample column")?.trim().to_string();
        let barcode = fields.next().context("missing barcode column")?.trim().to_string();
        let group = fields.next().context("missing group column")?.trim().to_string();

        let Some(&bam_idx) = label_to_bam.get(sample.as_str()) else {
            continue; // sample not in BAM list — skip
        };

        let group_idx = *group_index.entry(group.clone()).or_insert_with(|| {
            let idx = groups.len();
            groups.push(group);
            idx
        });

        cells.push((bam_idx, barcode, group_idx));
    }

    let n_groups = groups.len();
    let mut n_cells_per_group = vec![0usize; n_groups];
    for &(_, _, g) in &cells {
        n_cells_per_group[g] += 1;
    }

    Ok(ParsedGroups { groups, n_cells_per_group, cells })
}

// Core function

/// Compute pseudo-bulk coverage from one or more BAM files grouped by cell
/// identity, then write one bigWig or bedGraph per group.
///
/// `bam_paths` — list of `(path, sample_label)` pairs; labels must match the
/// `sample` column in `group_info_path`.
///
/// `group_info_path` — TSV with a header line followed by rows of
/// `sample\tbarcode\tgroup`.
///
/// Normalization denominators (CPM / RPKM) are computed from the total reads
/// over all bins, excluding chromosomes in `ignore_for_normalization`.
/// Chromosomes in `chr_to_skip` are not counted at all and are absent from the
/// output.
pub fn run_bulk_coverage(
    bam_paths: &[(&Path, &str)],
    group_info_path: &Path,
    output_prefix: &str,
    bin_size: usize,
    step_size: usize,
    bc_tag: &str,
    umi_tag: Option<&str>,
    min_mapq: Option<u8>,
    sam_flag_include: Option<u16>,
    sam_flag_exclude: Option<u16>,
    chr_to_skip: &[String],
    ignore_for_normalization: &[String],
    blacklist_path: Option<&Path>,
    extend_reads: Option<usize>,
    center_reads: bool,
    dup_method: Option<DupMethod>,
    genome_path: Option<&Path>,
    motifs: Option<&[(String, String)]>,
    qc_filter: Option<&QcFilter>,
    normalize_using: NormalizeMethod,
    scale_factor: f64,
    out_format: OutputFormat,
    read_mode: ReadMode,
    filter_rna_strand: Option<RnaStrand>,
    num_threads: usize,
    chunk_size: usize,
) -> Result<Vec<PathBuf>> {
    anyhow::ensure!(bin_size > 0, "bin_size must be greater than zero");
    anyhow::ensure!(step_size > 0, "step_size must be greater than zero");
    anyhow::ensure!(chunk_size > 0, "chunk_size must be greater than zero");
    anyhow::ensure!(!bam_paths.is_empty(), "at least one BAM file is required");

    let bam_labels: Vec<&str> = bam_paths.iter().map(|(_, l)| *l).collect();
    let parsed = parse_group_info(group_info_path, &bam_labels)?;

    let n_groups = parsed.groups.len();
    let n_cells = parsed.cells.len();

    if n_cells == 0 {
        anyhow::bail!("no cells matched between group_info and BAM labels");
    }

    // cell_global_idx → group_idx
    let cell_group: Vec<usize> = parsed.cells.iter().map(|(_, _, g)| *g).collect();

    // (bam_idx, barcode bytes) → cell_global_idx — keys borrow from parsed.cells.
    // Byte-keyed so the per-read lookup never allocates.
    let cell_index: AHashMap<(usize, &[u8]), usize> = parsed
        .cells
        .iter()
        .enumerate()
        .map(|(i, (bam_idx, bc, _))| ((*bam_idx, bc.as_bytes()), i))
        .collect();

    let has_motif = genome_path.is_some() && motifs.is_some();
    let record_opts = derive_record_opts(qc_filter, has_motif);
    let bc_tag_parsed = parse_tag(bc_tag)?;
    let umi_tag_parsed = umi_tag.map(parse_tag).transpose()?;

    let params = CountingParams {
        sam_flag_include,
        sam_flag_exclude,
        chr_to_skip: chr_to_skip.to_vec(),
        region: None,
        blacklist_path: blacklist_path.map(|p| p.to_path_buf()),
        extend_reads,
        center_reads,
        feature_type: None,
        exon_type: None,
        name_attr: None,
    };

    // Chromosomes to skip, as a byte-slice set to avoid per-contig allocation.
    let skip_set: AHashSet<&[u8]> = chr_to_skip.iter().map(|s| s.as_bytes()).collect();

    // Bin index from first BAM header.
    let chrom_sizes: Vec<(String, usize)> = {
        let (first, _) = bam_paths[0];
        let hdr = read_bam_header(first)?;
        hdr.reference_sequences()
            .iter()
            .filter(|(name, _)| {
                let bytes: &[u8] = name.as_ref();
                !skip_set.contains(bytes)
            })
            .map(|(name, seq)| (name.to_string(), seq.length().get()))
            .collect()
    };

    let (bin_index, _feature_names) = build_bin_index(&chrom_sizes, bin_size, step_size);

    let blacklist = params
        .blacklist_path
        .as_deref()
        .map(|p| parse_bed_file(p).map(|(idx, _)| idx))
        .transpose()?;

    // Build chunk work list sorted by descending size (LPT).
    let mut work: Vec<(usize, &Path, String, usize, usize)> = bam_paths
        .iter()
        .enumerate()
        .flat_map(|(bam_idx, &(bam_path, _))| {
            chrom_sizes.iter().flat_map(move |(chrom, chrom_len)| {
                (0..*chrom_len).step_by(chunk_size).map(move |start| {
                    (bam_idx, bam_path, chrom.clone(), start, (start + chunk_size).min(*chrom_len))
                })
            })
        })
        .collect();
    work.sort_unstable_by(|a, b| (b.4 - b.3).cmp(&(a.4 - a.3)));

    let n_threads = if num_threads == 0 { rayon::current_num_threads() } else { num_threads };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .context("failed to build thread pool")?;

    // Motif-filter ingredients, built once per worker rather than per chunk.
    let motif_ingredients = match (genome_path, motifs) {
        (Some(g), Some(m)) => Some((g, m)),
        _ => None,
    };

    // Parallel phase: accumulate per-(cell_idx, bin_idx) counts.
    let partial_accumulators: Vec<AHashMap<(usize, usize), u32>> = pool.install(|| {
        work.par_iter()
            .map_init(
                BamWorker::new,
                |worker, &(bam_idx, bam_path, ref chrom, chunk_start, chunk_end)| -> Result<AHashMap<(usize, usize), u32>> {
                    worker.prepare(bam_path, motif_ingredients)?;
                    let (reader, header, motif) = worker.parts();

                    // Each work item covers a single chromosome, so its bin
                    // geometry and blacklist segments are resolved once here
                    // rather than per read.
                    let Some(&(chrom_offset, n_bins)) = bin_index.chrom_bins.get(chrom.as_str())
                    else {
                        return Ok(AHashMap::new());
                    };
                    let chunk_blacklist = blacklist.as_ref().and_then(|bl| {
                        bl.get(chrom.as_str())
                            .or_else(|| chrom.strip_prefix("chr").and_then(|c| bl.get(c)))
                            .or_else(|| bl.get(&format!("chr{}", chrom)))
                    });

                    let region_str = format!("{}:{}-{}", chrom, chunk_start + 1, chunk_end);
                    let region: noodles::core::Region = region_str
                        .parse()
                        .with_context(|| format!("failed to parse region: {}", region_str))?;
                    let query = match reader.query(header, &region) {
                        Ok(q) => q,
                        Err(_) => return Ok(AHashMap::new()),
                    };

                    let mut dup_filter: Option<DuplicateFilter> = dup_method.map(DuplicateFilter::new);
                    let mut local_acc: AHashMap<(usize, usize), u32> = AHashMap::new();

                    for result in query.records() {
                        let record = result.context("failed to read BAM record")?;

                        let Some(sc_rec) = ScRecord::from_bam_record(
                            &record,
                            header,
                            &bc_tag_parsed,
                            umi_tag_parsed.as_ref(),
                            None,
                            min_mapq,
                            params.sam_flag_include,
                            params.sam_flag_exclude,
                            &record_opts,
                        )? else { continue };

                        if sc_rec.alignment_start < chunk_start {
                            continue;
                        }

                        let Some(barcode) = sc_rec.barcode else { continue };
                        let Some(&cell_idx) = cell_index.get(&(bam_idx, barcode)) else { continue };

                        if let Some(qc) = qc_filter {
                            if !qc.passes(&sc_rec) { continue; }
                        }
                        if let Some(ref mut dup) = dup_filter {
                            if !dup.passes(&sc_rec) { continue; }
                        }
                        if let Some(mf) = motif.as_mut() {
                            if !mf.passes(&sc_rec, chrom)? { continue; }
                        }
                        if let Some(strand) = filter_rna_strand {
                            if !passes_rna_strand_filter(&sc_rec, strand) { continue; }
                        }

                        if let Some(idx) = chunk_blacklist {
                            if idx.find(sc_rec.alignment_start, sc_rec.alignment_end).next().is_some() {
                                continue;
                            }
                        }

                        let Some((eff_start, eff_end)) = get_effective_interval(&sc_rec, &params, read_mode) else { continue };

                        if eff_end > 0 && n_bins > 0 {
                            let last_bin = ((eff_end - 1) / step_size).min(n_bins - 1);
                            let first_bin = if eff_start + 1 > bin_size {
                                (eff_start - bin_size + step_size) / step_size
                            } else {
                                0
                            };
                            for bin in first_bin..=last_bin {
                                *local_acc.entry((cell_idx, chrom_offset + bin)).or_insert(0) += sc_rec.count;
                            }
                        }
                    }

                    Ok(local_acc)
                },
            )
            .collect::<Result<Vec<_>>>()
    })?;

    // Sequential merge.
    let mut global_acc: AHashMap<(usize, usize), u32> = AHashMap::new();
    for partial in partial_accumulators {
        for ((cell, bin), val) in partial {
            *global_acc.entry((cell, bin)).or_insert(0) += val;
        }
    }

    // Aggregate per (group, bin)

    // sum of counts per (group, bin)
    let mut group_bin_sum: AHashMap<(usize, usize), f64> = AHashMap::new();
    // number of cells with ≥1 read in bin, per (group, bin) — needed for Frequency
    let mut group_bin_n_cells: AHashMap<(usize, usize), usize> = AHashMap::new();
    // total reads per group (excluding ignored chroms) — needed for CPM/RPKM denominators
    let mut group_total: Vec<f64> = vec![0.0; n_groups];

    // Build a fast set for ignore_for_normalization
    let ignore_set: std::collections::HashSet<&str> = ignore_for_normalization
        .iter()
        .map(String::as_str)
        .collect();

    // For each bin, know its chromosome so we can apply ignore_for_normalization.
    // Build bin_idx → chrom mapping lazily from chrom_sizes + bin_index.
    let bin_to_chrom: AHashMap<usize, &str> = chrom_sizes
        .iter()
        .filter_map(|(chrom, _)| {
            bin_index.chrom_bins.get(chrom).map(|&(offset, n_bins)| {
                (offset..offset + n_bins).map(move |bin_idx| (bin_idx, chrom.as_str()))
            })
        })
        .flatten()
        .collect();

    for (&(cell_idx, bin_idx), &count) in &global_acc {
        let group_idx = cell_group[cell_idx];
        *group_bin_sum.entry((group_idx, bin_idx)).or_insert(0.0) += count as f64;
        *group_bin_n_cells.entry((group_idx, bin_idx)).or_insert(0) += 1;
        // Exclude ignored chroms from normalization denominator.
        let is_ignored = bin_to_chrom
            .get(&bin_idx)
            .map_or(false, |c| ignore_set.contains(c));
        if !is_ignored {
            group_total[group_idx] += count as f64;
        }
    }

    // Write one file per group

    let ext = match out_format {
        OutputFormat::BigWig => "bw",
        OutputFormat::BedGraph => "bedgraph",
    };

    let mut output_files: Vec<PathBuf> = Vec::with_capacity(n_groups);

    // Sanitize group name for use as a file component.
    let safe_name = |s: &str| sanitize_filename::sanitize(s);

    for group_idx in 0..n_groups {
        let group_name = &parsed.groups[group_idx];
        let output_path = PathBuf::from(format!("{}.{}.{}", output_prefix, safe_name(group_name), ext));

        let total_reads = group_total[group_idx];
        let n_cells_in_group = parsed.n_cells_per_group[group_idx] as f64;

        // Build sorted (chrom, Value) list by iterating bins in chromosome order.
        let mut values: Vec<(String, Value)> = Vec::new();

        for (chrom_name, chrom_len) in &chrom_sizes {
            let Some(&(offset, n_chrom_bins)) = bin_index.chrom_bins.get(chrom_name) else { continue };

            for local_bin in 0..n_chrom_bins {
                let bin_idx = offset + local_bin;

                let raw_sum = match group_bin_sum.get(&(group_idx, bin_idx)) {
                    Some(&s) if s > 0.0 => s,
                    _ => continue,
                };

                let normalized = match normalize_using {
                    NormalizeMethod::None => raw_sum * scale_factor,
                    NormalizeMethod::Cpm => {
                        if total_reads > 0.0 { raw_sum / total_reads * 1e6 * scale_factor } else { 0.0 }
                    }
                    NormalizeMethod::Rpkm => {
                        if total_reads > 0.0 {
                            raw_sum / (total_reads * bin_size as f64 / 1e9) * scale_factor
                        } else {
                            0.0
                        }
                    }
                    NormalizeMethod::Mean => raw_sum / n_cells_in_group * scale_factor,
                    NormalizeMethod::Frequency => {
                        let n_cells_with = group_bin_n_cells
                            .get(&(group_idx, bin_idx))
                            .copied()
                            .unwrap_or(0) as f64;
                        n_cells_with / n_cells_in_group * scale_factor
                    }
                };

                if normalized == 0.0 {
                    continue;
                }

                let start = (local_bin * step_size) as u32;
                let end = ((local_bin * step_size + bin_size).min(*chrom_len)) as u32;
                values.push((chrom_name.clone(), Value { start, end, value: normalized as f32 }));
            }
        }

        match out_format {
            OutputFormat::BigWig => {
                write_bigwig(&output_path, &chrom_sizes, values)?;
            }
            OutputFormat::BedGraph => {
                write_bedgraph(&output_path, &values)?;
            }
        }

        output_files.push(output_path);
    }

    Ok(output_files)
}

// Output helpers

fn write_bigwig(
    path: &Path,
    chrom_sizes: &[(String, usize)],
    values: Vec<(String, Value)>,
) -> Result<()> {
    let chrom_map: HashMap<String, u32> = chrom_sizes
        .iter()
        .map(|(c, l)| (c.clone(), *l as u32))
        .collect();

    let writer = BigWigWrite::create_file(path, chrom_map)
        .with_context(|| format!("failed to create bigWig file: {}", path.display()))?;

    let iter = BedParserStreamingIterator::wrap_infallible_iter(values.into_iter(), false);

    let runtime = tokio::runtime::Builder::new_current_thread()
        .build()
        .context("failed to build tokio runtime for bigWig write")?;

    writer
        .write(iter, runtime)
        .with_context(|| format!("failed to write bigWig file: {}", path.display()))?;

    Ok(())
}

fn write_bedgraph(path: &Path, values: &[(String, Value)]) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("failed to create bedGraph file: {}", path.display()))?;
    let mut w = BufWriter::new(file);

    writeln!(w, "track type=bedGraph")?;
    for (chrom, v) in values {
        writeln!(w, "{}\t{}\t{}\t{}", chrom, v.start, v.end, v.value)?;
    }
    Ok(())
}

// Python entry point

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

/// Compute pseudo-bulk coverage tracks, one bigWig (or bedGraph) per cell group.
///
/// `bam_files` and `bam_labels` must be the same length; each label must match
/// the `sample` column in `group_info`.
///
/// `group_info` — path to a TSV file with a header line and columns:
/// `sample`, `barcode`, `group`.
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
    group_info: PathBuf,
    output_prefix: String,
    bin_size: usize,
    step_size: usize,
    bc_tag: &str,
    umi_tag: Option<String>,
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

    let normalize = parse_normalize_method(normalize_using)
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

    let format = match out_format {
        "bigwig" | "bw" => OutputFormat::BigWig,
        "bedgraph" | "bg" => OutputFormat::BedGraph,
        _ => return Err(PyRuntimeError::new_err(
            "out_format must be 'bigwig' or 'bedgraph'"
        )),
    };

    let dup = dup_method
        .as_deref()
        .map(parse_dup_method)
        .transpose()
        .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;

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

    let qc = {
        use crate::counting::filters::QcFilter;
        let needs = min_fragment_length.is_some()
            || max_fragment_length.is_some()
            || min_gc.is_some()
            || max_gc.is_some()
            || min_aligned_fraction.is_some();
        needs.then(|| QcFilter {
            min_fragment_length,
            max_fragment_length,
            min_gc,
            max_gc,
            min_aligned_fraction,
        })
    };

    // Validate and parse offset.
    let parsed_offset: Option<(i32, Option<i32>)> = match offset.as_deref() {
        None => None,
        Some([s]) => {
            if *s == 0 {
                return Err(PyRuntimeError::new_err("offset value 0 is not allowed (offsets are 1-based)"));
            }
            Some((*s, None))
        }
        Some([s, e]) => {
            if *s == 0 || *e == 0 {
                return Err(PyRuntimeError::new_err("offset value 0 is not allowed (offsets are 1-based)"));
            }
            if *e > 0 && *e < *s {
                return Err(PyRuntimeError::new_err("offset end must be >= offset start"));
            }
            Some((*s, Some(*e)))
        }
        Some(_) => return Err(PyRuntimeError::new_err("offset must have 1 or 2 elements")),
    };

    if mnase && parsed_offset.is_some() {
        return Err(PyRuntimeError::new_err("--MNase and --Offset are mutually exclusive"));
    }

    let read_mode = if mnase {
        ReadMode::MNase
    } else if let Some((s, e)) = parsed_offset {
        ReadMode::Offset(s, e)
    } else {
        ReadMode::Normal
    };

    let rna_strand = match filter_rna_strand.as_deref() {
        None => None,
        Some("forward") => Some(RnaStrand::Forward),
        Some("reverse") => Some(RnaStrand::Reverse),
        Some(other) => return Err(PyRuntimeError::new_err(
            format!("filter_rna_strand must be 'forward' or 'reverse', got {:?}", other)
        )),
    };

    let path_label: Vec<(PathBuf, String)> = bam_files
        .into_iter()
        .zip(bam_labels.into_iter())
        .collect();
    let bam_path_refs: Vec<(&Path, &str)> = path_label
        .iter()
        .map(|(p, l)| (p.as_path(), l.as_str()))
        .collect();

    run_bulk_coverage(
        &bam_path_refs,
        group_info.as_path(),
        &output_prefix,
        bin_size,
        step_size,
        bc_tag,
        umi_tag.as_deref(),
        min_mapq,
        sam_flag_include,
        sam_flag_exclude,
        &chr_to_skip,
        &ignore_for_normalization,
        blacklist_path.as_deref(),
        extend_reads,
        center_reads,
        dup,
        genome_2bit.as_deref(),
        motif_filter.as_deref(),
        qc.as_ref(),
        normalize,
        scale_factor,
        format,
        read_mode,
        rna_strand,
        num_threads,
        chunk_size,
    )
    .map_err(|e| PyRuntimeError::new_err(e.to_string()))
}
