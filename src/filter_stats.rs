use std::ops::AddAssign;
use std::path::{Path, PathBuf};

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::annotation::parse_annotation::parse_blacklist_bed;
use crate::annotation::region_index::GenomeIndex;
use crate::bam::bam_io::{BamWorker, ensure_barcode_tags_present, read_bam_header};
use crate::bam::filters::{DupMethod, DuplicateFilter, rna_strand_filter};
use crate::bam::sc_record::{ScRecord, ScRecordOptions, parse_tag};

#[derive(Default, Clone)]
struct BarcodeStat {
    total: u64,
    filtered: u64,
    blacklisted: u64,
    low_mapq: u64,
    missing_flags: u64,
    excluded_flags: u64,
    internal_dupes: u64,
    external_dupes: u64,
    singletons: u64,
    wrong_strand: u64,
    wrong_motif: u64,
    wrong_gc: u64,
    low_aligned_fraction: u64,
}

impl BarcodeStat {
    fn to_vec(&self) -> Vec<u64> {
        vec![
            self.total,
            self.filtered,
            self.blacklisted,
            self.low_mapq,
            self.missing_flags,
            self.excluded_flags,
            self.internal_dupes,
            self.external_dupes,
            self.singletons,
            self.wrong_strand,
            self.wrong_motif,
            self.wrong_gc,
            self.low_aligned_fraction,
        ]
    }
}

impl AddAssign for BarcodeStat {
    fn add_assign(&mut self, other: Self) {
        self.total += other.total;
        self.filtered += other.filtered;
        self.blacklisted += other.blacklisted;
        self.low_mapq += other.low_mapq;
        self.missing_flags += other.missing_flags;
        self.excluded_flags += other.excluded_flags;
        self.internal_dupes += other.internal_dupes;
        self.external_dupes += other.external_dupes;
        self.singletons += other.singletons;
        self.wrong_strand += other.wrong_strand;
        self.wrong_motif += other.wrong_motif;
        self.wrong_gc += other.wrong_gc;
        self.low_aligned_fraction += other.low_aligned_fraction;
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

fn is_blacklisted(bl: &GenomeIndex, chrom: &str, start: usize, end: usize) -> bool {
    bl.get(chrom)
        .or_else(|| chrom.strip_prefix("chr").and_then(|c| bl.get(c)))
        .or_else(|| bl.get(&format!("chr{chrom}")))
        .map(|idx| idx.find(start, end).next().is_some())
        .unwrap_or(false)
}

pub fn run_filter_stats(
    bam_path: &Path,
    barcodes: &[String],
    bc_tag: &str,
    umi_tag: Option<&str>,
    bin_size: usize,
    distance_between_bins: usize,
    min_mapq: Option<u8>,
    sam_flag_include: Option<u16>,
    sam_flag_exclude: Option<u16>,
    chr_to_skip: &[String],
    blacklist_path: Option<&Path>,
    dup_method: Option<DupMethod>,
    genome_path: Option<&Path>,
    motifs: Option<&[(String, String)]>,
    min_gc: Option<f32>,
    max_gc: Option<f32>,
    min_aligned_fraction: Option<f32>,
    filter_rna_strand: Option<&str>,
    num_threads: usize,
    chunk_size: usize,
) -> Result<(Vec<String>, Vec<Vec<u64>>)> {
    anyhow::ensure!(bin_size > 0, "bin_size must be greater than zero");
    anyhow::ensure!(chunk_size > 0, "chunk_size must be greater than zero");

    let bc_tag_parsed = parse_tag(bc_tag)?;
    let umi_tag_parsed = umi_tag.map(parse_tag).transpose()?;
    ensure_barcode_tags_present(&[bam_path], bc_tag_parsed, umi_tag_parsed)?;

    let record_opts = ScRecordOptions {
        compute_gc: min_gc.is_some() || max_gc.is_some(),
        compute_aligned_fraction: min_aligned_fraction.is_some(),
        store_sequence: motifs.is_some(),
    };

    let blacklist: Option<GenomeIndex> = blacklist_path.map(parse_blacklist_bed).transpose()?;

    // Keyed by raw bytes so the per-read barcode lookup never allocates.
    let barcode_idx: AHashMap<&[u8], usize> = barcodes
        .iter()
        .enumerate()
        .map(|(i, bc)| (bc.as_bytes(), i))
        .collect();

    let header = read_bam_header(bam_path)?;

    let skip_set: AHashSet<&[u8]> = chr_to_skip.iter().map(|s| s.as_bytes()).collect();
    let chrom_sizes: Vec<(String, usize)> = header
        .reference_sequences()
        .iter()
        .filter(|(name, _)| {
            let bytes: &[u8] = name.as_ref();
            !skip_set.contains(bytes)
        })
        .map(|(name, seq)| (name.to_string(), seq.length().get()))
        .collect();

    // stride = distance between consecutive sampling-bin starts on the chromosome grid
    let stride = bin_size.saturating_add(distance_between_bins).max(1);
    let n_barcodes = barcodes.len();

    // Build chunk work list sorted by descending size.
    // Each chunk will iterate the sampling bins whose start falls within it.
    let mut chunks: Vec<(String, usize, usize)> = chrom_sizes
        .iter()
        .flat_map(|(chrom, chrom_len)| {
            (0..*chrom_len)
                .step_by(chunk_size)
                .map(move |start| (chrom.clone(), start, (start + chunk_size).min(*chrom_len)))
        })
        .collect();
    chunks.sort_unstable_by_key(|b| std::cmp::Reverse(b.2 - b.1));

    let n_threads = if num_threads == 0 {
        rayon::current_num_threads()
    } else {
        num_threads
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .context("failed to build thread pool")?;

    // Motif-filter ingredients, built once per worker rather than per chunk.
    let motif_ingredients = match (genome_path, motifs) {
        (Some(g), Some(m)) => Some((g, m)),
        _ => None,
    };

    let partial_stats: Vec<Vec<BarcodeStat>> = pool.install(|| {
        chunks
            .par_iter()
            .map_init(
                BamWorker::new,
                |worker, &(ref chrom, chunk_start, chunk_end)| -> Result<Vec<BarcodeStat>> {
                    let (reader, header, motif) = worker.prepare(bam_path, motif_ingredients)?;

                    let mut dup_filter: Option<DuplicateFilter> =
                        dup_method.map(DuplicateFilter::new);

                    let mut local_stats: Vec<BarcodeStat> =
                        (0..n_barcodes).map(|_| BarcodeStat::default()).collect();

                    // Align to the global sampling grid so that bins are consistent
                    // across chunks. The grid starts at 0 on each chromosome with
                    // step = stride. First bin whose start falls in [chunk_start,
                    // chunk_end) is ceil(chunk_start / stride) * stride.
                    let first_bin = if chunk_start == 0 {
                        0
                    } else {
                        chunk_start.div_ceil(stride) * stride
                    };

                    let mut bin_start = first_bin;
                    while bin_start < chunk_end {
                        let bin_end = (bin_start + bin_size).min(chunk_end);

                        let region_str = format!("{}:{}-{}", chrom, bin_start + 1, bin_end);
                        let region: noodles::core::Region = match region_str.parse() {
                            Ok(r) => r,
                            Err(_) => {
                                bin_start += stride;
                                continue;
                            }
                        };
                        let query = match reader.query(header, &region) {
                            Ok(q) => q,
                            Err(_) => {
                                bin_start += stride;
                                continue;
                            }
                        };

                        for result in query.records() {
                            let record = result.context("failed to read BAM record")?;

                            let Some(sc_rec) = ScRecord::from_bam_record(
                                &record,
                                header,
                                &bc_tag_parsed,
                                umi_tag_parsed.as_ref(),
                                None,
                                &record_opts,
                            )?
                            else {
                                continue;
                            };

                            // Ownership: skip reads that started before this bin.
                            if sc_rec.alignment_start < bin_start {
                                continue;
                            }

                            let Some(barcode) = sc_rec.barcode else {
                                continue;
                            };
                            let Some(&cell_i) = barcode_idx.get(barcode) else {
                                continue;
                            };

                            let raw_flags = u16::from(record.flags());
                            let s = &mut local_stats[cell_i];
                            s.total += 1;
                            let mut fail = false;

                            if let Some(ref bl) = blacklist
                                && is_blacklisted(
                                    bl,
                                    chrom,
                                    sc_rec.alignment_start,
                                    sc_rec.alignment_end,
                                )
                            {
                                s.blacklisted += 1;
                                fail = true;
                            }

                            if let Some(min_q) = min_mapq
                                && record.mapping_quality().is_none_or(|q| q.get() < min_q)
                            {
                                s.low_mapq += 1;
                                fail = true;
                            }

                            if let Some(include) = sam_flag_include
                                && raw_flags & include != include
                            {
                                s.missing_flags += 1;
                                fail = true;
                            }

                            if let Some(exclude) = sam_flag_exclude
                                && raw_flags & exclude != 0
                            {
                                s.excluded_flags += 1;
                                fail = true;
                            }

                            if let Some(ref mut dup) = dup_filter
                                && !dup.passes(&sc_rec)
                            {
                                s.internal_dupes += 1;
                                fail = true;
                            }

                            if raw_flags & 0x400 != 0 {
                                s.external_dupes += 1;
                                fail = true;
                            }

                            if raw_flags & 0x1 != 0 && raw_flags & 0x8 != 0 {
                                s.singletons += 1;
                                fail = true;
                            }

                            if let Some(strand) = filter_rna_strand
                                && rna_strand_filter(raw_flags, strand)
                            {
                                s.wrong_strand += 1;
                                fail = true;
                            }

                            if (min_gc.is_some() || max_gc.is_some())
                                && let Some(gc) = sc_rec.gc_content
                            {
                                let gc_fail = min_gc.is_some_and(|lo| gc < lo)
                                    || max_gc.is_some_and(|hi| gc > hi);
                                if gc_fail {
                                    s.wrong_gc += 1;
                                    fail = true;
                                }
                            }

                            if let Some(mf) = motif.as_mut()
                                && !mf.passes(&sc_rec, chrom)?
                            {
                                s.wrong_motif += 1;
                                fail = true;
                            }

                            if let Some(min_af) = min_aligned_fraction
                                && let Some(af) = sc_rec.aligned_fraction
                                && af < min_af
                            {
                                s.low_aligned_fraction += 1;
                                fail = true;
                            }

                            if fail {
                                s.filtered += 1;
                            }
                        }

                        bin_start += stride;
                    }

                    Ok(local_stats)
                },
            )
            .collect::<Result<Vec<_>>>()
    })?;

    let mut stats: Vec<BarcodeStat> = (0..n_barcodes).map(|_| BarcodeStat::default()).collect();
    for partial in partial_stats {
        for (i, s) in partial.into_iter().enumerate() {
            stats[i] += s;
        }
    }

    let stat_vecs: Vec<Vec<u64>> = stats.iter().map(|s| s.to_vec()).collect();
    Ok((barcodes.to_vec(), stat_vecs))
}

#[pyfunction(signature = (
    bam_path,
    barcodes,
    bc_tag = "CB",
    umi_tag = None,
    bin_size = 100_000,
    distance_between_bins = 1_000_000,
    min_mapq = None,
    sam_flag_include = None,
    sam_flag_exclude = None,
    chr_to_skip = vec![],
    blacklist_path = None,
    dup_method = None,
    genome_2bit = None,
    motif_filter = None,
    min_gc = None,
    max_gc = None,
    min_aligned_fraction = None,
    filter_rna_strand = None,
    num_threads = 0,
    chunk_size = 1_000_000,
))]
pub fn filter_stats(
    bam_path: PathBuf,
    barcodes: Vec<String>,
    bc_tag: &str,
    umi_tag: Option<String>,
    bin_size: usize,
    distance_between_bins: usize,
    min_mapq: Option<u8>,
    sam_flag_include: Option<u16>,
    sam_flag_exclude: Option<u16>,
    chr_to_skip: Vec<String>,
    blacklist_path: Option<PathBuf>,
    dup_method: Option<String>,
    genome_2bit: Option<PathBuf>,
    motif_filter: Option<Vec<(String, String)>>,
    min_gc: Option<f32>,
    max_gc: Option<f32>,
    min_aligned_fraction: Option<f32>,
    filter_rna_strand: Option<String>,
    num_threads: usize,
    chunk_size: usize,
) -> PyResult<(Vec<String>, Vec<Vec<u64>>)> {
    let dup = dup_method
        .as_deref()
        .map(parse_dup_method)
        .transpose()
        .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?;

    run_filter_stats(
        bam_path.as_path(),
        &barcodes,
        bc_tag,
        umi_tag.as_deref(),
        bin_size,
        distance_between_bins,
        min_mapq,
        sam_flag_include,
        sam_flag_exclude,
        &chr_to_skip,
        blacklist_path.as_deref(),
        dup,
        genome_2bit.as_deref(),
        motif_filter.as_deref(),
        min_gc,
        max_gc,
        min_aligned_fraction,
        filter_rna_strand.as_deref(),
        num_threads,
        chunk_size,
    )
    .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))
}
