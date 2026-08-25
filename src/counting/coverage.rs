//! Pseudo-bulk coverage tracks from grouped single cells.
//!
//! Cells are pooled into groups by a group-info TSV, their reads counted into
//! genomic bins, and each group written out as its own bigWig or bedGraph.
//!
//! Unlike the count matrices, a read's signal interval can come from any
//! [`ReadMode`]: the extended/centered interval, the centre of an MNase
//! fragment, or a chosen offset within the read. Counts can then be normalized
//! ([`NormalizeMethod`]) against a denominator summed over the bins, optionally
//! ignoring some chromosomes.

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

use super::params::{CountingParams, parse_region};
use crate::annotation::parse_annotation::parse_blacklist_bed;
use crate::annotation::region_index::{
    bins_touched, build_bigwig_index, build_bigwig_index_in_window,
};
use crate::bam::bam_io::{
    BamWorker, ensure_barcode_tags_present, ensure_genome_matches_bams, read_bam_header,
    read_group_ids, warn_unknown_group,
};
use crate::bam::filters::{
    DupMethod, DuplicateFilter, QcFilter, RawRecordFilter, derive_record_opts,
};
use crate::bam::filters::{blacklist_chrom_index, read_is_blacklisted};
use crate::bam::fragment_length::{ensure_paired_end, resolve_extend_reads};
use crate::bam::sc_record::{AdjustRead, EffectiveIntervals, ScRecord, parse_tag};

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
    /// Standard behavior: apply extend_reads / center_reads from `AdjustRead`.
    Normal,
    /// MNase / CUT&RUN: use only the 2–3 central bp of the paired-end fragment.
    /// Reads that are not proper-pair forward reads return `None` and are skipped.
    MNase,
    /// Report only a specific offset position (or range) within the read.
    /// Values are 1-based; negative values count from the read end.
    Offset(i32, Option<i32>),
}

// Read-mode helpers

/// Returns the 2–3 central bp of a paired-end fragment (MNase / CUT&RUN mode).
/// Only proper-pair, non-reverse reads are considered; returns `None` otherwise.
fn apply_mnase(rec: &ScRecord) -> Option<(usize, usize)> {
    let tlen = rec.template_length.unsigned_abs() as usize;
    if !rec.is_proper_pair || rec.is_reverse || tlen <= 1 {
        return None;
    }

    let frag_start = rec.alignment_start + tlen / 2 - 1;
    let frag_end = frag_start + if tlen.is_multiple_of(2) { 2 } else { 3 };
    Some((frag_start, frag_end))
}

/// Returns the genomic intervals for an offset (or offset range) inside the
/// read. Offsets are 1-based; negative values index from the read's 3' end.
/// A single value gives one base, two give the range between them.
///
/// The offset counts **aligned** positions, walking the read's gapless blocks.
/// Counting along the alignment span instead would step through an intron or a
/// deletion and place the signal where the read has none. Counting along the
/// stored sequence would step through soft-clipped bases and run off the end.
/// Because a window can cross a gap, the result may be several intervals,
/// always in ascending genomic order.
fn apply_offset(rec: &ScRecord, start: i32, end: Option<i32>) -> Vec<(usize, usize)> {
    // Ungapped reads carry no block list, their one block being the alignment.
    let span = [(rec.alignment_start, rec.alignment_end)];
    let blocks: &[(usize, usize)] = rec.blocks.as_deref().unwrap_or(&span);

    let len: usize = blocks.iter().map(|(s, e)| e - s).sum();
    if len == 0 {
        return Vec::new();
    }
    let len_i = len as i32;

    // 1-based, possibly negative offset to a 0-based index into the aligned
    // positions, read 5' to 3'.
    let to_idx = |o: i32| -> i32 { if o > 0 { o - 1 } else { len_i + o } };
    let idx_s = to_idx(start);
    let idx_e = match end {
        None => idx_s + 1,
        Some(e) => to_idx(e) + 1,
    };
    if idx_s < 0 || idx_e <= idx_s || idx_s >= len_i {
        return Vec::new();
    }
    let (idx_s, idx_e) = (idx_s as usize, (idx_e.min(len_i)) as usize);

    // Walk the blocks in read order, taking each one's share of the window.
    // For a reverse read the 5' end is the genomic end, so the walk runs
    // backwards and each block is indexed from its own end.
    let reverse = rec.is_reverse;
    let mut out: Vec<(usize, usize)> = Vec::new();
    let mut consumed = 0usize;

    for i in 0..blocks.len() {
        let (block_start, block_end) = if reverse {
            blocks[blocks.len() - 1 - i]
        } else {
            blocks[i]
        };
        let block_len = block_end - block_start;
        let (from, to) = (consumed, consumed + block_len);
        consumed = to;

        let lo = idx_s.max(from);
        let hi = idx_e.min(to);
        if lo >= hi {
            continue;
        }
        out.push(if reverse {
            (block_end - (hi - from), block_end - (lo - from))
        } else {
            (block_start + (lo - from), block_start + (hi - from))
        });
    }

    if reverse {
        out.reverse(); // back into ascending genomic order
    }
    out
}

fn get_effective_intervals<'a>(
    rec: &'a ScRecord,
    adjust: &AdjustRead,
    mode: ReadMode,
) -> EffectiveIntervals<'a> {
    match mode {
        ReadMode::Normal => rec.effective_intervals(adjust),
        ReadMode::MNase => EffectiveIntervals::One(apply_mnase(rec)),
        ReadMode::Offset(start, end) => {
            EffectiveIntervals::Selected(apply_offset(rec, start, end).into_iter())
        }
    }
}

// Group-info parsing

#[derive(Debug)]
struct ParsedGroups {
    /// Ordered unique group names.
    groups: Vec<String>,
    /// Number of cells per group (indexed by group_idx).
    n_cells_per_group: Vec<usize>,
    /// All cells: (bam_idx, barcode, group_idx).
    cells: Vec<(usize, String, usize)>,
}

/// The two shapes a group-info file comes in, told apart by its column count.
///
/// Both name the same three things; the wide one is what `scClusterCells`
/// writes, where the cell is already one `sample::barcode` string and the two
/// UMAP coordinates sit between it and the group.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GroupInfoLayout {
    /// `sample`, `barcode`, `group`.
    Plain,
    /// `sample::barcode`, `UMAP1`, `UMAP2`, `group`.
    Clustered,
}

impl GroupInfoLayout {
    /// Read the layout off the header, which fixes the file's column count.
    fn from_header(header: &str, path: &Path) -> Result<Self> {
        match header.trim_end_matches(['\r', '\n']).split('\t').count() {
            3 => Ok(Self::Plain),
            4 => Ok(Self::Clustered),
            n => anyhow::bail!(
                "{} has {} columns; a group info file has 3 (sample, barcode, group) \
                 or 4 (sample::barcode, UMAP1, UMAP2, group)",
                path.display(),
                n
            ),
        }
    }

    fn n_columns(self) -> usize {
        match self {
            Self::Plain => 3,
            Self::Clustered => 4,
        }
    }

    /// Pull `(sample, barcode, group)` out of one data row.
    fn read_row<'a>(self, fields: &[&'a str]) -> Result<(&'a str, &'a str, &'a str)> {
        match self {
            Self::Plain => Ok((fields[0], fields[1], fields[2])),
            Self::Clustered => {
                let (sample, barcode) = fields[0].split_once("::").with_context(|| {
                    format!(
                        "cell id {:?} is not `sample::barcode`; a 4-column group info file \
                         names each cell that way",
                        fields[0]
                    )
                })?;
                Ok((sample, barcode, fields[3]))
            }
        }
    }
}

fn parse_group_info(path: &Path, bam_labels: &[&str]) -> Result<ParsedGroups> {
    let file = File::open(path)
        .with_context(|| format!("failed to open group info file: {}", path.display()))?;
    let reader = BufReader::new(file);

    // label -> bam index
    let label_to_bam: AHashMap<&str, usize> = bam_labels
        .iter()
        .enumerate()
        .map(|(i, l)| (*l, i))
        .collect();

    let mut group_index: AHashMap<String, usize> = AHashMap::new();
    let mut groups: Vec<String> = Vec::new();
    let mut cells: Vec<(usize, String, usize)> = Vec::new();

    let mut lines = reader.lines();
    // The header names the columns rather than holding data, and its width is
    // what says which of the two layouts this file uses.
    let header = lines
        .next()
        .transpose()
        .with_context(|| format!("failed to read {}", path.display()))?
        .with_context(|| format!("{} is empty", path.display()))?;
    let layout = GroupInfoLayout::from_header(&header, path)?;

    for (line_no, line) in lines.enumerate() {
        let line_no = line_no + 2;
        let line = line.with_context(|| format!("failed to read group info line {line_no}"))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').map(str::trim).collect();
        anyhow::ensure!(
            fields.len() == layout.n_columns(),
            "group info line {} has {} columns, but the header has {}",
            line_no,
            fields.len(),
            layout.n_columns()
        );
        let (sample, barcode, group) = layout
            .read_row(&fields)
            .with_context(|| format!("group info line {line_no}"))?;

        let Some(&bam_idx) = label_to_bam.get(sample) else {
            continue; // sample not in BAM list. Skip
        };

        let group_idx = *group_index.entry(group.to_string()).or_insert_with(|| {
            let idx = groups.len();
            groups.push(group.to_string());
            idx
        });

        cells.push((bam_idx, barcode.to_string(), group_idx));
    }

    let n_groups = groups.len();
    let mut n_cells_per_group = vec![0usize; n_groups];
    for &(_, _, g) in &cells {
        n_cells_per_group[g] += 1;
    }

    Ok(ParsedGroups {
        groups,
        n_cells_per_group,
        cells,
    })
}

// Core function

/// Compute pseudo-bulk coverage from one or more BAM files grouped by cell
/// identity, then write one bigWig or bedGraph per group.
///
/// `bam_paths`: list of `(path, sample_label)` pairs; labels must match the
/// `sample` column in `group_info_path`.
///
/// `group_info_path`: TSV with a header line followed by rows of
/// `sample\tbarcode\tgroup`.
///
/// Normalization denominators (CPM, RPKM) are computed from the total reads
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
    group_tag: Option<&str>,
    region: Option<&str>,
    chr_to_skip: &[String],
    ignore_for_normalization: &[String],
    blacklist_path: Option<&Path>,
    extend_reads: Option<usize>,
    center_reads: bool,
    dup_method: Option<DupMethod>,
    genome_path: Option<&Path>,
    motifs: Option<&[(String, String)]>,
    record_filter: Option<&RawRecordFilter>,
    qc_filter: Option<&QcFilter>,
    normalize_using: NormalizeMethod,
    scale_factor: f64,
    out_format: OutputFormat,
    read_mode: ReadMode,
    num_threads: usize,
    chunk_size: usize,
) -> Result<Vec<PathBuf>> {
    anyhow::ensure!(bin_size > 0, "bin_size must be greater than zero");
    anyhow::ensure!(step_size > 0, "step_size must be greater than zero");
    anyhow::ensure!(chunk_size > 0, "chunk_size must be greater than zero");
    anyhow::ensure!(!bam_paths.is_empty(), "at least one BAM file is required");

    // With --groupTag the samples come from the reads' group tag rather than
    // from separate files, so the group-info `sample` column names @RG IDs and
    // exactly one BAM is accepted.
    let group_ids: Option<Vec<Vec<u8>>> = match group_tag {
        Some(_) => {
            anyhow::ensure!(
                bam_paths.len() == 1,
                "--groupTag expects a single merged BAM, but {} were given",
                bam_paths.len()
            );
            let (path, _) = bam_paths[0];
            Some(read_group_ids(&read_bam_header(path)?, path)?)
        }
        None => None,
    };
    let group_names: Vec<String> = group_ids
        .iter()
        .flatten()
        .map(|id| String::from_utf8_lossy(id).into_owned())
        .collect();

    // The sample axis the group-info file is matched against.
    let sample_labels: Vec<&str> = match &group_ids {
        Some(_) => group_names.iter().map(String::as_str).collect(),
        None => bam_paths.iter().map(|(_, l)| *l).collect(),
    };
    let parsed = parse_group_info(group_info_path, &sample_labels)?;

    let group_index: AHashMap<&[u8], usize> = group_ids
        .iter()
        .flatten()
        .enumerate()
        .map(|(i, id)| (id.as_slice(), i))
        .collect();

    let n_groups = parsed.groups.len();
    let n_cells = parsed.cells.len();

    if n_cells == 0 {
        anyhow::bail!("no cells matched between group_info and BAM labels");
    }

    // cell_global_idx -> group_idx
    let cell_group: Vec<usize> = parsed.cells.iter().map(|(_, _, g)| *g).collect();

    // (bam_idx, barcode bytes) -> cell_global_idx: keys borrow from parsed.cells.
    // Byte-keyed so the per-read lookup never allocates.
    let cell_index: AHashMap<(usize, &[u8]), usize> = parsed
        .cells
        .iter()
        .enumerate()
        .map(|(i, (bam_idx, bc, _))| ((*bam_idx, bc.as_bytes()), i))
        .collect();

    // --mnase reads the fragment from the two mates, so a single-end library
    // has nothing to take a centre from. Checked before any counting, so the
    // run fails instead of writing an empty track.
    if matches!(read_mode, ReadMode::MNase) {
        ensure_paired_end(bam_paths, "--mnase")?;
    }

    // Resolved before the record options, which ask it whether the gapless
    // blocks are worth computing.
    let adjust = AdjustRead {
        extend_reads: resolve_extend_reads(extend_reads, bam_paths)?,
        center_reads,
        max_paired_fragment_length: qc_filter.and_then(|f| f.max_fragment_length),
    };

    let has_motif = genome_path.is_some() && motifs.is_some();
    let record_opts = derive_record_opts(qc_filter, has_motif, dup_method.is_some(), &adjust);
    let bc_tag_parsed = parse_tag(bc_tag)?;
    let umi_tag_parsed = umi_tag.map(parse_tag).transpose()?;
    let group_tag_parsed = group_tag.map(parse_tag).transpose()?;
    let all_bams: Vec<&Path> = bam_paths.iter().map(|(p, _)| *p).collect();
    ensure_barcode_tags_present(&all_bams, bc_tag_parsed, umi_tag_parsed)?;

    // A genome from the wrong assembly makes every motif lookup wrong, so this
    // is checked once here rather than a silent loss of reads.
    if has_motif && let Some(genome) = genome_path {
        ensure_genome_matches_bams(genome, &all_bams)?;
    }

    let params = CountingParams {
        chr_to_skip: chr_to_skip.to_vec(),
        region: region.map(String::from),
        blacklist_path: blacklist_path.map(|p| p.to_path_buf()),
        feature_type: None,
        exon_type: None,
        name_attr: None,
        metagene: false,
    };

    // Optional region restriction (chrom[:start-end]); reads outside it are
    // not counted, so the output covers only the requested region.
    let region_filter = params.region.as_deref().map(parse_region).transpose()?;

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
                if skip_set.contains(bytes) {
                    return false;
                }
                // A region confines the output to its own chromosome.
                match &region_filter {
                    Some((region_chrom, _, _)) => bytes == region_chrom.as_bytes(),
                    None => true,
                }
            })
            .map(|(name, seq)| (name.to_string(), seq.length().get()))
            .collect()
    };

    // What to tile on each chromosome: the whole thing, or just the region.
    // `parse_region` leaves an open end as usize::MAX, so clamp to the length.
    let windows: Vec<(String, usize, usize)> = chrom_sizes
        .iter()
        .map(|(chrom, chrom_len)| match &region_filter {
            Some((_, region_start, region_end)) => (
                chrom.clone(),
                (*region_start).min(*chrom_len),
                (*region_end).min(*chrom_len),
            ),
            None => (chrom.clone(), 0, *chrom_len),
        })
        .collect();

    let bin_index = match &region_filter {
        // No region: tile whole chromosomes, exactly as before.
        None => build_bigwig_index(&chrom_sizes, bin_size, step_size),
        // A region: tile only that window, so the track covers the region.
        Some(_) => build_bigwig_index_in_window(&windows, bin_size, step_size),
    };

    let blacklist = params
        .blacklist_path
        .as_deref()
        .map(parse_blacklist_bed)
        .transpose()?;

    // Build chunk work list sorted by descending size. When a region is
    // requested, only emit chunks on its chromosome that overlap it, so a
    // multi-chromosome BAM doesn't enumerate (and later skip) every chunk.
    let mut work: Vec<(usize, &Path, String, usize, usize)> = bam_paths
        .iter()
        .enumerate()
        .flat_map(|(bam_idx, &(bam_path, _))| {
            windows.iter().flat_map(move |(chrom, win_start, win_end)| {
                (*win_start..*win_end)
                    .step_by(chunk_size)
                    .map(move |start| {
                        (
                            bam_idx,
                            bam_path,
                            chrom.clone(),
                            start,
                            (start + chunk_size).min(*win_end),
                        )
                    })
            })
        })
        .filter(
            |(_, _, chrom, chunk_start, chunk_end)| match &region_filter {
                Some((region_chrom, region_start, region_end)) => {
                    chrom == region_chrom && chunk_start < region_end && chunk_end > region_start
                }
                None => true,
            },
        )
        .collect();
    work.sort_unstable_by_key(|b| std::cmp::Reverse(b.4 - b.3));

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

    // Parallel phase: accumulate per-(cell_idx, bin_idx) counts, then combine
    // the per-chunk maps with a parallel tree reduction. Collecting them all and
    // merging on one thread instead cost ~15% of wall-clock time in the bin
    // counter (see `bin_matrix.rs`), and held every partial map in memory at once.
    let global_acc: AHashMap<(usize, usize), u32> = pool.install(|| {
        work.par_iter()
            .map_init(
                BamWorker::new,
                |worker,
                 &(bam_idx, bam_path, ref chrom, chunk_start, chunk_end)|
                 -> Result<AHashMap<(usize, usize), u32>> {
                    let (reader, header, motif) = worker.prepare(bam_path, motif_ingredients)?;

                    // Each work chunk covers a single chromosome, so its bin
                    // geometry and blacklist are resolved once here
                    // rather than per read.
                    let Some(&(chrom_offset, n_bins)) = bin_index.chrom_bins.get(chrom.as_str())
                    else {
                        return Ok(AHashMap::new());
                    };
                    // Bin 0 of this chromosome starts here: 0 unless a region
                    // restricted the tiling to a window.
                    let window_start = bin_index.window_start(chrom.as_str());

                    // Hoisted region filter: a region on a different chromosome
                    // skips the whole chunk; otherwise carry its bounds.
                    let region_bounds = match &region_filter {
                        Some((region_chrom, region_start, region_end)) => {
                            if region_chrom.as_str() != chrom.as_str() {
                                return Ok(AHashMap::new());
                            }
                            Some((*region_start, *region_end))
                        }
                        None => None,
                    };

                    let chunk_blacklist = blacklist
                        .as_ref()
                        .and_then(|bl| blacklist_chrom_index(bl, chrom.as_str()));

                    let region_str = format!("{}:{}-{}", chrom, chunk_start + 1, chunk_end);
                    let region: noodles::core::Region = region_str
                        .parse()
                        .with_context(|| format!("failed to parse region: {}", region_str))?;
                    let query = match reader.query(header, &region) {
                        Ok(q) => q,
                        Err(_) => return Ok(AHashMap::new()),
                    };

                    let mut dup_filter: Option<DuplicateFilter> =
                        dup_method.map(DuplicateFilter::new);
                    let mut local_acc: AHashMap<(usize, usize), u32> = AHashMap::new();

                    for result in query.records() {
                        let record = result.context("failed to read BAM record")?;

                        if let Some(rf) = record_filter {
                            let flags = u16::from(record.flags());
                            let mapq = record.mapping_quality().map(|q| q.get());
                            if !rf.passes(flags, mapq) {
                                continue;
                            }
                        }

                        let Some(sc_rec) = ScRecord::from_bam_record(
                            &record,
                            header,
                            &bc_tag_parsed,
                            umi_tag_parsed.as_ref(),
                            None,
                            group_tag_parsed.as_ref(),
                            &record_opts,
                        )?
                        else {
                            continue;
                        };

                        if sc_rec.alignment_start < chunk_start {
                            continue;
                        }

                        if let Some((region_start, region_end)) = region_bounds
                            && (sc_rec.alignment_end <= region_start
                                || sc_rec.alignment_start >= region_end)
                        {
                            continue;
                        }

                        // Judged on the read's own alignment span, as the filter
                        // tools judge it, and before --extendReads / --centerReads
                        // move the interval that is covered.
                        if let Some(bl) = chunk_blacklist
                            && read_is_blacklisted(bl, sc_rec.alignment_start, sc_rec.alignment_end)
                        {
                            continue;
                        }

                        let Some(barcode) = sc_rec.barcode else {
                            continue;
                        };
                        // Under --groupTag the read's own group picks the sample,
                        // so a barcode shared across source samples stays two cells.
                        let sample_idx = if group_tag_parsed.is_some() {
                            let Some(group) = sc_rec.group else {
                                continue;
                            };
                            let Some(&group_i) = group_index.get(group) else {
                                warn_unknown_group(group);
                                continue;
                            };
                            group_i
                        } else {
                            bam_idx
                        };
                        let Some(&cell_idx) = cell_index.get(&(sample_idx, barcode)) else {
                            continue;
                        };

                        if let Some(qc) = qc_filter
                            && !qc.passes(&sc_rec)
                        {
                            continue;
                        }
                        if let Some(ref mut dup) = dup_filter
                            && !dup.passes(&sc_rec)
                        {
                            continue;
                        }
                        if let Some(mf) = motif.as_mut()
                            && !mf.passes(&sc_rec, chrom)?
                        {
                            continue;
                        }

                        // A spliced read covers its blocks and not the intron
                        // between them, and a bin reached twice over is still
                        // covered once.
                        for bin in bins_touched(
                            get_effective_intervals(&sc_rec, &adjust, read_mode),
                            window_start,
                            bin_size,
                            step_size,
                            n_bins,
                        ) {
                            *local_acc.entry((cell_idx, chrom_offset + bin)).or_insert(0) +=
                                sc_rec.count;
                        }
                    }

                    Ok(local_acc)
                },
            )
            .reduce(
                || Ok(AHashMap::new()),
                |a, b| {
                    let (a, b) = (a?, b?);
                    // Drain the smaller map into the larger. Merging costs one
                    // hash lookup per entry moved, so moving the shorter side
                    // does strictly less work.
                    let (mut keep, drain) = if a.len() >= b.len() { (a, b) } else { (b, a) };
                    for (key, val) in drain {
                        *keep.entry(key).or_insert(0) += val;
                    }
                    Ok(keep)
                },
            )
    })?;

    // Aggregate per (group, bin)

    let n_bins = bin_index.n_bins;

    // Sum of counts per (group, bin), as a dense row per group indexed
    // `group_idx * n_bins + bin_idx`.
    //
    // A pseudo-bulk track pools many cells, so most bins carry signal and a hash
    // map buys nothing for the sparsity. It costs, though: an entry keyed by
    // `(usize, usize)` spends 16 bytes on the key and a control byte on top of
    // the 8-byte value, against 8 bytes flat here. Dense is both the faster and
    // the smaller representation at this density.
    let mut group_bin_sum: Vec<f64> = vec![0.0; n_groups * n_bins];

    // Cells with >=1 read in a bin, per (group, bin). Only `Frequency` reads
    // this, so it stays unallocated for every other normalization.
    let mut group_bin_n_cells: Option<Vec<u32>> =
        matches!(normalize_using, NormalizeMethod::Frequency)
            .then(|| vec![0u32; n_groups * n_bins]);

    // total reads per group (excluding ignored chroms). Needed for CPM/RPKM denominators
    let mut group_total: Vec<f64> = vec![0.0; n_groups];

    // Bin ranges of the chromosomes left out of the normalization denominator.
    // Whether a bin is ignored varies only by chromosome, so this is a short
    // list of half-open ranges to scan.
    // A chromosome that is not in the index (it was skipped entirely) simply
    // contributes no range.
    let ignored_bin_ranges: Vec<(usize, usize)> = ignore_for_normalization
        .iter()
        .filter_map(|chrom| bin_index.chrom_bins.get(chrom.as_str()))
        .map(|&(offset, n_bins)| (offset, offset + n_bins))
        .collect();

    for (&(cell_idx, bin_idx), &count) in &global_acc {
        let group_idx = cell_group[cell_idx];
        let slot = group_idx * n_bins + bin_idx;
        group_bin_sum[slot] += count as f64;
        if let Some(n_cells) = group_bin_n_cells.as_mut() {
            n_cells[slot] += 1;
        }
        // Exclude ignored chroms from normalization denominator.
        let is_ignored = ignored_bin_ranges
            .iter()
            .any(|&(start, end)| bin_idx >= start && bin_idx < end);
        if !is_ignored {
            group_total[group_idx] += count as f64;
        }
    }

    // Write one file per group

    let ext = match out_format {
        OutputFormat::BigWig => "bw",
        OutputFormat::BedGraph => "bedgraph",
    };

    let output_files: Vec<PathBuf> = pool.install(|| {
        (0..n_groups)
            .into_par_iter()
            .map(|group_idx| -> Result<PathBuf> {
                // `{prefix}_{group}.{ext}`, the name the reference
                // implementation writes, so a pipeline that globs for these
                // files finds them either way. The group name is still made
                // safe to use as one path component, which the reference
                // implementation does not do.
                let group_name = &parsed.groups[group_idx];
                let output_path = PathBuf::from(format!(
                    "{}_{}.{}",
                    output_prefix,
                    sanitize_group_name(group_name),
                    ext
                ));

                let total_reads = group_total[group_idx];
                let n_cells_in_group = parsed.n_cells_per_group[group_idx] as f64;
                // Start of this group's dense row.
                let group_base = group_idx * n_bins;

                // Build sorted (chrom, Value) list by iterating bins in
                // chromosome order.
                let mut values: Vec<(String, Value)> = Vec::new();

                for (chrom_name, win_start, win_end) in &windows {
                    let Some(&(offset, n_chrom_bins)) = bin_index.chrom_bins.get(chrom_name) else {
                        continue;
                    };

                    for local_bin in 0..n_chrom_bins {
                        let bin_idx = offset + local_bin;

                        let raw_sum = group_bin_sum[group_base + bin_idx];
                        if raw_sum <= 0.0 {
                            continue;
                        }

                        let normalized = match normalize_using {
                            NormalizeMethod::None => raw_sum * scale_factor,
                            NormalizeMethod::Cpm => {
                                if total_reads > 0.0 {
                                    raw_sum / total_reads * 1e6 * scale_factor
                                } else {
                                    0.0
                                }
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
                                // Allocated exactly when this arm is reachable.
                                let n_cells_with = group_bin_n_cells
                                    .as_ref()
                                    .expect("Frequency normalization tracks per-bin cell counts")
                                    [group_base + bin_idx]
                                    as f64;
                                n_cells_with / n_cells_in_group * scale_factor
                            }
                        };

                        if normalized == 0.0 {
                            continue;
                        }

                        let start = (win_start + local_bin * step_size) as u32;
                        let end =
                            ((win_start + local_bin * step_size + bin_size).min(*win_end)) as u32;
                        values.push((
                            chrom_name.clone(),
                            Value {
                                start,
                                end,
                                value: normalized as f32,
                            },
                        ));
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

                Ok(output_path)
            })
            .collect::<Result<Vec<_>>>()
    })?;

    Ok(output_files)
}

// Output helpers

/// Make a cell-group label safe to use as a single path component.
///
/// Removes the characters that are illegal in filenames on common platforms
/// (path separators and Windows-reserved punctuation) plus control characters,
/// then trims trailing spaces and dots (also illegal on Windows). This covers
/// the subset of `sanitize_filename::sanitize` behavior we rely on for group
/// labels; Windows-reserved device names (e.g. `CON`) are not special-cased, as
/// they don't occur in practice for cell-group labels.
fn sanitize_group_name(s: &str) -> String {
    let cleaned: String = s
        .chars()
        .filter(|c| {
            !matches!(c, '/' | '\\' | ':' | '*' | '?' | '"' | '<' | '>' | '|') && !c.is_control()
        })
        .collect();
    cleaned.trim_end_matches([' ', '.']).to_string()
}

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
/// `group_info` is the path to a TSV file with a header line and columns:
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
    group_info: PathBuf,
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

    let normalize = parse_normalize_method(normalize_using)
        .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?;

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
        .map(parse_dup_method)
        .transpose()
        .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))?;

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

    if let Some(strand) = filter_rna_strand.as_deref()
        && strand != "forward"
        && strand != "reverse"
    {
        return Err(PyRuntimeError::new_err(format!(
            "filter_rna_strand must be 'forward' or 'reverse', got {:?}",
            strand
        )));
    }

    let qc = {
        let needs = min_fragment_length.is_some()
            || max_fragment_length.is_some()
            || min_gc.is_some()
            || max_gc.is_some()
            || min_aligned_fraction.is_some();
        needs.then_some(QcFilter {
            min_fragment_length,
            max_fragment_length,
            min_gc,
            max_gc,
            min_aligned_fraction,
        })
    };

    let record_filter = {
        let needs = min_mapq.is_some()
            || sam_flag_include.is_some()
            || sam_flag_exclude.is_some()
            || filter_rna_strand.is_some();
        needs.then_some(RawRecordFilter {
            min_mapq,
            sam_flag_include,
            sam_flag_exclude,
            filter_rna_strand,
        })
    };

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
        group_info.as_path(),
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
    .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::sc_record::test_record;
    use tempfile::TempDir;

    // Read-mode helpers: MNase

    #[test]
    fn mnase_takes_the_two_central_bases_of_an_even_fragment() {
        let mut rec = test_record(1000, 1050);
        rec.is_proper_pair = true;
        rec.template_length = 200;
        // Centre is 1000 + 200/2 - 1; an even fragment yields 2 bp.
        assert_eq!(apply_mnase(&rec), Some((1099, 1101)));
    }

    #[test]
    fn mnase_takes_three_bases_when_the_fragment_length_is_odd() {
        let mut rec = test_record(1000, 1050);
        rec.is_proper_pair = true;
        rec.template_length = 201;
        assert_eq!(apply_mnase(&rec), Some((1099, 1102)));
    }

    #[test]
    fn mnase_uses_the_absolute_template_length() {
        let mut fwd = test_record(1000, 1050);
        fwd.is_proper_pair = true;
        fwd.template_length = 200;

        let mut negative = test_record(1000, 1050);
        negative.is_proper_pair = true;
        negative.template_length = -200;

        assert_eq!(apply_mnase(&fwd), apply_mnase(&negative));
    }

    #[test]
    fn mnase_skips_anything_that_is_not_a_forward_proper_pair() {
        let mut unpaired = test_record(1000, 1050);
        unpaired.template_length = 200;
        assert_eq!(apply_mnase(&unpaired), None);

        let mut reverse = test_record(1000, 1050);
        reverse.is_proper_pair = true;
        reverse.is_reverse = true;
        reverse.template_length = 200;
        assert_eq!(apply_mnase(&reverse), None);

        // A fragment of one base or less has no centre to take.
        let mut degenerate = test_record(1000, 1050);
        degenerate.is_proper_pair = true;
        degenerate.template_length = 1;
        assert_eq!(apply_mnase(&degenerate), None);
    }

    // Read-mode helpers: offset

    #[test]
    fn a_positive_offset_is_one_based_from_the_read_start() {
        let rec = test_record(1000, 1050);
        assert_eq!(apply_offset(&rec, 1, None), [(1000, 1001)]);
        assert_eq!(apply_offset(&rec, 10, None), [(1009, 1010)]);
    }

    #[test]
    fn a_negative_offset_counts_back_from_the_read_end() {
        let rec = test_record(1000, 1050);
        assert_eq!(apply_offset(&rec, -1, None), [(1049, 1050)]);
        assert_eq!(apply_offset(&rec, -50, None), [(1000, 1001)]);
    }

    #[test]
    fn two_offsets_give_the_range_between_them() {
        let rec = test_record(1000, 1050);
        assert_eq!(apply_offset(&rec, 1, Some(5)), [(1000, 1005)]);
        assert_eq!(apply_offset(&rec, 3, Some(-1)), [(1002, 1050)]);
    }

    #[test]
    fn on_a_reverse_read_the_offset_is_measured_from_the_alignment_end() {
        let mut rec = test_record(1000, 1050);
        rec.is_reverse = true;
        // The first base of the read is the last base in genomic order.
        assert_eq!(apply_offset(&rec, 1, None), [(1049, 1050)]);
        assert_eq!(apply_offset(&rec, 1, Some(5)), [(1045, 1050)]);
    }

    #[test]
    fn offsets_outside_the_read_are_rejected() {
        let rec = test_record(1000, 1050);
        assert!(apply_offset(&rec, 51, None).is_empty(), "past the read end");
        assert!(
            apply_offset(&rec, 5, Some(2)).is_empty(),
            "end before start"
        );

        let empty = test_record(1000, 1000);
        assert!(apply_offset(&empty, 1, None).is_empty(), "zero-length read");
    }

    #[test]
    fn a_range_that_runs_past_the_read_end_is_clamped_to_it() {
        let rec = test_record(1000, 1050);
        assert_eq!(apply_offset(&rec, 48, Some(100)), [(1047, 1050)]);
    }

    // Read-mode dispatch

    /// A spliced read: two 10 bp exons either side of a 100 bp intron.
    fn spliced_record<'a>() -> ScRecord<'a> {
        let mut rec = test_record(1000, 1120);
        rec.blocks = Some(vec![(1000, 1010), (1110, 1120)]);
        rec
    }

    #[test]
    fn an_offset_counts_aligned_bases_and_steps_over_an_intron() {
        let rec = spliced_record();

        // Bases 1-10 are the first exon; the 11th aligned base is the first of
        // the second exon, not the 11th base of the intron.
        assert_eq!(apply_offset(&rec, 10, None), [(1009, 1010)]);
        assert_eq!(apply_offset(&rec, 11, None), [(1110, 1111)]);

        // The RiboSeq P-site offset, which is what this option exists for.
        assert_eq!(apply_offset(&rec, 12, None), [(1111, 1112)]);

        // The last aligned base is the end of the second exon.
        assert_eq!(apply_offset(&rec, -1, None), [(1119, 1120)]);
    }

    #[test]
    fn an_offset_range_crossing_a_gap_comes_back_as_two_intervals() {
        let rec = spliced_record();

        // Bases 8-13 straddle the junction: three in each exon, none between.
        assert_eq!(
            apply_offset(&rec, 8, Some(13)),
            [(1007, 1010), (1110, 1113)]
        );
        // The whole read is its two exons, and nothing of the intron.
        assert_eq!(
            apply_offset(&rec, 1, Some(-1)),
            [(1000, 1010), (1110, 1120)]
        );
    }

    #[test]
    fn an_offset_past_the_aligned_length_yields_nothing() {
        // 20 aligned bases over a 120 bp footprint: asking for the 21st must
        // report nothing rather than a position inside the intron.
        let rec = spliced_record();
        assert!(apply_offset(&rec, 21, None).is_empty());
        assert!(apply_offset(&rec, -21, None).is_empty());
    }

    #[test]
    fn on_a_reverse_read_the_offset_counts_from_the_far_exon() {
        // The 5' end of a reverse read is the genomic end, so counting starts
        // at the last base of the second exon and runs back through it.
        let mut rec = spliced_record();
        rec.is_reverse = true;

        assert_eq!(apply_offset(&rec, 1, None), [(1119, 1120)]);
        assert_eq!(apply_offset(&rec, 10, None), [(1110, 1111)]);
        assert_eq!(apply_offset(&rec, 11, None), [(1009, 1010)]);

        // A range crossing the junction still comes back ascending. Aligned
        // bases 9-12 of this read are 1111 and 1110 in the far exon, then 1009
        // and 1008 in the near one.
        assert_eq!(
            apply_offset(&rec, 9, Some(12)),
            [(1008, 1010), (1110, 1112)]
        );
    }

    #[test]
    fn an_offset_ignores_soft_clipped_bases() {
        // `5S10M5S`: the stored sequence is 20 bases but only 10 align, so the
        // 12th aligned base does not exist. Counting the clips would have put
        // it two bases past the alignment.
        let mut rec = test_record(1000, 1010);
        rec.read_length = 20;

        assert!(apply_offset(&rec, 12, None).is_empty());
        assert_eq!(apply_offset(&rec, -1, None), [(1009, 1010)]);
    }

    fn intervals_of(rec: &ScRecord, adjust: &AdjustRead, mode: ReadMode) -> Vec<(usize, usize)> {
        get_effective_intervals(rec, adjust, mode).collect()
    }

    #[test]
    fn the_read_mode_selects_which_interval_rule_applies() {
        let mut rec = test_record(1000, 1050);
        rec.is_proper_pair = true;
        rec.template_length = 200;
        let adjust = AdjustRead::default();

        assert_eq!(
            intervals_of(&rec, &adjust, ReadMode::Normal),
            [(1000, 1050)]
        );
        assert_eq!(intervals_of(&rec, &adjust, ReadMode::MNase), [(1099, 1101)]);
        assert_eq!(
            intervals_of(&rec, &adjust, ReadMode::Offset(1, None)),
            [(1000, 1001)]
        );
    }

    #[test]
    fn a_read_the_mode_cannot_place_yields_no_interval() {
        let rec = test_record(1000, 1050);
        // Not a proper pair, so MNase has no fragment to centre on.
        let intervals = intervals_of(&rec, &AdjustRead::default(), ReadMode::MNase);
        assert!(intervals.is_empty(), "{intervals:?}");
    }

    #[test]
    fn only_the_normal_mode_follows_a_spliced_read() {
        // A synthetic window inside the read is one interval whatever the read's
        // own shape, so only `Normal` splits.
        let mut rec = test_record(1000, 1120);
        rec.is_proper_pair = true;
        rec.template_length = 200;
        rec.blocks = Some(vec![(1000, 1010), (1110, 1120)]);
        let adjust = AdjustRead::default();

        assert_eq!(
            intervals_of(&rec, &adjust, ReadMode::Normal),
            [(1000, 1010), (1110, 1120)]
        );
        assert_eq!(intervals_of(&rec, &adjust, ReadMode::MNase).len(), 1);
        assert_eq!(
            intervals_of(&rec, &adjust, ReadMode::Offset(1, None)).len(),
            1
        );
    }

    // Group-name sanitizing

    #[test]
    fn sanitizing_removes_path_and_wildcard_characters() {
        assert_eq!(sanitize_group_name("T cells/CD8"), "T cellsCD8");
        assert_eq!(sanitize_group_name("a\\b:c*d?e\"f<g>h|i"), "abcdefghi");
    }

    #[test]
    fn sanitizing_drops_control_characters_and_trailing_dots_and_spaces() {
        assert_eq!(sanitize_group_name("group\tone"), "groupone");
        assert_eq!(sanitize_group_name("cluster1..  "), "cluster1");
    }

    #[test]
    fn sanitizing_keeps_interior_dots_and_spaces() {
        assert_eq!(sanitize_group_name("cluster 1.2"), "cluster 1.2");
    }

    // Enum parsing

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

    #[test]
    fn every_documented_dup_method_parses() {
        assert!(matches!(
            parse_dup_method("barcode_start").unwrap(),
            DupMethod::BarcodeStart
        ));
        assert!(matches!(
            parse_dup_method("barcode_start_end").unwrap(),
            DupMethod::BarcodeStartEnd
        ));
        assert!(matches!(
            parse_dup_method("barcode_umi_start").unwrap(),
            DupMethod::BarcodeUmiStart
        ));
        assert!(matches!(
            parse_dup_method("barcode_umi_start_end").unwrap(),
            DupMethod::BarcodeUmiStartEnd
        ));
    }

    #[test]
    fn an_unknown_dup_method_names_the_valid_ones() {
        let err = parse_dup_method("start_umi").unwrap_err().to_string();
        assert!(err.contains("barcode_umi_start"), "{err}");
    }

    // Group-info parsing

    fn write_group_info(dir: &TempDir, body: &str) -> PathBuf {
        let path = dir.path().join("groups.tsv");
        let mut file = File::create(&path).unwrap();
        file.write_all(body.as_bytes()).unwrap();
        path
    }

    #[test]
    fn group_info_skips_the_header_and_numbers_groups_in_first_seen_order() {
        let dir = TempDir::new().unwrap();
        let path = write_group_info(
            &dir,
            "sample\tbarcode\tgroup\ns1\tAAA\tB\ns1\tCCC\tA\ns2\tGGG\tB\n",
        );

        let parsed = parse_group_info(&path, &["s1", "s2"]).unwrap();

        assert_eq!(parsed.groups, vec!["B", "A"]);
        assert_eq!(parsed.n_cells_per_group, vec![2, 1]);
        assert_eq!(
            parsed.cells,
            vec![
                (0, "AAA".to_string(), 0),
                (0, "CCC".to_string(), 1),
                (1, "GGG".to_string(), 0),
            ]
        );
    }

    #[test]
    fn a_four_column_group_info_names_each_cell_as_sample_and_barcode() {
        // The shape `scClusterCells` writes: the cell is one `sample::barcode`
        // string with the two UMAP coordinates between it and the group.
        let dir = TempDir::new().unwrap();
        let path = write_group_info(
            &dir,
            "Cell_ID\tUMAP1\tUMAP2\tcluster\n\
             s1::AAA\t0.1\t-2.3\tg1\n\
             s1::CCC\t1.4\t0.7\tg2\n\
             s2::GGG\t-0.8\t1.1\tg1\n",
        );

        let parsed = parse_group_info(&path, &["s1", "s2"]).unwrap();

        assert_eq!(parsed.groups, vec!["g1", "g2"]);
        assert_eq!(
            parsed.cells,
            vec![
                (0, "AAA".to_string(), 0),
                (0, "CCC".to_string(), 1),
                (1, "GGG".to_string(), 0),
            ]
        );
        assert_eq!(parsed.n_cells_per_group, vec![2, 1]);
    }

    #[test]
    fn a_group_info_with_another_column_count_is_rejected() {
        let dir = TempDir::new().unwrap();

        for header in ["sample\tgroup", "a\tb\tc\td\te"] {
            let path = write_group_info(&dir, &format!("{header}\n"));
            let err = parse_group_info(&path, &["s1"]).unwrap_err().to_string();
            assert!(err.contains("columns"), "{err}");
        }
    }

    #[test]
    fn a_four_column_cell_id_must_name_its_sample() {
        // Without the `::` there is no sample to match against a BAM, which is
        // an error rather than a row to skip: every row would be dropped and
        // the run would report an empty track.
        let dir = TempDir::new().unwrap();
        let path = write_group_info(
            &dir,
            "Cell_ID\tUMAP1\tUMAP2\tcluster\ns1_AAA\t0.1\t-2.3\tg1\n",
        );

        let err = format!("{:#}", parse_group_info(&path, &["s1"]).unwrap_err());
        assert!(err.contains("sample::barcode"), "{err}");
    }

    #[test]
    fn a_row_that_does_not_match_the_header_width_is_rejected() {
        let dir = TempDir::new().unwrap();
        let path = write_group_info(&dir, "sample\tbarcode\tgroup\ns1\tAAA\tg1\ns1\tCCC\n");

        let err = parse_group_info(&path, &["s1"]).unwrap_err().to_string();
        assert!(err.contains("line 3"), "{err}");
    }

    #[test]
    fn an_empty_group_info_is_rejected() {
        let dir = TempDir::new().unwrap();
        let path = write_group_info(&dir, "");

        let err = parse_group_info(&path, &["s1"]).unwrap_err().to_string();
        assert!(err.contains("empty"), "{err}");
    }

    #[test]
    fn group_info_ignores_blank_lines_comments_and_unknown_samples() {
        let dir = TempDir::new().unwrap();
        let path = write_group_info(
            &dir,
            "sample\tbarcode\tgroup\n\n# a comment\ns1\tAAA\tA\nabsent\tTTT\tA\n",
        );

        let parsed = parse_group_info(&path, &["s1"]).unwrap();

        assert_eq!(parsed.groups, vec!["A"]);
        assert_eq!(parsed.n_cells_per_group, vec![1]);
        assert_eq!(parsed.cells, vec![(0, "AAA".to_string(), 0)]);
    }

    #[test]
    fn group_info_trims_whitespace_around_every_field() {
        let dir = TempDir::new().unwrap();
        let path = write_group_info(&dir, "sample\tbarcode\tgroup\ns1 \t AAA \t A \n");

        let parsed = parse_group_info(&path, &["s1"]).unwrap();

        assert_eq!(parsed.groups, vec!["A"]);
        assert_eq!(parsed.cells, vec![(0, "AAA".to_string(), 0)]);
    }

    #[test]
    fn a_missing_group_info_file_is_reported_with_its_path() {
        // ParsedGroups is not Debug, so unwrap_err() is not available here.
        let err = match parse_group_info(Path::new("/nonexistent/groups.tsv"), &["s1"]) {
            Ok(_) => panic!("expected an error for a missing file"),
            Err(e) => e.to_string(),
        };
        assert!(err.contains("group info file"), "{err}");
    }

    // bedGraph writing

    #[test]
    fn bedgraph_writes_a_track_line_then_one_row_per_interval() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.bedgraph");
        let values = vec![
            (
                "chr1".to_string(),
                Value {
                    start: 0,
                    end: 10,
                    value: 1.5,
                },
            ),
            (
                "chr2".to_string(),
                Value {
                    start: 10,
                    end: 20,
                    value: 0.0,
                },
            ),
        ];

        write_bedgraph(&path, &values).unwrap();

        let text = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines,
            ["track type=bedGraph", "chr1\t0\t10\t1.5", "chr2\t10\t20\t0"]
        );
    }

    // Pseudo-bulk coverage, end to end

    /// `test_group_info.tsv` puts each of the 8 cells of `test_i1.bam` in a
    /// group of its own and pools all 7 cells of `test_i2.bam` into `i2_pool`,
    /// so a run produces 9 groups. `GCGAGCAT` occurs in both BAMs and must be
    /// treated as two different cells.
    fn testdata() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/testdata")
    }

    const N_GROUPS: usize = 9;

    #[allow(clippy::too_many_arguments)]
    fn run_coverage(
        prefix: &Path,
        labels: (&str, &str),
        format: OutputFormat,
        normalize: NormalizeMethod,
        scale_factor: f64,
        region: Option<&str>,
        mode: ReadMode,
        bin_size: usize,
        step_size: usize,
        chunk_size: usize,
    ) -> Result<Vec<PathBuf>> {
        let i1 = testdata().join("test_i1.bam");
        let i2 = testdata().join("test_i2.bam");
        run_bulk_coverage(
            &[(i1.as_path(), labels.0), (i2.as_path(), labels.1)],
            &testdata().join("test_group_info.tsv"),
            prefix.to_str().unwrap(),
            bin_size,
            step_size,
            "BC",
            None,
            None,
            region,
            &[],
            &[],
            None,
            None,
            false,
            None,
            None,
            None,
            None,
            None,
            normalize,
            scale_factor,
            format,
            mode,
            1,
            chunk_size,
        )
    }

    /// The defaults every test below starts from: bedGraph, unnormalized.
    fn run_default(prefix: &Path) -> Result<Vec<PathBuf>> {
        run_coverage(
            prefix,
            ("test_i1", "test_i2"),
            OutputFormat::BedGraph,
            NormalizeMethod::None,
            1.0,
            None,
            ReadMode::Normal,
            100_000,
            100_000,
            1_000_000,
        )
    }

    /// `run_default`, with a blacklist and its own output prefix.
    fn run_with_blacklist(prefix: &Path, blacklist: Option<&Path>) -> Vec<PathBuf> {
        let i1 = testdata().join("test_i1.bam");
        let i2 = testdata().join("test_i2.bam");
        run_bulk_coverage(
            &[(i1.as_path(), "test_i1"), (i2.as_path(), "test_i2")],
            &testdata().join("test_group_info.tsv"),
            prefix.to_str().unwrap(),
            100_000,
            100_000,
            "BC",
            None,
            None,
            None,
            &[],
            &[],
            blacklist,
            None,
            false,
            None,
            None,
            None,
            None,
            None,
            NormalizeMethod::None,
            1.0,
            OutputFormat::BedGraph,
            ReadMode::Normal,
            1,
            1_000_000,
        )
        .unwrap()
    }

    /// Total signal over every group's track.
    fn total_signal(files: &[PathBuf]) -> f64 {
        files
            .iter()
            .flat_map(|f| bedgraph_rows(f))
            .map(|(_, _, _, value)| value)
            .sum()
    }

    #[test]
    fn a_blacklist_removes_reads_from_the_coverage_tracks() {
        let dir = TempDir::new().unwrap();

        let base = total_signal(&run_with_blacklist(&dir.path().join("base"), None));
        assert!(base > 0.0, "the unfiltered run produced no signal");

        // 20 bp, shorter than any read here, so no read reaches the threshold
        // and the tracks are untouched.
        let narrow = dir.path().join("narrow.bed");
        std::fs::write(&narrow, "5\t65966820\t65966840\tregion\t0\t.\n").unwrap();
        let clipped = total_signal(&run_with_blacklist(
            &dir.path().join("clipped"),
            Some(&narrow),
        ));
        assert_eq!(clipped, base, "a clipping blacklist dropped a read");

        // Both BED regions cover the whole locus these BAMs were cut from, so
        // every read is at least half blacklisted.
        let filtered = total_signal(&run_with_blacklist(
            &dir.path().join("filtered"),
            Some(&testdata().join("Chrna9_regions.bed")),
        ));
        assert!(
            filtered < base,
            "the blacklist removed no read at all ({base})"
        );
    }

    /// bedGraph rows, without the leading track line.
    fn bedgraph_rows(path: &Path) -> Vec<(String, u32, u32, f64)> {
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<&str> = line.split('\t').collect();
                (
                    f[0].to_string(),
                    f[1].parse().unwrap(),
                    f[2].parse().unwrap(),
                    f[3].parse().unwrap(),
                )
            })
            .collect()
    }

    #[test]
    fn one_output_file_is_written_per_group() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");

        let files = run_default(&prefix).unwrap();

        assert_eq!(files.len(), N_GROUPS);
        for file in &files {
            assert!(file.exists(), "{} was not written", file.display());
        }
    }

    #[test]
    fn each_output_is_named_after_its_group() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");

        let files = run_default(&prefix).unwrap();
        let names: Vec<String> = files
            .iter()
            .map(|f| f.file_name().unwrap().to_string_lossy().into_owned())
            .collect();

        // The 8 singleton groups of test_i1 plus the one pooled group of test_i2.
        assert!(
            names.contains(&"cov_i2_pool.bedgraph".to_string()),
            "{names:?}"
        );
        for bc in ["ACGGTAAT", "ATATAACT", "GCGAGCAT", "TAGACTTG"] {
            assert!(
                names.contains(&format!("cov_i1_{bc}.bedgraph")),
                "missing i1_{bc} in {names:?}"
            );
        }
    }

    #[test]
    fn a_barcode_present_in_both_bams_is_two_separate_cells() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");

        let files = run_default(&prefix).unwrap();
        let names: Vec<String> = files
            .iter()
            .map(|f| f.file_name().unwrap().to_string_lossy().into_owned())
            .collect();

        // GCGAGCAT is in both BAMs: once as its own test_i1 group, once inside
        // the pooled test_i2 group. It must not collapse into one cell.
        assert!(names.contains(&"cov_i1_GCGAGCAT.bedgraph".to_string()));
        assert!(names.contains(&"cov_i2_pool.bedgraph".to_string()));
        assert_eq!(files.len(), N_GROUPS);
    }

    #[test]
    fn every_output_is_a_bedgraph_with_a_track_line() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");

        for file in run_default(&prefix).unwrap() {
            let text = std::fs::read_to_string(&file).unwrap();
            let first = text.lines().next().unwrap();
            assert_eq!(first, "track type=bedGraph", "in {}", file.display());

            for (chrom, start, end, _) in bedgraph_rows(&file) {
                assert_eq!(chrom, "5", "the test BAMs only hold chromosome 5");
                assert!(end > start, "empty interval {start}-{end}");
            }
        }
    }

    #[test]
    fn the_pooled_group_covers_at_least_as_much_as_one_of_its_own_cells() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");
        run_default(&prefix).unwrap();

        // i2_pool holds all 7 test_i2 cells, so its signal cannot be less than
        // the single test_i1 cell that shares a barcode with one of them.
        let pooled: f64 = bedgraph_rows(&dir.path().join("cov_i2_pool.bedgraph"))
            .iter()
            .map(|(_, _, _, v)| v)
            .sum();
        assert!(pooled > 0.0, "the pooled group counted nothing");
    }

    #[test]
    fn mnase_refuses_a_single_end_library() {
        // The centre of a fragment is defined by its two mates, so there is
        // nothing to compute from a single-end file. Better to say so than to
        // write an empty track, which is what skipping every read would give.
        let dir = TempDir::new().unwrap();
        let single_end = testdata().join("test_spliced.bam");
        let group_info =
            write_group_info(&dir, "sample\tbarcode\tgroup\ntest_spliced\tATATAACT\tg1\n");

        let err = run_bulk_coverage(
            &[(single_end.as_path(), "test_spliced")],
            &group_info,
            dir.path().join("cov").to_str().unwrap(),
            10_000,
            10_000,
            "BC",
            None,
            None,
            None,
            &[],
            &[],
            None,
            None,
            false,
            None,
            None,
            None,
            None,
            None,
            NormalizeMethod::None,
            1.0,
            OutputFormat::BedGraph,
            ReadMode::MNase,
            1,
            1_000_000,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("--mnase"), "{err}");
        assert!(err.contains("single-end"), "{err}");
    }

    #[test]
    fn bigwig_output_is_written_with_the_bw_extension() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");

        let files = run_coverage(
            &prefix,
            ("test_i1", "test_i2"),
            OutputFormat::BigWig,
            NormalizeMethod::None,
            1.0,
            None,
            ReadMode::Normal,
            100_000,
            100_000,
            1_000_000,
        )
        .unwrap();

        assert_eq!(files.len(), N_GROUPS);
        for file in &files {
            assert_eq!(file.extension().unwrap(), "bw");
            assert!(file.metadata().unwrap().len() > 0, "empty bigWig");
        }

        // `{prefix}_{group}.bw`, as the reference implementation names it.
        let names: Vec<String> = files
            .iter()
            .map(|f| f.file_name().unwrap().to_string_lossy().into_owned())
            .collect();
        assert!(names.contains(&"cov_i2_pool.bw".to_string()), "{names:?}");
        assert!(
            names
                .iter()
                .all(|n| n.starts_with("cov_") && n.ends_with(".bw")),
            "{names:?}"
        );
    }

    #[test]
    fn a_scale_factor_multiplies_the_values_but_keeps_the_intervals() {
        let dir = TempDir::new().unwrap();

        let plain = dir.path().join("plain");
        run_coverage(
            &plain,
            ("test_i1", "test_i2"),
            OutputFormat::BedGraph,
            NormalizeMethod::None,
            1.0,
            None,
            ReadMode::Normal,
            100_000,
            100_000,
            1_000_000,
        )
        .unwrap();

        let scaled = dir.path().join("scaled");
        run_coverage(
            &scaled,
            ("test_i1", "test_i2"),
            OutputFormat::BedGraph,
            NormalizeMethod::None,
            10.0,
            None,
            ReadMode::Normal,
            100_000,
            100_000,
            1_000_000,
        )
        .unwrap();

        let a = bedgraph_rows(&dir.path().join("plain_i2_pool.bedgraph"));
        let b = bedgraph_rows(&dir.path().join("scaled_i2_pool.bedgraph"));
        assert_eq!(a.len(), b.len(), "scaling changed the interval layout");

        for ((c1, s1, e1, v1), (c2, s2, e2, v2)) in a.iter().zip(b.iter()) {
            assert_eq!((c1, s1, e1), (c2, s2, e2));
            assert!((v2 - v1 * 10.0).abs() < 1e-3, "{v1} scaled to {v2}");
        }
    }

    #[test]
    fn cpm_normalization_changes_the_values_but_keeps_the_intervals() {
        let dir = TempDir::new().unwrap();

        let raw = dir.path().join("raw");
        run_default(&raw).unwrap();

        let cpm = dir.path().join("cpm");
        run_coverage(
            &cpm,
            ("test_i1", "test_i2"),
            OutputFormat::BedGraph,
            NormalizeMethod::Cpm,
            1.0,
            None,
            ReadMode::Normal,
            100_000,
            100_000,
            1_000_000,
        )
        .unwrap();

        let a = bedgraph_rows(&dir.path().join("raw_i2_pool.bedgraph"));
        let b = bedgraph_rows(&dir.path().join("cpm_i2_pool.bedgraph"));

        assert_eq!(a.len(), b.len());
        let intervals = |rows: &[(String, u32, u32, f64)]| -> Vec<(String, u32, u32)> {
            rows.iter()
                .map(|(c, s, e, _)| (c.clone(), *s, *e))
                .collect()
        };
        assert_eq!(intervals(&a), intervals(&b));
        assert_ne!(
            a.iter().map(|r| r.3).collect::<Vec<_>>(),
            b.iter().map(|r| r.3).collect::<Vec<_>>(),
            "CPM left the values unchanged"
        );
    }

    #[test]
    fn an_offset_read_mode_still_produces_every_group() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");

        let files = run_coverage(
            &prefix,
            ("test_i1", "test_i2"),
            OutputFormat::BedGraph,
            NormalizeMethod::None,
            1.0,
            None,
            ReadMode::Offset(1, None),
            100_000,
            100_000,
            1_000_000,
        )
        .unwrap();

        assert_eq!(files.len(), N_GROUPS);
    }

    #[test]
    fn restricting_to_an_empty_window_leaves_no_intervals() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");

        // The reads sit around 65.97 Mb, so this window holds none of them.
        let files = run_coverage(
            &prefix,
            ("test_i1", "test_i2"),
            OutputFormat::BedGraph,
            NormalizeMethod::None,
            1.0,
            Some("5:0-100000"),
            ReadMode::Normal,
            100_000,
            100_000,
            1_000_000,
        )
        .unwrap();

        assert_eq!(
            files.len(),
            N_GROUPS,
            "a group with no signal still gets a file"
        );
        let total: f64 = files
            .iter()
            .flat_map(|f| bedgraph_rows(f))
            .map(|(_, _, _, v)| v)
            .sum();
        assert_eq!(total, 0.0, "an empty window should carry no signal");
    }

    // Argument validation

    #[test]
    fn labels_that_match_no_group_info_row_are_rejected() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");

        let err = run_coverage(
            &prefix,
            ("other1", "other2"),
            OutputFormat::BedGraph,
            NormalizeMethod::None,
            1.0,
            None,
            ReadMode::Normal,
            100_000,
            100_000,
            1_000_000,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("no cells matched"), "{err}");
    }

    #[test]
    fn zero_sized_bins_steps_and_chunks_are_rejected() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");

        let run = |bin, step, chunk| {
            run_coverage(
                &prefix,
                ("test_i1", "test_i2"),
                OutputFormat::BedGraph,
                NormalizeMethod::None,
                1.0,
                None,
                ReadMode::Normal,
                bin,
                step,
                chunk,
            )
            .unwrap_err()
            .to_string()
        };

        assert!(run(0, 100, 1_000).contains("bin_size"));
        assert!(run(100, 0, 1_000).contains("step_size"));
        assert!(run(100, 100, 0).contains("chunk_size"));
    }

    #[test]
    fn running_with_no_bam_files_is_rejected() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");

        let err = run_bulk_coverage(
            &[],
            &testdata().join("test_group_info.tsv"),
            prefix.to_str().unwrap(),
            100_000,
            100_000,
            "BC",
            None,
            None,
            None,
            &[],
            &[],
            None,
            None,
            false,
            None,
            None,
            None,
            None,
            None,
            NormalizeMethod::None,
            1.0,
            OutputFormat::BedGraph,
            ReadMode::Normal,
            1,
            1_000_000,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("at least one BAM"), "{err}");
    }

    // Region-restricted tiling

    /// Coverage keyed by genomic interval, so runs with different bin layouts
    /// can be compared interval for interval.
    fn signal_by_interval(path: &Path) -> std::collections::BTreeMap<(u32, u32), f64> {
        bedgraph_rows(path)
            .into_iter()
            .map(|(_, start, end, value)| ((start, end), value))
            .collect()
    }

    #[test]
    fn a_region_confines_the_track_to_its_own_window() {
        let dir = TempDir::new().unwrap();
        let prefix = dir.path().join("cov");

        // The reads sit around 65.97 Mb; tile just the 100 kb around them.
        let files = run_coverage(
            &prefix,
            ("test_i1", "test_i2"),
            OutputFormat::BedGraph,
            NormalizeMethod::None,
            1.0,
            Some("5:65900000-66000000"),
            ReadMode::Normal,
            10_000,
            10_000,
            1_000_000,
        )
        .unwrap();

        let mut saw_signal = false;
        for file in &files {
            for (chrom, start, end, _) in bedgraph_rows(file) {
                saw_signal = true;
                assert_eq!(chrom, "5");
                assert!(
                    start >= 65_900_000 && end <= 66_000_000,
                    "interval {start}-{end} escaped the requested window"
                );
            }
        }
        assert!(saw_signal, "the window should contain the test reads");
    }

    #[test]
    fn a_read_lands_at_the_same_coordinates_with_or_without_a_region() {
        // Windowing shifts every bin index. If the offset arithmetic were
        // wrong, the same read would be written at a different coordinate.
        let dir = TempDir::new().unwrap();

        let whole = dir.path().join("whole");
        run_coverage(
            &whole,
            ("test_i1", "test_i2"),
            OutputFormat::BedGraph,
            NormalizeMethod::None,
            1.0,
            None,
            ReadMode::Normal,
            10_000,
            10_000,
            1_000_000,
        )
        .unwrap();

        let region = dir.path().join("region");
        run_coverage(
            &region,
            ("test_i1", "test_i2"),
            OutputFormat::BedGraph,
            NormalizeMethod::None,
            1.0,
            Some("5:65900000-66000000"),
            ReadMode::Normal,
            10_000,
            10_000,
            1_000_000,
        )
        .unwrap();

        let whole_signal = signal_by_interval(&dir.path().join("whole_i2_pool.bedgraph"));
        let region_signal = signal_by_interval(&dir.path().join("region_i2_pool.bedgraph"));

        assert!(!region_signal.is_empty(), "the region carried no signal");
        for (interval, value) in &region_signal {
            let expected = whole_signal.get(interval).copied().unwrap_or(0.0);
            assert!(
                (value - expected).abs() < 1e-6,
                "interval {interval:?} held {value} inside the region but {expected} without it"
            );
        }
    }

    #[test]
    fn a_bare_chromosome_region_covers_the_whole_chromosome() {
        // `parse_region` leaves an open end as usize::MAX; clamping it to the
        // chromosome length must give back the unrestricted track.
        let dir = TempDir::new().unwrap();

        let whole = dir.path().join("whole");
        run_default(&whole).unwrap();

        let bare = dir.path().join("bare");
        run_coverage(
            &bare,
            ("test_i1", "test_i2"),
            OutputFormat::BedGraph,
            NormalizeMethod::None,
            1.0,
            Some("5"),
            ReadMode::Normal,
            100_000,
            100_000,
            1_000_000,
        )
        .unwrap();

        assert_eq!(
            signal_by_interval(&dir.path().join("bare_i2_pool.bedgraph")),
            signal_by_interval(&dir.path().join("whole_i2_pool.bedgraph"))
        );
    }
}
