use std::ops::AddAssign;
use std::path::{Path, PathBuf};

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::annotation::parse_annotation::parse_blacklist_bed;
use crate::annotation::region_index::GenomeIndex;
use crate::bam::bam_io::{
    BamWorker, ensure_barcode_tags_present, read_bam_header, read_group_ids, warn_unknown_group,
};
use crate::bam::filters::{DupMethod, DuplicateFilter, is_blacklisted, rna_strand_filter};
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

pub fn run_filter_stats(
    bam_path: &Path,
    barcodes: &[String],
    bc_tag: &str,
    umi_tag: Option<&str>,
    group_tag: Option<&str>,
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
        compute_covered_span: dup_method.is_some(),
        // This tool counts reads, not intervals, so an alignment's shape does
        // not matter to it.
        compute_blocks: false,
    };

    let blacklist: Option<GenomeIndex> = blacklist_path.map(parse_blacklist_bed).transpose()?;

    // Keyed by raw bytes so the per-read barcode lookup never allocates.
    let barcode_idx: AHashMap<&[u8], usize> = barcodes
        .iter()
        .enumerate()
        .map(|(i, bc)| (bc.as_bytes(), i))
        .collect();

    let header = read_bam_header(bam_path)?;

    // With --groupTag the row unit is `group::barcode`: a merged BAM's reads
    // carry their sample of origin, and the valid values are its @RG IDs.
    let group_tag_parsed = group_tag.map(parse_tag).transpose()?;
    let group_ids: Option<Vec<Vec<u8>>> = match group_tag {
        Some(_) => Some(read_group_ids(&header, bam_path)?),
        None => None,
    };
    let group_index: AHashMap<&[u8], usize> = group_ids
        .iter()
        .flatten()
        .enumerate()
        .map(|(i, id)| (id.as_slice(), i))
        .collect();

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
    // One row block per group when grouping, otherwise a single block.
    let n_rows = group_ids.as_ref().map_or(1, Vec::len) * n_barcodes;

    // Build chunk work list sorted by descending size.
    // Each chunk iterates the sampling bins whose start falls within it, and
    // carries its chromosome's length to bound the bins.
    let mut chunks: Vec<(String, usize, usize, usize)> = chrom_sizes
        .iter()
        .flat_map(|(chrom, chrom_len)| {
            (0..*chrom_len).step_by(chunk_size).map(move |start| {
                (
                    chrom.clone(),
                    start,
                    (start + chunk_size).min(*chrom_len),
                    *chrom_len,
                )
            })
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
                |worker,
                 &(ref chrom, chunk_start, chunk_end, chrom_len)|
                 -> Result<Vec<BarcodeStat>> {
                    let (reader, header, motif) = worker.prepare(bam_path, motif_ingredients)?;

                    let mut dup_filter: Option<DuplicateFilter> =
                        dup_method.map(DuplicateFilter::new);

                    let mut local_stats: Vec<BarcodeStat> =
                        (0..n_rows).map(|_| BarcodeStat::default()).collect();

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
                        // A bin is bounded by its chromosome, not by the chunk
                        // that owns it. Cutting it at `chunk_end` would drop the
                        // rest of a bin that reaches past the boundary. Ownership
                        // goes by the bin's *start*, which lies back here. Reading
                        // past the edge cannot double count for the same reason.
                        let bin_end = (bin_start + bin_size).min(chrom_len);

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
                                group_tag_parsed.as_ref(),
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
                            let Some(&local_bc) = barcode_idx.get(barcode) else {
                                continue;
                            };
                            // Under --groupTag the read's own group picks the row
                            // block, so a barcode shared by two source samples
                            // stays two rows.
                            let cell_i = if group_tag_parsed.is_some() {
                                let Some(group) = sc_rec.group else {
                                    continue;
                                };
                                let Some(&group_i) = group_index.get(group) else {
                                    warn_unknown_group(group);
                                    continue;
                                };
                                group_i * n_barcodes + local_bc
                            } else {
                                local_bc
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

    let mut stats: Vec<BarcodeStat> = (0..n_rows).map(|_| BarcodeStat::default()).collect();
    for partial in partial_stats {
        for (i, s) in partial.into_iter().enumerate() {
            stats[i] += s;
        }
    }

    let stat_vecs: Vec<Vec<u64>> = stats.iter().map(|s| s.to_vec()).collect();

    // Row labels: bare barcodes normally, `group::barcode` when grouping, so the
    // caller can use them as Cell_IDs without knowing which mode ran.
    let row_labels: Vec<String> = match &group_ids {
        None => barcodes.to_vec(),
        Some(ids) => ids
            .iter()
            .flat_map(|id| {
                let group = String::from_utf8_lossy(id).into_owned();
                barcodes.iter().map(move |bc| format!("{group}::{bc}"))
            })
            .collect(),
    };
    Ok((row_labels, stat_vecs))
}

#[pyfunction(signature = (
    bam_path,
    barcodes,
    bc_tag = "CB",
    umi_tag = None,
    group_tag = None,
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
    group_tag: Option<String>,
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
        group_tag.as_deref(),
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

#[cfg(test)]
mod tests {
    use super::*;

    // Column order of BarcodeStat::to_vec, used to index the returned rows.
    const TOTAL: usize = 0;
    const FILTERED: usize = 1;
    const BLACKLISTED: usize = 2;
    const LOW_MAPQ: usize = 3;
    const N_STATS: usize = 13;

    fn testdata() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/testdata")
    }

    fn test_barcodes() -> Vec<String> {
        std::fs::read_to_string(testdata().join("test_barcodes.txt"))
            .unwrap()
            .lines()
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty())
            .collect()
    }

    // Stat accumulation

    #[test]
    fn the_stat_vector_lists_every_counter_in_a_fixed_order() {
        let stat = BarcodeStat {
            total: 1,
            filtered: 2,
            blacklisted: 3,
            low_mapq: 4,
            missing_flags: 5,
            excluded_flags: 6,
            internal_dupes: 7,
            external_dupes: 8,
            singletons: 9,
            wrong_strand: 10,
            wrong_motif: 11,
            wrong_gc: 12,
            low_aligned_fraction: 13,
        };

        assert_eq!(stat.to_vec(), (1..=13).collect::<Vec<u64>>());
        assert_eq!(stat.to_vec().len(), N_STATS);
    }

    #[test]
    fn adding_two_stat_sets_sums_every_counter() {
        let mut a = BarcodeStat {
            total: 1,
            filtered: 1,
            blacklisted: 1,
            low_mapq: 1,
            missing_flags: 1,
            excluded_flags: 1,
            internal_dupes: 1,
            external_dupes: 1,
            singletons: 1,
            wrong_strand: 1,
            wrong_motif: 1,
            wrong_gc: 1,
            low_aligned_fraction: 1,
        };
        let b = BarcodeStat {
            total: 2,
            filtered: 2,
            blacklisted: 2,
            low_mapq: 2,
            missing_flags: 2,
            excluded_flags: 2,
            internal_dupes: 2,
            external_dupes: 2,
            singletons: 2,
            wrong_strand: 2,
            wrong_motif: 2,
            wrong_gc: 2,
            low_aligned_fraction: 2,
        };

        a += b;
        assert_eq!(a.to_vec(), vec![3u64; N_STATS]);
    }

    // Whole run against the test BAM

    #[allow(clippy::too_many_arguments)]
    fn stats(
        barcodes: &[String],
        bc_tag: &str,
        min_mapq: Option<u8>,
        chr_to_skip: &[String],
        bin_size: usize,
        chunk_size: usize,
    ) -> Result<(Vec<String>, Vec<Vec<u64>>)> {
        run_filter_stats(
            &testdata().join("test_i1.bam"),
            barcodes,
            bc_tag,
            None,
            None,
            bin_size,
            0,
            min_mapq,
            None,
            None,
            chr_to_skip,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            1,
            chunk_size,
        )
    }

    #[test]
    fn one_row_of_counters_comes_back_per_barcode() {
        let barcodes = test_barcodes();
        let (returned, rows) = stats(&barcodes, "BC", None, &[], 2_000, 1_000_000).unwrap();

        assert_eq!(returned, barcodes, "the barcodes come back in input order");
        assert_eq!(rows.len(), barcodes.len());
        for row in &rows {
            assert_eq!(row.len(), N_STATS);
        }
    }

    #[test]
    fn barcodes_in_the_file_accumulate_reads_and_absent_ones_stay_at_zero() {
        // ATATAACT is in test_i1.bam; the other is not a barcode at all.
        let barcodes = vec!["ATATAACT".to_string(), "TTTTTTTT".to_string()];
        let (_, rows) = stats(&barcodes, "BC", None, &[], 2_000, 1_000_000).unwrap();

        assert!(rows[0][TOTAL] > 0, "a real barcode saw no reads");
        assert_eq!(
            rows[1],
            vec![0u64; N_STATS],
            "an absent barcode counted reads"
        );
    }

    #[test]
    fn no_counter_can_exceed_the_total_it_is_drawn_from() {
        let (_, rows) = stats(&test_barcodes(), "BC", None, &[], 2_000, 1_000_000).unwrap();

        for row in &rows {
            assert!(row[FILTERED] <= row[TOTAL], "filtered {row:?}");
            assert!(row[BLACKLISTED] <= row[TOTAL], "blacklisted {row:?}");
            assert!(row[LOW_MAPQ] <= row[TOTAL], "low_mapq {row:?}");
        }
    }

    #[test]
    fn an_impossible_mapping_quality_threshold_rejects_every_read() {
        let barcodes = test_barcodes();

        let (_, lenient) = stats(&barcodes, "BC", None, &[], 2_000, 1_000_000).unwrap();
        let (_, strict) = stats(&barcodes, "BC", Some(255), &[], 2_000, 1_000_000).unwrap();

        let sum = |rows: &[Vec<u64>], col: usize| -> u64 { rows.iter().map(|r| r[col]).sum() };

        assert_eq!(
            sum(&lenient, LOW_MAPQ),
            0,
            "no threshold, no low-mapq reads"
        );
        assert!(
            sum(&strict, LOW_MAPQ) > 0,
            "a 255 threshold rejected nothing"
        );
        assert_eq!(
            sum(&strict, TOTAL),
            sum(&lenient, TOTAL),
            "the total is counted before filtering, so it must not move"
        );
    }

    #[test]
    fn skipping_the_only_chromosome_leaves_nothing_to_count() {
        let barcodes = test_barcodes();
        let (_, rows) = stats(&barcodes, "BC", None, &["5".to_string()], 2_000, 1_000_000).unwrap();

        let total: u64 = rows.iter().map(|r| r[TOTAL]).sum();
        assert_eq!(total, 0, "chromosome 5 was skipped but reads were counted");
    }

    #[test]
    fn the_work_chunk_size_does_not_change_what_is_sampled() {
        // Chunking is how the work is spread over threads, nothing more. A bin
        // reaching past a chunk boundary must still be sampled whole, so the
        // answer cannot depend on where those boundaries fall.
        //
        // Every read of the file sits between 65.96 and 65.98 Mb, so at a
        // 200 kb sampling bin they all fall in the bin starting at 65.8 Mb.
        // 65,900,000 is the chunk size that matters here: it puts a boundary
        // inside that bin, which is the case that used to lose them. The other
        // sizes leave the bin whole and are here to keep the property honest.
        let barcodes = test_barcodes();

        let reference = stats(&barcodes, "BC", None, &[], 200_000, 4_000_000).unwrap();
        for chunk_size in [65_900_000, 300_000, 700_000, 1_000_001] {
            let other = stats(&barcodes, "BC", None, &[], 200_000, chunk_size).unwrap();
            assert_eq!(
                other, reference,
                "a chunk size of {chunk_size} changed the sampling"
            );
        }
    }

    #[test]
    fn a_bin_wider_than_a_work_chunk_is_still_sampled_whole() {
        // The reads of test_i1.bam sit at 65.97 Mb, inside the bin that starts
        // at 64 Mb when the bin is 2 Mb wide. Cutting that bin at its chunk's
        // end used to hide every one of them, and the run reported nothing.
        let barcodes = test_barcodes();

        let (_, rows) = stats(&barcodes, "BC", None, &[], 2_000_000, 1_000_000).unwrap();
        let sampled: u64 = rows.iter().map(|r| r[TOTAL]).sum();

        assert_eq!(sampled, 13, "the BAM's 13 reads were not all sampled");
    }

    // Argument validation

    #[test]
    fn a_zero_bin_size_is_rejected() {
        let err = stats(&test_barcodes(), "BC", None, &[], 0, 1_000_000)
            .unwrap_err()
            .to_string();
        assert!(err.contains("bin_size"), "{err}");
    }

    #[test]
    fn a_barcode_tag_the_file_does_not_carry_is_rejected() {
        let err = stats(&test_barcodes(), "ZZ", None, &[], 2_000, 1_000_000)
            .unwrap_err()
            .to_string();
        assert!(err.contains("ZZ"), "{err}");
    }

    #[test]
    fn a_malformed_tag_name_is_rejected() {
        assert!(stats(&test_barcodes(), "B", None, &[], 2_000, 1_000_000).is_err());
        assert!(stats(&test_barcodes(), "BCX", None, &[], 2_000, 1_000_000).is_err());
    }
}
