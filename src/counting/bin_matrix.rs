//! Counts reads into a cell × genomic-bin matrix.
//!
//! Bin geometry is taken from the first BAM's header, so no annotation file is
//! involved. Work is parallel over sub-chromosome chunks, each an independent
//! BAI query.
//!
//! A read is credited to exactly one bin — the one its effective interval
//! overlaps most — and is owned by the chunk holding its alignment start, so a
//! read whose interval spills into the next chunk is still counted once.

use std::path::Path;

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result};
use rayon::prelude::*;

use super::count_utils::{build_csr, write_counts_anndata};
use super::params::{CountingParams, parse_region};
use crate::annotation::parse_annotation::parse_blacklist_bed;
use crate::annotation::region_index::{BinIndex, GenomeIndex, build_bin_index};
use crate::bam::bam_io::{BamWorker, ensure_barcode_tags_present, read_bam_header};
use crate::bam::filters::{
    DupMethod, DuplicateFilter, QcFilter, RawRecordFilter, derive_record_opts,
};
use crate::bam::sc_record::{AdjustRead, ScRecord, parse_tag};

/// Count reads from one or more BAM files into a cell × genomic-bin matrix,
/// then write the result as an AnnData HDF5 file.
///
/// `bam_paths` is a list of `(path, sample_name)` pairs. Bin geometry is
/// derived from the first BAM's header. Parallelism is over sub-chromosome
/// chunks of `chunk_size` bp; each chunk is an independent BAI query. Reads
/// are assigned to chunks by `alignment_start` to avoid double-counting. A
/// read's effective interval can still overlap bins in adjacent chunks.
pub fn count_bam_bins(
    bam_paths: &[(&Path, &str)],
    bin_size: usize,
    step_size: usize,
    barcodes: &[String],
    bc_tag: &str,
    umi_tag: Option<&str>,
    count_tag: Option<&str>,
    params: &CountingParams,
    adjust: &AdjustRead,
    record_filter: Option<&RawRecordFilter>,
    qc_filter: Option<&QcFilter>,
    dup_method: Option<DupMethod>,
    genome_path: Option<&Path>,
    motifs: Option<&[(String, String)]>,
    output_path: &Path,
    compression: &str,
    compression_level: u8,
    num_threads: usize,
    chunk_size: usize,
) -> Result<()> {
    anyhow::ensure!(bin_size > 0, "bin_size must be greater than zero");
    anyhow::ensure!(step_size > 0, "step_size must be greater than zero");
    anyhow::ensure!(!bam_paths.is_empty(), "at least one BAM file is required");
    anyhow::ensure!(chunk_size > 0, "chunk_size must be greater than zero");

    let has_motif = genome_path.is_some() && motifs.is_some();
    let record_opts = derive_record_opts(qc_filter, has_motif);
    let bc_tag_parsed = parse_tag(bc_tag)?;
    let umi_tag_parsed = umi_tag.map(parse_tag).transpose()?;
    let all_bams: Vec<&Path> = bam_paths.iter().map(|(p, _)| *p).collect();
    ensure_barcode_tags_present(&all_bams, bc_tag_parsed, umi_tag_parsed)?;
    let count_tag_parsed = count_tag.map(parse_tag).transpose()?;
    let region_filter = params.region.as_deref().map(parse_region).transpose()?;

    let n_barcodes = barcodes.len();
    // Keyed by raw bytes so the per-read barcode lookup never allocates.
    let barcode_index: AHashMap<&[u8], usize> = barcodes
        .iter()
        .enumerate()
        .map(|(i, bc)| (bc.as_bytes(), i))
        .collect();
    let n_cells = bam_paths.len() * n_barcodes;

    // Chromosomes to skip, as a byte-slice set so membership tests over the
    // header don't allocate a `String` per reference sequence.
    let skip_set: AHashSet<&[u8]> = params.chr_to_skip.iter().map(|s| s.as_bytes()).collect();

    // Bin geometry from first BAM's header.
    let chrom_sizes: Vec<(String, usize)> = {
        let (first_path, _) = bam_paths[0];
        let hdr = read_bam_header(first_path)?;
        hdr.reference_sequences()
            .iter()
            .filter(|(name, _)| {
                let bytes: &[u8] = name.as_ref();
                !skip_set.contains(bytes)
            })
            .map(|(name, seq)| (name.to_string(), seq.length().get()))
            .collect()
    };

    let (bin_index, var_meta) = build_bin_index(&chrom_sizes, bin_size, step_size);
    let n_features = var_meta.len();

    let blacklist = params
        .blacklist_path
        .as_deref()
        .map(parse_blacklist_bed)
        .transpose()?;
    let blacklisted_bins = blacklisted_bin_indices(&bin_index, blacklist.as_ref());

    // Build chunk work list: (bam_idx, bam_path, chrom, chunk_start, chunk_end).
    // Use chrom_sizes (from first BAM) as the master chromosome set.
    // Sort by descending chunk size.
    let mut work: Vec<(usize, &Path, String, usize, usize)> = bam_paths
        .iter()
        .enumerate()
        .flat_map(|(bam_idx, &(bam_path, _))| {
            chrom_sizes.iter().flat_map(move |(chrom, chrom_len)| {
                (0..*chrom_len).step_by(chunk_size).map(move |start| {
                    (
                        bam_idx,
                        bam_path,
                        chrom.clone(),
                        start,
                        (start + chunk_size).min(*chrom_len),
                    )
                })
            })
        })
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

    // Motif-filter ingredients, passed to each worker so it can build its own
    // filter once per thread (rather than once per chunk).
    let motif_ingredients = match (genome_path, motifs) {
        (Some(g), Some(m)) => Some((g, m)),
        _ => None,
    };

    // Count each chunk into its own map, then combine them with a parallel
    // tree reduction.
    let global_acc: AHashMap<(usize, usize), u32> = pool.install(|| {
        work.par_iter()
            .map_init(
                BamWorker::new,
                |worker,
                 &(bam_idx, bam_path, ref chrom, chunk_start, chunk_end)|
                 -> Result<AHashMap<(usize, usize), u32>> {
                    let (reader, header, motif) = worker.prepare(bam_path, motif_ingredients)?;

                    // Each work chunk covers a single chromosome, so its bin
                    // geometry is looked up once here rather than per read.
                    let Some(&(chrom_offset, n_bins)) = bin_index.chrom_bins.get(chrom.as_str())
                    else {
                        return Ok(AHashMap::new());
                    };

                    // Hoisted region filter: if a region was requested on a
                    // different chromosome, the whole chunk is skipped.
                    let region_bounds = match &region_filter {
                        Some((region_chrom, region_start, region_end)) => {
                            if region_chrom.as_str() != chrom.as_str() {
                                return Ok(AHashMap::new());
                            }
                            Some((*region_start, *region_end))
                        }
                        None => None,
                    };

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
                    let cell_offset = bam_idx * n_barcodes;

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
                            count_tag_parsed.as_ref(),
                            &record_opts,
                        )?
                        else {
                            continue;
                        };

                        // Ownership: a read belongs to the chunk containing its
                        // alignment_start. Skip reads owned by a previous chunk.
                        if sc_rec.alignment_start < chunk_start {
                            continue;
                        }

                        if let Some((region_start, region_end)) = region_bounds
                            && (sc_rec.alignment_end <= region_start
                                || sc_rec.alignment_start >= region_end)
                        {
                            continue;
                        }

                        let Some(barcode) = sc_rec.barcode else {
                            continue;
                        };
                        let Some(&local_bc_idx) = barcode_index.get(barcode) else {
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

                        let (eff_start, eff_end) = sc_rec.effective_interval(adjust);
                        let cell_idx = cell_offset + local_bc_idx;

                        // Largest-overlap-wins: assign the read to exactly one
                        // bin, the one whose [bin_start, bin_start+bin_size)
                        // interval overlaps the effective read interval the most.
                        // For non-overlapping bins this is equivalent to placing
                        // the read at the center of its effective interval.
                        if eff_end > 0 && n_bins > 0 {
                            let last_bin = ((eff_end - 1) / step_size).min(n_bins - 1);
                            let first_bin = if eff_start + 1 > bin_size {
                                (eff_start - bin_size + step_size) / step_size
                            } else {
                                0
                            };
                            let mut best_bin = first_bin;
                            let mut best_overlap = 0usize;
                            for bin in first_bin..=last_bin {
                                let bin_start = bin * step_size;
                                let bin_end = bin_start + bin_size;
                                let overlap = eff_end.min(bin_end) - eff_start.max(bin_start);
                                if overlap > best_overlap {
                                    best_overlap = overlap;
                                    best_bin = bin;
                                }
                            }
                            if best_overlap > 0 {
                                let feature_idx = chrom_offset + best_bin;
                                if blacklisted_bins.is_empty()
                                    || !blacklisted_bins.contains(&feature_idx)
                                {
                                    *local_acc.entry((cell_idx, feature_idx)).or_insert(0) +=
                                        sc_rec.count;
                                }
                            }
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

    let matrix = build_csr(&global_acc, n_cells, n_features)?;

    write_counts_anndata(
        output_path,
        matrix,
        bam_paths,
        barcodes,
        &var_meta,
        compression,
        compression_level,
    )?;

    Ok(())
}

fn blacklisted_bin_indices(
    bin_index: &BinIndex,
    blacklist: Option<&GenomeIndex>,
) -> AHashSet<usize> {
    let Some(bl) = blacklist else {
        return AHashSet::new();
    };

    let bin_size = bin_index.bin_size;
    let step_size = bin_index.step_size;
    let mut set = AHashSet::new();

    for (chrom, bl_chrom_idx) in bl {
        let Some(&(chrom_offset, n_bins)) = bin_index.chrom_bins.get(chrom) else {
            continue;
        };
        for bl_iv in bl_chrom_idx.iter() {
            if bl_iv.end == 0 {
                continue;
            }
            let last_bl_bin = ((bl_iv.end - 1) / step_size).min(n_bins - 1);
            let first_bl_bin = if bl_iv.start + 1 > bin_size {
                (bl_iv.start - bin_size + step_size) / step_size
            } else {
                0
            };
            for bin in first_bl_bin..=last_bl_bin {
                set.insert(chrom_offset + bin);
            }
        }
    }
    set
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::region_index::{ChromIndex, Interval};

    fn blacklist(chrom: &str, spans: &[(usize, usize)]) -> GenomeIndex {
        let intervals = spans
            .iter()
            .map(|&(start, end)| Interval {
                start,
                end,
                var_idx: 0,
            })
            .collect();
        [(chrom.to_string(), ChromIndex::build(intervals))]
            .into_iter()
            .collect()
    }

    /// chr1 (300 bp) and chr2 (100 bp) tiled into 100 bp bins:
    /// feature indices 0,1,2 = chr1 bins; 3 = chr2 bin.
    fn two_chrom_index() -> BinIndex {
        let (index, _) = build_bin_index(
            &[("chr1".to_string(), 300), ("chr2".to_string(), 100)],
            100,
            100,
        );
        index
    }

    fn sorted(set: AHashSet<usize>) -> Vec<usize> {
        let mut v: Vec<usize> = set.into_iter().collect();
        v.sort_unstable();
        v
    }

    #[test]
    fn no_blacklist_means_no_excluded_bins() {
        assert!(blacklisted_bin_indices(&two_chrom_index(), None).is_empty());
    }

    #[test]
    fn a_blacklist_region_excludes_every_bin_it_touches() {
        // [150, 160) lies inside the second chr1 bin.
        let excluded =
            blacklisted_bin_indices(&two_chrom_index(), Some(&blacklist("chr1", &[(150, 160)])));
        assert_eq!(sorted(excluded), vec![1]);

        // A region spanning a bin boundary excludes both bins.
        let excluded =
            blacklisted_bin_indices(&two_chrom_index(), Some(&blacklist("chr1", &[(90, 110)])));
        assert_eq!(sorted(excluded), vec![0, 1]);

        // A region covering the whole chromosome excludes all of its bins.
        let excluded =
            blacklisted_bin_indices(&two_chrom_index(), Some(&blacklist("chr1", &[(0, 300)])));
        assert_eq!(sorted(excluded), vec![0, 1, 2]);
    }

    #[test]
    fn excluded_bins_are_offset_by_chromosome() {
        // chr2's only bin is feature index 3, not 0.
        let excluded =
            blacklisted_bin_indices(&two_chrom_index(), Some(&blacklist("chr2", &[(10, 20)])));
        assert_eq!(sorted(excluded), vec![3]);
    }

    #[test]
    fn blacklist_regions_on_unknown_chromosomes_are_ignored() {
        let excluded = blacklisted_bin_indices(
            &two_chrom_index(),
            Some(&blacklist("chrUnplaced", &[(0, 50)])),
        );
        assert!(excluded.is_empty());
    }

    #[test]
    fn a_blacklist_region_past_the_chromosome_end_clamps_to_the_last_bin() {
        let excluded =
            blacklisted_bin_indices(&two_chrom_index(), Some(&blacklist("chr1", &[(250, 5000)])));
        assert_eq!(sorted(excluded), vec![2]);
    }

    #[test]
    fn sliding_bins_that_overlap_a_blacklist_region_are_all_excluded() {
        // 100 bp bins every 50 bp on a 300 bp chromosome: bins start at
        // 0, 50, 100, 150, 200, 250. A region at [120, 130) is covered by the
        // bins starting at 50 and 100.
        let (index, _) = build_bin_index(&[("chr1".to_string(), 300)], 100, 50);
        let excluded = blacklisted_bin_indices(&index, Some(&blacklist("chr1", &[(120, 130)])));
        assert_eq!(sorted(excluded), vec![1, 2]);
    }

    // Counting a real BAM end to end

    use anndata::{AnnData, AnnDataOp, Backend};
    use anndata_hdf5::H5;
    use std::path::PathBuf;
    use tempfile::TempDir;

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

    /// `count_bam_bins` with everything optional switched off.
    #[allow(clippy::too_many_arguments)]
    fn count_into(
        out: &Path,
        bin_size: usize,
        step_size: usize,
        barcodes: &[String],
        params: &CountingParams,
        chunk_size: usize,
    ) -> Result<()> {
        let bam = testdata().join("test_i1.bam");
        count_bam_bins(
            &[(bam.as_path(), "s1")],
            bin_size,
            step_size,
            barcodes,
            "BC",
            None,
            None,
            params,
            &AdjustRead::default(),
            None,
            None,
            None,
            None,
            None,
            out,
            "none",
            0,
            1,
            chunk_size,
        )
    }

    #[test]
    fn counting_bins_writes_one_matrix_row_per_cell_and_one_column_per_bin() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("bins.h5ad");
        let barcodes = test_barcodes();

        count_into(
            &out,
            100_000,
            100_000,
            &barcodes,
            &CountingParams::default(),
            1_000_000,
        )
        .unwrap();

        assert!(out.exists(), "no output file was written");
        let adata = AnnData::<H5>::open(H5::open(&out).unwrap()).unwrap();
        assert_eq!(adata.n_obs(), barcodes.len());
        assert!(adata.n_vars() > 0, "expected at least one bin");
    }

    #[test]
    fn restricting_to_a_region_keeps_the_bin_geometry_but_counts_fewer_reads() {
        let barcodes = test_barcodes();
        let dir = TempDir::new().unwrap();

        let whole = dir.path().join("whole.h5ad");
        count_into(
            &whole,
            100_000,
            100_000,
            &barcodes,
            &CountingParams::default(),
            1_000_000,
        )
        .unwrap();

        // The test BAM's reads sit around 65.97 Mb, so this window holds none.
        let region_params = CountingParams {
            region: Some("5:1-100000".to_string()),
            ..CountingParams::default()
        };
        let region = dir.path().join("region.h5ad");
        count_into(
            &region,
            100_000,
            100_000,
            &barcodes,
            &region_params,
            1_000_000,
        )
        .unwrap();

        let open = |p: &Path| AnnData::<H5>::open(H5::open(p).unwrap()).unwrap();
        let total = |p: &Path| -> f64 {
            super::super::count_utils::read_x_f64(&open(p))
                .unwrap()
                .triplet_iter()
                .map(|(_, _, &v)| v)
                .sum()
        };

        // `region` filters reads, not the bin layout: the var axis still spans
        // the whole chromosome, so both runs are the same width.
        assert_eq!(open(&region).n_vars(), open(&whole).n_vars());

        let whole_total = total(&whole);
        assert!(whole_total > 0.0, "the unrestricted run counted nothing");
        assert!(
            total(&region) < whole_total,
            "restricting to an empty window did not reduce the counts"
        );
    }

    #[test]
    fn skipping_the_only_chromosome_leaves_no_bins_to_count() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("skipped.h5ad");
        let params = CountingParams {
            chr_to_skip: vec!["5".to_string()],
            ..CountingParams::default()
        };

        // The test BAM has only chromosome 5, so skipping it removes everything.
        let result = count_into(&out, 100_000, 100_000, &test_barcodes(), &params, 1_000_000);
        match result {
            Ok(()) => {
                let adata = AnnData::<H5>::open(H5::open(&out).unwrap()).unwrap();
                assert_eq!(adata.n_vars(), 0);
            }
            // Refusing outright is also a defensible answer for an empty genome.
            Err(e) => {
                let msg = e.to_string();
                assert!(!msg.is_empty());
            }
        }
    }

    #[test]
    fn a_smaller_chunk_size_does_not_change_the_result() {
        // Chunking is a parallelism detail: reads are owned by the chunk holding
        // their alignment start, so the totals must not depend on chunk size.
        let barcodes = test_barcodes();
        let dir = TempDir::new().unwrap();

        let big = dir.path().join("big.h5ad");
        count_into(
            &big,
            50_000,
            50_000,
            &barcodes,
            &CountingParams::default(),
            10_000_000,
        )
        .unwrap();

        let small = dir.path().join("small.h5ad");
        count_into(
            &small,
            50_000,
            50_000,
            &barcodes,
            &CountingParams::default(),
            250_000,
        )
        .unwrap();

        let sum = |p: &Path| -> f64 {
            let adata = AnnData::<H5>::open(H5::open(p).unwrap()).unwrap();
            super::super::count_utils::read_x_f64(&adata)
                .unwrap()
                .triplet_iter()
                .map(|(_, _, &v)| v)
                .sum()
        };
        assert_eq!(sum(&big), sum(&small));
    }

    // Argument validation

    #[test]
    fn zero_sized_bins_and_steps_are_rejected() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("unused.h5ad");
        let barcodes = test_barcodes();
        let params = CountingParams::default();

        let err = count_into(&out, 0, 100, &barcodes, &params, 1_000)
            .unwrap_err()
            .to_string();
        assert!(err.contains("bin_size"), "{err}");

        let err = count_into(&out, 100, 0, &barcodes, &params, 1_000)
            .unwrap_err()
            .to_string();
        assert!(err.contains("step_size"), "{err}");
    }

    #[test]
    fn a_zero_chunk_size_is_rejected() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("unused.h5ad");
        let err = count_into(
            &out,
            100,
            100,
            &test_barcodes(),
            &CountingParams::default(),
            0,
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("chunk_size"), "{err}");
    }

    #[test]
    fn counting_with_no_bam_files_is_rejected() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("unused.h5ad");
        let err = count_bam_bins(
            &[],
            100,
            100,
            &test_barcodes(),
            "BC",
            None,
            None,
            &CountingParams::default(),
            &AdjustRead::default(),
            None,
            None,
            None,
            None,
            None,
            &out,
            "none",
            0,
            1,
            1_000,
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("at least one BAM"), "{err}");
    }

    #[test]
    fn a_barcode_tag_the_bam_does_not_carry_is_rejected_before_any_counting() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("unused.h5ad");
        let bam = testdata().join("test_i1.bam");

        let err = count_bam_bins(
            &[(bam.as_path(), "s1")],
            100_000,
            100_000,
            &test_barcodes(),
            "ZZ",
            None,
            None,
            &CountingParams::default(),
            &AdjustRead::default(),
            None,
            None,
            None,
            None,
            None,
            &out,
            "none",
            0,
            1,
            1_000_000,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("ZZ"), "{err}");
        assert!(
            !out.exists(),
            "nothing should be written when the tag is wrong"
        );
    }
}
