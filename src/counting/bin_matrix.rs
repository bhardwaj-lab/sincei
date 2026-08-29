//! Counts reads into a cell × genomic-bin matrix.
//!
//! Bin geometry is taken from the first BAM's header, so no annotation file is
//! involved. Work is parallel over sub-chromosome chunks, each an independent
//! BAI query.
//!
//! A read is credited to every bin its effective interval overlaps, and is
//! owned by the chunk holding its alignment start, so a read whose interval
//! spills into the next chunk is still counted exactly once per bin.

use std::path::Path;

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result};
use rayon::prelude::*;

use super::count_utils::{build_csr, write_counts_anndata};
use super::params::{CountingParams, parse_region};
use crate::annotation::parse_annotation::parse_blacklist_bed;
use crate::annotation::region_index::{bins_touched, build_bin_index, build_bin_index_in_window};
use crate::bam::bam_io::{
    BamWorker, ensure_barcode_tags_present, ensure_genome_matches_bams, read_bam_header,
    read_group_ids, warn_unknown_group,
};
use crate::bam::filters::{
    DupMethod, DuplicateFilter, QcFilter, RawRecordFilter, blacklist_chrom_index,
    derive_record_opts, read_is_blacklisted,
};
use crate::bam::sc_record::{AdjustRead, ScRecord, parse_tag};

/// Count reads from one or more BAM files into a cell × genomic-bin matrix,
/// then write the result as an AnnData HDF5 file.
///
/// `bam_paths` is a list of `(path, sample_name)` pairs. Bin geometry is
/// derived from the first BAM's header. Parallelism is over sub-chromosome
/// chunks of `chunk_size` bp; each chunk is an independent BAI query. Reads
/// are assigned to chunks by `alignment_start` to avoid double-counting. A
/// read's effective interval can still overlap bins in adjacent chunks, and it
/// is counted in each of them.
pub fn count_bam_bins(
    bam_paths: &[(&Path, &str)],
    bin_size: usize,
    step_size: usize,
    barcodes: &[String],
    bc_tag: &str,
    umi_tag: Option<&str>,
    count_tag: Option<&str>,
    group_tag: Option<&str>,
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
    let record_opts = derive_record_opts(qc_filter, has_motif, dup_method.is_some(), adjust);
    let bc_tag_parsed = parse_tag(bc_tag)?;
    let umi_tag_parsed = umi_tag.map(parse_tag).transpose()?;
    let all_bams: Vec<&Path> = bam_paths.iter().map(|(p, _)| *p).collect();
    ensure_barcode_tags_present(&all_bams, bc_tag_parsed, umi_tag_parsed)?;

    // A genome from the wrong assembly makes every motif lookup wrong, so this
    // is checked once here rather than a silent loss of reads.
    if has_motif && let Some(genome) = genome_path {
        ensure_genome_matches_bams(genome, &all_bams)?;
    }

    let count_tag_parsed = count_tag.map(parse_tag).transpose()?;
    let group_tag_parsed = group_tag.map(parse_tag).transpose()?;
    let region_filter = params.region.as_deref().map(parse_region).transpose()?;

    let n_barcodes = barcodes.len();
    // Keyed by raw bytes so the per-read barcode lookup never allocates.
    let barcode_index: AHashMap<&[u8], usize> = barcodes
        .iter()
        .enumerate()
        .map(|(i, bc)| (bc.as_bytes(), i))
        .collect();

    // With --groupTag the sample axis comes from the reads' group tag rather
    // than from their BAM of origin, so exactly one input is allowed and the
    // row space is the header's @RG IDs x barcodes.
    let group_ids: Option<Vec<Vec<u8>>> = if group_tag.is_some() {
        anyhow::ensure!(
            bam_paths.len() == 1,
            "--groupTag expects a single merged BAM, but {} were given",
            bam_paths.len()
        );
        let (path, _) = bam_paths[0];
        Some(read_group_ids(&read_bam_header(path)?, path)?)
    } else {
        None
    };
    let group_index: AHashMap<&[u8], usize> = group_ids
        .iter()
        .flatten()
        .enumerate()
        .map(|(i, id)| (id.as_slice(), i))
        .collect();

    // Row labels, and hence the row count: group IDs when grouping, otherwise
    // one block per input BAM.
    let sample_labels: Vec<String> = match &group_ids {
        Some(ids) => ids
            .iter()
            .map(|id| String::from_utf8_lossy(id).into_owned())
            .collect(),
        None => bam_paths.iter().map(|(_, l)| (*l).to_string()).collect(),
    };
    let n_cells = sample_labels.len() * n_barcodes;

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
                if skip_set.contains(bytes) {
                    return false;
                }
                // A region confines the output to its own chromosome, so the
                // other chromosomes contribute no bins at all.
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

    let (bin_index, var_meta) = match &region_filter {
        // No region: tile whole chromosomes, exactly as before.
        None => build_bin_index(&chrom_sizes, bin_size, step_size),
        // A region: tile only that window, so `var` holds the region's bins.
        Some(_) => build_bin_index_in_window(&windows, bin_size, step_size),
    };
    let n_features = var_meta.len();

    let blacklist = params
        .blacklist_path
        .as_deref()
        .map(parse_blacklist_bed)
        .transpose()?;

    // Build chunk work list: (bam_idx, bam_path, chrom, chunk_start, chunk_end).
    // Use chrom_sizes (from first BAM) as the master chromosome set.
    // Sort by descending chunk size.
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
                    // Bin 0 of this chromosome starts here: 0 unless a region
                    // restricted the tiling to a window.
                    let window_start = bin_index.window_start(chrom.as_str());
                    // Likewise the blacklist: one chromosome per chunk, so a
                    // read costs an interval query and no name hashing.
                    let chunk_blacklist = blacklist
                        .as_ref()
                        .and_then(|bl| blacklist_chrom_index(bl, chrom.as_str()));

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
                    // Without --groupTag the sample is fixed for the whole chunk, so
                    // the row offset is hoisted; with it the sample varies per read
                    // and the offset is resolved below.
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
                            group_tag_parsed.as_ref(),
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

                        // Judged on the read's own alignment span, as the filter
                        // tools judge it, and before --extendReads / --centerReads
                        // move the interval that is counted.
                        if let Some(bl) = chunk_blacklist
                            && read_is_blacklisted(bl, sc_rec.alignment_start, sc_rec.alignment_end)
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

                        // Under --groupTag the read's own group tag picks the row
                        // block, so two reads sharing a barcode but coming from
                        // different source samples stay separate cells.
                        let cell_idx = if group_tag_parsed.is_some() {
                            let Some(group) = sc_rec.group else {
                                continue;
                            };
                            let Some(&group_idx) = group_index.get(group) else {
                                warn_unknown_group(group);
                                continue;
                            };
                            group_idx * n_barcodes + local_bc_idx
                        } else {
                            cell_offset + local_bc_idx
                        };

                        // A read counts once in every bin it overlaps, so a read
                        // straddling a boundary adds to both sides. The column
                        // sums therefore exceed the read count, which is what the
                        // reference implementation reports. A spliced read
                        // contributes its blocks, not the intron between them,
                        // and a bin reached by two of its blocks is still counted
                        // once.
                        for bin in bins_touched(
                            sc_rec.effective_intervals(adjust),
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

    let matrix = build_csr(&global_acc, n_cells, n_features)?;

    write_counts_anndata(
        output_path,
        matrix,
        &sample_labels,
        barcodes,
        &var_meta,
        compression,
        compression_level,
    )?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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

    /// Column sums keyed by bin name, so two runs with different bin layouts
    /// can still be compared bin for bin.
    fn column_sums(path: &Path) -> std::collections::BTreeMap<String, f64> {
        let adata = AnnData::<H5>::open(H5::open(path).unwrap()).unwrap();
        let names = adata.var_names().into_vec();
        let x = super::super::count_utils::read_x_f64(&adata).unwrap();
        let mut sums = std::collections::BTreeMap::new();
        for (_, col, &v) in x.triplet_iter() {
            *sums.entry(names[col].clone()).or_insert(0.0) += v;
        }
        sums
    }

    fn n_vars(path: &Path) -> usize {
        AnnData::<H5>::open(H5::open(path).unwrap())
            .unwrap()
            .n_vars()
    }

    /// Every matrix entry, keyed by cell row and bin name, so two runs can be
    /// compared entry for entry.
    fn entries(path: &Path) -> std::collections::BTreeMap<(usize, String), f64> {
        let adata = AnnData::<H5>::open(H5::open(path).unwrap()).unwrap();
        let names = adata.var_names().into_vec();
        let x = super::super::count_utils::read_x_f64(&adata).unwrap();
        let mut out = std::collections::BTreeMap::new();
        for (row, col, &v) in x.triplet_iter() {
            *out.entry((row, names[col].clone())).or_insert(0.0) += v;
        }
        out
    }

    /// Every value in the matrix, summed.
    fn total(path: &Path) -> f64 {
        column_sums(path).values().sum()
    }

    #[test]
    fn a_read_counts_in_every_bin_it_overlaps() {
        // test_i1.bam holds 13 reads of 47-62 bp, all inside one 100 kb bin, so
        // a coarse run counts each of them exactly once.
        //
        // At 100 bp six of those reads straddle a bin boundary. Each of the six
        // adds one to the bin on either side, which takes the total from 13 to
        // 19. Crediting a read to its best bin alone would leave it at 13.
        let barcodes = test_barcodes();
        let dir = TempDir::new().unwrap();

        let coarse = dir.path().join("coarse.h5ad");
        count_into(
            &coarse,
            100_000,
            100_000,
            &barcodes,
            &CountingParams::default(),
            1_000_000,
        )
        .unwrap();
        assert_eq!(total(&coarse), 13.0, "one count per read was expected");

        let fine = dir.path().join("fine.h5ad");
        count_into(
            &fine,
            100,
            100,
            &barcodes,
            &CountingParams::default(),
            1_000_000,
        )
        .unwrap();
        assert_eq!(total(&fine), 19.0, "a straddling read counts on both sides");
    }

    #[test]
    fn a_read_lands_in_two_neighbouring_bins_at_once() {
        // The two reads at 65,966,811 are 61 bp long. A 50 bp tiling puts a
        // boundary at 65,966,850, inside them, so both bins must hold them.
        let barcodes = test_barcodes();
        let dir = TempDir::new().unwrap();

        let out = dir.path().join("split.h5ad");
        count_into(
            &out,
            50,
            50,
            &barcodes,
            &CountingParams::default(),
            1_000_000,
        )
        .unwrap();

        let sums = column_sums(&out);
        assert_eq!(
            sums.get("5:65966800-65966850").copied(),
            Some(2.0),
            "the bin holding the read starts counted {sums:?}"
        );
        assert_eq!(
            sums.get("5:65966850-65966900").copied(),
            Some(2.0),
            "the bin holding the read ends counted {sums:?}"
        );
    }

    /// `count_bam_bins` over a named BAM, everything optional off.
    fn count_bam_into(
        out: &Path,
        bam: &Path,
        bin_size: usize,
        barcodes: &[String],
        adjust: &AdjustRead,
    ) -> Result<()> {
        count_bam_bins(
            &[(bam, "s1")],
            bin_size,
            bin_size,
            barcodes,
            "BC",
            None,
            None,
            None,
            &CountingParams::default(),
            adjust,
            None,
            None,
            None,
            None,
            None,
            out,
            "none",
            0,
            1,
            1_000_000,
        )
    }

    #[test]
    fn a_spliced_read_counts_its_exons_and_not_its_intron() {
        // test_spliced.bam holds two reads with a 20 kb intron whose exons land
        // in bins 65,960,000 and 65,980,000, one plain read in the first bin,
        // and one read whose 10 bp deletion leaves both its blocks in the
        // middle bin. At a 10 kb tiling the intron covers that middle bin
        // entirely.
        let barcodes = test_barcodes();
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("spliced.h5ad");

        count_bam_into(
            &out,
            &testdata().join("test_spliced.bam"),
            10_000,
            &barcodes,
            &AdjustRead::default(),
        )
        .unwrap();

        let sums = column_sums(&out);
        assert_eq!(
            sums.get("5:65960000-65970000").copied(),
            Some(3.0),
            "two spliced first exons and the plain read {sums:?}"
        );
        assert_eq!(
            sums.get("5:65980000-65990000").copied(),
            Some(2.0),
            "two spliced second exons {sums:?}"
        );
        // The middle bin is pure intron for the spliced reads. Only the
        // deletion read, which really covers it, is counted there.
        assert_eq!(
            sums.get("5:65970000-65980000").copied(),
            Some(1.0),
            "the intron was counted {sums:?}"
        );
        assert_eq!(total(&out), 6.0, "{sums:?}");
    }

    #[test]
    fn two_blocks_of_one_read_in_one_bin_count_once() {
        // The deletion read's blocks are 10 bp apart, so a 10 kb bin holds both.
        // Counting per block rather than per read would make it 2.
        let barcodes = test_barcodes();
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("deleted.h5ad");

        count_bam_into(
            &out,
            &testdata().join("test_spliced.bam"),
            10_000,
            &barcodes,
            &AdjustRead::default(),
        )
        .unwrap();

        let per_cell: f64 = entries(&out)
            .iter()
            .filter(|((_, bin), _)| bin == "5:65970000-65980000")
            .map(|(_, v)| *v)
            .sum();
        assert_eq!(per_cell, 1.0);
    }

    #[test]
    fn extending_a_spliced_read_covers_its_intron_again() {
        // Extension replaces the alignment with one interval, so the blocks go
        // with it. The reference implementation does the same, and warns that
        // --extendReads is not for spliced data.
        let barcodes = test_barcodes();
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("extended.h5ad");

        count_bam_into(
            &out,
            &testdata().join("test_spliced.bam"),
            10_000,
            &barcodes,
            &AdjustRead {
                extend_reads: Some(200),
                center_reads: false,
                ..AdjustRead::default()
            },
        )
        .unwrap();

        let sums = column_sums(&out);
        assert_eq!(
            sums.get("5:65970000-65980000").copied(),
            Some(3.0),
            "both spliced reads should now span the middle bin {sums:?}"
        );
    }

    /// A BED holding one region, written where the test wants it.
    fn blacklist_bed(dir: &Path, name: &str, start: usize, end: usize) -> PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, format!("5\t{start}\t{end}\tregion\t0\t.\n")).unwrap();
        path
    }

    fn count_with_blacklist(out: &Path, blacklist: Option<PathBuf>, barcodes: &[String]) {
        count_into(
            out,
            100_000,
            100_000,
            barcodes,
            &CountingParams {
                blacklist_path: blacklist,
                ..CountingParams::default()
            },
            1_000_000,
        )
        .unwrap();
    }

    #[test]
    fn a_blacklist_removes_reads_and_never_moves_them() {
        let barcodes = test_barcodes();
        let dir = TempDir::new().unwrap();

        let base = dir.path().join("base.h5ad");
        count_with_blacklist(&base, None, &barcodes);

        // Both BED regions cover the whole locus this BAM was cut from, so
        // every read is at least half blacklisted.
        let filtered = dir.path().join("filtered.h5ad");
        count_with_blacklist(
            &filtered,
            Some(testdata().join("Chrna9_regions.bed")),
            &barcodes,
        );

        // The bins themselves survive: a blacklist filters reads, it does not
        // delete columns.
        assert_eq!(n_vars(&filtered), n_vars(&base));

        let base = entries(&base);
        let filtered = entries(&filtered);
        assert!(!base.is_empty(), "the unfiltered run counted nothing");

        // No read is moved: an entry can only shrink, and no entry appears that
        // the unfiltered run did not have.
        for (key, value) in &filtered {
            let before = base.get(key).copied().unwrap_or(0.0);
            assert!(
                *value <= before,
                "{key:?} grew from {before} to {value} under a blacklist"
            );
        }
        let total: f64 = filtered.values().sum();
        assert!(total < base.values().sum::<f64>(), "no read was removed");
    }

    #[test]
    fn a_blacklist_region_that_only_clips_reads_removes_none() {
        let barcodes = test_barcodes();
        let dir = TempDir::new().unwrap();

        let base = dir.path().join("base.h5ad");
        count_with_blacklist(&base, None, &barcodes);

        // 20 bp inside the leftmost reads, which are 61 bp long: under the
        // threshold, so nothing is dropped. Excluding whole bins instead would
        // empty the 100 kb bin that holds them.
        let bed = blacklist_bed(dir.path(), "narrow.bed", 65_966_820, 65_966_840);
        let clipped = dir.path().join("clipped.h5ad");
        count_with_blacklist(&clipped, Some(bed), &barcodes);

        assert_eq!(entries(&clipped), entries(&base));
    }

    #[test]
    fn restricting_to_a_region_tiles_only_that_region() {
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

        let region = dir.path().join("region.h5ad");
        count_into(
            &region,
            100_000,
            100_000,
            &barcodes,
            &CountingParams {
                region: Some("5:0-100000".to_string()),
                ..CountingParams::default()
            },
            1_000_000,
        )
        .unwrap();

        let open = |p: &Path| AnnData::<H5>::open(H5::open(p).unwrap()).unwrap();
        let n_whole = open(&whole).n_vars();
        let n_region = open(&region).n_vars();

        // One 100 kb bin for the window, against the whole chromosome before.
        assert_eq!(n_region, 1, "expected a single bin for a 100 kb window");
        assert!(n_region < n_whole, "{n_region} is not fewer than {n_whole}");
    }

    #[test]
    fn a_region_that_does_not_start_at_zero_tiles_from_its_own_start() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("offset.h5ad");

        // 65,900,000..66,000,000 in 10 kb bins: ten bins, the first starting at
        // the region start rather than at the chromosome start.
        count_into(
            &out,
            10_000,
            10_000,
            &test_barcodes(),
            &CountingParams {
                region: Some("5:65900000-66000000".to_string()),
                ..CountingParams::default()
            },
            1_000_000,
        )
        .unwrap();

        let adata = AnnData::<H5>::open(H5::open(&out).unwrap()).unwrap();
        let names = adata.var_names().into_vec();

        assert_eq!(names.len(), 10, "{names:?}");
        assert_eq!(names[0], "5:65900000-65910000");
        assert_eq!(names[9], "5:65990000-66000000");
    }

    #[test]
    fn a_read_lands_in_the_same_genomic_bin_with_or_without_a_region() {
        // The window shifts every bin index; if the offset arithmetic were
        // wrong, reads would be credited to a different genomic interval.
        let barcodes = test_barcodes();
        let dir = TempDir::new().unwrap();

        let whole = dir.path().join("whole.h5ad");
        count_into(
            &whole,
            10_000,
            10_000,
            &barcodes,
            &CountingParams::default(),
            1_000_000,
        )
        .unwrap();

        let region = dir.path().join("region.h5ad");
        count_into(
            &region,
            10_000,
            10_000,
            &barcodes,
            &CountingParams {
                region: Some("5:65900000-66000000".to_string()),
                ..CountingParams::default()
            },
            1_000_000,
        )
        .unwrap();

        let whole_sums = column_sums(&whole);
        let region_sums = column_sums(&region);

        assert!(!region_sums.is_empty(), "the region counted nothing");
        for (bin_name, region_count) in &region_sums {
            let whole_count = whole_sums.get(bin_name).copied().unwrap_or(0.0);
            assert_eq!(
                *region_count, whole_count,
                "bin {bin_name} counted {region_count} inside the region but {whole_count} without it"
            );
        }
    }

    #[test]
    fn a_region_narrower_than_one_bin_still_produces_a_bin() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("narrow.h5ad");

        count_into(
            &out,
            10_000,
            10_000,
            &test_barcodes(),
            &CountingParams {
                region: Some("5:65970000-65971000".to_string()),
                ..CountingParams::default()
            },
            1_000_000,
        )
        .unwrap();

        let adata = AnnData::<H5>::open(H5::open(&out).unwrap()).unwrap();
        let names = adata.var_names().into_vec();
        // The last bin is clamped to the window end, not to the bin size.
        assert_eq!(names, vec!["5:65970000-65971000"]);
    }

    #[test]
    fn a_region_start_tiles_from_the_window_start() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("offgrid.h5ad");

        count_into(
            &out,
            10_000,
            10_000,
            &test_barcodes(),
            &CountingParams {
                region: Some("5:65905000-65965000".to_string()),
                ..CountingParams::default()
            },
            1_000_000,
        )
        .unwrap();

        let adata = AnnData::<H5>::open(H5::open(&out).unwrap()).unwrap();
        let names = adata.var_names().into_vec();

        // 65,905,000..65,965,000 is 60 kb, so exactly six 10 kb bins, the
        // first starting at the requested start.
        assert_eq!(
            names,
            vec![
                "5:65905000-65915000",
                "5:65915000-65925000",
                "5:65925000-65935000",
                "5:65935000-65945000",
                "5:65945000-65955000",
                "5:65955000-65965000",
            ]
        );
    }

    #[test]
    fn a_region_length_that_is_not_a_whole_number_of_bins_ends_short() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("shorttail.h5ad");

        count_into(
            &out,
            10_000,
            10_000,
            &test_barcodes(),
            &CountingParams {
                region: Some("5:65900000-65925000".to_string()),
                ..CountingParams::default()
            },
            1_000_000,
        )
        .unwrap();

        let adata = AnnData::<H5>::open(H5::open(&out).unwrap()).unwrap();
        let names = adata.var_names().into_vec();

        assert_eq!(
            names,
            vec![
                "5:65900000-65910000",
                "5:65910000-65920000",
                "5:65920000-65925000",
            ]
        );
    }

    #[test]
    fn a_bare_chromosome_region_still_tiles_the_whole_chromosome() {
        // parse_region leaves an open end as usize::MAX; it must be clamped to
        // the chromosome length rather than overflowing the bin count.
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

        let bare = dir.path().join("bare.h5ad");
        count_into(
            &bare,
            100_000,
            100_000,
            &barcodes,
            &CountingParams {
                region: Some("5".to_string()),
                ..CountingParams::default()
            },
            1_000_000,
        )
        .unwrap();

        let open = |p: &Path| AnnData::<H5>::open(H5::open(p).unwrap()).unwrap();
        assert_eq!(open(&bare).n_vars(), open(&whole).n_vars());
        assert_eq!(column_sums(&bare), column_sums(&whole));
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

    // --groupTag: recovering per-sample identity from a merged BAM

    fn merged_bam() -> PathBuf {
        testdata().join("test_i1_i2.bam")
    }

    #[allow(clippy::too_many_arguments)]
    fn count_grouped(out: &Path, bam: &Path, group_tag: Option<&str>) -> Result<()> {
        count_bam_bins(
            &[(bam, "merged")],
            100_000,
            100_000,
            &test_barcodes(),
            "BC",
            None,
            None,
            group_tag,
            &CountingParams::default(),
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
            1_000_000,
        )
    }

    #[test]
    fn a_merged_bam_counted_by_group_matches_counting_its_sources_separately() {
        // The whole point of --groupTag: one merged file must give the same
        // matrix as the two files it was made from.
        let dir = TempDir::new().unwrap();

        let separate = dir.path().join("separate.h5ad");
        let i1 = testdata().join("test_i1.bam");
        let i2 = testdata().join("test_i2.bam");
        count_bam_bins(
            &[(i1.as_path(), "test_i1"), (i2.as_path(), "test_i2")],
            100_000,
            100_000,
            &test_barcodes(),
            "BC",
            None,
            None,
            None,
            &CountingParams::default(),
            &AdjustRead::default(),
            None,
            None,
            None,
            None,
            None,
            &separate,
            "none",
            0,
            1,
            1_000_000,
        )
        .unwrap();

        let grouped = dir.path().join("grouped.h5ad");
        count_grouped(&grouped, &merged_bam(), Some("RG")).unwrap();

        let open = |p: &Path| AnnData::<H5>::open(H5::open(p).unwrap()).unwrap();
        assert_eq!(
            open(&grouped).obs_names().into_vec(),
            open(&separate).obs_names().into_vec()
        );

        let total = |p: &Path| -> f64 {
            super::super::count_utils::read_x_f64(&open(p))
                .unwrap()
                .triplet_iter()
                .map(|(_, _, &v)| v)
                .sum()
        };
        assert_eq!(total(&grouped), total(&separate));
    }

    #[test]
    fn without_a_group_tag_a_shared_barcode_collapses_two_cells_into_one() {
        // GCGAGCAT occurs in both source samples. Counted per-barcode the two
        // cells merge; counted per-group they stay apart. This is the failure
        // --groupTag exists to prevent.
        let dir = TempDir::new().unwrap();

        let flat = dir.path().join("flat.h5ad");
        count_grouped(&flat, &merged_bam(), None).unwrap();
        let grouped = dir.path().join("grouped.h5ad");
        count_grouped(&grouped, &merged_bam(), Some("RG")).unwrap();

        let open = |p: &Path| AnnData::<H5>::open(H5::open(p).unwrap()).unwrap();
        let rows = |p: &Path| open(p).obs_names().into_vec();

        let flat_rows = rows(&flat);
        let grouped_rows = rows(&grouped);
        assert_eq!(flat_rows.len(), test_barcodes().len());
        assert_eq!(grouped_rows.len(), 2 * test_barcodes().len());

        // One row for the shared barcode without grouping, two with it.
        let shared = |rs: &[String]| rs.iter().filter(|r| r.ends_with("GCGAGCAT")).count();
        assert_eq!(shared(&flat_rows), 1);
        assert_eq!(shared(&grouped_rows), 2);
    }

    #[test]
    fn a_group_tag_with_several_bams_is_rejected() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("unused.h5ad");
        let i1 = testdata().join("test_i1.bam");
        let i2 = testdata().join("test_i2.bam");

        let err = count_bam_bins(
            &[(i1.as_path(), "a"), (i2.as_path(), "b")],
            100_000,
            100_000,
            &test_barcodes(),
            "BC",
            None,
            None,
            Some("RG"),
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

        assert!(err.contains("single merged BAM"), "{err}");
    }

    #[test]
    fn a_group_tag_on_a_bam_without_read_groups_is_rejected() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("unused.h5ad");
        let err = count_grouped(&out, &testdata().join("test_i1.bam"), Some("RG"))
            .unwrap_err()
            .to_string();
        assert!(err.contains("no @RG"), "{err}");
    }
}
