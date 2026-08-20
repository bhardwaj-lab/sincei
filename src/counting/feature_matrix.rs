//! Counts reads into a cell × annotation-feature matrix.
//!
//! Regions come from a BED/GTF/GFF3 file. Work is parallel over sub-chromosome
//! chunks, each an independent BAI query, with reads owned by the chunk holding
//! their alignment start so none is counted twice.
//!
//! A read is credited to exactly one feature: the overlaps of a feature's
//! sub-intervals, metagene exons, or the pieces left by blacklist subtraction
//! are summed first, and the feature with the largest total wins. So a read
//! spanning two exons of one gene counts once for that gene.

use ahash::{AHashMap, AHashSet};
use std::path::Path;

use anyhow::{Context, Result};
use rayon::prelude::*;

use super::count_utils::{build_csr, write_counts_anndata};
use super::params::{CountingParams, parse_region};
use crate::annotation::parse_annotation::{
    build_counting_index, parse_annotation_files, parse_blacklist_bed,
};
use crate::bam::bam_io::{
    BamWorker, ensure_barcode_tags_present, read_bam_header, read_group_ids, warn_unknown_group,
};
use crate::bam::filters::{
    DupMethod, DuplicateFilter, QcFilter, RawRecordFilter, derive_record_opts,
};

// Maximum number of distinct regions a single read can overlap and still be
// assigned correctly. In practice reads are short and features rarely pile up
// this densely, so 16 is a safe ceiling.
const MAX_REGION_HITS: usize = 16;
use crate::bam::sc_record::{AdjustRead, ScRecord, parse_tag};

/// Count reads from one or more BAM files into a cell × feature matrix, then
/// write the result as an AnnData HDF5 file.
///
/// `bam_paths` is a list of `(path, sample_name)` pairs. Parallelism is over
/// sub-chromosome chunks of `chunk_size` bp; each chunk is an independent BAI
/// query. Reads are assigned to chunks by `alignment_start` to avoid
/// double-counting. A read's effective interval can still extend into any
/// feature regardless of chunk boundaries.
pub fn count_bam_features(
    bam_paths: &[(&Path, &str)],
    annotation_path: &Path,
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
    anyhow::ensure!(!bam_paths.is_empty(), "at least one BAM file is required");
    anyhow::ensure!(chunk_size > 0, "chunk_size must be greater than zero");

    let has_motif = genome_path.is_some() && motifs.is_some();
    let record_opts = derive_record_opts(qc_filter, has_motif, dup_method.is_some());
    let bc_tag_parsed = parse_tag(bc_tag)?;
    let umi_tag_parsed = umi_tag.map(parse_tag).transpose()?;
    let all_bams: Vec<&Path> = bam_paths.iter().map(|(p, _)| *p).collect();
    ensure_barcode_tags_present(&all_bams, bc_tag_parsed, umi_tag_parsed)?;
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
    // than from which BAM they were read out of, so exactly one input is
    // allowed and the row space is the header's @RG IDs x barcodes.
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

    // `Vec<String>` -> `&[&str]` for the parser, which borrows the type names.
    let feature_types: Option<Vec<&str>> = params
        .feature_type
        .as_ref()
        .map(|t| t.iter().map(String::as_str).collect());
    let exon_types: Option<Vec<&str>> = params
        .exon_type
        .as_ref()
        .map(|t| t.iter().map(String::as_str).collect());

    let (mut feature_index, var_meta) = parse_annotation_files(
        [annotation_path],
        feature_types.as_deref(),
        exon_types.as_deref(),
        params.name_attr.as_deref(),
        params.metagene,
    )?;
    let n_features = var_meta.len();
    feature_index.retain(|chrom, _| !params.chr_to_skip.contains(chrom));

    let blacklist = params
        .blacklist_path
        .as_deref()
        .map(parse_blacklist_bed)
        .transpose()?;
    let counting_index = build_counting_index(&feature_index, blacklist.as_ref());

    // Chromosomes to skip, as a byte-slice set for O(1) membership tests.
    let skip_set: AHashSet<&[u8]> = params.chr_to_skip.iter().map(|s| s.as_bytes()).collect();

    // Build chunk work list: (bam_idx, bam_path, chrom, chunk_start, chunk_end).
    // Only include chromosomes that appear in counting_index.
    // Sort by descending chunk size.
    let mut work: Vec<(usize, &Path, String, usize, usize)> = Vec::new();
    for (bam_idx, &(bam_path, _)) in bam_paths.iter().enumerate() {
        let hdr = read_bam_header(bam_path)?;
        for (name, seq) in hdr.reference_sequences().iter() {
            let name_bytes: &[u8] = name.as_ref();
            if skip_set.contains(name_bytes) {
                continue;
            }
            let chrom = name.to_string();
            if !counting_index.contains_key(&chrom) {
                continue;
            }
            let chrom_len = seq.length().get();
            let mut start = 0;
            while start < chrom_len {
                let end = (start + chunk_size).min(chrom_len);
                work.push((bam_idx, bam_path, chrom.clone(), start, end));
                start += chunk_size;
            }
        }
    }
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
    // tree reduction rather than a single-threaded merge.
    let global_acc: AHashMap<(usize, usize), u32> = pool.install(|| {
        work.par_iter()
            .map_init(
                BamWorker::new,
                |worker,
                 &(bam_idx, bam_path, ref chrom, chunk_start, chunk_end)|
                 -> Result<AHashMap<(usize, usize), u32>> {
                    let (reader, header, motif) = worker.prepare(bam_path, motif_ingredients)?;

                    // Each work chunk covers a single chromosome, so its feature
                    // index is looked up once here rather than per read.
                    let Some(chrom_index) = counting_index.get(chrom.as_str()) else {
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
                        // alignment_start. Reads that started before this chunk are
                        // handled by the previous work chunk.
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

                        // Largest-overlap-wins: each read contributes to exactly
                        // one region. Sub-intervals of the same region (from
                        // blacklist splitting or metagene exons) accumulate their
                        // overlaps before the comparison, so a read spanning two
                        // exons of gene A still counts once for gene A.
                        let mut hits = [(0usize, 0usize); MAX_REGION_HITS];
                        let mut n_hits = 0usize;
                        for sub in chrom_index.find(eff_start, eff_end) {
                            let overlap = eff_end.min(sub.end) - eff_start.max(sub.start);
                            let mut merged = false;
                            for entry in hits[..n_hits].iter_mut() {
                                if entry.0 == sub.var_idx {
                                    entry.1 += overlap;
                                    merged = true;
                                    break;
                                }
                            }
                            if !merged && n_hits < MAX_REGION_HITS {
                                hits[n_hits] = (sub.var_idx, overlap);
                                n_hits += 1;
                            }
                        }
                        if let Some(&(best_val, _)) =
                            hits[..n_hits].iter().max_by_key(|&&(_, ov)| ov)
                        {
                            *local_acc.entry((cell_idx, best_val)).or_insert(0) += sc_rec.count;
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

    /// `count_bam_features` against `Chrna9.gtf`, everything optional off.
    fn count_into(out: &Path, annotation: &Path, params: &CountingParams) -> Result<()> {
        let bam = testdata().join("test_i1.bam");
        count_bam_features(
            &[(bam.as_path(), "s1")],
            annotation,
            &test_barcodes(),
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
            1_000_000,
        )
    }

    fn open(path: &Path) -> AnnData<H5> {
        AnnData::<H5>::open(H5::open(path).unwrap()).unwrap()
    }

    #[test]
    fn counting_features_gives_one_row_per_cell_and_one_column_per_feature() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("features.h5ad");

        count_into(
            &out,
            &testdata().join("Chrna9.gtf"),
            &CountingParams::default(),
        )
        .unwrap();

        assert!(out.exists(), "no output file was written");
        let adata = open(&out);
        assert_eq!(adata.n_obs(), test_barcodes().len());
        assert!(
            adata.n_vars() > 0,
            "the GTF should yield at least one feature"
        );
    }

    #[test]
    fn a_bed_annotation_is_accepted_as_well_as_a_gtf() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("bed.h5ad");

        count_into(
            &out,
            &testdata().join("Chrna9_regions.bed"),
            &CountingParams::default(),
        )
        .unwrap();

        assert!(open(&out).n_vars() > 0);
    }

    #[test]
    fn metagene_mode_collapses_transcripts_onto_their_genes() {
        let dir = TempDir::new().unwrap();
        let gtf = testdata().join("Chrna9.gtf");

        let per_transcript = dir.path().join("transcripts.h5ad");
        count_into(&per_transcript, &gtf, &CountingParams::default()).unwrap();

        let metagene_params = CountingParams {
            metagene: true,
            ..CountingParams::default()
        };
        let per_gene = dir.path().join("genes.h5ad");
        count_into(&per_gene, &gtf, &metagene_params).unwrap();

        // One gene can have many transcripts, never the other way round.
        assert!(
            open(&per_gene).n_vars() <= open(&per_transcript).n_vars(),
            "metagene mode produced more columns than transcript mode"
        );
    }

    #[test]
    fn metagene_counting_is_exonic_so_it_never_exceeds_whole_transcript_counting() {
        let dir = TempDir::new().unwrap();
        let gtf = testdata().join("Chrna9.gtf");

        let per_gene = dir.path().join("genes.h5ad");
        count_into(
            &per_gene,
            &gtf,
            &CountingParams {
                metagene: true,
                ..CountingParams::default()
            },
        )
        .unwrap();

        let per_transcript = dir.path().join("transcripts.h5ad");
        count_into(&per_transcript, &gtf, &CountingParams::default()).unwrap();

        let total = |p: &Path| -> f64 {
            super::super::count_utils::read_x_f64(&open(p))
                .unwrap()
                .triplet_iter()
                .map(|(_, _, &v)| v)
                .sum()
        };

        let gene_total = total(&per_gene);
        let transcript_total = total(&per_transcript);

        assert!(gene_total > 0.0, "nothing was counted against the exons");
        // A transcript feature spans its introns too, so an intronic read is
        // counted there but has no exon to fall in under metagene grouping.
        assert!(
            gene_total <= transcript_total,
            "exonic total {gene_total} exceeded whole-transcript total {transcript_total}"
        );
    }

    // Argument validation

    #[test]
    fn a_missing_annotation_file_is_reported() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("unused.h5ad");
        let missing = dir.path().join("no_such.gtf");

        assert!(count_into(&out, &missing, &CountingParams::default()).is_err());
    }

    #[test]
    fn counting_with_no_bam_files_is_rejected() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("unused.h5ad");

        let err = count_bam_features(
            &[],
            &testdata().join("Chrna9.gtf"),
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
    fn a_zero_chunk_size_is_rejected() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("unused.h5ad");
        let bam = testdata().join("test_i1.bam");

        let err = count_bam_features(
            &[(bam.as_path(), "s1")],
            &testdata().join("Chrna9.gtf"),
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
            0,
        )
        .unwrap_err()
        .to_string();

        assert!(err.contains("chunk_size"), "{err}");
    }

    #[test]
    fn a_barcode_tag_the_bam_does_not_carry_is_rejected_before_any_counting() {
        let dir = TempDir::new().unwrap();
        let out = dir.path().join("unused.h5ad");
        let bam = testdata().join("test_i1.bam");

        let err = count_bam_features(
            &[(bam.as_path(), "s1")],
            &testdata().join("Chrna9.gtf"),
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
        assert!(!out.exists());
    }
}
