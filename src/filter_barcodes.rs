use std::path::{Path, PathBuf};

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result};
use noodles::bam;
use noodles::sam::alignment::Record as AlignmentRecord;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rayon::prelude::*;
use triple_accel::hamming::hamming;

use crate::annotation::parse_annotation::parse_blacklist_bed;
use crate::annotation::region_index::GenomeIndex;
use crate::bam::bam_io::{
    BamWorker, ensure_barcode_tags_present, read_bam_header, read_group_ids, warn_unknown_group,
};
use crate::bam::filters::is_blacklisted;

/// A map of barcodes stored as bytes in a `Vec<u8>` (directly read from the BAM
/// record) to the bins is was detected in, stored as their index.
/// Bin indices are packed as `chromosome index : bin index` in a `u64`.`
type BinsByBarcode = AHashMap<Vec<u8>, AHashSet<u64>>;

/// Get a bin index from its chromosome index and bin position.
#[inline]
fn pack_bin(chrom_idx: usize, bin_pos: usize) -> u64 {
    ((chrom_idx as u64) << 32) | bin_pos as u64
}

fn run_filter_barcodes(
    bamfile: &Path,
    whitelist: Option<Vec<String>>,
    blacklist_file_name: Option<&Path>,
    cell_tag: &str,
    group_tag: Option<&str>,
    min_hamming_dist: usize,
    min_mapping_quality: Option<u8>,
    bin_size: usize,
    chr_to_skip: &[String],
    num_threads: usize,
    chunk_size: usize,
) -> Result<Vec<(String, usize)>> {
    if bin_size == 0 {
        anyhow::bail!("bin_size must be greater than zero");
    }
    anyhow::ensure!(chunk_size > 0, "chunk_size must be greater than zero");

    let whitelist = whitelist.unwrap_or_default();
    let whitelist_is_active = !whitelist.is_empty();
    let whitelist_matcher = WhitelistMatcher::build(&whitelist, min_hamming_dist);
    let blacklist_index = if let Some(p) = blacklist_file_name {
        parse_blacklist_bed(p)?
    } else {
        GenomeIndex::new()
    };

    let tag = parse_tag(cell_tag)?;
    ensure_barcode_tags_present(&[bamfile], tag, None)?;

    let header = read_bam_header(bamfile)?;

    // With --groupTag a barcode alone no longer names a cell, so the reported
    // unit becomes `group::barcode`. The valid groups are the BAM's own @RG IDs.
    let group_tag_parsed = group_tag.map(parse_tag).transpose()?;
    let group_ids: Option<Vec<Vec<u8>>> = match group_tag {
        Some(_) => Some(read_group_ids(&header, bamfile)?),
        None => None,
    };
    let known_groups: AHashSet<&[u8]> = group_ids.iter().flatten().map(Vec::as_slice).collect();

    let chrom_sizes: Vec<(String, usize)> = header
        .reference_sequences()
        .iter()
        .filter(|(name, _)| !chr_to_skip.contains(&name.to_string()))
        .map(|(name, seq)| (name.to_string(), seq.length().get()))
        .collect();

    // Build chunk work list sorted by descending size. Each chunk carries its
    // chromosome's index so the per-read bin key needs no chromosome name.
    let mut chunks: Vec<(usize, String, usize, usize)> = chrom_sizes
        .iter()
        .enumerate()
        .flat_map(|(chrom_idx, (chrom, chrom_len))| {
            (0..*chrom_len).step_by(chunk_size).map(move |start| {
                (
                    chrom_idx,
                    chrom.clone(),
                    start,
                    (start + chunk_size).min(*chrom_len),
                )
            })
        })
        .collect();
    chunks.sort_unstable_by_key(|b| std::cmp::Reverse(b.3 - b.2));

    let n_threads = if num_threads == 0 {
        rayon::current_num_threads()
    } else {
        num_threads
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .context("failed to build thread pool")?;

    let partial_maps: Vec<BinsByBarcode> = pool.install(|| {
        chunks
            .par_iter()
            .map_init(
                BamWorker::new,
                |worker, (chrom_idx, chrom, chunk_start, chunk_end)| -> Result<BinsByBarcode> {
                    // One reader per rayon thread rather than per chunk.
                    let (reader, header, _motif) = worker.prepare(bamfile, None)?;

                    let region_str = format!("{}:{}-{}", chrom, chunk_start + 1, chunk_end);
                    let region: noodles::core::Region = region_str
                        .parse()
                        .with_context(|| format!("failed to parse region: {}", region_str))?;

                    let query = match reader.query(header, &region) {
                        Ok(q) => q,
                        Err(_) => return Ok(AHashMap::new()),
                    };

                    let mut local_bins: BinsByBarcode = AHashMap::new();
                    // Reused across records so building the composite key
                    // does not allocate per read.
                    let mut composite: Vec<u8> = Vec::new();

                    for result in query.records() {
                        let record = result.context("failed to read BAM record")?;

                        let flags = record.flags();
                        if flags.is_unmapped() {
                            continue;
                        }

                        if let Some(min_mq) = min_mapping_quality {
                            match record.mapping_quality() {
                                Some(mq) if mq.get() >= min_mq => {}
                                _ => continue,
                            }
                        }

                        let Some(aln_start) = record
                            .alignment_start()
                            .transpose()
                            .context("failed to read alignment start")?
                        else {
                            continue;
                        };
                        let Some(aln_end) = record
                            .alignment_end()
                            .transpose()
                            .context("failed to read alignment end")?
                        else {
                            continue;
                        };

                        let start = aln_start.get().saturating_sub(1);
                        // Ownership: a read belongs to the chunk that contains its
                        // alignment_start. The BAI query returns overlapping reads,
                        // so skip anything that started before this chunk.
                        if start < *chunk_start {
                            continue;
                        }

                        let end: usize = aln_end.get();
                        if end <= start {
                            continue;
                        }

                        if is_blacklisted(&blacklist_index, chrom, start, end) {
                            continue;
                        }

                        let Some(barcode) = read_cell_barcode(&record, &tag)? else {
                            continue;
                        };

                        if whitelist_is_active && !whitelist_matcher.matches(barcode) {
                            continue;
                        }

                        // The counted unit: the barcode alone, or `group::barcode`
                        // when the reads carry their sample of origin.
                        let key: &[u8] = match &group_tag_parsed {
                            Some(gtag) => {
                                let Some(group) = read_cell_barcode(&record, gtag)? else {
                                    continue;
                                };
                                if !known_groups.contains(group) {
                                    warn_unknown_group(group);
                                    continue;
                                }
                                composite.clear();
                                composite.extend_from_slice(group);
                                composite.extend_from_slice(b"::");
                                composite.extend_from_slice(barcode);
                                &composite
                            }
                            None => barcode,
                        };

                        // A read counts once, in the bin holding its start.
                        let bin_idx = start / bin_size;
                        let bin_key = pack_bin(*chrom_idx, bin_idx);
                        // Look up by the borrowed bytes first and only copy the
                        // key the first time this chunk sees it.
                        match local_bins.get_mut(key) {
                            Some(bins) => {
                                bins.insert(bin_key);
                            }
                            None => {
                                let mut bins = AHashSet::new();
                                bins.insert(bin_key);
                                local_bins.insert(key.to_vec(), bins);
                            }
                        }
                    }

                    Ok(local_bins)
                },
            )
            .collect::<Result<Vec<_>>>()
    })?;

    let mut bins_by_barcode: BinsByBarcode = AHashMap::new();
    for partial in partial_maps {
        for (barcode, bins) in partial {
            bins_by_barcode.entry(barcode).or_default().extend(bins);
        }
    }

    let mut barcode_counts: Vec<(String, usize)> = bins_by_barcode
        .into_iter()
        .map(|(barcode, bins)| (String::from_utf8_lossy(&barcode).into_owned(), bins.len()))
        .collect();

    barcode_counts.sort_by(|l, r| r.1.cmp(&l.1).then_with(|| l.0.cmp(&r.0)));

    Ok(barcode_counts)
}

fn parse_tag(cell_tag: &str) -> Result<Tag> {
    let bytes = cell_tag.as_bytes();
    if bytes.len() != 2 {
        anyhow::bail!("barcode tag must be exactly two characters");
    }
    Ok(Tag::new(bytes[0], bytes[1]))
}

/// Borrow the barcode tag's bytes, tied to the record's lifetime, so the hot
/// path never allocates. Returns `None` if the tag is absent or not a string.
fn read_cell_barcode<'a>(record: &'a bam::Record, tag: &Tag) -> Result<Option<&'a [u8]>> {
    match record.data().get(tag) {
        Some(Ok(Value::String(value))) => Ok(Some(value.as_ref())),
        Some(Ok(_)) => Ok(None),
        Some(Err(error)) => Err(error.into()),
        None => Ok(None),
    }
}

/// Half-open bounds of block `b` when a string of `len` bytes is cut into
/// `n_blocks` near-equal pieces.
fn block_bounds(len: usize, n_blocks: usize, b: usize) -> (usize, usize) {
    (len * b / n_blocks, len * (b + 1) / n_blocks)
}

/// Decides whether a read's barcode matches the whitelist.
///
/// Exact matching is a single hash lookup. Fuzzy matching uses a pigeonhole
/// index to match the barcode to the whitelist.
enum WhitelistMatcher<'a> {
    Exact(AHashSet<&'a [u8]>),
    Fuzzy(FuzzyWhitelist<'a>),
}

/// Pigeonhole ("partition") index over the whitelist.
///
/// Two equal-length strings at Hamming distance `<= d` must agree on at least
/// one of `d + 1` disjoint blocks, because `d` substitutions cannot touch every
/// block. Indexing each entry under all of its blocks yields a candidate set
/// that contains every true match, and the exact `hamming` check then runs on
/// those candidates alone.
///
/// Entries are bucketed by length as well, since only equal-length entries are
/// comparable and the block boundaries depend on the length.
/// `(length, block position) -> block bytes -> entry indices`.
///
/// Nested rather than keyed by one `(usize, usize, &[u8])` tuple so a lookup
/// can borrow the block straight out of the barcode.
type BlockBuckets<'a> = AHashMap<(usize, usize), AHashMap<&'a [u8], Vec<u32>>>;

struct FuzzyWhitelist<'a> {
    max_dist: usize,
    n_blocks: usize,
    entries: &'a [String],
    buckets: BlockBuckets<'a>,
}

impl<'a> FuzzyWhitelist<'a> {
    fn build(entries: &'a [String], max_dist: usize) -> Self {
        let n_blocks = max_dist + 1;
        let mut buckets: BlockBuckets<'a> = AHashMap::new();

        for (idx, entry) in entries.iter().enumerate() {
            let bytes = entry.as_bytes();
            for b in 0..n_blocks {
                let (start, end) = block_bounds(bytes.len(), n_blocks, b);
                buckets
                    .entry((bytes.len(), b))
                    .or_default()
                    .entry(&bytes[start..end])
                    .or_default()
                    .push(idx as u32);
            }
        }

        Self {
            max_dist,
            n_blocks,
            entries,
            buckets,
        }
    }

    fn matches(&self, barcode: &[u8]) -> bool {
        let len = barcode.len();
        for b in 0..self.n_blocks {
            let (start, end) = block_bounds(len, self.n_blocks, b);
            let Some(by_block) = self.buckets.get(&(len, b)) else {
                continue;
            };
            let Some(candidates) = by_block.get(&barcode[start..end]) else {
                continue;
            };
            for &idx in candidates {
                // Same length by construction, which `hamming` requires.
                if hamming(barcode, self.entries[idx as usize].as_bytes()) as usize <= self.max_dist
                {
                    return true;
                }
            }
        }
        false
    }
}

impl<'a> WhitelistMatcher<'a> {
    fn build(whitelist: &'a [String], min_hamming_dist: usize) -> Self {
        if min_hamming_dist == 0 {
            // Keyed by bytes so the per-read lookup needs no UTF-8 validation either.
            Self::Exact(whitelist.iter().map(|bc| bc.as_bytes()).collect())
        } else {
            Self::Fuzzy(FuzzyWhitelist::build(whitelist, min_hamming_dist))
        }
    }

    fn matches(&self, barcode: &[u8]) -> bool {
        match self {
            Self::Exact(set) => set.contains(barcode),
            Self::Fuzzy(index) => index.matches(barcode),
        }
    }
}

/// Detect the cell barcodes in a BAM file and count the bins each one occupies.
///
/// Every mapped read carrying a ``cell_tag`` tag is placed in the ``bin_size``
/// bin holding its alignment start. A barcode's count is the number of
/// *distinct* bins it reaches, not its number of reads: a real cell spreads
/// signal over many bins, while a PCR pile-up puts many reads in few bins. This
/// is the statistic behind the barcode rank ("knee") plot.
///
/// Parameters
/// ----------
/// bamfile
///     Path to a coordinate-sorted BAM file. A BAI index must sit beside it.
/// whitelist
///     Barcodes to keep. ``None`` keeps every barcode found.
/// blacklist_file_name
///     BED file of regions to exclude; overlapping reads are ignored.
/// cell_tag
///     Two-character BAM auxiliary tag holding the cell barcode.
/// min_hamming_dist
///     Substitutions allowed when matching ``whitelist``. ``0`` demands an exact
///     match. Higher values cost noticeably more.
/// min_mapping_quality
///     Drop reads below this MAPQ. ``None`` applies no threshold.
/// bin_size
///     Width in bp of the bins a barcode is counted across.
/// chr_to_skip
///     Chromosome names to exclude entirely.
/// num_threads
///     Worker threads. ``0`` uses every available core.
/// chunk_size
///     Width in bp of one unit of parallel work.
///
/// Returns
/// -------
/// list of (str, int)
///     One ``(barcode, count)`` pair per detected barcode, sorted by descending
///     count and then by barcode. Every barcode seen is returned: applying a
///     minimum count and assigning ranks is left to the caller.
///
/// Raises
/// ------
/// RuntimeError
///     If the BAM or its index cannot be read, if ``cell_tag`` is not two
///     characters, or if ``bin_size`` or ``chunk_size`` is zero.
#[pyfunction(signature = (
    bamfile,
    whitelist = None,
    blacklist_file_name = None,
    cell_tag = "CB",
    group_tag = None,
    min_hamming_dist = 0,
    min_mapping_quality = None,
    bin_size = 100_000,
    chr_to_skip = vec![],
    num_threads = 0,
    chunk_size = 1_000_000,
))]
pub fn filter_barcodes(
    bamfile: PathBuf,
    whitelist: Option<Vec<String>>,
    blacklist_file_name: Option<PathBuf>,
    cell_tag: &str,
    group_tag: Option<String>,
    min_hamming_dist: usize,
    min_mapping_quality: Option<u8>,
    bin_size: usize,
    chr_to_skip: Vec<String>,
    num_threads: usize,
    chunk_size: usize,
) -> PyResult<Vec<(String, usize)>> {
    run_filter_barcodes(
        bamfile.as_path(),
        whitelist,
        blacklist_file_name.as_deref(),
        cell_tag,
        group_tag.as_deref(),
        min_hamming_dist,
        min_mapping_quality,
        bin_size,
        &chr_to_skip,
        num_threads,
        chunk_size,
    )
    .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wl(entries: &[&str]) -> Vec<String> {
        entries.iter().map(|s| s.to_string()).collect()
    }

    fn brute_force_match(whitelist: &[String], barcode: &[u8], max_dist: usize) -> bool {
        whitelist.iter().any(|entry| {
            entry.len() == barcode.len() && hamming(barcode, entry.as_bytes()) as usize <= max_dist
        })
    }

    #[test]
    fn exact_matcher_accepts_only_whole_barcodes() {
        let whitelist = wl(&["AAAACCCC", "GGGGTTTT"]);
        let matcher = WhitelistMatcher::build(&whitelist, 0);

        assert!(matcher.matches(b"AAAACCCC"));
        assert!(matcher.matches(b"GGGGTTTT"));
        // One substitution is not an exact match.
        assert!(!matcher.matches(b"AAAACCCG"));
        assert!(!matcher.matches(b"AAAACCC"));
    }

    #[test]
    fn fuzzy_matcher_accepts_within_distance_and_rejects_beyond_it() {
        let whitelist = wl(&["AAAACCCC", "GGGGTTTT"]);
        let matcher = WhitelistMatcher::build(&whitelist, 1);

        assert!(matcher.matches(b"AAAACCCC")); // distance 0
        assert!(matcher.matches(b"TAAACCCC")); // distance 1, first block
        assert!(matcher.matches(b"AAAACCCT")); // distance 1, second block
        assert!(!matcher.matches(b"TTAACCCC")); // distance 2
    }

    #[test]
    fn fuzzy_matcher_only_compares_equal_lengths() {
        // `hamming` requires equal lengths, so entries of another length must
        // never reach it.
        let whitelist = wl(&["AAAA", "AAAAAAAA"]);
        let matcher = WhitelistMatcher::build(&whitelist, 1);

        assert!(matcher.matches(b"AAAT"));
        assert!(matcher.matches(b"AAAAAAAT"));
        assert!(!matcher.matches(b"AAAAAA"));
    }

    #[test]
    fn fuzzy_matcher_agrees_with_whole_whitelist_scan() {
        // The index is an optimization, so its answers must be identical to the
        // scan it replaced for every query, not merely similar.
        let bases = *b"ACGT";
        let whitelist: Vec<String> = (0..256u32)
            .map(|i| {
                (0..4)
                    .map(|p| bases[((i >> (2 * p)) & 0b11) as usize] as char)
                    .collect()
            })
            .collect();

        for max_dist in 1..=2 {
            let matcher = WhitelistMatcher::build(&whitelist, max_dist);
            // Every 4-mer over a 5-letter alphabet, so queries include the
            // off-alphabet 'N' the index must still handle.
            let alphabet = *b"ACGTN";
            for a in alphabet {
                for b in alphabet {
                    for c in alphabet {
                        for d in alphabet {
                            let query = [a, b, c, d];
                            assert_eq!(
                                matcher.matches(&query),
                                brute_force_match(&whitelist, &query, max_dist),
                                "disagreement at {:?} with max_dist {}",
                                std::str::from_utf8(&query).unwrap(),
                                max_dist
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn pack_bin_is_unique_per_chromosome_and_bin() {
        assert_eq!(pack_bin(0, 0), 0);
        assert_ne!(pack_bin(0, 1), pack_bin(1, 0));
        assert_ne!(pack_bin(1, 0), pack_bin(0, u32::MAX as usize));
        // A bin index at the edge of its 32-bit half must not bleed into the
        // chromosome half.
        assert_eq!(pack_bin(2, u32::MAX as usize) >> 32, 2);
    }

    // Tag parsing and barcode extraction

    #[test]
    fn a_barcode_tag_must_be_exactly_two_characters() {
        assert_eq!(parse_tag("BC").unwrap(), Tag::new(b'B', b'C'));
        assert!(parse_tag("B").is_err());
        assert!(parse_tag("BCX").is_err());
        assert!(parse_tag("").is_err());
    }

    fn testdata() -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/testdata")
    }

    fn first_record() -> (bam::Record, ()) {
        let path = testdata().join("test_i1.bam");
        let mut reader = bam::io::reader::Builder.build_from_path(&path).unwrap();
        reader.read_header().unwrap();
        let record = reader.records().next().unwrap().unwrap();
        (record, ())
    }

    #[test]
    fn the_barcode_is_borrowed_from_the_record_when_the_tag_holds_a_string() {
        let (record, _) = first_record();
        let barcode = read_cell_barcode(&record, &Tag::new(b'B', b'C')).unwrap();
        assert_eq!(barcode, Some(b"ATATAACT".as_slice()));
    }

    #[test]
    fn a_tag_that_is_absent_or_not_a_string_yields_no_barcode() {
        let (record, _) = first_record();
        // NM is an integer, ZZ is not present at all.
        assert_eq!(
            read_cell_barcode(&record, &Tag::new(b'N', b'M')).unwrap(),
            None
        );
        assert_eq!(
            read_cell_barcode(&record, &Tag::new(b'Z', b'Z')).unwrap(),
            None
        );
    }

    // Block splitting

    #[test]
    fn blocks_tile_the_whole_string_without_gaps_or_overlap() {
        for (len, n_blocks) in [(8usize, 3usize), (10, 4), (7, 7), (5, 1)] {
            let mut previous_end = 0;
            for b in 0..n_blocks {
                let (start, end) = block_bounds(len, n_blocks, b);
                assert_eq!(start, previous_end, "gap or overlap at block {b}");
                assert!(end >= start);
                previous_end = end;
            }
            assert_eq!(previous_end, len, "blocks did not reach the end");
        }
    }

    #[test]
    fn blocks_of_an_evenly_divisible_string_are_equal() {
        assert_eq!(block_bounds(9, 3, 0), (0, 3));
        assert_eq!(block_bounds(9, 3, 1), (3, 6));
        assert_eq!(block_bounds(9, 3, 2), (6, 9));
    }

    #[test]
    fn a_string_shorter_than_its_block_count_leaves_empty_blocks() {
        // Two bytes over three blocks: one block has to be empty.
        let bounds: Vec<(usize, usize)> = (0..3).map(|b| block_bounds(2, 3, b)).collect();
        assert_eq!(bounds.last().unwrap().1, 2);
        assert!(bounds.iter().any(|(s, e)| s == e), "{bounds:?}");
    }

    // Whole run against the test BAM

    fn barcodes_of(bam: &str) -> Vec<String> {
        std::fs::read_to_string(testdata().join(bam))
            .unwrap()
            .lines()
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty())
            .collect()
    }

    fn run(
        whitelist: Option<Vec<String>>,
        min_hamming: usize,
        bin_size: usize,
    ) -> Result<Vec<(String, usize)>> {
        run_filter_barcodes(
            &testdata().join("test_i1.bam"),
            whitelist,
            None,
            "BC",
            None,
            min_hamming,
            None,
            bin_size,
            &[],
            1,
            1_000_000,
        )
    }

    #[test]
    fn every_detected_barcode_occupies_at_least_one_bin() {
        let counts = run(None, 0, 2_000).unwrap();

        assert!(!counts.is_empty(), "no barcodes were detected");
        for (barcode, bins) in &counts {
            assert!(*bins > 0, "{barcode} was reported with no bins");
            assert!(!barcode.is_empty());
        }
    }

    #[test]
    fn the_detected_barcodes_are_the_ones_the_file_carries() {
        let counts = run(None, 0, 2_000).unwrap();
        let found: AHashSet<&str> = counts.iter().map(|(bc, _)| bc.as_str()).collect();

        // These are the eight barcodes present in test_i1.bam.
        for expected in ["ATATAACT", "ACGGTAAT", "GTCAAGCA", "TAGACTTG"] {
            assert!(found.contains(expected), "missing {expected} in {found:?}");
        }
    }

    #[test]
    fn a_whitelist_restricts_the_output_to_its_own_barcodes() {
        let whitelist = vec!["ATATAACT".to_string()];
        let counts = run(Some(whitelist), 0, 2_000).unwrap();

        assert_eq!(counts.len(), 1);
        assert_eq!(counts[0].0, "ATATAACT");
    }

    #[test]
    fn a_whitelist_entry_the_file_never_uses_reports_nothing() {
        let counts = run(Some(vec!["TTTTTTTT".to_string()]), 0, 2_000).unwrap();
        assert!(counts.is_empty(), "{counts:?}");
    }

    #[test]
    fn fuzzy_matching_recovers_a_barcode_with_one_substitution() {
        // ATATAACT with its last base changed; a Hamming distance of 1 should
        // still match it, while exact matching should not.
        let near_miss = vec!["ATATAACG".to_string()];

        let exact = run(Some(near_miss.clone()), 0, 2_000).unwrap();
        assert!(exact.is_empty(), "exact matching should not reach it");

        let fuzzy = run(Some(near_miss), 1, 2_000).unwrap();
        assert_eq!(fuzzy.len(), 1, "fuzzy matching missed the neighbour");
    }

    #[test]
    fn a_smaller_bin_size_never_reports_fewer_bins() {
        // Bins only subdivide, so the same reads can only spread over more.
        let coarse = run(None, 0, 1_000_000).unwrap();
        let fine = run(None, 0, 1_000).unwrap();

        let lookup: AHashMap<&str, usize> =
            coarse.iter().map(|(bc, n)| (bc.as_str(), *n)).collect();
        for (barcode, fine_bins) in &fine {
            let coarse_bins = lookup.get(barcode.as_str()).copied().unwrap_or(0);
            assert!(
                *fine_bins >= coarse_bins,
                "{barcode}: {fine_bins} fine bins is fewer than {coarse_bins} coarse"
            );
        }
    }

    #[test]
    fn the_whole_barcode_file_can_be_used_as_a_whitelist() {
        let counts = run(Some(barcodes_of("test_barcodes.txt")), 0, 2_000).unwrap();
        // Eight of the fourteen listed barcodes occur in test_i1.bam.
        assert_eq!(counts.len(), 8, "{counts:?}");
    }

    // Argument validation

    #[test]
    fn a_zero_bin_size_is_rejected() {
        let err = run(None, 0, 0).unwrap_err().to_string();
        assert!(err.contains("bin_size"), "{err}");
    }

    #[test]
    fn a_zero_chunk_size_is_rejected() {
        let err = run_filter_barcodes(
            &testdata().join("test_i1.bam"),
            None,
            None,
            "BC",
            None,
            0,
            None,
            2_000,
            &[],
            1,
            0,
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("chunk_size"), "{err}");
    }
}
