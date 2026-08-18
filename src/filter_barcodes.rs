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
use crate::annotation::region_index::{ChromIndex, GenomeIndex};
use crate::bam::bam_io::{BamWorker, ensure_barcode_tags_present, read_bam_header};

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

                        // A read counts once, in the bin holding its start.
                        let bin_idx = start / bin_size;
                        let bin_key = pack_bin(*chrom_idx, bin_idx);
                        // Look up by the borrowed bytes first and only copy the
                        // barcode the first time this chunk sees it.
                        match local_bins.get_mut(barcode) {
                            Some(bins) => {
                                bins.insert(bin_key);
                            }
                            None => {
                                let mut bins = AHashSet::new();
                                bins.insert(bin_key);
                                local_bins.insert(barcode.to_vec(), bins);
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
struct FuzzyWhitelist<'a> {
    max_dist: usize,
    n_blocks: usize,
    entries: &'a [String],
    /// `(length, block position) -> block bytes -> entry indices`.
    ///
    /// Nested rather than keyed by one `(usize, usize, &[u8])` tuple so a
    /// lookup can borrow the block straight out of the barcode.
    buckets: AHashMap<(usize, usize), AHashMap<&'a [u8], Vec<u32>>>,
}

impl<'a> FuzzyWhitelist<'a> {
    fn build(entries: &'a [String], max_dist: usize) -> Self {
        let n_blocks = max_dist + 1;
        let mut buckets: AHashMap<(usize, usize), AHashMap<&'a [u8], Vec<u32>>> = AHashMap::new();

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

fn is_blacklisted(
    blacklist_index: &GenomeIndex,
    chromosome: &str,
    start: usize,
    end: usize,
) -> bool {
    blacklist_chrom_index(blacklist_index, chromosome)
        .map(|idx| idx.find(start, end).next().is_some())
        .unwrap_or(false)
}

fn blacklist_chrom_index<'a>(
    blacklist_index: &'a GenomeIndex,
    chromosome: &str,
) -> Option<&'a ChromIndex> {
    blacklist_index
        .get(chromosome)
        .or_else(|| {
            chromosome
                .strip_prefix("chr")
                .and_then(|c| blacklist_index.get(c))
        })
        .or_else(|| blacklist_index.get(&format!("chr{chromosome}")))
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
        let bases = [b'A', b'C', b'G', b'T'];
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
            let alphabet = [b'A', b'C', b'G', b'T', b'N'];
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
}
