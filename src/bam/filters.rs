//! Filter BAM records to be counted.
//!
//! The filters are split by how much of a record they need, so the counting
//! loop can reject a read as cheaply as possible:
//! [`RawRecordFilter`] sees only a raw record's flags and mapping quality,
//! before any tag parsing, while [`QcFilter`], [`DuplicateFilter`] and
//! [`MotifFilter`] work on a parsed [`ScRecord`].
//!
//! `derive_record_opts` reports which optional fields (GC content, aligned
//! fraction, read sequence) the active filters actually need, so none of them
//! are computed for a run that does not require them.

use ahash::AHashSet;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use anyhow::Result;
use twobit::TwoBitFile;

use super::sc_record::{AdjustRead, ScRecord, ScRecordOptions};
use crate::annotation::region_index::{ChromIndex, GenomeIndex};

/// Cheap per-record filter evaluated directly on a raw BAM record's flags and
/// mapping quality, before any tag parsing or [`ScRecord`] construction.
///
/// Keeping this separate from [`QcFilter`] lets callers reject reads before
/// paying for barcode/UMI tag lookups or sequence/CIGAR processing.
pub struct RawRecordFilter {
    /// Minimum mapping quality (MAPQ).
    pub min_mapq: Option<u8>,
    /// Only keep reads for which `flags & include == include`.
    pub sam_flag_include: Option<u16>,
    /// Drop reads for which `flags & exclude != 0`.
    pub sam_flag_exclude: Option<u16>,
    /// Strand filter for RNA-seq reads (dUTP library protocol); `"forward"`
    /// (keeps minus-strand reads) or `"reverse"` (keeps plus-strand reads).
    pub filter_rna_strand: Option<String>,
}

impl RawRecordFilter {
    pub fn new() -> Self {
        Self {
            min_mapq: None,
            sam_flag_include: None,
            sam_flag_exclude: None,
            filter_rna_strand: None,
        }
    }

    /// Returns `true` if the record passes all active thresholds.
    #[inline]
    pub fn passes(&self, flags: u16, mapq: Option<u8>) -> bool {
        if let Some(min_q) = self.min_mapq {
            match mapq {
                Some(q) if q >= min_q => {}
                _ => return false,
            }
        }

        if let Some(include) = self.sam_flag_include
            && flags & include != include
        {
            return false;
        }

        if let Some(exclude) = self.sam_flag_exclude
            && flags & exclude != 0
        {
            return false;
        }

        if let Some(ref strand) = self.filter_rna_strand
            && rna_strand_filter(flags, strand)
        {
            return false;
        }

        true
    }
}

impl Default for RawRecordFilter {
    fn default() -> Self {
        Self::new()
    }
}

/// Per-record quality-control filter.
///
/// All active thresholds must pass; an unset threshold is always satisfied.
pub struct QcFilter {
    /// Minimum fragment length (uses |TLEN| for paired-end, alignment span for single-end).
    pub min_fragment_length: Option<usize>,
    /// Maximum fragment length.
    pub max_fragment_length: Option<usize>,
    /// Minimum GC content in `[0, 1]`. Requires `ScRecordOptions::compute_gc`.
    pub min_gc: Option<f32>,
    /// Maximum GC content in `[0, 1]`. Requires `ScRecordOptions::compute_gc`.
    pub max_gc: Option<f32>,
    /// Minimum fraction of read bases in aligned CIGAR ops (`M`, `=`, `X`).
    /// Requires `ScRecordOptions::compute_aligned_fraction`.
    pub min_aligned_fraction: Option<f32>,
}

impl QcFilter {
    pub fn new() -> Self {
        Self {
            min_fragment_length: None,
            max_fragment_length: None,
            min_gc: None,
            max_gc: None,
            min_aligned_fraction: None,
        }
    }

    /// Returns `true` if the record passes all active thresholds.
    pub fn passes(&self, rec: &ScRecord<'_>) -> bool {
        // Fragment length: |TLEN| for paired-end reads, and otherwise the
        // reference bases the read covers.
        //
        // Deliberately `covered_span` rather than `alignment_end -
        // alignment_start`: the span counts a skipped region, so a spliced read
        // would be measured as long as its intron and fail the maximum length.
        if self.min_fragment_length.is_some() || self.max_fragment_length.is_some() {
            let frag_len = if rec.template_length != 0 {
                rec.template_length.unsigned_abs() as usize
            } else {
                rec.covered_span
                    .unwrap_or(rec.alignment_end - rec.alignment_start)
            };
            if let Some(min) = self.min_fragment_length
                && frag_len < min
            {
                return false;
            }
            if let Some(max) = self.max_fragment_length
                && frag_len > max
            {
                return false;
            }
        }

        // GC content.
        if (self.min_gc.is_some() || self.max_gc.is_some())
            && let Some(gc) = rec.gc_content
        {
            if let Some(min) = self.min_gc
                && gc < min
            {
                return false;
            }
            if let Some(max) = self.max_gc
                && gc > max
            {
                return false;
            }
        }

        // Aligned fraction.
        if let Some(min_af) = self.min_aligned_fraction
            && let Some(af) = rec.aligned_fraction
            && af < min_af
        {
            return false;
        }

        true
    }
}

impl Default for QcFilter {
    fn default() -> Self {
        Self::new()
    }
}

impl QcFilter {
    /// Returns whether GC content must be computed to evaluate this filter.
    pub(super) fn needs_gc(&self) -> bool {
        self.min_gc.is_some() || self.max_gc.is_some()
    }

    /// Returns whether the aligned fraction must be computed to evaluate this filter.
    pub(super) fn needs_aligned_fraction(&self) -> bool {
        self.min_aligned_fraction.is_some()
    }

    /// Returns whether the covered span must be computed to evaluate this
    /// filter, which the fragment-length bounds measure a single-end read with.
    pub(super) fn needs_covered_span(&self) -> bool {
        self.min_fragment_length.is_some() || self.max_fragment_length.is_some()
    }
}

/// Derive the [`ScRecordOptions`] needed to evaluate `qc` and, when
/// `has_motif` is set, the motif filter (which needs the raw read sequence).
///
/// `adjust` decides whether the gapless blocks shoudl be computed.
pub(crate) fn derive_record_opts(
    qc: Option<&QcFilter>,
    has_motif: bool,
    dedup: bool,
    adjust: &AdjustRead,
) -> ScRecordOptions {
    ScRecordOptions {
        compute_gc: qc.is_some_and(|f| f.needs_gc()),
        compute_aligned_fraction: qc.is_some_and(|f| f.needs_aligned_fraction()),
        store_sequence: has_motif,
        compute_covered_span: dedup || qc.is_some_and(|f| f.needs_covered_span()),
        compute_blocks: adjust.extend_reads.is_none() && !adjust.center_reads,
    }
}

/// Strategy for identifying duplicate reads.
#[derive(Clone, Copy, Debug)]
pub enum DupMethod {
    /// Duplicate key = barcode + alignment start.
    BarcodeStart,
    /// Duplicate key = barcode + alignment start + alignment end.
    BarcodeStartEnd,
    /// Duplicate key = barcode + UMI + alignment start.
    BarcodeUmiStart,
    /// Duplicate key = barcode + UMI + alignment start + alignment end.
    BarcodeUmiStartEnd,
}

// Key tuple: (barcode, umi, fragment_start, fragment_end, mate_reference, strand)
//
// Deduplication is per *fragment*, not per read: two reads that begin at the
// same base can belong to different templates, so the key is built from TLEN
// and the mate position rather than from this read's own alignment span.
//
// `fragment_start` and `fragment_end` are `Option` because the start-only
// methods key on the 5' end alone (the fragment start for a forward read,
// the fragment end for a reverse one) leaving the other side out of the key.
//
// The barcode/UMI bytes are copied into the key only here, when a record is
// actually deduplicated. The hot path borrows them from the record.
//
// This read's own chromosome is deliberately *not* part of the key: each
// `DuplicateFilter` lives for exactly one work chunk, and every chunk covers a
// single chromosome, so it is constant for the filter's whole lifetime. Reads
// sharing an alignment start always land in the same chunk, so duplicates are
// never split across filters. The *mate's* reference does vary, and is keyed.
type DupKey = (
    Option<Vec<u8>>,
    Option<Vec<u8>>,
    Option<usize>,
    Option<usize>,
    Option<usize>,
    bool,
);

/// Signed observed template length, matching the reference implementation's
/// `getTLen(read, notAbs=True)`: TLEN when it is set, and otherwise the read's
/// own span on the reference, so single-end reads still describe a fragment.
fn signed_template_length(rec: &ScRecord<'_>) -> i64 {
    if rec.template_length != 0 {
        rec.template_length
    } else {
        // The reference bases actually covered, so a spliced read describes its
        // exons rather than its whole genomic footprint. `covered_span` is
        // always populated when a duplicate filter is in use, because
        // `derive_record_opts` requests it; the fallback only guards misuse.
        rec.covered_span
            .unwrap_or_else(|| rec.alignment_end.saturating_sub(rec.alignment_start)) as i64
    }
}

/// The fragment `[start, end)` this read belongs to.
///
/// A negative TLEN means this read is the rightmost of the pair, so the
/// fragment begins at the mate and runs back to here.
///
/// The mate position is therefore always present when TLEN is negative: the
/// spec sets TLEN to 0 for a single-segment template or when the placement is
/// unknown, so a negative value implies a mapped mate. The fallback below is
/// only reached on a malformed record, where any stable key will do.
fn fragment_span(rec: &ScRecord<'_>) -> (usize, usize) {
    let tlen = signed_template_length(rec);
    if tlen >= 0 {
        let start = rec.alignment_start;
        (start, start + tlen as usize)
    } else {
        let start = rec.mate_alignment_start.unwrap_or(rec.alignment_start);
        (start, start + tlen.unsigned_abs() as usize)
    }
}

/// Streaming duplicate filter backed by an in-memory fingerprint set.
///
/// A record is a duplicate if its fingerprint (determined by [`DupMethod`]) has
/// already been seen. The first occurrence is kept; subsequent ones are dropped.
pub struct DuplicateFilter {
    pub method: DupMethod,
    seen: AHashSet<DupKey>,
    /// Alignment start of the previous record, which bounds the comparison
    /// window (see `passes`).
    window_start: Option<usize>,
}

impl DuplicateFilter {
    pub fn new(method: DupMethod) -> Self {
        Self {
            method,
            seen: AHashSet::new(),
            window_start: None,
        }
    }

    /// Returns `true` if the record is the **first** occurrence of its fingerprint
    /// and `false` if it is a duplicate.
    pub fn passes(&mut self, rec: &ScRecord<'_>) -> bool {
        // Reads are only ever compared against others starting at the very same
        // base: the window resets whenever the alignment start moves on.  Two
        // reads at different starts are never duplicates of each other, however
        // alike their fragments, so nothing is inferred about their origin.
        // Input is coordinate-sorted, so a start is visited exactly once.
        if self.window_start != Some(rec.alignment_start) {
            self.seen.clear();
            self.window_start = Some(rec.alignment_start);
        }

        let uses_end = matches!(
            self.method,
            DupMethod::BarcodeStartEnd | DupMethod::BarcodeUmiStartEnd
        );
        let uses_umi = matches!(
            self.method,
            DupMethod::BarcodeUmiStart | DupMethod::BarcodeUmiStartEnd
        );

        let (fragment_start, fragment_end) = fragment_span(rec);

        let (start, end, mate_reference) = if uses_end {
            // A chimeric pair spans two references, so it has no fragment end
            // to speak of; the mate's start stands in for one.
            let end = if rec.mate_reference_id == rec.reference_id {
                Some(fragment_end)
            } else {
                rec.mate_alignment_start
            };
            (Some(fragment_start), end, rec.mate_reference_id)
        } else {
            // Key on the 5' end only, so a read and its mate are not mistaken
            // for one another.
            if rec.is_reverse {
                (None, Some(fragment_end), rec.reference_id)
            } else {
                (Some(fragment_start), None, rec.reference_id)
            }
        };

        let key: DupKey = (
            rec.barcode.map(<[u8]>::to_vec),
            if uses_umi {
                rec.umi.map(<[u8]>::to_vec)
            } else {
                None
            },
            start,
            end,
            mate_reference,
            rec.is_reverse,
        );
        // insert returns true when a new key is inserted.
        self.seen.insert(key)
    }
}

/// Filter that checks for a nucleotide motif at the 5′ end of the read and the
/// corresponding genomic overhang.
///
/// A record passes if **any** of the supplied `(read_motif, ref_motif)` pairs
/// match it. For forward reads the read motif is compared to the first N bases
/// of the forward read sequence and the reference motif to the genomic bases
/// immediately upstream of the alignment start. For reverse reads both windows
/// are mirrored to the 3′ end of the alignment on the reference.
///
/// Requires `ScRecordOptions::store_sequence = true` on the records.
pub struct MotifFilter {
    genome: TwoBitFile<BufReader<File>>,
    /// Each entry is (read_motif, ref_motif); the filter passes on any match.
    pub motifs: Vec<(String, String)>,
}

impl MotifFilter {
    pub fn new(genome_path: &Path, motifs: Vec<(String, String)>) -> Result<Self> {
        let genome = TwoBitFile::open(genome_path)?;
        Ok(Self { genome, motifs })
    }

    /// Returns `true` if the record passes the motif filter.
    ///
    /// `chrom` is the record's chromosome (known by the caller from the work
    /// chunk). If `read_sequence` is `None` on the record (sequence was not
    /// stored), the filter is vacuously passed.
    pub fn passes(&mut self, rec: &ScRecord<'_>, chrom: &str) -> Result<bool> {
        let Some(stored_seq) = &rec.read_sequence else {
            return Ok(true);
        };

        let seq_len = stored_seq.len();

        for (read_motif, ref_motif) in &self.motifs {
            let r_len = read_motif.len();
            let g_len = ref_motif.len();

            if seq_len < r_len {
                continue;
            }

            // Compare the read motif against the 5′-most `r_len` bases of the
            // forward-strand read. For forward reads the BAM sequence is
            // already in that orientation; for reverse reads the forward-strand
            // prefix is the reverse-complement of the sequence's `r_len`-base
            // suffix. Both are computed without allocating the full sequence.
            let read_ok = if rec.is_reverse {
                stored_seq[seq_len - r_len..]
                    .iter()
                    .rev()
                    .zip(read_motif.as_bytes())
                    .all(|(&s, &m)| complement(s).eq_ignore_ascii_case(&m))
            } else {
                stored_seq[..r_len].eq_ignore_ascii_case(read_motif.as_bytes())
            };
            if !read_ok {
                continue;
            }

            // Fetch the genomic reference motif.
            let ref_seq = if rec.is_reverse {
                // Reverse read: overhang begins at alignment_end on the reference.
                let g_start = rec.alignment_end.saturating_sub(1);
                let g_end = rec.alignment_end + g_len - 1;
                self.genome.read_sequence(chrom, g_start..g_end)
            } else {
                // Forward read: overhang ends at alignment_start (inclusive).
                let g_end = rec.alignment_start + 1;
                let g_start = g_end.saturating_sub(g_len);
                self.genome.read_sequence(chrom, g_start..g_end)
            };

            let ref_seq = match ref_seq {
                Ok(s) => s,
                Err(_) => continue,
            };

            if ref_seq
                .as_bytes()
                .eq_ignore_ascii_case(ref_motif.as_bytes())
            {
                return Ok(true);
            }
        }

        Ok(false)
    }
}

/// Overlap threshold for a read to be considered blacklisted.
pub const BLACKLIST_MIN_OVERLAP_PERCENT: usize = 50;

/// Returns `true` if at least [`BLACKLIST_MIN_OVERLAP_PERCENT`] of the read
/// interval `[start, end)` is covered by blacklisted regions on `chromosome`.
///
/// Coordinates are 0-based half-open on both sides. The read interval is the
/// full alignment span, CIGAR is not consulted.
pub fn is_blacklisted(
    blacklist_index: &GenomeIndex,
    chromosome: &str,
    start: usize,
    end: usize,
) -> bool {
    let Some(idx) = blacklist_chrom_index(blacklist_index, chromosome) else {
        return false;
    };
    read_is_blacklisted(idx, start, end)
}

/// [`is_blacklisted`] against a chromosome index the caller already holds.
///
/// A work chunk covers one chromosome, so the counting loops resolve the index
/// once per chunk and no read hashes a chromosome name.
pub fn read_is_blacklisted(idx: &ChromIndex, start: usize, end: usize) -> bool {
    let read_len = end.saturating_sub(start);
    if read_len == 0 {
        return false;
    }
    // Rounded up.
    let required = (read_len * BLACKLIST_MIN_OVERLAP_PERCENT).div_ceil(100);

    // Most reads meet no blacklist region at all, and most of the rest meet
    // exactly one, which the early return settles without allocating.
    let mut spans: Vec<(usize, usize)> = Vec::new();
    for iv in idx.find(start, end) {
        let s = iv.start.max(start);
        let e = iv.end.min(end);
        if e <= s {
            continue;
        }
        if e - s >= required {
            return true;
        }
        spans.push((s, e));
    }
    if spans.len() < 2 {
        return false;
    }

    spans.sort_unstable();
    let mut covered = 0usize;
    let (mut cur_start, mut cur_end) = spans[0];
    for &(s, e) in &spans[1..] {
        if s <= cur_end {
            cur_end = cur_end.max(e);
        } else {
            covered += cur_end - cur_start;
            if covered >= required {
                return true;
            }
            (cur_start, cur_end) = (s, e);
        }
    }
    covered += cur_end - cur_start;
    covered >= required
}

/// Look up a chromosome in a blacklist index, bridging the `chr` prefix.
///
/// BAM headers and BED files sometimes disagree about the prefix, so `chr1` and
/// `1` are treated as the same chromosome.
pub fn blacklist_chrom_index<'a>(
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

/// Returns `true` if the read should be **excluded** based on the RNA strand filter.
///
/// This assumes a dUTP-based paired-end library: "forward" keeps reads from genes
/// on the forward strand (= read2-forward or read1-with-forward-mate). Single-end
/// logic inverts the strand relative to paired-end since there is only one read
/// per fragment.
///
/// `flags` is the raw SAM flag field (u16 little-endian bit pattern).
#[inline]
pub fn rna_strand_filter(flags: u16, strand: &str) -> bool {
    let paired = flags & 0x1 != 0;
    if paired {
        match strand {
            // Keep read2-on-forward (0x80 set, 0x10 clear) OR read1-with-forward-mate (0x40 set, 0x20 clear).
            "forward" => !((flags & 0x90 == 0x80) || (flags & 0x60 == 0x40)),
            // Keep read2-on-reverse (0x80 | 0x10 both set) OR read1-with-reverse-mate (0x40 | 0x20 both set).
            "reverse" => !((flags & 0x90 == 0x90) || (flags & 0x60 == 0x60)),
            _ => false,
        }
    } else {
        match strand {
            // dUTP single-end forward: keep reads on reverse strand (0x10 set).
            "forward" => flags & 0x10 == 0,
            // dUTP single-end reverse: keep reads on forward strand (0x10 clear).
            "reverse" => flags & 0x10 != 0,
            _ => false,
        }
    }
}

#[inline]
fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        _ => b'N',
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::region_index::Interval;
    use crate::bam::sc_record::test_record;

    // Blacklist lookup

    fn genome_index(chrom: &str, spans: &[(usize, usize)]) -> GenomeIndex {
        let intervals = spans
            .iter()
            .enumerate()
            .map(|(i, &(start, end))| Interval {
                start,
                end,
                var_idx: i,
            })
            .collect();
        let mut index = GenomeIndex::new();
        index.insert(chrom.to_string(), ChromIndex::build(intervals));
        index
    }

    #[test]
    fn a_read_inside_a_blacklisted_span_is_blacklisted() {
        let bl = genome_index("chr1", &[(100, 200)]);
        assert!(is_blacklisted(&bl, "chr1", 150, 160));
        assert!(
            !is_blacklisted(&bl, "chr1", 200, 300),
            "half-open at the end"
        );
        assert!(
            !is_blacklisted(&bl, "chr1", 0, 100),
            "half-open at the start"
        );
    }

    #[test]
    fn a_read_is_blacklisted_only_once_half_of_it_is_covered() {
        let bl = genome_index("chr1", &[(100, 200)]);
        // 100 bp read, 50 bp inside: exactly at the threshold.
        assert!(is_blacklisted(&bl, "chr1", 150, 250));
        // One base less is one base short.
        assert!(!is_blacklisted(&bl, "chr1", 151, 251));
        // A read that only clips the edge stays.
        assert!(!is_blacklisted(&bl, "chr1", 199, 299));
    }

    #[test]
    fn an_odd_length_read_needs_more_than_half_its_bases() {
        // 101 bp read: 50 bp is under half, 51 bp is over it.
        let bl = genome_index("chr1", &[(100, 150)]);
        assert!(!is_blacklisted(&bl, "chr1", 100, 201));
        let bl = genome_index("chr1", &[(100, 151)]);
        assert!(is_blacklisted(&bl, "chr1", 100, 201));
    }

    #[test]
    fn several_blacklist_regions_add_up_without_double_counting() {
        // Two disjoint 30 bp spans cover 60 of the read's 100 bases.
        let bl = genome_index("chr1", &[(100, 130), (160, 190)]);
        assert!(is_blacklisted(&bl, "chr1", 100, 200));

        // Two 20 bp spans cover only 40, which is under the threshold.
        let bl = genome_index("chr1", &[(100, 120), (160, 180)]);
        assert!(!is_blacklisted(&bl, "chr1", 100, 200));

        // Overlapping spans reach 60 bases together, but a naive sum would
        // count the shared 40 twice and wrongly pass on 100 bases of read.
        let bl = genome_index("chr1", &[(100, 140), (110, 145), (120, 160)]);
        assert!(!is_blacklisted(&bl, "chr1", 100, 300));
    }

    #[test]
    fn the_chr_prefix_is_bridged_in_both_directions() {
        // BAMs and BED files disagree about the "chr" prefix all the time.
        let with_prefix = genome_index("chr1", &[(100, 200)]);
        assert!(is_blacklisted(&with_prefix, "1", 150, 160));
        assert!(blacklist_chrom_index(&with_prefix, "1").is_some());

        let without_prefix = genome_index("1", &[(100, 200)]);
        assert!(is_blacklisted(&without_prefix, "chr1", 150, 160));
        assert!(blacklist_chrom_index(&without_prefix, "chr1").is_some());
    }

    #[test]
    fn a_chromosome_absent_from_the_blacklist_is_never_blacklisted() {
        let bl = genome_index("chr1", &[(100, 200)]);
        assert!(blacklist_chrom_index(&bl, "chr9").is_none());
        assert!(!is_blacklisted(&bl, "chr9", 150, 160));
    }

    #[test]
    fn an_empty_read_interval_is_never_blacklisted() {
        let bl = genome_index("chr1", &[(100, 200)]);
        assert!(!is_blacklisted(&bl, "chr1", 150, 150));
    }

    #[test]
    fn the_hoisted_call_answers_exactly_as_the_lookup_one() {
        // What the counting loops call once the chromosome is resolved.
        let bl = genome_index("chr1", &[(100, 200)]);
        let idx = blacklist_chrom_index(&bl, "chr1").unwrap();
        for (start, end) in [(150, 160), (150, 250), (151, 251), (199, 299), (150, 150)] {
            assert_eq!(
                read_is_blacklisted(idx, start, end),
                is_blacklisted(&bl, "chr1", start, end),
                "[{start}, {end})"
            );
        }
    }

    // SAM flag bits used below.
    const PAIRED: u16 = 0x1;
    const PROPER_PAIR: u16 = 0x2;
    const REVERSE: u16 = 0x10;
    const MATE_REVERSE: u16 = 0x20;
    const READ1: u16 = 0x40;
    const READ2: u16 = 0x80;
    const DUPLICATE: u16 = 0x400;

    #[test]
    fn a_filter_with_no_thresholds_keeps_everything() {
        let f = RawRecordFilter::new();
        assert!(f.passes(0, None));
        assert!(f.passes(DUPLICATE, Some(0)));
    }

    #[test]
    fn min_mapq_rejects_low_and_missing_mapping_qualities() {
        let f = RawRecordFilter {
            min_mapq: Some(30),
            ..RawRecordFilter::new()
        };

        assert!(f.passes(0, Some(30)));
        assert!(f.passes(0, Some(60)));
        assert!(!f.passes(0, Some(29)));
        // An absent MAPQ cannot clear the bar.
        assert!(!f.passes(0, None));
    }

    #[test]
    fn sam_flag_include_requires_every_requested_bit() {
        let f = RawRecordFilter {
            sam_flag_include: Some(PAIRED | PROPER_PAIR),
            ..RawRecordFilter::new()
        };

        assert!(f.passes(PAIRED | PROPER_PAIR, None));
        assert!(f.passes(PAIRED | PROPER_PAIR | REVERSE, None));
        // Only one of the two required bits is set.
        assert!(!f.passes(PAIRED, None));
        assert!(!f.passes(0, None));
    }

    #[test]
    fn sam_flag_exclude_rejects_any_forbidden_bit() {
        let f = RawRecordFilter {
            sam_flag_exclude: Some(DUPLICATE | 0x200),
            ..RawRecordFilter::new()
        };

        assert!(f.passes(PAIRED, None));
        assert!(!f.passes(DUPLICATE, None));
        assert!(!f.passes(0x200, None));
    }

    #[test]
    fn paired_end_rna_strand_filter_follows_the_dutp_convention() {
        // "forward" keeps read2-on-forward and read1-with-forward-mate.
        assert!(!rna_strand_filter(PAIRED | READ2, "forward"));
        assert!(!rna_strand_filter(PAIRED | READ1, "forward"));
        assert!(rna_strand_filter(PAIRED | READ2 | REVERSE, "forward"));
        assert!(rna_strand_filter(PAIRED | READ1 | MATE_REVERSE, "forward"));

        // "reverse" is the mirror image.
        assert!(!rna_strand_filter(PAIRED | READ2 | REVERSE, "reverse"));
        assert!(!rna_strand_filter(PAIRED | READ1 | MATE_REVERSE, "reverse"));
        assert!(rna_strand_filter(PAIRED | READ2, "reverse"));
        assert!(rna_strand_filter(PAIRED | READ1, "reverse"));
    }

    #[test]
    fn single_end_rna_strand_filter_inverts_the_paired_end_logic() {
        // dUTP single-end: a "forward" gene yields reads on the reverse strand.
        assert!(!rna_strand_filter(REVERSE, "forward"));
        assert!(rna_strand_filter(0, "forward"));

        assert!(!rna_strand_filter(0, "reverse"));
        assert!(rna_strand_filter(REVERSE, "reverse"));
    }

    #[test]
    fn an_unrecognized_strand_name_excludes_nothing() {
        assert!(!rna_strand_filter(PAIRED | READ2, "sideways"));
        assert!(!rna_strand_filter(0, ""));
    }

    #[test]
    fn record_filter_applies_the_rna_strand_filter() {
        let f = RawRecordFilter {
            filter_rna_strand: Some("forward".to_string()),
            ..RawRecordFilter::new()
        };

        assert!(f.passes(PAIRED | READ2, None));
        assert!(!f.passes(PAIRED | READ2 | REVERSE, None));
    }

    #[test]
    fn qc_filter_with_no_thresholds_keeps_everything() {
        assert!(QcFilter::new().passes(&test_record(100, 200)));
    }

    #[test]
    fn fragment_length_uses_the_alignment_span_for_single_end_reads() {
        let f = QcFilter {
            min_fragment_length: Some(50),
            max_fragment_length: Some(150),
            ..QcFilter::new()
        };

        assert!(f.passes(&test_record(1000, 1100))); // span 100
        assert!(f.passes(&test_record(1000, 1050))); // span 50, at the minimum
        assert!(f.passes(&test_record(1000, 1150))); // span 150, at the maximum
        assert!(!f.passes(&test_record(1000, 1049))); // span 49
        assert!(!f.passes(&test_record(1000, 1151))); // span 151
    }

    #[test]
    fn fragment_length_uses_the_insert_size_for_paired_end_reads() {
        let f = QcFilter {
            max_fragment_length: Some(150),
            ..QcFilter::new()
        };

        // A 100 bp alignment on a 500 bp fragment is judged by |TLEN|, not the span.
        let mut rec = test_record(1000, 1100);
        rec.template_length = 500;
        assert!(!f.passes(&rec));

        // The sign of TLEN is irrelevant.
        rec.template_length = -500;
        assert!(!f.passes(&rec));

        rec.template_length = 120;
        assert!(f.passes(&rec));

        // A covered span is not consulted when the read carries an insert size.
        rec.covered_span = Some(100);
        rec.template_length = 500;
        assert!(!f.passes(&rec));
    }

    #[test]
    fn gc_thresholds_are_inclusive_at_the_bounds() {
        let f = QcFilter {
            min_gc: Some(0.3),
            max_gc: Some(0.7),
            ..QcFilter::new()
        };

        let mut rec = test_record(0, 100);
        for (gc, expected) in [
            (0.3, true),
            (0.5, true),
            (0.7, true),
            (0.29, false),
            (0.71, false),
        ] {
            rec.gc_content = Some(gc);
            assert_eq!(f.passes(&rec), expected, "gc = {gc}");
        }
    }

    #[test]
    fn gc_and_aligned_fraction_thresholds_pass_when_the_value_was_not_computed() {
        // The counting loop only computes these when a filter asks for them, so
        // an uncomputed value must never silently drop a read.
        let f = QcFilter {
            min_gc: Some(0.9),
            min_aligned_fraction: Some(0.9),
            ..QcFilter::new()
        };

        let rec = test_record(0, 100);
        assert!(rec.gc_content.is_none());
        assert!(rec.aligned_fraction.is_none());
        assert!(f.passes(&rec));
    }

    #[test]
    fn min_aligned_fraction_rejects_poorly_matched_reads() {
        let f = QcFilter {
            min_aligned_fraction: Some(0.8),
            ..QcFilter::new()
        };

        let mut rec = test_record(0, 100);
        rec.aligned_fraction = Some(0.8);
        assert!(f.passes(&rec));
        rec.aligned_fraction = Some(0.79);
        assert!(!f.passes(&rec));
    }

    #[test]
    fn record_options_are_derived_from_the_active_filters() {
        let none = derive_record_opts(None, false, false, &AdjustRead::default());
        assert!(!none.compute_gc);
        assert!(!none.compute_aligned_fraction);
        assert!(!none.store_sequence);

        let gc_only = QcFilter {
            max_gc: Some(0.6),
            ..QcFilter::new()
        };
        let opts = derive_record_opts(Some(&gc_only), false, false, &AdjustRead::default());
        assert!(opts.compute_gc);
        assert!(!opts.compute_aligned_fraction);

        let af_only = QcFilter {
            min_aligned_fraction: Some(0.5),
            ..QcFilter::new()
        };
        let opts = derive_record_opts(Some(&af_only), true, false, &AdjustRead::default());
        assert!(!opts.compute_gc);
        assert!(opts.compute_aligned_fraction);
        // The motif filter is the only consumer of the stored sequence.
        assert!(opts.store_sequence);
    }

    #[test]
    fn the_covered_span_is_computed_for_the_fragment_length_bounds() {
        let opts = |qc: Option<&QcFilter>, dedup: bool| {
            derive_record_opts(qc, false, dedup, &AdjustRead::default()).compute_covered_span
        };

        assert!(!opts(None, false));
        assert!(opts(None, true), "deduplication still needs it");

        let bounded = QcFilter {
            max_fragment_length: Some(500),
            ..QcFilter::new()
        };
        assert!(opts(Some(&bounded), false));

        // A filter that does not measure a fragment that does not need it.
        let gc_only = QcFilter {
            max_gc: Some(0.6),
            ..QcFilter::new()
        };
        assert!(!opts(Some(&gc_only), false));
    }

    #[test]
    fn the_fragment_length_bounds_ignore_a_spliced_read_intron() {
        // A 100 bp read over a 20 kb intron is a 100 bp fragment.
        let filter = QcFilter {
            max_fragment_length: Some(500),
            ..QcFilter::new()
        };

        let mut rec = test_record(1000, 21_100);
        rec.covered_span = Some(100);
        assert!(filter.passes(&rec));

        // The bound still filters a read that really is that long.
        let mut long = test_record(1000, 21_100);
        long.covered_span = Some(20_100);
        assert!(!filter.passes(&long));
    }

    #[test]
    fn the_gapless_blocks_are_computed_only_when_the_read_keeps_its_own_shape() {
        let blocks_of =
            |adjust: &AdjustRead| derive_record_opts(None, false, false, adjust).compute_blocks;

        // Counted on the alignment itself: the blocks are the interval.
        assert!(blocks_of(&AdjustRead::default()));

        // Both adjustments replace the alignment with a synthetic interval, so
        // there is no point splitting it.
        assert!(!blocks_of(&AdjustRead {
            extend_reads: Some(300),
            center_reads: false,
            ..AdjustRead::default()
        }));
        assert!(!blocks_of(&AdjustRead {
            extend_reads: None,
            center_reads: true,
            ..AdjustRead::default()
        }));
    }

    /// A record carrying a barcode/UMI, since the duplicate key is built from them.
    fn dup_record<'a>(barcode: &'a [u8], umi: &'a [u8], start: usize, end: usize) -> ScRecord<'a> {
        let mut rec = test_record(start, end);
        rec.barcode = Some(barcode);
        rec.umi = Some(umi);
        rec
    }

    #[test]
    fn barcode_start_dedup_keeps_only_the_first_read_at_a_position() {
        let mut f = DuplicateFilter::new(DupMethod::BarcodeStart);

        assert!(f.passes(&dup_record(b"AAA", b"U1", 100, 200)));
        // Same barcode and start: a duplicate, even with a different end and UMI.
        assert!(!f.passes(&dup_record(b"AAA", b"U2", 100, 250)));
        // A different barcode at the same start is a different cell.
        assert!(f.passes(&dup_record(b"CCC", b"U1", 100, 200)));
        // A different start is a different fragment.
        assert!(f.passes(&dup_record(b"AAA", b"U1", 101, 200)));
    }

    #[test]
    fn barcode_start_end_dedup_distinguishes_reads_by_their_end() {
        let mut f = DuplicateFilter::new(DupMethod::BarcodeStartEnd);

        assert!(f.passes(&dup_record(b"AAA", b"U1", 100, 200)));
        assert!(!f.passes(&dup_record(b"AAA", b"U1", 100, 200)));
        // Same start, different end: kept.
        assert!(f.passes(&dup_record(b"AAA", b"U1", 100, 250)));
    }

    #[test]
    fn umi_aware_dedup_keeps_distinct_umis_at_the_same_position() {
        let mut f = DuplicateFilter::new(DupMethod::BarcodeUmiStart);

        assert!(f.passes(&dup_record(b"AAA", b"U1", 100, 200)));
        assert!(!f.passes(&dup_record(b"AAA", b"U1", 100, 200)));
        // A different UMI at the same position is an independent molecule.
        assert!(f.passes(&dup_record(b"AAA", b"U2", 100, 200)));
    }

    #[test]
    fn barcode_umi_start_end_dedup_uses_the_full_fingerprint() {
        let mut f = DuplicateFilter::new(DupMethod::BarcodeUmiStartEnd);

        assert!(f.passes(&dup_record(b"AAA", b"U1", 100, 200)));
        assert!(!f.passes(&dup_record(b"AAA", b"U1", 100, 200)));
        assert!(f.passes(&dup_record(b"AAA", b"U1", 100, 250)));
        assert!(f.passes(&dup_record(b"AAA", b"U2", 100, 200)));
    }

    #[test]
    fn reads_without_a_barcode_still_deduplicate_against_each_other() {
        let mut f = DuplicateFilter::new(DupMethod::BarcodeStart);

        assert!(f.passes(&test_record(100, 200)));
        assert!(!f.passes(&test_record(100, 200)));
        assert!(f.passes(&test_record(300, 400)));
    }

    // Fragment-aware deduplication

    #[test]
    fn two_reads_sharing_a_start_but_not_a_fragment_are_not_duplicates() {
        // Taken from tests/testdata/test_i1.bam: both ATATAACT reads align at
        // 65966812 with 61M, so their read spans are identical, but one has
        // TLEN 0 and the other TLEN -612. They are different templates.
        let mut f = DuplicateFilter::new(DupMethod::BarcodeStartEnd);

        let mut lone = dup_record(b"AAA", b"U1", 65_966_811, 65_966_872);
        lone.is_reverse = true;
        lone.template_length = 0;

        let mut paired = dup_record(b"AAA", b"U1", 65_966_811, 65_966_872);
        paired.is_reverse = true;
        paired.template_length = -612;
        paired.mate_alignment_start = Some(65_966_260);

        assert!(f.passes(&lone));
        assert!(
            f.passes(&paired),
            "reads from different fragments must not deduplicate against each other"
        );
    }

    #[test]
    fn a_forward_and_a_reverse_read_at_one_position_stay_distinct() {
        let mut f = DuplicateFilter::new(DupMethod::BarcodeStart);

        let forward = dup_record(b"AAA", b"U1", 100, 200);
        let mut reverse = dup_record(b"AAA", b"U1", 100, 200);
        reverse.is_reverse = true;

        assert!(f.passes(&forward));
        assert!(f.passes(&reverse), "strand is part of the fingerprint");
    }

    #[test]
    fn reads_at_different_starts_are_never_duplicates() {
        // The GTCAAGCA pair in tests/testdata/test_i1.bam: same barcode, same
        // mate, same TLEN, ends flush -- but starts 4 bp apart because of 5'
        // trimming.  Nothing is inferred from that similarity; a differing
        // start is enough to keep both.
        let mut f = DuplicateFilter::new(DupMethod::BarcodeStart);

        let mut first = dup_record(b"AAA", b"U1", 65_977_534, 65_977_596);
        first.is_reverse = true;
        first.template_length = -270;
        first.mate_alignment_start = Some(65_977_326);

        let mut second = dup_record(b"AAA", b"U1", 65_977_538, 65_977_596);
        second.is_reverse = true;
        second.template_length = -270;
        second.mate_alignment_start = Some(65_977_326);

        assert!(f.passes(&first));
        assert!(
            f.passes(&second),
            "a differing alignment start is never a duplicate"
        );
    }

    #[test]
    fn a_spliced_read_is_sized_by_its_exons_not_its_intron() {
        // A read like 50M1000N50M covers 100 reference bases across a 1100 bp
        // footprint.  `alignment_end - alignment_start` would say 1100 and
        // swallow the intron; the fragment is the 100 bases actually covered.
        let mut f = DuplicateFilter::new(DupMethod::BarcodeStartEnd);

        let mut spliced = dup_record(b"AAA", b"U1", 1_000, 2_100);
        spliced.template_length = 0;
        spliced.covered_span = Some(100);

        // A second read covering the same 100 bases from the same start is a
        // duplicate, however its intron is placed.
        let mut other_intron = dup_record(b"AAA", b"U1", 1_000, 3_500);
        other_intron.template_length = 0;
        other_intron.covered_span = Some(100);

        assert!(f.passes(&spliced));
        assert!(
            !f.passes(&other_intron),
            "the intron must not enter the fragment size"
        );
    }

    #[test]
    fn a_spliced_read_is_distinct_from_an_unspliced_one_of_the_same_footprint() {
        let mut f = DuplicateFilter::new(DupMethod::BarcodeStartEnd);

        let mut spliced = dup_record(b"AAA", b"U1", 1_000, 2_100);
        spliced.template_length = 0;
        spliced.covered_span = Some(100);

        let mut contiguous = dup_record(b"AAA", b"U1", 1_000, 2_100);
        contiguous.template_length = 0;
        contiguous.covered_span = Some(1_100);

        assert!(f.passes(&spliced));
        assert!(
            f.passes(&contiguous),
            "different covered spans, different fragments"
        );
    }

    #[test]
    fn the_comparison_window_spans_only_one_alignment_start() {
        let mut f = DuplicateFilter::new(DupMethod::BarcodeStart);
        let here = dup_record(b"AAA", b"U1", 100, 200);

        assert!(f.passes(&here));
        assert!(!f.passes(&here), "same start and same key is a duplicate");

        // Moving to a new start clears the window ...
        assert!(f.passes(&dup_record(b"AAA", b"U1", 300, 400)));
        // ... so the earlier fingerprint is no longer remembered.
        assert!(f.passes(&here));
    }

    #[test]
    fn a_chimeric_pair_is_kept_apart_from_a_normal_one() {
        let mut f = DuplicateFilter::new(DupMethod::BarcodeStartEnd);

        let mut normal = dup_record(b"AAA", b"U1", 100, 200);
        normal.reference_id = Some(0);
        normal.mate_reference_id = Some(0);

        let mut chimeric = dup_record(b"AAA", b"U1", 100, 200);
        chimeric.reference_id = Some(0);
        chimeric.mate_reference_id = Some(1);
        chimeric.mate_alignment_start = Some(5_000);

        assert!(f.passes(&normal));
        assert!(
            f.passes(&chimeric),
            "a mate on another reference is another fragment"
        );
    }

    #[test]
    fn a_read_without_a_template_length_uses_its_own_span_as_the_fragment() {
        // TLEN 0 means single-end, so the fragment is the read's own span.
        let mut f = DuplicateFilter::new(DupMethod::BarcodeStartEnd);

        let implied = dup_record(b"AAA", b"U1", 100, 200);
        let mut stated = dup_record(b"AAA", b"U1", 100, 200);
        stated.template_length = 100;

        assert!(f.passes(&implied));
        assert!(
            !f.passes(&stated),
            "an explicit TLEN of 100 describes the same fragment"
        );
    }

    #[test]
    fn complement_maps_bases_and_treats_anything_else_as_n() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'C'), b'G');
        // Lower case is handled, and unknown bases collapse to N.
        assert_eq!(complement(b'a'), b'T');
        assert_eq!(complement(b'N'), b'N');
        assert_eq!(complement(b'X'), b'N');
    }
}
