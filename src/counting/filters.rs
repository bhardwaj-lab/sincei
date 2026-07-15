use ahash::AHashSet;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use anyhow::Result;
use twobit::TwoBitFile;

use crate::bam::sc_record::{ScRecord, ScRecordOptions};

// Raw-record filter

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

// QC filter

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
    /// Minimum fraction of read bases in M-type CIGAR ops.
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
        // Fragment length: use |TLEN| for paired-end reads, alignment span otherwise.
        if self.min_fragment_length.is_some() || self.max_fragment_length.is_some() {
            let frag_len = if rec.template_length != 0 {
                rec.template_length.unsigned_abs() as usize
            } else {
                rec.alignment_end - rec.alignment_start
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
}

/// Derive the [`ScRecordOptions`] needed to evaluate `qc` and, when
/// `has_motif` is set, the motif filter (which needs the raw read sequence).
pub(crate) fn derive_record_opts(qc: Option<&QcFilter>, has_motif: bool) -> ScRecordOptions {
    ScRecordOptions {
        compute_gc: qc.is_some_and(|f| f.needs_gc()),
        compute_aligned_fraction: qc.is_some_and(|f| f.needs_aligned_fraction()),
        store_sequence: has_motif,
    }
}

// Duplicate filter

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

// Key tuple: (barcode, start, end_or_zero, umi_or_none)
//
// The barcode/UMI bytes are copied into the key only here, when a record is
// actually deduplicated. The hot path borrows them from the record.
//
// The chromosome is deliberately *not* part of the key: each `DuplicateFilter`
// lives for exactly one work chunk, and every chunk covers a single
// chromosome, so the chromosome is constant for the filter's whole lifetime.
// Reads sharing an alignment start always land in the same chunk, so
// duplicates are never split across filters.
type DupKey = (Option<Vec<u8>>, usize, usize, Option<Vec<u8>>);

/// Streaming duplicate filter backed by an in-memory fingerprint set.
///
/// A record is a duplicate if its fingerprint (determined by [`DupMethod`]) has
/// already been seen. The first occurrence is kept; subsequent ones are dropped.
pub struct DuplicateFilter {
    pub method: DupMethod,
    seen: AHashSet<DupKey>,
}

impl DuplicateFilter {
    pub fn new(method: DupMethod) -> Self {
        Self {
            method,
            seen: AHashSet::new(),
        }
    }

    /// Returns `true` if the record is the **first** occurrence of its fingerprint
    /// and should be kept; `false` if it is a duplicate and should be dropped.
    pub fn passes(&mut self, rec: &ScRecord<'_>) -> bool {
        let barcode = || rec.barcode.map(<[u8]>::to_vec);
        let umi = || rec.umi.map(<[u8]>::to_vec);
        let key: DupKey = match self.method {
            DupMethod::BarcodeStart => (barcode(), rec.alignment_start, 0, None),
            DupMethod::BarcodeStartEnd => (barcode(), rec.alignment_start, rec.alignment_end, None),
            DupMethod::BarcodeUmiStart => (barcode(), rec.alignment_start, 0, umi()),
            DupMethod::BarcodeUmiStartEnd => {
                (barcode(), rec.alignment_start, rec.alignment_end, umi())
            }
        };
        // insert returns true when a new key is inserted.
        self.seen.insert(key)
    }
}

// Motif filter

/// Filter that checks for a nucleotide motif at the 5′ end of the read and the
/// corresponding genomic overhang.
///
/// A record passes if **any** of the supplied `(read_motif, ref_motif)` pairs
/// matches. For forward reads the read motif is compared to the first N bases
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

/// Returns `true` if the read should be **excluded** based on the RNA strand filter.
///
/// This assumes a dUTP-based paired-end library: "forward" keeps reads from genes
/// on the forward strand (= read2-forward or read1-with-forward-mate). Single-end
/// logic inverts the strand relative to paired-end since there is only one read
/// per fragment.
///
/// `flags` is the raw SAM flag field (u16 little-endian bit pattern).
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
    use crate::bam::sc_record::test_record;

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
        let none = derive_record_opts(None, false);
        assert!(!none.compute_gc);
        assert!(!none.compute_aligned_fraction);
        assert!(!none.store_sequence);

        let gc_only = QcFilter {
            max_gc: Some(0.6),
            ..QcFilter::new()
        };
        let opts = derive_record_opts(Some(&gc_only), false);
        assert!(opts.compute_gc);
        assert!(!opts.compute_aligned_fraction);

        let af_only = QcFilter {
            min_aligned_fraction: Some(0.5),
            ..QcFilter::new()
        };
        let opts = derive_record_opts(Some(&af_only), true);
        assert!(!opts.compute_gc);
        assert!(opts.compute_aligned_fraction);
        // The motif filter is the only consumer of the stored sequence.
        assert!(opts.store_sequence);
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
