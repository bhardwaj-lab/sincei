use ahash::AHashSet;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use anyhow::Result;
use twobit::TwoBitFile;

use super::sc_record::{ScRecord, ScRecordOptions};

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
// actually deduplicated — the hot path borrows them from the record.
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
/// already been seen.  The first occurrence is kept; subsequent ones are dropped.
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
/// matches.  For forward reads the read motif is compared to the first N bases
/// of the forward read sequence and the reference motif to the genomic bases
/// immediately upstream of the alignment start.  For reverse reads both windows
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
    /// chunk).  If `read_sequence` is `None` on the record (sequence was not
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
            // forward-strand read.  For forward reads the BAM sequence is
            // already in that orientation; for reverse reads the forward-strand
            // prefix is the reverse-complement of the sequence's `r_len`-base
            // suffix.  Both are computed without allocating the full sequence.
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
/// on the forward strand (= read2-forward or read1-with-forward-mate).  Single-end
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
