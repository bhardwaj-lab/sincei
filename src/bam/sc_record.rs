//! Abstraction of a single-cell sequencing read for filtering and counting.
//!
//! [`ScRecord`] is a raw record parsed once into the fields the filters and
//! counting loops use: barcode/UMI (borrowed, not copied), coordinates, flags,
//! and the optional QC values that `ScRecordOptions` switches on. The optional
//! work is opt-in because this runs on every record in the file.
//!
//! [`ScRecord::effective_interval`] then turns an alignment into the interval
//! the read is credited to, applying the extension and centering in
//! [`AdjustRead`].

use anyhow::{Context, Result};
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::Record as _;
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles::sam::alignment::record::data::field::{Tag, Value};

/// Controls which optional fields are computed when building an [`ScRecord`].
pub(crate) struct ScRecordOptions {
    pub(crate) compute_gc: bool,
    pub(crate) compute_aligned_fraction: bool,
    /// Needed by the motif filter.
    pub(crate) store_sequence: bool,
    /// Needed by the duplicate filter to size a fragment when TLEN is 0.
    pub(crate) compute_covered_span: bool,
}

/// How a read's position can be adjusted.
#[derive(Clone, Copy, Debug, Default)]
pub struct AdjustRead {
    /// Extend each read to this fragment length (bp).  For a properly paired
    /// read the observed insert size (TLEN) is used, up to `4 × extend_reads`.
    /// `None` leaves the alignment span unchanged.
    pub extend_reads: Option<usize>,
    /// After extension, replace the interval with a `read_length` window
    /// centered on its midpoint.
    pub center_reads: bool,
}

/// A parsed, filter-passing BAM record carrying the fields needed for
/// single-cell feature counting.
///
/// Coordinates are 0-based half-open `[start, end)`.
///
/// The chromosome is intentionally **not** stored: every caller queries one
/// chromosome at a time, so the chromosome name is already known from the
/// work chunk and need not be re-derived (and re-allocated) per record.
pub struct ScRecord<'a> {
    pub alignment_start: usize,
    pub alignment_end: usize,
    /// Whether the read is on the reverse strand.
    pub is_reverse: bool,
    /// True if the read is part of a paired-end library (FLAG 0x1).
    pub is_paired: bool,
    /// SAM proper-pair flag (FLAG 0x2).
    pub is_proper_pair: bool,
    /// True if this is the first read in the pair (FLAG 0x40).
    pub is_read1: bool,
    /// True if the mate maps to the reverse strand (FLAG 0x20).
    pub mate_is_reverse: bool,
    /// Signed template length (TLEN); 0 for single-end reads.
    pub template_length: i64,
    /// Mate's 0-based alignment start, if available.
    pub mate_alignment_start: Option<usize>,
    /// Reference bases the alignment covers, counting M/D/=/X but **not**
    /// skipped (N) regions. This is the span a single-end read contributes as a
    /// fragment.
    pub covered_span: Option<usize>,
    /// Reference this read aligns to.
    pub reference_id: Option<usize>,
    /// Reference the mate aligns to. Differs from `reference_id` for a chimeric
    /// pair, which the duplicate filter has to treat as a distinct fragment.
    pub mate_reference_id: Option<usize>,
    /// Number of bases in the read sequence (used by center_reads).
    pub read_length: usize,
    /// Barcode tag bytes, borrowed from the BAM record (no allocation).
    pub barcode: Option<&'a [u8]>,
    /// UMI tag bytes, borrowed from the BAM record (no allocation).
    pub umi: Option<&'a [u8]>,
    /// Value to add to the count matrix (defaults to 1).
    pub count: u32,
    /// GC fraction in `[0, 1]`. `None` if not requested.
    pub gc_content: Option<f32>,
    /// Fraction of read bases in M-type CIGAR ops. `None` if not requested.
    pub aligned_fraction: Option<f32>,
    /// Raw BAM read sequence (read orientation). `None` if not requested.
    pub read_sequence: Option<Vec<u8>>,
}

impl<'a> ScRecord<'a> {
    /// Try to build a `ScRecord` from a raw BAM record.
    ///
    /// Returns `Ok(None)` when the record is filtered out.
    /// Returns `Err` only on I/O / decoding failures.
    pub(crate) fn from_bam_record(
        record: &'a bam::Record,
        _header: &Header,
        bc_tag: &Tag,
        umi_tag: Option<&Tag>,
        count_tag: Option<&Tag>,
        opts: &ScRecordOptions,
    ) -> Result<Option<Self>> {
        let flags = record.flags();

        if flags.is_unmapped() {
            return Ok(None);
        }

        let Some(aln_start) = record
            .alignment_start()
            .transpose()
            .context("failed to read alignment start")?
        else {
            return Ok(None);
        };
        let Some(aln_end) = record
            .alignment_end()
            .transpose()
            .context("failed to read alignment end")?
        else {
            return Ok(None);
        };

        // noodles returns 1-based inclusive; convert to 0-based half-open.
        let start = aln_start.get().saturating_sub(1);
        let end = aln_end.get();
        if end <= start {
            return Ok(None);
        }

        let is_reverse = flags.is_reverse_complemented();
        let is_proper_pair = flags.is_properly_segmented();
        let is_paired = flags.is_segmented();
        let is_read1 = flags.is_first_segment();
        let mate_is_reverse = flags.is_mate_reverse_complemented();
        let template_length = record.template_length() as i64;

        // Mate position (1-based -> 0-based).
        let mate_alignment_start = record
            .mate_alignment_start()
            .map(|res| res.map(|pos| pos.get().saturating_sub(1)))
            .transpose()
            .context("failed to read mate alignment start")?;

        let reference_id = record
            .reference_sequence_id()
            .transpose()
            .context("failed to read reference sequence id")?;
        let mate_reference_id = record
            .mate_reference_sequence_id()
            .transpose()
            .context("failed to read mate reference sequence id")?;

        // Borrow the barcode/UMI bytes directly from the record with no per-read
        // allocation. Callers look these up against a byte-keyed whitelist and
        // only the duplicate filter ever copies them.
        let barcode = get_tag_bytes(record, bc_tag)?;

        let umi = match umi_tag {
            Some(utag) => get_tag_bytes(record, utag)?,
            None => None,
        };

        let count = match count_tag {
            Some(ctag) => get_count_tag(record, ctag)?,
            None => Some(1u32),
        };
        let count = count.unwrap_or(1);

        // Sequence handling. Only materialize a `Vec` when the motif filter
        // needs the full read sequence (`store_sequence`). When only GC is
        // required we stream the sequence once and avoid the allocation.
        let (read_sequence, read_length, gc_content) = if opts.store_sequence {
            let seq: Vec<u8> = record.sequence().iter().collect();
            let len = seq.len();
            let gc = if opts.compute_gc {
                gc_fraction(&seq)
            } else {
                None
            };
            (Some(seq), len, gc)
        } else if opts.compute_gc {
            let mut gc_count = 0usize;
            let mut len = 0usize;
            for b in record.sequence().iter() {
                if b == b'G' || b == b'C' {
                    gc_count += 1;
                }
                len += 1;
            }
            let gc = if len > 0 {
                Some(gc_count as f32 / len as f32)
            } else {
                None
            };
            (None, len, gc)
        } else {
            (None, record.sequence().len(), None)
        };

        // Reference span excluding skips, for the duplicate filter's fragment
        // sizing. Deliberately not `alignment_end - alignment_start`: that
        // counts N, so a spliced read would look as long as its intron.
        let covered_span = if opts.compute_covered_span {
            let mut span = 0usize;
            for result in record.cigar().iter() {
                let op = result.context("failed to decode CIGAR op")?;
                if matches!(
                    op.kind(),
                    CigarKind::Match
                        | CigarKind::Deletion
                        | CigarKind::SequenceMatch
                        | CigarKind::SequenceMismatch
                ) {
                    span += op.len();
                }
            }
            Some(span)
        } else {
            None
        };

        // Aligned fraction: M operations / read-consuming length.
        let aligned_fraction = if opts.compute_aligned_fraction {
            let mut match_len: usize = 0;
            let mut read_consuming: usize = 0;
            for result in record.cigar().iter() {
                let op = result.context("failed to decode CIGAR op")?;
                if matches!(op.kind(), CigarKind::Match) {
                    match_len += op.len();
                }
                if op.kind().consumes_read() {
                    read_consuming += op.len();
                }
            }
            if read_consuming > 0 {
                Some(match_len as f32 / read_consuming as f32)
            } else {
                None
            }
        } else {
            None
        };

        Ok(Some(Self {
            alignment_start: start,
            alignment_end: end,
            is_reverse,
            is_proper_pair,
            is_paired,
            is_read1,
            mate_is_reverse,
            template_length,
            mate_alignment_start,
            covered_span,
            reference_id,
            mate_reference_id,
            read_length,
            barcode,
            umi,
            count,
            gc_content,
            aligned_fraction,
            read_sequence,
        }))
    }

    /// Genomic interval to credit this read to, after applying the extension
    /// and centering in `adjust`.
    pub(crate) fn effective_interval(&self, adjust: &AdjustRead) -> (usize, usize) {
        let (start, end) = match adjust.extend_reads {
            None => (self.alignment_start, self.alignment_end),
            Some(frag_len) => {
                let max_paired = 4 * frag_len;
                let tlen_abs = self.template_length.unsigned_abs() as usize;
                let use_tlen = self.is_proper_pair && tlen_abs > 0 && tlen_abs <= max_paired;

                if use_tlen {
                    if self.is_reverse {
                        let mate = self.mate_alignment_start.unwrap_or(self.alignment_start);
                        (mate, self.alignment_end)
                    } else {
                        (self.alignment_start, self.alignment_start + tlen_abs)
                    }
                } else if self.is_reverse {
                    // Extend leftward to `frag_len`, but never shrink a read that
                    // already spans more than that (min keeps the original start).
                    (
                        self.alignment_start
                            .min(self.alignment_end.saturating_sub(frag_len)),
                        self.alignment_end,
                    )
                } else {
                    // Extend rightward, without shrinking a longer read.
                    (
                        self.alignment_start,
                        self.alignment_end.max(self.alignment_start + frag_len),
                    )
                }
            }
        };

        if adjust.center_reads && self.read_length > 0 {
            let center = (start + end) / 2;
            let half = self.read_length / 2;
            let cs = center.saturating_sub(half);
            (cs, cs + self.read_length)
        } else {
            (start, end)
        }
    }
}

/// GC fraction of a sequence in `[0, 1]`, or `None` for an empty sequence.
/// Counts uppercase `G`/`C` only (BAM sequences are upper-cased on decode).
#[inline]
fn gc_fraction(seq: &[u8]) -> Option<f32> {
    if seq.is_empty() {
        return None;
    }
    let gc = seq.iter().filter(|&&b| b == b'G' || b == b'C').count();
    Some(gc as f32 / seq.len() as f32)
}

/// Parse a two-character SAM auxiliary tag string into a [`Tag`].
pub fn parse_tag(tag_str: &str) -> Result<Tag> {
    let bytes = tag_str.as_bytes();
    anyhow::ensure!(
        bytes.len() == 2,
        "tag must be exactly 2 ASCII characters, got: {:?}",
        tag_str
    );
    Ok(Tag::new(bytes[0], bytes[1]))
}

/// Borrow a string-valued auxiliary tag as raw bytes, tied to the record's
/// lifetime. Returns `None` if the tag is absent or not a string.
fn get_tag_bytes<'a>(record: &'a bam::Record, tag: &Tag) -> Result<Option<&'a [u8]>> {
    match record.data().get(tag) {
        Some(Ok(Value::String(v))) => Ok(Some(v.as_ref())),
        Some(Ok(_)) => Ok(None),
        Some(Err(e)) => Err(e).context("failed to decode auxiliary tag"),
        None => Ok(None),
    }
}

/// A minimal single-end, forward-strand `ScRecord` spanning `[start, end)`.
///
/// Test-only: the counting filters and `effective_interval` operate on parsed
/// records, so tests build them directly rather than round-tripping a BAM.
#[cfg(test)]
pub(crate) fn test_record<'a>(start: usize, end: usize) -> ScRecord<'a> {
    ScRecord {
        alignment_start: start,
        alignment_end: end,
        is_reverse: false,
        is_proper_pair: false,
        is_paired: false,
        is_read1: false,
        mate_is_reverse: false,
        template_length: 0,
        mate_alignment_start: None,
        covered_span: None,
        reference_id: None,
        mate_reference_id: None,
        read_length: end - start,
        barcode: None,
        umi: None,
        count: 1,
        gc_content: None,
        aligned_fraction: None,
        read_sequence: None,
    }
}

fn get_count_tag(record: &bam::Record, tag: &Tag) -> Result<Option<u32>> {
    match record.data().get(tag) {
        Some(Ok(value)) => match value {
            Value::String(v) => {
                let s = String::from_utf8_lossy(v.as_ref()).into_owned();
                Ok(Some(s.parse::<u32>().with_context(|| {
                    format!("failed to parse count tag as u32: {:?}", s)
                })?))
            }
            Value::Int8(v) => Ok(Some(v.unsigned_abs() as u32)),
            Value::UInt8(v) => Ok(Some(v as u32)),
            Value::Int16(v) => Ok(Some(v.unsigned_abs() as u32)),
            Value::UInt16(v) => Ok(Some(v as u32)),
            Value::Int32(v) => Ok(Some(v.unsigned_abs())),
            Value::UInt32(v) => Ok(Some(v)),
            _ => Ok(None),
        },
        Some(Err(e)) => Err(e).context("failed to decode count tag"),
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_tag_requires_exactly_two_characters() {
        assert_eq!(parse_tag("CB").unwrap(), Tag::new(b'C', b'B'));
        assert!(parse_tag("C").is_err());
        assert!(parse_tag("CBX").is_err());
        assert!(parse_tag("").is_err());
    }

    #[test]
    fn gc_fraction_counts_g_and_c() {
        assert_eq!(gc_fraction(b"GCGC"), Some(1.0));
        assert_eq!(gc_fraction(b"ATAT"), Some(0.0));
        assert_eq!(gc_fraction(b"GATC"), Some(0.5));
        // N is neither GC nor a reason to change the denominator.
        assert_eq!(gc_fraction(b"GN"), Some(0.5));
        assert_eq!(gc_fraction(b""), None);
    }

    #[test]
    fn without_extension_the_alignment_span_is_used_verbatim() {
        let rec = test_record(100, 150);
        assert_eq!(rec.effective_interval(&AdjustRead::default()), (100, 150));
    }

    #[test]
    fn extend_reads_lengthens_single_end_reads_in_the_read_direction() {
        let adjust = AdjustRead {
            extend_reads: Some(200),
            ..Default::default()
        };

        // Forward: extend downstream from the start.
        let fwd = test_record(1000, 1050);
        assert_eq!(fwd.effective_interval(&adjust), (1000, 1200));

        // Reverse: extend upstream from the end.
        let mut rev = test_record(1000, 1050);
        rev.is_reverse = true;
        assert_eq!(rev.effective_interval(&adjust), (850, 1050));
    }

    #[test]
    fn extend_reads_never_shrinks_a_read_longer_than_the_fragment() {
        // A read already longer than the requested fragment is left unchanged,
        // not truncated to `frag_len` — matching "reads that already exceed the
        // fragment length are not extended".
        let adjust = AdjustRead {
            extend_reads: Some(200),
            ..Default::default()
        };

        // Forward span 300 > 200: keep [1000, 1300), not [1000, 1200).
        let fwd = test_record(1000, 1300);
        assert_eq!(fwd.effective_interval(&adjust), (1000, 1300));

        // Reverse span 300 > 200: keep [1000, 1300), not [1100, 1300).
        let mut rev = test_record(1000, 1300);
        rev.is_reverse = true;
        assert_eq!(rev.effective_interval(&adjust), (1000, 1300));
    }

    #[test]
    fn extend_reads_prefers_the_observed_insert_size_for_proper_pairs() {
        let adjust = AdjustRead {
            extend_reads: Some(200),
            ..Default::default()
        };

        // Forward mate: TLEN wins over the requested fragment length.
        let mut fwd = test_record(1000, 1050);
        fwd.is_proper_pair = true;
        fwd.template_length = 300;
        assert_eq!(fwd.effective_interval(&adjust), (1000, 1300));

        // Reverse mate: the fragment runs from the mate's start to this read's end.
        let mut rev = test_record(1250, 1300);
        rev.is_proper_pair = true;
        rev.is_reverse = true;
        rev.template_length = -300;
        rev.mate_alignment_start = Some(1000);
        assert_eq!(rev.effective_interval(&adjust), (1000, 1300));
    }

    #[test]
    fn implausible_insert_sizes_fall_back_to_the_requested_length() {
        let adjust = AdjustRead {
            extend_reads: Some(200),
            ..Default::default()
        };

        // |TLEN| > 4 × extend_reads is treated as a mapping artifact.
        let mut rec = test_record(1000, 1050);
        rec.is_proper_pair = true;
        rec.template_length = 801;
        assert_eq!(rec.effective_interval(&adjust), (1000, 1200));

        // Exactly 4 × extend_reads is still plausible.
        rec.template_length = 800;
        assert_eq!(rec.effective_interval(&adjust), (1000, 1800));
    }

    #[test]
    fn reverse_read_extension_saturates_at_the_chromosome_start() {
        let adjust = AdjustRead {
            extend_reads: Some(200),
            ..Default::default()
        };
        let mut rec = test_record(10, 60);
        rec.is_reverse = true;
        assert_eq!(rec.effective_interval(&adjust), (0, 60));
    }

    #[test]
    fn center_reads_replaces_the_fragment_with_a_read_length_window() {
        let adjust = AdjustRead {
            center_reads: true,
            ..Default::default()
        };

        // Fragment [100, 300) centered at 200; a 50 bp read spans [175, 225).
        let mut rec = test_record(100, 300);
        rec.read_length = 50;
        assert_eq!(rec.effective_interval(&adjust), (175, 225));
    }

    #[test]
    fn center_reads_applies_after_extension() {
        let adjust = AdjustRead {
            extend_reads: Some(200),
            center_reads: true,
        };

        // Extended to [1000, 1200), centered at 1100; a 50 bp read spans [1075, 1125).
        let mut rec = test_record(1000, 1050);
        rec.read_length = 50;
        assert_eq!(rec.effective_interval(&adjust), (1075, 1125));
    }

    #[test]
    fn center_reads_is_a_no_op_for_zero_length_reads() {
        let adjust = AdjustRead {
            center_reads: true,
            ..Default::default()
        };
        let mut rec = test_record(100, 300);
        rec.read_length = 0;
        assert_eq!(rec.effective_interval(&adjust), (100, 300));
    }

    // Parsing real BAM records

    /// The first records of `test_i1.bam`. Its reads carry `BC` (barcode),
    /// `RX` (UMI) and `NM` (an integer edit distance, reused here as a stand-in
    /// for a numeric count tag).
    fn read_test_bam(n: usize) -> (Header, Vec<bam::Record>) {
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/testdata/test_i1.bam");
        let mut reader = bam::io::reader::Builder.build_from_path(&path).unwrap();
        let header = reader.read_header().unwrap();
        let records = reader.records().take(n).map(|r| r.unwrap()).collect();
        (header, records)
    }

    fn plain_opts() -> ScRecordOptions {
        ScRecordOptions {
            compute_gc: false,
            compute_aligned_fraction: false,
            store_sequence: false,
            compute_covered_span: false,
        }
    }

    #[test]
    fn a_mapped_record_parses_into_coordinates_and_a_barcode() {
        let (header, records) = read_test_bam(1);
        let bc = Tag::new(b'B', b'C');

        let rec = ScRecord::from_bam_record(&records[0], &header, &bc, None, None, &plain_opts())
            .unwrap()
            .expect("the first record of the test BAM is mapped");

        assert!(
            rec.alignment_end > rec.alignment_start,
            "half-open interval"
        );
        assert_eq!(rec.barcode, Some(b"ATATAACT".as_slice()));
        assert_eq!(rec.umi, None, "no UMI tag was requested");
        assert_eq!(rec.count, 1, "no count tag means one read");
        assert!(rec.read_length > 0);

        // Nothing optional was switched on.
        assert_eq!(rec.gc_content, None);
        assert_eq!(rec.aligned_fraction, None);
        assert_eq!(rec.read_sequence, None);
    }

    #[test]
    fn the_umi_is_read_only_when_a_umi_tag_is_asked_for() {
        let (header, records) = read_test_bam(1);
        let bc = Tag::new(b'B', b'C');
        let rx = Tag::new(b'R', b'X');

        let rec =
            ScRecord::from_bam_record(&records[0], &header, &bc, Some(&rx), None, &plain_opts())
                .unwrap()
                .unwrap();
        assert_eq!(rec.umi, Some(b"AGC".as_slice()));
    }

    #[test]
    fn a_tag_the_record_lacks_leaves_the_field_empty_rather_than_failing() {
        let (header, records) = read_test_bam(1);
        let absent = Tag::new(b'Z', b'Z');

        let rec =
            ScRecord::from_bam_record(&records[0], &header, &absent, None, None, &plain_opts())
                .unwrap()
                .unwrap();
        assert_eq!(rec.barcode, None);
    }

    #[test]
    fn gc_content_is_computed_without_keeping_the_sequence() {
        let (header, records) = read_test_bam(1);
        let bc = Tag::new(b'B', b'C');
        let opts = ScRecordOptions {
            compute_gc: true,
            compute_aligned_fraction: false,
            store_sequence: false,
            compute_covered_span: false,
        };

        let rec = ScRecord::from_bam_record(&records[0], &header, &bc, None, None, &opts)
            .unwrap()
            .unwrap();

        let gc = rec.gc_content.expect("gc was requested");
        assert!((0.0..=1.0).contains(&gc), "gc out of range: {gc}");
        assert_eq!(rec.read_sequence, None, "the sequence was not asked for");
    }

    #[test]
    fn storing_the_sequence_yields_one_base_per_read_position() {
        let (header, records) = read_test_bam(1);
        let bc = Tag::new(b'B', b'C');
        let opts = ScRecordOptions {
            compute_gc: true,
            compute_aligned_fraction: false,
            store_sequence: true,
            compute_covered_span: false,
        };

        let rec = ScRecord::from_bam_record(&records[0], &header, &bc, None, None, &opts)
            .unwrap()
            .unwrap();

        let seq = rec
            .read_sequence
            .as_ref()
            .expect("the sequence was requested");
        assert_eq!(seq.len(), rec.read_length);
        assert!(
            rec.gc_content.is_some(),
            "gc is derived from the stored sequence"
        );
    }

    #[test]
    fn the_aligned_fraction_is_matches_over_the_read_consuming_length() {
        let (header, records) = read_test_bam(1);
        let bc = Tag::new(b'B', b'C');
        let opts = ScRecordOptions {
            compute_gc: false,
            compute_aligned_fraction: true,
            store_sequence: false,
            compute_covered_span: false,
        };

        let rec = ScRecord::from_bam_record(&records[0], &header, &bc, None, None, &opts)
            .unwrap()
            .unwrap();

        let fraction = rec
            .aligned_fraction
            .expect("the aligned fraction was requested");
        assert!((0.0..=1.0).contains(&fraction), "out of range: {fraction}");
    }

    #[test]
    fn every_record_in_the_test_file_either_parses_or_is_filtered_out() {
        let (header, records) = read_test_bam(200);
        let bc = Tag::new(b'B', b'C');
        assert!(!records.is_empty());

        for record in &records {
            // Never an Err: the file is well formed, so each record is either
            // a usable ScRecord or a deliberate skip.
            let parsed =
                ScRecord::from_bam_record(record, &header, &bc, None, None, &plain_opts()).unwrap();
            if let Some(rec) = parsed {
                assert!(rec.alignment_end > rec.alignment_start);
            }
        }
    }

    // Auxiliary tag helpers

    #[test]
    fn string_tags_are_borrowed_and_other_types_are_ignored() {
        let (_header, records) = read_test_bam(1);
        let record = &records[0];

        assert_eq!(
            get_tag_bytes(record, &Tag::new(b'B', b'C')).unwrap(),
            Some(b"ATATAACT".as_slice())
        );
        // NM is an integer, so it is not a byte string.
        assert_eq!(get_tag_bytes(record, &Tag::new(b'N', b'M')).unwrap(), None);
        // A tag the record does not carry at all.
        assert_eq!(get_tag_bytes(record, &Tag::new(b'Z', b'Z')).unwrap(), None);
    }

    #[test]
    fn an_integer_count_tag_is_read_as_a_number() {
        let (_header, records) = read_test_bam(1);
        // NM:i:1 on the first record of the test file.
        let count = get_count_tag(&records[0], &Tag::new(b'N', b'M')).unwrap();
        assert_eq!(count, Some(1));
    }

    #[test]
    fn a_count_tag_that_is_not_a_number_is_an_error_naming_the_value() {
        let (_header, records) = read_test_bam(1);
        // BC holds a barcode, which cannot be parsed as a count.
        let err = get_count_tag(&records[0], &Tag::new(b'B', b'C'))
            .unwrap_err()
            .to_string();
        assert!(err.contains("count tag"), "{err}");
    }

    #[test]
    fn a_missing_count_tag_is_absent_rather_than_an_error() {
        let (_header, records) = read_test_bam(1);
        assert_eq!(
            get_count_tag(&records[0], &Tag::new(b'Z', b'Z')).unwrap(),
            None
        );
    }
}
