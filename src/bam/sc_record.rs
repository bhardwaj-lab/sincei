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
    /// Split the alignment into its gapless blocks, so a spliced read credits
    /// its exons and not the intron between them.
    pub(crate) compute_blocks: bool,
}

/// How a read's position can be adjusted.
#[derive(Clone, Copy, Debug, Default)]
pub struct AdjustRead {
    /// Extend each read to this fragment length (bp). A pair that is properly
    /// paired, on one reference, facing inward and no longer than
    /// [`Self::max_paired_fragment_length`] is extended to its observed insert
    /// size (TLEN) instead; any other read is extended to this length.
    /// `None` leaves the alignment span unchanged.
    pub extend_reads: Option<usize>,
    /// After extension, replace the interval with a `read_length` window
    /// centered on its midpoint.
    pub center_reads: bool,
    /// Largest insert size that may be accepted as a fragment. `None` means
    /// `4 × extend_reads`. If`--maxFragmentLength` is given, that value
    /// replaces it.
    pub max_paired_fragment_length: Option<usize>,
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
    /// Group tag bytes (the read's sample of origin in a merged BAM), borrowed
    /// from the record. `None` unless a group tag was asked for.
    pub group: Option<&'a [u8]>,
    /// Value to add to the count matrix (defaults to 1).
    pub count: u32,
    /// GC fraction in `[0, 1]`. `None` if not requested.
    pub gc_content: Option<f32>,
    /// Fraction of read bases in aligned CIGAR ops (`M`, `=`, `X`).
    /// `None` if not requested.
    pub aligned_fraction: Option<f32>,
    /// Raw BAM read sequence (read orientation). `None` if not requested.
    pub read_sequence: Option<Vec<u8>>,
    /// The alignment's gapless blocks, `None` unless it has more than one.
    ///
    /// A read without a `N` or `D` in its CIGAR covers one unbroken interval,
    /// which `alignment_start .. alignment_end` already describes, so the
    /// common case stores nothing and allocates nothing.
    pub blocks: Option<Vec<(usize, usize)>>,
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
        group_tag: Option<&Tag>,
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

        let group = match group_tag {
            Some(gtag) => get_tag_bytes(record, gtag)?,
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

        // Gapless blocks, the reference intervals the read actually covers.
        // `D` breaks a block as `N` does, matching pysam's `get_blocks`. An
        // insertion consumes no reference, so it leaves the block it sits in
        // unbroken.
        //
        // Almost every read is one block, and that block is the alignment span
        // the record already carries. So the first block is held in a local and
        // `later` stays empty (and unallocated) until a second one turns up:
        // an ungapped read walks its CIGAR and allocates nothing.
        //
        // A single-operation CIGAR cannot hold a gap *between* aligned bases, so
        // it is not walked at all. That is 90% of the reads in a typical file.
        let blocks = if opts.compute_blocks && record.cigar().len() > 1 {
            let mut first: Option<(usize, usize)> = None;
            let mut later: Vec<(usize, usize)> = Vec::new();
            let mut pos = start;
            let mut open: Option<usize> = None;

            for result in record.cigar().iter() {
                let op = result.context("failed to decode CIGAR op")?;
                match op.kind() {
                    CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch => {
                        open.get_or_insert(pos);
                        pos += op.len();
                    }
                    CigarKind::Deletion | CigarKind::Skip => {
                        if let Some(block_start) = open.take() {
                            match first {
                                None => first = Some((block_start, pos)),
                                Some(_) => later.push((block_start, pos)),
                            }
                        }
                        pos += op.len();
                    }
                    // Consumes no reference, so it cannot end a block.
                    _ => {}
                }
            }
            if let Some(block_start) = open {
                match first {
                    None => first = Some((block_start, pos)),
                    Some(_) => later.push((block_start, pos)),
                }
            }

            first.filter(|_| !later.is_empty()).map(|first| {
                let mut all = Vec::with_capacity(later.len() + 1);
                all.push(first);
                all.append(&mut later);
                all
            })
        } else {
            None
        };

        // Aligned fraction: the aligned operations over the length of the read
        // as it was sequenced.
        //
        // `=` and `X` are counted alongside `M`: all three consume a reference
        // base for a read base, and an aligner that spells its matches out
        // (minimap2 --eqx, bwa-mem2 -M alternatives) describes the same
        // alignment as one that emits `M`. Counting `M` alone would score every
        // read of such a file at 0. The rest of this module already treats the
        // three the same way.
        //
        // The denominator counts hard-clipped bases as well as the ones the
        // record still carries. Those bases were part of the read, so leaving
        // them out would overstate how much of it aligned. Soft clips are
        // already counted by `consumes_read`.
        let aligned_fraction = if opts.compute_aligned_fraction {
            let mut aligned_len: usize = 0;
            let mut read_consuming: usize = 0;
            for result in record.cigar().iter() {
                let op = result.context("failed to decode CIGAR op")?;
                if matches!(
                    op.kind(),
                    CigarKind::Match | CigarKind::SequenceMatch | CigarKind::SequenceMismatch
                ) {
                    aligned_len += op.len();
                }
                if op.kind().consumes_read() || matches!(op.kind(), CigarKind::HardClip) {
                    read_consuming += op.len();
                }
            }
            if read_consuming > 0 {
                Some(aligned_len as f32 / read_consuming as f32)
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
            group,
            count,
            gc_content,
            aligned_fraction,
            read_sequence,
            blocks,
        }))
    }

    /// The genomic intervals to credit this read to.
    ///
    /// Normally the one interval [`Self::effective_interval`] returns. A spliced
    /// read counted on its own alignment yields one interval per gapless block
    /// instead, so its intron collects no signal. Extension and centering both
    /// replace the alignment with a synthetic interval that has no blocks to
    /// speak of, so they always yield a single interval.
    ///
    /// Callers must credit a bin or feature **once per read**, not once per
    /// interval: two blocks of one read can reach the same bin.
    pub(crate) fn effective_intervals(&self, adjust: &AdjustRead) -> EffectiveIntervals<'_> {
        match &self.blocks {
            Some(blocks) if adjust.extend_reads.is_none() && !adjust.center_reads => {
                EffectiveIntervals::Blocks(blocks.iter())
            }
            _ => EffectiveIntervals::One(Some(self.effective_interval(adjust))),
        }
    }

    /// Genomic interval to credit this read to, after applying the extension
    /// and centering in `adjust`.
    pub(crate) fn effective_interval(&self, adjust: &AdjustRead) -> (usize, usize) {
        let (start, end) = match adjust.extend_reads {
            None => (self.alignment_start, self.alignment_end),
            Some(frag_len) => {
                // Four times the requested fragment length, unless a maximum was
                // given.
                let max_paired = adjust.max_paired_fragment_length.unwrap_or(4 * frag_len);
                let tlen_abs = self.template_length.unsigned_abs() as usize;
                // The proper-pair flag alone cannot be trusted to mean that TLEN
                // describes a fragment, so the pair is checked for: one reference,
                // opposite strands, facing each other, and an insert size within
                // `max_paired`. A pair failing any of those is treated as
                // single-ended and extended by a fixed length.
                //
                // The orientation test is what keeps the interval sound. An
                // outward-facing pair puts a reverse read's mate *after* the
                // read, and taking the fragment from the mate would then produce
                // an interval that starts after it ends.
                let fragment_start = self.mate_alignment_start.filter(|&mate| {
                    self.is_proper_pair
                        && self.reference_id == self.mate_reference_id
                        && self.is_reverse != self.mate_is_reverse
                        && tlen_abs > 0
                        && tlen_abs <= max_paired
                        && if self.is_reverse {
                            mate <= self.alignment_start
                        } else {
                            self.alignment_start <= mate
                        }
                });

                if let Some(mate) = fragment_start {
                    if self.is_reverse {
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

        // Every branch above keeps the bounds in order, as the reference
        // implementation asserts of its own fragment. An inverted interval would
        // not be caught downstream: the bin loop would silently drop the read
        // and the feature loop would underflow its overlap subtraction.
        debug_assert!(start < end, "inverted interval {start}..{end}");

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

/// The intervals of one read: either a single span or its gapless blocks.
///
/// An enum rather than a boxed iterator so the common single-interval case
/// costs no allocation and no indirect call.
pub(crate) enum EffectiveIntervals<'a> {
    One(Option<(usize, usize)>),
    Blocks(std::slice::Iter<'a, (usize, usize)>),
    /// A window selected inside the read, which the record itself does not
    /// hold: `--Offset` trims the blocks at both ends, and a window crossing a
    /// gap comes back as several intervals. Built per read, so it owns them.
    Selected(std::vec::IntoIter<(usize, usize)>),
}

impl Iterator for EffectiveIntervals<'_> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::One(interval) => interval.take(),
            Self::Blocks(blocks) => blocks.next().copied(),
            Self::Selected(intervals) => intervals.next(),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let n = self.len();
        (n, Some(n))
    }
}

impl ExactSizeIterator for EffectiveIntervals<'_> {
    fn len(&self) -> usize {
        match self {
            Self::One(interval) => usize::from(interval.is_some()),
            Self::Blocks(blocks) => blocks.len(),
            Self::Selected(intervals) => intervals.len(),
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
        group: None,
        count: 1,
        gc_content: None,
        aligned_fraction: None,
        read_sequence: None,
        blocks: None,
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
        // not truncated to `frag_len`. Matching "reads that already exceed the
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

    /// A record whose pair passes every proper-pair test: one reference,
    /// opposite strands and facing inward.
    ///
    /// The sign of `tlen` says which end of the fragment this read is. A
    /// positive one makes it the leftmost, forward read of the pair, and a
    /// negative one the rightmost, reverse read, which is the orientation a
    /// library is expected to have.
    fn proper_pair<'a>(start: usize, end: usize, tlen: i64, mate_start: usize) -> ScRecord<'a> {
        let mut rec = test_record(start, end);
        rec.is_paired = true;
        rec.is_proper_pair = true;
        rec.template_length = tlen;
        rec.mate_alignment_start = Some(mate_start);
        rec.reference_id = Some(0);
        rec.mate_reference_id = Some(0);
        rec.is_reverse = tlen < 0;
        rec.mate_is_reverse = tlen >= 0;
        rec
    }

    /// `--extendReads 200`, the setting every extension test below uses.
    fn extend_to_200() -> AdjustRead {
        AdjustRead {
            extend_reads: Some(200),
            ..Default::default()
        }
    }

    #[test]
    fn extend_reads_prefers_the_observed_insert_size_for_proper_pairs() {
        let adjust = extend_to_200();

        // Forward mate: TLEN wins over the requested fragment length.
        let fwd = proper_pair(1000, 1050, 300, 1250);
        assert_eq!(fwd.effective_interval(&adjust), (1000, 1300));

        // Reverse mate: the fragment runs from the mate's start to this read's end.
        let rev = proper_pair(1250, 1300, -300, 1000);
        assert_eq!(rev.effective_interval(&adjust), (1000, 1300));
    }

    #[test]
    fn implausible_insert_sizes_fall_back_to_the_requested_length() {
        let adjust = extend_to_200();

        // |TLEN| > 4 × extend_reads is treated as a mapping artifact.
        let mut rec = proper_pair(1000, 1050, 801, 1750);
        assert_eq!(rec.effective_interval(&adjust), (1000, 1200));

        // Exactly 4 × extend_reads is still plausible.
        rec.template_length = 800;
        assert_eq!(rec.effective_interval(&adjust), (1000, 1800));
    }

    #[test]
    fn a_maximum_fragment_length_replaces_the_insert_size_bound() {
        // Without a maximum the bound is 4 x the requested fragment length, so a
        // 700 bp insert is believed at `--extendReads 200`.
        let rec = proper_pair(1000, 1050, 700, 1650);
        assert_eq!(rec.effective_interval(&extend_to_200()), (1000, 1700));

        // `--maxFragmentLength` replaces that bound outright, as it does in the
        // reference implementation, so the same pair is now too long to believe
        // and falls back to the fixed extension.
        let bounded = AdjustRead {
            extend_reads: Some(200),
            center_reads: false,
            max_paired_fragment_length: Some(500),
        };
        assert_eq!(rec.effective_interval(&bounded), (1000, 1200));

        // It replaces rather than tightens: a bound above 4x widens what counts.
        let widened = AdjustRead {
            extend_reads: Some(200),
            center_reads: false,
            max_paired_fragment_length: Some(2000),
        };
        let long = proper_pair(1000, 1050, 1500, 2450);
        assert_eq!(long.effective_interval(&widened), (1000, 2500));
        assert_eq!(long.effective_interval(&extend_to_200()), (1000, 1200));
    }

    // The proper-pair flag is not enough on its own: a pair failing any of the
    // tests below is extended as if it were single-ended (1000..1050 + 200).

    #[test]
    fn a_mate_on_another_reference_is_not_a_fragment() {
        let mut rec = proper_pair(1000, 1050, 300, 1250);
        rec.mate_reference_id = Some(1);

        assert_eq!(rec.effective_interval(&extend_to_200()), (1000, 1200));
    }

    #[test]
    fn a_mate_on_the_same_strand_is_not_a_fragment() {
        let mut rec = proper_pair(1000, 1050, 300, 1250);
        rec.mate_is_reverse = false;

        assert_eq!(rec.effective_interval(&extend_to_200()), (1000, 1200));
    }

    #[test]
    fn a_mate_without_a_recorded_position_is_not_a_fragment() {
        let mut rec = proper_pair(1000, 1050, 300, 1250);
        rec.mate_alignment_start = None;

        assert_eq!(rec.effective_interval(&extend_to_200()), (1000, 1200));
    }

    #[test]
    fn an_outward_facing_pair_is_not_a_fragment() {
        let adjust = extend_to_200();

        // A forward read whose mate lies behind it: the fragment would run
        // backwards from the read.
        let fwd = proper_pair(2000, 2050, 300, 1000);
        assert_eq!(fwd.effective_interval(&adjust), (2000, 2200));

        // The reverse read of such a pair is the dangerous one. Its fragment
        // start would be the mate at 2000, past its own end at 1300, giving an
        // interval that starts after it ends.
        let rev = proper_pair(1250, 1300, -300, 2000);
        let (start, end) = rev.effective_interval(&adjust);
        assert!(start < end, "inverted interval {start}..{end}");
        assert_eq!((start, end), (1100, 1300));
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
            ..AdjustRead::default()
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
        read_bam("test_i1.bam", n)
    }

    fn read_bam(name: &str, n: usize) -> (Header, Vec<bam::Record>) {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/testdata")
            .join(name);
        let mut reader = bam::io::reader::Builder.build_from_path(&path).unwrap();
        let header = reader.read_header().unwrap();
        let records = reader.records().take(n).map(|r| r.unwrap()).collect();
        (header, records)
    }

    /// Every read of `test_cigars.bam`, keyed by the CIGAR shape its name gives.
    fn cigar_shapes(opts: &ScRecordOptions) -> std::collections::BTreeMap<String, ScRecord<'_>> {
        let (header, records) = read_bam("test_cigars.bam", 100);
        let (header, records) = (Box::leak(Box::new(header)), Box::leak(Box::new(records)));
        let bc = Tag::new(b'B', b'C');

        records
            .iter()
            .map(|record| {
                let name = String::from_utf8_lossy(record.name().unwrap().as_ref()).into_owned();
                let rec = ScRecord::from_bam_record(record, header, &bc, None, None, None, opts)
                    .unwrap()
                    .unwrap();
                (name, rec)
            })
            .collect()
    }

    fn plain_opts() -> ScRecordOptions {
        ScRecordOptions {
            compute_gc: false,
            compute_aligned_fraction: false,
            store_sequence: false,
            compute_covered_span: false,
            compute_blocks: false,
        }
    }

    #[test]
    fn a_mapped_record_parses_into_coordinates_and_a_barcode() {
        let (header, records) = read_test_bam(1);
        let bc = Tag::new(b'B', b'C');

        let rec =
            ScRecord::from_bam_record(&records[0], &header, &bc, None, None, None, &plain_opts())
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

        let rec = ScRecord::from_bam_record(
            &records[0],
            &header,
            &bc,
            Some(&rx),
            None,
            None,
            &plain_opts(),
        )
        .unwrap()
        .unwrap();
        assert_eq!(rec.umi, Some(b"AGC".as_slice()));
    }

    #[test]
    fn a_tag_the_record_lacks_leaves_the_field_empty_rather_than_failing() {
        let (header, records) = read_test_bam(1);
        let absent = Tag::new(b'Z', b'Z');

        let rec = ScRecord::from_bam_record(
            &records[0],
            &header,
            &absent,
            None,
            None,
            None,
            &plain_opts(),
        )
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
            compute_blocks: false,
        };

        let rec = ScRecord::from_bam_record(&records[0], &header, &bc, None, None, None, &opts)
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
            compute_blocks: false,
        };

        let rec = ScRecord::from_bam_record(&records[0], &header, &bc, None, None, None, &opts)
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
            compute_blocks: false,
        };

        let rec = ScRecord::from_bam_record(&records[0], &header, &bc, None, None, None, &opts)
            .unwrap()
            .unwrap();

        let fraction = rec
            .aligned_fraction
            .expect("the aligned fraction was requested");
        assert!((0.0..=1.0).contains(&fraction), "out of range: {fraction}");
    }

    #[test]
    fn the_aligned_fraction_divides_by_the_length_of_the_sequenced_read() {
        // The denominator is the read as it was sequenced, so hard-clipped
        // bases count towards it exactly as soft-clipped ones do. Leaving them
        // out would report a hard-clipped read as fully aligned.
        let opts = ScRecordOptions {
            compute_gc: false,
            compute_aligned_fraction: true,
            store_sequence: false,
            compute_covered_span: false,
            compute_blocks: false,
        };
        let shapes = cigar_shapes(&opts);
        let fraction = |name: &str| shapes[name].aligned_fraction.unwrap();

        // 10H10M and 5S10M5S: 10 aligned bases of a 20-base read.
        assert_eq!(fraction("hardclip"), 0.5);
        assert_eq!(fraction("softclip"), 0.5);

        // 10M5I10M: the insertion is part of the read and is not an M.
        assert_eq!(fraction("ins"), 20.0 / 25.0);

        // Deletions and introns consume no read base, so they do not dilute it.
        assert_eq!(fraction("del2"), 1.0);
        assert_eq!(fraction("intron"), 1.0);

        // `=` and `X` are aligned ops like `M`, so 5=5X10N5= is 15 aligned
        // bases of a 15-base read. The intron consumes none of it.
        assert_eq!(fraction("eqx"), 1.0);
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
                ScRecord::from_bam_record(record, &header, &bc, None, None, None, &plain_opts())
                    .unwrap();
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
