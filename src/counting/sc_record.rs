use anyhow::{Context, Result};
use noodles::bam;
use noodles::sam::Header;
use noodles::sam::alignment::Record as _;
use noodles::sam::alignment::record::cigar::op::Kind as CigarKind;
use noodles::sam::alignment::record::data::field::{Tag, Value};

use super::params::CountingParams;

/// Controls which optional fields are computed when building an [`ScRecord`].
pub(crate) struct ScRecordOptions {
    pub(crate) compute_gc: bool,
    pub(crate) compute_aligned_fraction: bool,
    pub(crate) store_sequence: bool,
}

/// A parsed, filter-passing BAM record carrying the fields needed for
/// single-cell feature counting.
///
/// Coordinates are 0-based half-open `[start, end)`.
///
/// The chromosome is intentionally *not* stored: every caller queries one
/// chromosome at a time, so the chromosome name is already known from the
/// work chunk and need not be re-derived (and re-allocated) per record.
pub struct ScRecord<'a> {
    pub alignment_start: usize,
    pub alignment_end: usize,
    /// Whether the read is on the reverse strand.
    pub is_reverse: bool,
    /// SAM proper-pair flag (FLAG 0x2).
    pub is_proper_pair: bool,
    /// True if the read is part of a paired-end library (FLAG 0x1).
    pub is_paired: bool,
    /// True if this is the first read in the pair (FLAG 0x40).
    pub is_read1: bool,
    /// True if the mate maps to the reverse strand (FLAG 0x20).
    pub mate_is_reverse: bool,
    /// Signed template length (TLEN); 0 for single-end reads.
    pub template_length: i64,
    /// Mate's 0-based alignment start, if available.
    pub next_alignment_start: Option<usize>,
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

        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
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

        // Mate position (1-based → 0-based).
        let next_alignment_start = record
            .mate_alignment_start()
            .map(|res| res.map(|pos| pos.get().saturating_sub(1)))
            .transpose()
            .context("failed to read mate alignment start")?;

        // Borrow the barcode/UMI bytes directly from the record — no per-read
        // allocation.  Callers look these up against a byte-keyed whitelist and
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

        // Sequence handling.  Only materialize a `Vec` when the motif filter
        // needs the full read sequence (`store_sequence`).  When only GC is
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
            next_alignment_start,
            read_length,
            barcode,
            umi,
            count,
            gc_content,
            aligned_fraction,
            read_sequence,
        }))
    }

    /// Genomic interval to credit this read to, after applying `extend_reads`
    /// and `center_reads` from `params`.
    pub(crate) fn effective_interval(&self, params: &CountingParams) -> (usize, usize) {
        let (start, end) = match params.extend_reads {
            None => (self.alignment_start, self.alignment_end),
            Some(frag_len) => {
                let max_paired = 4 * frag_len;
                let tlen_abs = self.template_length.unsigned_abs() as usize;
                let use_tlen = self.is_proper_pair && tlen_abs > 0 && tlen_abs <= max_paired;

                if use_tlen {
                    if self.is_reverse {
                        let mate = self.next_alignment_start.unwrap_or(self.alignment_start);
                        (mate, self.alignment_end)
                    } else {
                        (self.alignment_start, self.alignment_start + tlen_abs)
                    }
                } else if self.is_reverse {
                    (
                        self.alignment_end.saturating_sub(frag_len),
                        self.alignment_end,
                    )
                } else {
                    (self.alignment_start, self.alignment_start + frag_len)
                }
            }
        };

        if params.center_reads && self.read_length > 0 {
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
/// lifetime.  Returns `None` if the tag is absent or not a string.
fn get_tag_bytes<'a>(record: &'a bam::Record, tag: &Tag) -> Result<Option<&'a [u8]>> {
    match record.data().get(tag) {
        Some(Ok(Value::String(v))) => Ok(Some(v.as_ref())),
        Some(Ok(_)) => Ok(None),
        Some(Err(e)) => Err(e).context("failed to decode auxiliary tag"),
        None => Ok(None),
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
