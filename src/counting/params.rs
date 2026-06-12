use std::path::PathBuf;

/// Options that are common to all BAM counting functions.
///
/// These control pre-counting filtering and fragment-interval adjustments that
/// are independent of the per-record QC/duplicate/motif filters.
pub struct CountingParams {
    /// Only count reads for which `flags & include == include`.
    /// Handled inside [`ScRecord::from_bam_record`]; does not require extra
    /// fields on the record.
    pub sam_flag_include: Option<u16>,

    /// Skip reads for which `flags & exclude != 0`.
    pub sam_flag_exclude: Option<u16>,

    /// Chromosomes to exclude entirely (e.g. `["chrM", "chrUn"]`).
    pub chr_to_skip: Vec<String>,

    /// Restrict counting to a single genomic region.
    /// Format: `"chrom:start-end"` (1-based inclusive, samtools style) or
    /// just `"chrom"` for a whole chromosome.
    pub region: Option<String>,

    /// BED file whose regions define a blacklist.  Reads whose effective
    /// interval overlaps any blacklisted region are discarded.
    pub blacklist_path: Option<PathBuf>,

    /// Extend reads to this fragment length.
    ///
    /// For proper paired-end reads the insert size from TLEN is used instead
    /// (up to `4 × extend_reads`); the supplied value is only applied to
    /// single-end or improper-pair reads.
    ///
    /// Mirrors `extendReads` / `defaultFragmentLength` in `ReadCounter.py`.
    pub extend_reads: Option<usize>,

    /// After computing the extended interval, replace it with a window of
    /// `read_length` bp centered on the fragment midpoint.
    ///
    /// Mirrors `center_read` in `ReadCounter.py`.
    pub center_reads: bool,

    /// For GTF / GFF3 annotation files: the column-3 feature type that defines
    /// a counting region.  `None` selects the per-format default
    /// (`"transcript"`).  Ignored for BED files.
    pub feature_type: Option<String>,

    /// For GTF / GFF3 annotation files: the column-3 feature type that defines
    /// an exon.  `None` selects the per-format default (`"exon"`).  Only used
    /// in metagene mode (not yet implemented); accepted and stored now so the
    /// surface is stable.  Ignored for BED files.
    pub exon_type: Option<String>,

    /// For GTF / GFF3 annotation files: the column-9 attribute key whose value
    /// names the feature.  `None` selects the per-format default
    /// (`"transcript_id"` for GTF, `"ID"` for GFF3).  Falls back to
    /// `"chrom:start-end"` when the attribute is absent. Ignored for BED files
    /// (uses the 4th column instead).
    pub name_attr: Option<String>,
}

impl CountingParams {
    pub fn new() -> Self {
        Self {
            sam_flag_include: None,
            sam_flag_exclude: None,
            chr_to_skip: Vec::new(),
            region: None,
            blacklist_path: None,
            extend_reads: None,
            center_reads: false,
            feature_type: None,
            exon_type: None,
            name_attr: None,
        }
    }
}

impl Default for CountingParams {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse a region string (`"chrom:start-end"` or `"chrom"`) into
/// `(chrom, start, end)` using 0-based half-open coordinates.
///
/// Input coordinates are treated as 1-based inclusive (samtools style).
pub(super) fn parse_region(region: &str) -> anyhow::Result<(String, usize, usize)> {
    if let Some((chrom, rest)) = region.split_once(':') {
        if let Some((start_str, end_str)) = rest.split_once('-') {
            let start: usize = start_str
                .parse()
                .map_err(|_| anyhow::anyhow!("invalid region start in {:?}", region))?;
            let end: usize = end_str
                .parse()
                .map_err(|_| anyhow::anyhow!("invalid region end in {:?}", region))?;
            // Convert 1-based inclusive → 0-based half-open.
            Ok((chrom.to_string(), start.saturating_sub(1), end))
        } else {
            Err(anyhow::anyhow!(
                "region {:?} must be 'chrom:start-end' or 'chrom'",
                region
            ))
        }
    } else {
        Ok((region.to_string(), 0, usize::MAX))
    }
}
