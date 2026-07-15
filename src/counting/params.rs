use std::path::PathBuf;

/// Options that are common to all BAM counting functions.
///
/// These control pre-counting filtering and fragment-interval adjustments that
/// are independent of the per-record QC/duplicate/motif filters.
pub struct CountingParams {
    /// Chromosomes to exclude entirely (e.g. `["chrM", "chrUn"]`).
    pub chr_to_skip: Vec<String>,

    /// Restrict counting to a single genomic region.
    /// Format: `"chrom:start-end"` (1-based inclusive, samtools style) or
    /// just `"chrom"` for a whole chromosome.
    pub region: Option<String>,

    /// BED file whose regions define a blacklist. Reads whose effective
    /// interval overlaps any blacklisted region are discarded.
    pub blacklist_path: Option<PathBuf>,

    /// Extend reads to this fragment length.
    ///
    /// For properly paired reads the insert size from TLEN is used instead
    /// (up to `4 × extend_reads`); the supplied value is only applied to
    /// single-end or improperly paired reads.
    pub extend_reads: Option<usize>,

    /// After computing the extended interval, replace it with a window of
    /// `read_length` bp centered on the fragment midpoint. This can be useful
    /// to get sharper signal around enriched regions.
    pub center_reads: bool,

    /// For GTF / GFF3 annotation files: the column-3 feature types that define
    /// a counting region. A record is a region when its type is any of them,
    /// so an Ensembl-style GFF3 can contribute several transcript biotypes
    /// (`mRNA`, `lnc_RNA`, `tRNA`, …). `None` means `"transcript"` for GTF and,
    /// for GFF3, the transcript types read out of the file itself, so every
    /// biotype is kept. Ignored for BED files.
    pub feature_type: Option<Vec<String>>,

    /// For GTF / GFF3 annotation files: the column-3 feature types that define
    /// an exon (e.g. `["exon"]`, or `["CDS"]` to count coding sequence only).
    /// `None` selects `DEFAULT_EXON_TYPES`. Only used when `metagene` is true.
    /// Ignored for BED files.
    pub exon_type: Option<Vec<String>>,

    /// For GTF / GFF3 annotation files: the column-9 attribute key whose value
    /// names the feature. `None` selects the per-format default
    /// (`"transcript_id"` for GTF, `"ID"` for GFF3, or `"gene_id"` for both
    /// when `metagene` is set). Ignored for BED files (uses the 4th column
    /// instead). Falls back to `"chrom:start-end"` when the attribute is absent.
    pub name_attr: Option<String>,

    /// When `true`, count reads per gene rather than per transcript.
    ///
    /// Exon intervals (type determined by `exon_type`, default `"exon"`) are
    /// grouped per gene: a read overlapping several exons of one gene is
    /// counted once for that gene. The `feature_type` field is ignored in
    /// metagene mode.
    ///
    /// An exon's gene is taken from its `gene_id` attribute or, for the
    /// Ensembl / RefSeq GFF3 style, whose exons name only their transcript in
    /// `Parent`, by following that `Parent` to the transcript record and on to
    /// its gene. Setting `name_attr` overrides this and groups on that
    /// attribute alone.
    pub metagene: bool,
}

impl CountingParams {
    pub fn new() -> Self {
        Self {
            chr_to_skip: Vec::new(),
            region: None,
            blacklist_path: None,
            extend_reads: None,
            center_reads: false,
            feature_type: None,
            exon_type: None,
            name_attr: None,
            metagene: false,
        }
    }
}

impl Default for CountingParams {
    fn default() -> Self {
        Self::new()
    }
}

/// Parse a region string into `(chrom, start, end)` using 0-based half-open
/// coordinates.
///
/// Accepts `"chrom"`, `"chrom:start"`, or `"chrom:start:end"`. The coordinate
/// separator may be `':'` (the form the CLI normalizes to) or `'-'` (samtools
/// style). Input coordinates are 1-based inclusive; the returned start is
/// 0-based and the end is exclusive. A missing start defaults to 0 and a
/// missing end to `usize::MAX`.
pub(super) fn parse_region(region: &str) -> anyhow::Result<(String, usize, usize)> {
    let Some((chrom, rest)) = region.split_once(':') else {
        // Bare chromosome (note: chromosome names may themselves contain '-',
        // so we only split on the first ':').
        return Ok((region.to_string(), 0, usize::MAX));
    };

    let mut coords = rest.splitn(2, [':', '-']);
    let start: usize = coords
        .next()
        .unwrap_or("")
        .parse::<usize>()
        .map_err(|_| anyhow::anyhow!("invalid region start in {:?}", region))?;
    let end: usize = match coords.next() {
        Some(end_str) => end_str
            .parse()
            .map_err(|_| anyhow::anyhow!("invalid region end in {:?}", region))?,
        None => usize::MAX,
    };
    // Convert 1-based inclusive -> 0-based half-open.
    Ok((chrom.to_string(), start.saturating_sub(1), end))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bare_chromosome_spans_everything() {
        assert_eq!(
            parse_region("chr1").unwrap(),
            ("chr1".to_string(), 0, usize::MAX)
        );
    }

    #[test]
    fn start_only_is_open_ended() {
        assert_eq!(
            parse_region("chr1:100").unwrap(),
            ("chr1".to_string(), 99, usize::MAX)
        );
    }

    #[test]
    fn both_coordinate_separators_are_accepted() {
        // samtools style and the CLI-normalized style must agree.
        let dash = parse_region("chr1:100-200").unwrap();
        let colon = parse_region("chr1:100:200").unwrap();
        assert_eq!(dash, ("chr1".to_string(), 99, 200));
        assert_eq!(dash, colon);
    }

    #[test]
    fn chromosome_names_may_contain_dashes() {
        // Only the first ':' splits the chromosome off, so a '-' in the name
        // must not be mistaken for a coordinate separator.
        assert_eq!(
            parse_region("chr1-alt").unwrap(),
            ("chr1-alt".to_string(), 0, usize::MAX)
        );
        assert_eq!(
            parse_region("chr1-alt:5-10").unwrap(),
            ("chr1-alt".to_string(), 4, 10)
        );
    }

    #[test]
    fn start_of_zero_does_not_underflow() {
        assert_eq!(parse_region("chr1:0-10").unwrap().1, 0);
    }

    #[test]
    fn non_numeric_coordinates_are_rejected() {
        assert!(parse_region("chr1:abc").is_err());
        assert!(parse_region("chr1:100-xyz").is_err());
        assert!(parse_region("chr1:").is_err());
    }
}
