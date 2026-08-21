//! Read BED / GTF / GFF3 annotations into searchable regions.
//!
//! [`parse_annotation_files`] is the entry point: it picks a parser from each
//! file's extension, decompresses gzip/BGZF, and merges every file into a
//! [`GenomeIndex`] to search the ordered [`Feature`]s that correspond to each
//! of a cell × region count matrix's columns.
//! [`parse_blacklist_bed`] is the same machinery for genomic regions that are
//! only ever excluded from counting.

use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use ahash::AHashMap;
use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use noodles::{bed, gff, gtf};

use super::region_index::{ChromIndex, Feature, GenomeIndex, Interval};

/// Open an annotation file, decompressing it when it is gzip- or BGZF-compressed.
///
/// Compression is detected from the gzip magic bytes (`1f 8b`).
/// `MultiGzDecoder` is used (not `GzDecoder`) because it works with both gzip
/// and BGZF compression.
fn open_annotation(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open annotation file: {}", path.display()))?;
    let mut reader = BufReader::new(file);

    // Check magic bytes.
    let is_gzip = {
        let header = reader
            .fill_buf()
            .with_context(|| format!("failed to read: {}", path.display()))?;
        header.starts_with(&[0x1f, 0x8b])
    };

    if is_gzip {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
}

/// The column-3 type GTF and GENCODE GFF3 give every transcript, whatever its
/// biotype.
pub const GENCODE_TRANSCRIPT_TYPE: &str = "transcript";

/// The `ID` prefix an Ensembl-style GFF3 gives its transcript records.
///
/// GENCODE GFF3 does not namespace ids this way, so the presence of this string
/// distinguishes the two styles.
const ENSEMBL_TRANSCRIPT_ID_PREFIX: &str = "ID=transcript:";

/// Parent transcript ID of an exon in a GTF file.
pub const GTF_TRANSCRIPT_ATTR: &str = "transcript_id";

/// Parent ID of a GFF3 feature.
pub const GFF_TRANSCRIPT_ATTR: &str = "Parent";

/// Default exon type.
pub const DEFAULT_EXON_TYPES: &[&str] = &["exon"];

/// Whether a record's column-3 type is one of `types`.
fn type_matches(record_type: &[u8], types: &[&str]) -> bool {
    types.iter().any(|t| t.as_bytes() == record_type)
}

/// Read which column-3 types a GFF3 file uses for its transcripts.
///
/// GTF and GENCODE GFF3 type every transcript as `transcript`, whatever its
/// biotype. Ensembl-style GFF3 instead uses Sequence Ontology terms for the
/// biotypes (`mRNA`, `lnc_RNA`, `snoRNA`, `tRNA`, etc.).
///
/// Ensembl namespaces a transcript's id as `ID=transcript:...`, GENCODE as
/// `ID=transcript...`, so we can distinguish between them.
///
/// Once we know the style of the GFF3, we can read the transcript types from
/// the file itself. Every column-3 value on a line containing `ID=transcript:`
/// is a transcript type.
///
/// This is a plain line scan; records are not parsed.
fn detect_gff_transcript_types(path: &Path) -> Result<Vec<String>> {
    let mut types: Vec<String> = Vec::new();

    for line in open_annotation(path)?.lines() {
        let line = line.context("failed to read annotation line")?;
        if line.starts_with('#') || !line.contains(ENSEMBL_TRANSCRIPT_ID_PREFIX) {
            continue;
        }
        if let Some(col3) = line.split('\t').nth(2)
            && !types.iter().any(|t| t == col3)
        {
            types.push(col3.to_string());
        }
    }

    if types.is_empty() {
        types.push(GENCODE_TRANSCRIPT_TYPE.to_string());
    }
    Ok(types)
}

/// Strip Ensembl's GFF3 type prefix from an identifier (`"gene:ENSMUSG…"` to
/// `"ENSMUSG…"`), so gene names match the bare ids GENCODE and GTF files use.
fn strip_gff_id_prefix(id: &str) -> String {
    for prefix in ["gene:", "transcript:"] {
        if let Some(bare) = id.strip_prefix(prefix) {
            return bare.to_string();
        }
    }
    id.to_string()
}

/// Read a column-9 attribute as a string.
///
/// A `Value` is either a single string or an array of them (e.g. a GFF3 `tag`
/// with several comma-separated values); the first is taken.
fn attr_string(record: &gff::feature::RecordBuf, key: &str) -> Option<String> {
    record
        .attributes()
        .get(key.as_bytes())
        .and_then(|v| {
            v.as_string()
                .or_else(|| v.as_array().and_then(|a| a.first()).map(|s| s.as_ref()))
        })
        .map(|s| s.to_string())
}

/// Map a noodles GFF/GTF strand to a single character (`'+'`, `'-'`, or `'*'`).
fn gff_strand_to_char(s: gff::feature::record::Strand) -> char {
    use gff::feature::record::Strand;
    match s {
        Strand::Forward => '+',
        Strand::Reverse => '-',
        _ => '*',
    }
}

/// What an annotation file parses into, before any index is built.
///
/// The two halves are the [`Interval`] / [`Feature`] split: `intervals` is the
/// geometry to search (many per feature such as the exons of a transcript, or the
/// pieces a blacklist splits a feature into), and `features` is the ordered list
/// of matrix columns those intervals point at, via their `var_idx`.
///
/// The intervals are left **unsorted**, per chromosome, rather than built into a
/// [`GenomeIndex`] here: [`parse_annotation_files`] may merge several files, and
/// shifting each file's `var_idx` by a running offset means the buckets are
/// still growing. Building the [`ChromIndex`] once at the end sorts each
/// chromosome exactly once, rather than sorting per file and re-sorting on merge.
struct ParsedAnnotation {
    intervals: AHashMap<String, Vec<Interval>>,
    features: Vec<Feature>,
}

/// Load a BED file of blacklisted regions into a searchable [`GenomeIndex`].
///
/// Blacklisted regions are only ever *excluded*, never counted, so unlike an
/// annotation this drops the [`Feature`]s. Counting *into* the regions of a
/// BED file goes through [`parse_annotation_files`] instead.
pub fn parse_blacklist_bed(path: &Path) -> Result<GenomeIndex> {
    Ok(build_genome_index(parse_bed_file(path)?.intervals))
}

/// Parse a BED file.
///
/// Feature names are returned in BED-file order. The 4th BED column is used
/// as the feature name when present and non-empty (and not `.`); otherwise
/// `"chrom:start-end"`. Lines beginning with `#` are skipped by the noodles
/// parser.
///
/// Uses `bed::io::Reader<3, _>` so extra columns (5+) are silently accepted.
/// Coordinate conversion: noodles adds 1 to BED's 0-based start when storing
/// it as `Position`, so `start_0 = feature_start().get() - 1`; the end is
/// stored verbatim, so `end_0 = feature_end().get()`.
fn parse_bed_file(path: &Path) -> Result<ParsedAnnotation> {
    let mut reader = bed::io::Reader::<3, _>::new(open_annotation(path)?);
    let mut record = bed::Record::<3>::default();

    let mut intervals_by_chrom: AHashMap<String, Vec<Interval>> = AHashMap::new();
    let mut var: Vec<Feature> = Vec::new();
    let mut region_idx: usize = 0;

    loop {
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => return Err(anyhow::Error::from(e)).context("failed to read BED record"),
        }

        let chrom = record.reference_sequence_name().to_string();

        let start = record
            .feature_start()
            .context("invalid BED start coordinate")?
            .get()
            - 1;

        let end = match record
            .feature_end()
            .transpose()
            .context("invalid BED end coordinate")?
        {
            Some(p) => p.get(),
            None => continue,
        };

        if end <= start {
            continue;
        }

        // 4th column lives in other_fields[0] for a Reader<3, _>.
        let name = record
            .other_fields()
            .get(0)
            .and_then(|n| {
                let s = n.to_string();
                if s.is_empty() || s == "." {
                    None
                } else {
                    Some(s)
                }
            })
            .unwrap_or_else(|| format!("{}:{}-{}", chrom, start, end));

        // 6th column (strand) lives in other_fields[2] for a Reader<3, _>.
        let strand = record
            .other_fields()
            .get(2)
            .and_then(|s| s.to_string().chars().next())
            .map(|c| if c == '+' || c == '-' { c } else { '*' })
            .unwrap_or('*');

        var.push(Feature {
            chrom: chrom.clone(),
            start,
            end,
            name: name.clone(),
            strand,
        });
        intervals_by_chrom.entry(chrom).or_default().push(Interval {
            start,
            end,
            var_idx: region_idx,
        });
        region_idx += 1;
    }

    Ok(ParsedAnnotation {
        intervals: intervals_by_chrom,
        features: var,
    })
}

/// Build a [`GenomeIndex`] from per-chromosome unsorted interval lists,
/// sorting each chromosome's intervals exactly once.
fn build_genome_index(intervals_by_chrom: AHashMap<String, Vec<Interval>>) -> GenomeIndex {
    intervals_by_chrom
        .into_iter()
        .map(|(c, ivs)| (c, ChromIndex::build(ivs)))
        .collect()
}

// GTF / GFF3 parsers
// Both `gtf::io::Reader::record_bufs()` and `gff::io::Reader::record_bufs()`
// yield `gff::feature::RecordBuf`, so a single helper covers both.
fn build_annotation_index<I>(
    iter: I,
    feature_types: Option<&[&str]>,
    name_attr: &str,
) -> Result<ParsedAnnotation>
where
    I: Iterator<Item = io::Result<gff::feature::RecordBuf>>,
{
    let mut intervals_by_chrom: AHashMap<String, Vec<Interval>> = AHashMap::new();
    let mut var: Vec<Feature> = Vec::new();
    let mut region_idx: usize = 0;

    for result in iter {
        let record = result.context("failed to read annotation record")?;

        if let Some(types) = feature_types
            && !type_matches(record.ty(), types)
        {
            continue;
        }

        let chrom = record.reference_sequence_name().to_string();
        // GTF and GFF3 use 1-based inclusive coordinates.
        let start = record.start().get() - 1;
        let end = record.end().get();

        if end <= start {
            continue;
        }

        let name = attr_string(&record, name_attr)
            .unwrap_or_else(|| format!("{}:{}-{}", chrom, start, end));

        let strand = gff_strand_to_char(record.strand());

        var.push(Feature {
            chrom: chrom.clone(),
            start,
            end,
            name: name.clone(),
            strand,
        });
        intervals_by_chrom.entry(chrom).or_default().push(Interval {
            start,
            end,
            var_idx: region_idx,
        });
        region_idx += 1;
    }

    Ok(ParsedAnnotation {
        intervals: intervals_by_chrom,
        features: var,
    })
}

/// Parse a GTF file.
///
/// `feature_types` filters by the feature-type column (e.g. `["gene"]`,
/// `["exon", "CDS"]`); a record is kept when its type is any of them, and
/// `None` includes every type. `name_attr` is the attribute key used as the
/// feature name (e.g. `"gene_id"`, `"gene_name"`); falls back to
/// `"chrom:start-end"` when absent. Coordinates are converted from 1-based
/// inclusive to 0-based half-open.
fn parse_gtf_file(
    path: &Path,
    feature_types: Option<&[&str]>,
    name_attr: &str,
) -> Result<ParsedAnnotation> {
    let mut reader = gtf::io::Reader::new(open_annotation(path)?);
    build_annotation_index(reader.record_bufs(), feature_types, name_attr)
}

/// Parse a GFF3 file.
///
/// `feature_types` filters by the type column (e.g. `["mRNA", "lnc_RNA"]`); a
/// record is kept when its type is any of them, and `None` includes every type.
/// `name_attr` is the attribute tag used as the feature name (e.g. `"ID"`,
/// `"Name"`); falls back to `"chrom:start-end"`. Coordinates are converted
/// from 1-based inclusive to 0-based half-open.
fn parse_gff_file(
    path: &Path,
    feature_types: Option<&[&str]>,
    name_attr: &str,
) -> Result<ParsedAnnotation> {
    let mut reader = gff::io::Reader::new(open_annotation(path)?);
    build_annotation_index(reader.record_bufs(), feature_types, name_attr)
}

/// Parse one or more annotation files (BED, GTF, GFF3, or any mix) and merge
/// their intervals into a single [`GenomeIndex`].
///
/// Format is detected from each file's extension:
///
/// | Extension(s)                | Parser | Default region types    | Default `name_attr` |
/// |-----------------------------|--------|-------------------------|---------------------|
/// | `.gtf`, `.gtf.gz`           | GTF    | `"transcript"`          | `"transcript_id"`   |
/// | `.gff`, `.gff3`, `.gff3.gz` | GFF3   | detected from the file  | `"ID"`              |
/// | anything else               | BED    | -                       | -                   |
///
/// A record becomes a region when its column-3 type is **any** of
/// `feature_types` (several may be given) which is what lets an Ensembl-style
/// GFF3 contribute all of its transcript biotypes (`mRNA`, `lnc_RNA`, `tRNA`, …)
/// rather than just one. `None` reads the types out of the file itself (see
/// [`detect_gff_transcript_types`]), so every transcript is kept regardless of
/// biotype. In metagene mode `exon_types` plays the same role for exons,
/// defaulting to [`DEFAULT_EXON_TYPES`], and `name_attr` names the *grouping*
/// attribute instead, defaulting to `"gene_id"` for both formats.
///
/// Feature names are concatenated in file order. Each file's `var_idx` values are
/// shifted by the cumulative feature count so every feature in the combined
/// index maps to a unique slot in the returned name vector.
///
/// `feature_types`, `exon_types` and `name_attr` only affect GTF / GFF3 files;
/// BED files always use the 4th column (or `"chrom:start-end"`) as the name.
///
/// When `metagene = true`, exon intervals are grouped per gene, so that every
/// exon of a gene shares one `var_idx` (= gene index). This lets the counting loop
/// count a read for a gene regardless of how many of its exons it overlaps.
///
/// Exons are grouped **per transcript** by default, which every style
/// supports: a GTF exon names its transcript in [`GTF_TRANSCRIPT_ATTR`], a GFF3
/// exon in [`GFF_TRANSCRIPT_ATTR`]. Ensembl's `gene:` / `transcript:` id
/// prefixes are stripped, so names come out as the bare ids the other formats
/// use.
///
/// Passing `name_attr` groups on that attribute instead (`Some("gene_id")` for
/// per-gene features). That works for GTF and GENCODE GFF3, whose exons carry a
/// gene id. An Ensembl-style GFF3 does not: its exons name only their
/// transcript, so grouping it by gene is an error rather than a silently wrong
/// per-exon result.
/// Fold one file's [`ParsedAnnotation`] into the running merge.
///
/// Each file's features are appended to `all_features`, and its intervals'
/// `var_idx` shifted by `offset` (the feature count before this file) so every
/// feature keeps a unique column in the combined output.
fn merge(
    parsed: ParsedAnnotation,
    offset: usize,
    all_features: &mut Vec<Feature>,
    all_intervals: &mut AHashMap<String, Vec<Interval>>,
) {
    all_features.extend(parsed.features);
    for (chrom, intervals) in parsed.intervals {
        let bucket = all_intervals.entry(chrom).or_default();
        for mut interval in intervals {
            interval.var_idx += offset;
            bucket.push(interval);
        }
    }
}

pub fn parse_annotation_files<P: AsRef<Path>>(
    paths: impl IntoIterator<Item = P>,
    feature_types: Option<&[&str]>,
    exon_types: Option<&[&str]>,
    name_attr: Option<&str>,
    metagene: bool,
) -> Result<(GenomeIndex, Vec<Feature>)> {
    let mut all_var: Vec<Feature> = Vec::new();
    let mut intervals_by_chrom: AHashMap<String, Vec<Interval>> = AHashMap::new();

    for path in paths {
        let path = path.as_ref();

        // Strip any compression suffix so the format is detected from the
        // underlying extension (e.g. `genes.gtf.gz` -> `genes.gtf`). The
        // decompression itself is decided by magic bytes, not by this name.
        let filename = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let stem = [".gz", ".bgz", ".bgzf"]
            .iter()
            .find_map(|suffix| filename.strip_suffix(suffix))
            .unwrap_or(filename);
        let is_gtf = stem.ends_with(".gtf");
        let is_gff = stem.ends_with(".gff") || stem.ends_with(".gff3");

        // Metagene mode groups exon records so all exons of one feature share a
        // value. The offset shift below still applies so multi-file merges work.
        //
        // Grouped per transcript by default, which every style supports: a GTF
        // exon names its transcript in `transcript_id`, a GFF3 exon in `Parent`.
        // Naming `name_attr` groups on that attribute instead: `"gene_id"` for
        // per-gene features, which works for GTF and GENCODE GFF3. An
        // Ensembl-style GFF3 carries no gene id on its exons at all, so asking
        // it to group by gene is an error rather than a silent per-exon result.
        if metagene && (is_gtf || is_gff) {
            let exon_types = exon_types.unwrap_or(DEFAULT_EXON_TYPES);
            let offset = all_var.len();
            let parsed = if is_gtf {
                let group_attr = name_attr.unwrap_or(GTF_TRANSCRIPT_ATTR);
                let mut reader = gtf::io::Reader::new(open_annotation(path)?);
                build_metagene_index(reader.record_bufs(), exon_types, group_attr)
            } else {
                let group_attr = name_attr.unwrap_or(GFF_TRANSCRIPT_ATTR);
                let mut reader = gff::io::Reader::new(open_annotation(path)?);
                build_metagene_index(reader.record_bufs(), exon_types, group_attr)
            }
            .with_context(|| format!("failed to parse annotation file: {}", path.display()))?;

            merge(parsed, offset, &mut all_var, &mut intervals_by_chrom);
            continue;
        }

        let offset = all_var.len();

        // Region types, when the caller named none: GTF always types its
        // transcripts `transcript`, while a GFF3 has to be asked which style
        // it is (see `detect_gff_transcript_types`). The two formats also
        // differ in how a region names itself: `transcript_id` vs `ID`.
        //
        // `detected` owns the type strings for as long as the borrowed slice
        // handed to the parser is alive.
        let detected: Vec<String>;
        let parsed = if is_gtf {
            let region_types = feature_types.unwrap_or(&[GENCODE_TRANSCRIPT_TYPE]);
            parse_gtf_file(
                path,
                Some(region_types),
                name_attr.unwrap_or("transcript_id"),
            )
        } else if is_gff {
            let borrowed: Vec<&str>;
            let region_types = match feature_types {
                Some(types) => types,
                None => {
                    detected = detect_gff_transcript_types(path)?;
                    borrowed = detected.iter().map(String::as_str).collect();
                    &borrowed
                }
            };
            parse_gff_file(path, Some(region_types), name_attr.unwrap_or("ID"))
        } else {
            parse_bed_file(path)
        }
        .with_context(|| format!("failed to parse annotation file: {}", path.display()))?;

        merge(parsed, offset, &mut all_var, &mut intervals_by_chrom);
    }

    Ok((build_genome_index(intervals_by_chrom), all_var))
}

/// Group exon records by the feature they belong to, so all exons of one
/// transcript (or gene) share a single region index (`var_idx`).
///
/// This is the metagene parsing path. The resulting intervals have `var_idx` =
/// group index (into the returned `Feature` vec) rather than a per-exon index,
/// so the counting inner loop naturally treats all exons of one group as one
/// region. A group's span runs from its first exon's start to its last exon's
/// end.
///
/// Exons are grouped by the value of their `group_attr` attribute,
/// `transcript_id` or `Parent` for per-transcript features, `gene_id` for
/// per-gene ones. Ensembl's `gene:` / `transcript:` id prefixes are
/// stripped, so names come out as the bare ids the other formats use.
///
/// Errors when no exon in the file carries `group_attr`: that annotation simply
/// cannot be grouped that way (an Ensembl-style GFF3 asked to group by
/// `gene_id`, say, whose exons name only their transcript), and grouping every
/// exon into a feature of its own would be a silently wrong answer.
fn build_metagene_index<I>(
    iter: I,
    exon_types: &[&str],
    group_attr: &str,
) -> Result<ParsedAnnotation>
where
    I: Iterator<Item = io::Result<gff::feature::RecordBuf>>,
{
    let mut intervals_by_chrom: AHashMap<String, Vec<Interval>> = AHashMap::new();
    let mut var: Vec<Feature> = Vec::new();
    let mut group_to_idx: AHashMap<String, usize> = AHashMap::new();
    let mut n_exons = 0usize;

    for result in iter {
        let record = result.context("failed to read annotation record")?;

        if !type_matches(record.ty(), exon_types) {
            continue;
        }

        let chrom = record.reference_sequence_name().to_string();
        let start = record.start().get() - 1;
        let end = record.end().get();

        if end <= start {
            continue;
        }
        n_exons += 1;

        let Some(group_key) = attr_string(&record, group_attr) else {
            continue;
        };
        let group_name = strip_gff_id_prefix(&group_key);

        let group_idx = match group_to_idx.get(&group_name).copied() {
            Some(idx) => {
                // Extend the group's genomic span to cover this exon.
                var[idx].start = var[idx].start.min(start);
                var[idx].end = var[idx].end.max(end);
                idx
            }
            None => {
                let idx = var.len();
                group_to_idx.insert(group_name.clone(), idx);
                var.push(Feature {
                    chrom: chrom.clone(),
                    start,
                    end,
                    name: group_name.clone(),
                    strand: gff_strand_to_char(record.strand()),
                });
                idx
            }
        };

        intervals_by_chrom.entry(chrom).or_default().push(Interval {
            start,
            end,
            var_idx: group_idx,
        });
    }

    if n_exons > 0 && var.is_empty() {
        bail!(
            "cannot group exons by {group_attr:?}: no exon record carries that attribute. \
             An Ensembl-style GFF3 names only its transcript, in \"Parent\", so it cannot be \
             grouped by \"gene_id\"; drop --transcriptIDtag to group per transcript."
        );
    }

    Ok(ParsedAnnotation {
        intervals: intervals_by_chrom,
        features: var,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;

    use tempfile::TempDir;

    // The same mm39 locus in five files: three genes Prim2, Rab23, and the
    // lncRNA (A930006A01Rik),.
    // The styles differ in ways the parsers have to absorb:
    //
    //   * GENCODE names chromosomes "chr1" and versions its ids
    //     ("ENSMUSG00000026134.13"); Ensembl uses "1" and bare ids.
    //   * GTF and GENCODE GFF3 type every transcript "transcript", whatever its
    //     biotype. Ensembl GFF3 uses the biotype itself: "mRNA" for the coding
    //     genes, "lnc_RNA" for the lncRNA.
    //   * Ensembl GFF3 namespaces ids ("gene:ENSMUSG…"), its exons carry no gene
    //     id (only a Parent naming their transcript), and it types the lncRNA's
    //     gene record "ncRNA_gene" rather than "gene".
    const BED: &str = "mm39_genes.bed";
    const GTF_GENCODE: &str = "mm39_gencode.gtf.gz";
    const GTF_ENSEMBL: &str = "mm39_ensembl.gtf.gz";
    const GFF_GENCODE: &str = "mm39_gencode.gff3.gz";
    const GFF_ENSEMBL: &str = "mm39_ensembl.gff3.gz";

    // The three genes, in the order they appear in the annotation files, as
    // 0-based half-open spans.
    const LNCRNA: (usize, usize) = (4_855_978, 4_874_316);
    const PRIM2: (usize, usize) = (33_492_887, 33_709_466);
    const RAB23: (usize, usize) = (33_758_944, 33_781_645);

    // The lncRNA's BED entry covers only its transcript, so it stops short of
    // the gene's annotated end.
    const LNCRNA_BED: (usize, usize) = (4_855_978, 4_857_360);

    const GENCODE_IDS: [&str; 3] = [
        "ENSMUSG00000120403.2",
        "ENSMUSG00000026134.13",
        "ENSMUSG00000004768.16",
    ];
    const ENSEMBL_IDS: [&str; 3] = [
        "ENSMUSG00000120403",
        "ENSMUSG00000026134",
        "ENSMUSG00000004768",
    ];

    // Every GTF/GFF here holds 48 transcripts of these 3 genes. Ensembl GFF3
    // annotates one exon more for the lncRNA than the other files do.
    const N_TRANSCRIPTS: usize = 48;

    fn data(name: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/testdata/annotation_files")
            .join(name)
    }

    fn span(v: &Feature) -> (usize, usize) {
        (v.start, v.end)
    }

    fn names(var: &[Feature]) -> Vec<String> {
        var.iter().map(|v| v.name.clone()).collect()
    }

    fn spans(var: &[Feature]) -> Vec<(usize, usize)> {
        var.iter().map(span).collect()
    }

    /// Feature indices the index reports as overlapping `[start, end)`.
    fn hits(index: &GenomeIndex, chrom: &str, start: usize, end: usize) -> Vec<usize> {
        let mut v: Vec<usize> = index
            .get(chrom)
            .unwrap_or_else(|| panic!("{chrom} missing from index"))
            .find(start, end)
            .map(|iv| iv.var_idx)
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    }

    // BED

    #[test]
    fn bed_features_are_parsed_with_names_and_strands() {
        let var = parse_bed_file(&data(BED)).unwrap().features;

        // Columns 7-9 are present and must be tolerated, not choke the reader.
        assert_eq!(names(&var), ["A930006A01Rik", "Prim2", "Rab23"]);

        // BED starts are already 0-based and ends exclusive, so they carry over.
        assert_eq!(spans(&var), [LNCRNA_BED, PRIM2, RAB23]);
        assert!(var.iter().all(|v| v.chrom == "chr1"));
        assert_eq!(
            var.iter().map(|v| v.strand).collect::<Vec<_>>(),
            ['+', '-', '+']
        );
    }

    #[test]
    fn bed_index_finds_features_by_position() {
        let index = build_genome_index(parse_bed_file(&data(BED)).unwrap().intervals);

        assert_eq!(hits(&index, "chr1", PRIM2.0 + 100, PRIM2.0 + 200), vec![1]);
        // The intergenic gap between Prim2 and Rab23.
        assert_eq!(hits(&index, "chr1", PRIM2.1, RAB23.0), Vec::<usize>::new());
        // Half-open: the base at `end` is outside the feature.
        assert_eq!(hits(&index, "chr1", PRIM2.1 - 1, PRIM2.1), vec![1]);
        assert_eq!(
            hits(&index, "chr1", PRIM2.1, PRIM2.1 + 1),
            Vec::<usize>::new()
        );
    }

    #[test]
    fn bed_falls_back_to_a_positional_name_and_skips_degenerate_lines() {
        // Every shipped BED feature is named, so the fallback paths (missing
        // name, "." name, comments, end <= start) need a purpose-built file.
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("edge.bed");
        let mut f = File::create(&path).unwrap();
        writeln!(f, "# a comment line").unwrap();
        writeln!(f, "chr1\t100\t200\tnamed").unwrap();
        writeln!(f, "chr1\t300\t400\t.").unwrap();
        writeln!(f, "chr1\t500\t600").unwrap();
        writeln!(f, "chr1\t900\t900\tempty").unwrap();
        f.flush().unwrap();

        let var = parse_bed_file(&path).unwrap().features;

        // The zero-length feature is dropped; "." and a missing column both
        // fall back to "chrom:start-end".
        assert_eq!(names(&var), ["named", "chr1:300-400", "chr1:500-600"]);
        // No strand column means unstranded.
        assert!(var.iter().all(|v| v.strand == '*'));
    }

    // Transcript types: every biotype must survive, in every style

    #[test]
    fn transcript_types_are_read_out_of_an_ensembl_gff() {
        // Ensembl marks a transcript by namespacing its id; the column-3 values
        // on those lines are exactly the biotypes this file uses, in the order
        // they first appear (the lncRNA locus comes before the coding genes).
        let types = detect_gff_transcript_types(&data(GFF_ENSEMBL)).unwrap();
        assert_eq!(types, ["lnc_RNA", "mRNA", "transcript"]);
    }

    #[test]
    fn a_gff_without_namespaced_ids_is_gencode_style() {
        // GENCODE never writes "ID=transcript:", and types every transcript
        // "transcript" whatever its biotype.
        let types = detect_gff_transcript_types(&data(GFF_GENCODE)).unwrap();
        assert_eq!(types, [GENCODE_TRANSCRIPT_TYPE]);
    }

    #[test]
    fn the_detected_types_keep_every_transcript_of_every_style() {
        // The point of detecting rather than fixing the types: the lncRNA is
        // typed "transcript" by GTF/GENCODE and "lnc_RNA" by Ensembl GFF3, and
        // all 48 transcripts have to come through either way.
        for file in [GTF_GENCODE, GTF_ENSEMBL, GFF_GENCODE, GFF_ENSEMBL] {
            let var = parse_annotation_files([data(file)], None, None, None, false)
                .unwrap()
                .1;
            assert_eq!(var.len(), N_TRANSCRIPTS, "{file}");
        }
    }

    #[test]
    fn an_unknown_biotype_is_still_detected() {
        // Detection reads the file rather than a fixed list, so a biotype no
        // list would have anticipated is picked up all the same.
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("exotic.gff3");
        let mut f = File::create(&path).unwrap();
        for row in [
            "1\tensembl\tncRNA_gene\t101\t900\t.\t+\t.\tID=gene:G1;Name=Odd",
            "1\tensembl\tvault_RNA\t101\t900\t.\t+\t.\tID=transcript:T1;Parent=gene:G1",
            "1\tensembl\texon\t101\t900\t.\t+\t.\tParent=transcript:T1",
        ] {
            writeln!(f, "{row}").unwrap();
        }
        f.flush().unwrap();

        assert_eq!(detect_gff_transcript_types(&path).unwrap(), ["vault_RNA"]);
        let var = parse_annotation_files([&path], None, None, None, false)
            .unwrap()
            .1;
        assert_eq!(names(&var), ["transcript:T1"]);
    }

    #[test]
    fn the_lncrna_transcript_is_kept_whatever_its_biotype_is_called() {
        // Ensembl types this transcript "lnc_RNA", so a "transcript"-only
        // default would have dropped it. It must appear, on the lncRNA locus.
        let ensembl = parse_annotation_files([data(GFF_ENSEMBL)], None, None, None, false)
            .unwrap()
            .1;
        let lnc: Vec<&Feature> = ensembl
            .iter()
            .filter(|v| v.start >= LNCRNA.0 && v.end <= LNCRNA.1)
            .collect();
        assert!(!lnc.is_empty());
        assert!(lnc.iter().all(|v| v.chrom == "1" && v.strand == '+'));

        let coding_only =
            parse_annotation_files([data(GFF_ENSEMBL)], Some(&["mRNA"]), None, None, false)
                .unwrap()
                .1;
        assert_eq!(coding_only.len(), 37);
        assert!(coding_only.iter().all(|v| v.start >= PRIM2.0));
    }

    #[test]
    fn several_feature_types_can_be_requested_at_once() {
        // Ensembl GFF3 splits the transcripts across mRNA (37), lnc_RNA (10)
        // and transcript (1); asking for all three recovers the full set.
        let combined = parse_gff_file(
            &data(GFF_ENSEMBL),
            Some(&["mRNA", "lnc_RNA", "transcript"]),
            "ID",
        )
        .unwrap()
        .features;
        assert_eq!(combined.len(), N_TRANSCRIPTS);

        // The same mechanism selects unrelated types together.
        let genes = parse_gff_file(&data(GFF_ENSEMBL), Some(&["gene", "ncRNA_gene"]), "gene_id")
            .unwrap()
            .features;
        assert_eq!(genes.len(), 3);
    }

    #[test]
    fn ensembl_types_the_lncrna_gene_record_as_ncrna_gene() {
        // Asking only for "gene" misses it.
        let only_gene = parse_gff_file(&data(GFF_ENSEMBL), Some(&["gene"]), "gene_id")
            .unwrap()
            .features;
        assert_eq!(only_gene.len(), 2);
        assert_eq!(spans(&only_gene), [PRIM2, RAB23]);

        // The GTFs and the GENCODE GFF3 type all three the same way.
        for (file, ids) in [(GTF_GENCODE, GENCODE_IDS), (GTF_ENSEMBL, ENSEMBL_IDS)] {
            let genes = parse_gtf_file(&data(file), Some(&["gene"]), "gene_id")
                .unwrap()
                .features;
            assert_eq!(names(&genes), ids, "{file}");
            assert_eq!(spans(&genes), [LNCRNA, PRIM2, RAB23], "{file}");
        }
    }

    // GTF

    #[test]
    fn gtf_transcripts_are_parsed_from_both_styles() {
        for (file, chrom) in [(GTF_GENCODE, "chr1"), (GTF_ENSEMBL, "1")] {
            let var = parse_gtf_file(&data(file), Some(&["transcript"]), "transcript_id")
                .unwrap()
                .features;

            assert_eq!(var.len(), N_TRANSCRIPTS, "{file}");
            assert_eq!(var[0].chrom, chrom, "{file}");
            // GTF is 1-based inclusive, so the start comes back one lower.
            assert_eq!(var[0].start, LNCRNA.0, "{file}");
            assert!(var.iter().all(|v| v.name.starts_with("ENSMUST")), "{file}");
        }
    }

    #[test]
    fn gtf_gene_ids_are_versioned_in_gencode_but_not_in_ensembl() {
        let gencode = parse_gtf_file(&data(GTF_GENCODE), Some(&["gene"]), "gene_id")
            .unwrap()
            .features;
        let ensembl = parse_gtf_file(&data(GTF_ENSEMBL), Some(&["gene"]), "gene_id")
            .unwrap()
            .features;

        assert_eq!(names(&gencode), GENCODE_IDS);
        assert_eq!(names(&ensembl), ENSEMBL_IDS);

        // The ids and chromosome names differ, but the loci must not.
        assert_eq!(spans(&gencode), spans(&ensembl));
        assert_eq!(gencode[0].chrom, "chr1");
        assert_eq!(ensembl[0].chrom, "1");
    }

    #[test]
    fn gtf_name_attr_chooses_the_feature_name() {
        let by_name = parse_gtf_file(&data(GTF_GENCODE), Some(&["gene"]), "gene_name")
            .unwrap()
            .features;
        assert_eq!(names(&by_name), ["A930006A01Rik", "Prim2", "Rab23"]);

        // An attribute no record carries falls back to the position.
        let fallback = parse_gtf_file(&data(GTF_GENCODE), Some(&["gene"]), "no_such_attr")
            .unwrap()
            .features;
        assert_eq!(fallback[0].name, format!("chr1:{}-{}", LNCRNA.0, LNCRNA.1));
    }

    #[test]
    fn no_feature_type_filter_keeps_every_record() {
        let all = parse_gtf_file(&data(GTF_GENCODE), None, "gene_id")
            .unwrap()
            .features;
        // 391 exon + 269 CDS + 162 UTR + 48 transcript + 37 start + 37 stop + 3 gene
        assert_eq!(all.len(), 391 + 269 + 162 + 48 + 37 + 37 + 3);
    }

    // GFF3

    #[test]
    fn gff_regions_are_named_by_id_and_agree_across_styles() {
        let gencode = parse_gff_file(&data(GFF_GENCODE), Some(&["transcript"]), "ID")
            .unwrap()
            .features;
        assert_eq!(gencode.len(), N_TRANSCRIPTS);
        assert!(gencode.iter().all(|v| v.name.starts_with("ENSMUST")));

        // Ensembl namespaces the id; only metagene grouping strips the prefix.
        let ensembl = parse_gff_file(
            &data(GFF_ENSEMBL),
            Some(&["mRNA", "lnc_RNA", "transcript"]),
            "ID",
        )
        .unwrap()
        .features;
        assert_eq!(ensembl.len(), N_TRANSCRIPTS);
        assert!(
            ensembl
                .iter()
                .all(|v| v.name.starts_with("transcript:ENSMUST"))
        );
    }

    #[test]
    fn every_format_agrees_on_where_the_coding_genes_are() {
        // The invariant tying the files together: whatever the style, Prim2
        // and Rab23 come out on identical coordinates. (The lncRNA is excluded
        // only because its BED entry covers the transcript, not the gene.)
        let bed = parse_bed_file(&data(BED)).unwrap().features;
        let gencode_gtf = parse_gtf_file(&data(GTF_GENCODE), Some(&["gene"]), "gene_id")
            .unwrap()
            .features;
        let ensembl_gtf = parse_gtf_file(&data(GTF_ENSEMBL), Some(&["gene"]), "gene_id")
            .unwrap()
            .features;
        let ensembl_gff = parse_gff_file(&data(GFF_ENSEMBL), Some(&["gene"]), "gene_id")
            .unwrap()
            .features;

        assert_eq!(spans(&bed)[1..], [PRIM2, RAB23]);
        assert_eq!(spans(&gencode_gtf)[1..], [PRIM2, RAB23]);
        assert_eq!(spans(&ensembl_gtf)[1..], [PRIM2, RAB23]);
        assert_eq!(spans(&ensembl_gff), [PRIM2, RAB23]);
    }

    // Compression

    #[test]
    fn a_plain_annotation_parses_the_same_as_its_bgzf_copy() {
        // The shipped GTFs are BGZF; compression is detected from the magic
        // bytes, so an uncompressed copy under any name must parse identically.
        let dir = TempDir::new().unwrap();
        let plain_path = dir.path().join("genes.gtf");
        let mut decoder = MultiGzDecoder::new(File::open(data(GTF_GENCODE)).unwrap());
        let mut plain = Vec::new();
        std::io::Read::read_to_end(&mut decoder, &mut plain).unwrap();
        File::create(&plain_path)
            .unwrap()
            .write_all(&plain)
            .unwrap();

        let from_plain = parse_gtf_file(&plain_path, Some(&["transcript"]), "transcript_id")
            .unwrap()
            .features;
        let from_bgzf = parse_gtf_file(&data(GTF_GENCODE), Some(&["transcript"]), "transcript_id")
            .unwrap()
            .features;

        assert_eq!(from_plain.len(), N_TRANSCRIPTS);
        assert_eq!(names(&from_plain), names(&from_bgzf));
        assert_eq!(spans(&from_plain), spans(&from_bgzf));
    }

    // parse_annotation_files: format dispatch, defaults, merging

    #[test]
    fn format_is_detected_from_the_extension_beneath_any_compression_suffix() {
        let bed = parse_annotation_files([data(BED)], None, None, None, false)
            .unwrap()
            .1;
        assert_eq!(names(&bed), ["A930006A01Rik", "Prim2", "Rab23"]);

        // ".gtf.gz" / ".gff3.gz" are dispatched on the extension under the
        // compression suffix, and named by transcript_id / ID respectively.
        for file in [GTF_GENCODE, GFF_GENCODE] {
            let var = parse_annotation_files([data(file)], None, None, None, false)
                .unwrap()
                .1;
            assert_eq!(var.len(), N_TRANSCRIPTS, "{file}");
            assert!(var.iter().all(|v| v.name.starts_with("ENSMUST")), "{file}");
        }
    }

    #[test]
    fn merging_several_files_shifts_feature_indices_so_they_stay_unique() {
        let (index, var) =
            parse_annotation_files([data(BED), data(GTF_GENCODE)], None, None, None, false)
                .unwrap();

        // 3 BED genes followed by 48 GTF transcripts.
        assert_eq!(var.len(), 3 + N_TRANSCRIPTS);
        assert_eq!(var[0].name, "A930006A01Rik");
        assert!(var[3].name.starts_with("ENSMUST"));

        // Inside Prim2 both files contribute: BED feature 1, plus transcripts
        // whose indices were shifted past the BED's 3 features.
        let overlapping = hits(&index, "chr1", PRIM2.0 + 100, PRIM2.0 + 200);
        assert!(overlapping.contains(&1));
        assert!(overlapping.iter().any(|&v| v >= 3));
        assert!(overlapping.iter().all(|&v| v < 3 + N_TRANSCRIPTS));
    }

    // Metagene grouping

    #[test]
    fn metagene_groups_per_transcript_by_default() {
        // Every style names an exon's transcript on the exon itself: GTF in
        // transcript_id, GFF3 in Parent, so the default grouping works on all
        // four, and gives one feature per transcript.
        for (file, chrom) in [
            (GTF_GENCODE, "chr1"),
            (GTF_ENSEMBL, "1"),
            (GFF_GENCODE, "chr1"),
            (GFF_ENSEMBL, "1"),
        ] {
            let (index, var) =
                parse_annotation_files([data(file)], None, None, None, true).unwrap();

            assert_eq!(var.len(), N_TRANSCRIPTS, "{file}");
            // Ensembl's "transcript:" id prefix is stripped, so the names match
            // the bare ids the other styles use.
            assert!(var.iter().all(|v| v.name.starts_with("ENSMUST")), "{file}");
            assert!(var.iter().all(|v| v.chrom == chrom), "{file}");

            // Every exon is indexed, and points at one of the transcripts.
            let ivs: Vec<&Interval> = index.get(chrom).unwrap().iter().collect();
            assert!(ivs.len() >= 391, "{file}");
            assert!(ivs.iter().all(|iv| iv.var_idx < N_TRANSCRIPTS), "{file}");
        }
    }

    #[test]
    fn metagene_groups_per_gene_when_the_gene_attribute_is_named() {
        // Grouping on gene_id collapses the exons onto the 3 genes. The
        // Ensembl GFF3 is absent here: its exons carry no gene id (see below).
        for (file, chrom, ids) in [
            (GTF_GENCODE, "chr1", GENCODE_IDS),
            (GTF_ENSEMBL, "1", ENSEMBL_IDS),
            (GFF_GENCODE, "chr1", GENCODE_IDS),
        ] {
            let (index, var) =
                parse_annotation_files([data(file)], None, None, Some("gene_id"), true).unwrap();

            assert_eq!(names(&var), ids, "{file}");
            assert_eq!(spans(&var), [LNCRNA, PRIM2, RAB23], "{file}");
            assert_eq!(
                var.iter().map(|v| v.strand).collect::<Vec<_>>(),
                ['+', '-', '+'],
                "{file}"
            );

            // Every exon of a gene carries that gene's index, so a read spanning
            // several exons of one gene still resolves to a single feature.
            let ivs: Vec<&Interval> = index.get(chrom).unwrap().iter().collect();
            assert_eq!(ivs.iter().filter(|iv| iv.var_idx == 0).count(), 3, "{file}");
            assert_eq!(
                ivs.iter().filter(|iv| iv.var_idx == 1).count(),
                172,
                "{file}"
            );
            assert_eq!(
                ivs.iter().filter(|iv| iv.var_idx == 2).count(),
                216,
                "{file}"
            );
        }
    }

    #[test]
    fn grouping_an_ensembl_gff_by_gene_is_an_error() {
        // Its exons name only their transcript, so there is no gene id to group
        // on.
        let err = parse_annotation_files([data(GFF_ENSEMBL)], None, None, Some("gene_id"), true)
            .err()
            .expect("grouping by gene_id should fail on an Ensembl GFF3");
        let msg = format!("{err:#}");
        assert!(msg.contains("gene_id"), "{msg}");
        assert!(
            msg.contains("no exon record carries that attribute"),
            "{msg}"
        );
    }

    #[test]
    fn metagene_grouping_attribute_can_be_overridden() {
        // GENCODE GFF3 exons carry transcript_id as well as Parent, so grouping
        // on it is equivalent to the default.
        let by_tx_id =
            parse_annotation_files([data(GFF_GENCODE)], None, None, Some("transcript_id"), true)
                .unwrap()
                .1;
        let by_default = parse_annotation_files([data(GFF_GENCODE)], None, None, None, true)
            .unwrap()
            .1;
        assert_eq!(names(&by_tx_id), names(&by_default));
        assert_eq!(by_tx_id.len(), N_TRANSCRIPTS);
    }

    #[test]
    fn metagene_exon_types_are_configurable() {
        // --exonID reaches the parser: grouping CDS records per gene keeps the 2
        // coding genes and drops the lncRNA, which has no CDS.
        let (index, var) = parse_annotation_files(
            [data(GTF_GENCODE)],
            None,
            Some(&["CDS"]),
            Some("gene_id"),
            true,
        )
        .unwrap();

        assert_eq!(names(&var), [GENCODE_IDS[1], GENCODE_IDS[2]]);
        assert_eq!(index.get("chr1").unwrap().len(), 269);

        // Several exon types can be combined.
        let both = parse_annotation_files(
            [data(GTF_GENCODE)],
            None,
            Some(&["exon", "CDS"]),
            Some("gene_id"),
            true,
        )
        .unwrap();
        assert_eq!(both.1.len(), 3);
        assert_eq!(both.0.get("chr1").unwrap().len(), 391 + 269);
    }

    #[test]
    fn metagene_gene_ids_agree_across_styles_once_versions_are_dropped() {
        // GENCODE versions its ids and Ensembl does not, but strip the suffix
        // and each style that can group by gene names the same three genes.
        let strip_version = |var: &[Feature]| -> Vec<String> {
            var.iter()
                .map(|v| v.name.split('.').next().unwrap().to_string())
                .collect()
        };

        for file in [GTF_GENCODE, GTF_ENSEMBL, GFF_GENCODE] {
            let var = parse_annotation_files([data(file)], None, None, Some("gene_id"), true)
                .unwrap()
                .1;
            assert_eq!(strip_version(&var), ENSEMBL_IDS, "{file}");
        }
    }
}
