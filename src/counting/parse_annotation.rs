use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result};
use noodles::{bed, gff, gtf};

use super::region_index::{ChromIndex, Interval, RegionIndex, VarMeta};

/// Map a noodles GFF/GTF strand to a single character (`'+'`, `'-'`, or `'*'`).
fn gff_strand_to_char(s: gff::feature::record::Strand) -> char {
    use gff::feature::record::Strand;
    match s {
        Strand::Forward => '+',
        Strand::Reverse => '-',
        _ => '*',
    }
}

/// Per-chromosome unsorted interval lists plus the ordered per-feature metadata,
/// as produced by the raw (pre-index) parsers below.
type RawIntervals = (AHashMap<String, Vec<Interval>>, Vec<VarMeta>);

/// Parse a BED file and build a per-chromosome interval index.
///
/// Feature names are returned in BED-file order.  The 4th BED column is used
/// as the feature name when present and non-empty (and not `.`); otherwise
/// `"chrom:start-end"`.  Lines beginning with `#` are skipped by the noodles
/// parser.
///
/// Uses `bed::io::Reader<3, _>` so extra columns (5+) are silently accepted.
/// Coordinate conversion: noodles adds 1 to BED's 0-based start when storing
/// it as `Position`, so `start_0 = feature_start().get() - 1`; the end is
/// stored verbatim, so `end_0 = feature_end().get()`.
pub fn parse_bed_file(path: &Path) -> Result<(RegionIndex, Vec<VarMeta>)> {
    let (intervals_by_chrom, var) = parse_bed_raw(path)?;
    Ok((build_region_index(intervals_by_chrom), var))
}

/// Parse a BED file into per-chromosome **unsorted** interval lists plus the
/// ordered per-feature [`VarMeta`].
///
/// Returned separately from index construction so callers that merge several
/// files can build each [`ChromIndex`] exactly once instead of building and
/// then immediately rebuilding it.
fn parse_bed_raw(path: &Path) -> Result<RawIntervals> {
    let file =
        File::open(path).with_context(|| format!("failed to open BED file: {}", path.display()))?;

    let mut reader = bed::io::Reader::<3, _>::new(BufReader::new(file));
    let mut record = bed::Record::<3>::default();

    let mut intervals_by_chrom: AHashMap<String, Vec<Interval>> = AHashMap::new();
    let mut var: Vec<VarMeta> = Vec::new();
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

        var.push(VarMeta {
            chrom: chrom.clone(),
            start,
            end,
            name: name.clone(),
            strand,
        });
        intervals_by_chrom.entry(chrom).or_default().push(Interval {
            start,
            end,
            val: region_idx,
            name,
        });
        region_idx += 1;
    }

    Ok((intervals_by_chrom, var))
}

/// Build a [`RegionIndex`] from per-chromosome unsorted interval lists,
/// sorting each chromosome's intervals exactly once.
fn build_region_index(intervals_by_chrom: AHashMap<String, Vec<Interval>>) -> RegionIndex {
    intervals_by_chrom
        .into_iter()
        .map(|(c, ivs)| (c, ChromIndex::build(ivs)))
        .collect()
}

// GTF / GFF3 parsers
// Both `gtf::io::Reader::record_bufs()` and `gff::io::Reader::record_bufs()`
// yield `gff::feature::RecordBuf`, so a single generic helper covers both.
fn build_annotation_index_raw<I>(
    iter: I,
    feature_type: Option<&str>,
    name_attr: &str,
) -> Result<RawIntervals>
where
    I: Iterator<Item = io::Result<gff::feature::RecordBuf>>,
{
    let mut intervals_by_chrom: AHashMap<String, Vec<Interval>> = AHashMap::new();
    let mut var: Vec<VarMeta> = Vec::new();
    let mut region_idx: usize = 0;

    for result in iter {
        let record = result.context("failed to read annotation record")?;

        if let Some(ft) = feature_type
            && record.ty() != ft.as_bytes()
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

        // attributes().get() uses the inherent method on record_buf::Attributes
        // (takes &[u8], returns Option<&Value>).  Value is String(BString) or
        // Array(Vec<BString>); we take the first available string.
        let name = record
            .attributes()
            .get(name_attr.as_bytes())
            .and_then(|v| {
                v.as_string()
                    .or_else(|| v.as_array().and_then(|a| a.first()).map(|s| s.as_ref()))
            })
            .map(|s| s.to_string())
            .unwrap_or_else(|| format!("{}:{}-{}", chrom, start, end));

        let strand = gff_strand_to_char(record.strand());

        var.push(VarMeta {
            chrom: chrom.clone(),
            start,
            end,
            name: name.clone(),
            strand,
        });
        intervals_by_chrom.entry(chrom).or_default().push(Interval {
            start,
            end,
            val: region_idx,
            name,
        });
        region_idx += 1;
    }

    Ok((intervals_by_chrom, var))
}

/// Parse a GTF file and build a per-chromosome interval index.
///
/// `feature_type` filters by the feature-type column (e.g. `"gene"`, `"exon"`);
/// `None` includes all types.  `name_attr` is the attribute key used as the
/// feature name (e.g. `"gene_id"`, `"gene_name"`); falls back to
/// `"chrom:start-end"` when absent.  Coordinates are converted from 1-based
/// inclusive to 0-based half-open.
pub fn parse_gtf_file(
    path: &Path,
    feature_type: Option<&str>,
    name_attr: &str,
) -> Result<(RegionIndex, Vec<VarMeta>)> {
    let (intervals_by_chrom, var) = parse_gtf_raw(path, feature_type, name_attr)?;
    Ok((build_region_index(intervals_by_chrom), var))
}

/// Parse a GTF file into per-chromosome **unsorted** interval lists plus var metadata.
fn parse_gtf_raw(path: &Path, feature_type: Option<&str>, name_attr: &str) -> Result<RawIntervals> {
    let file =
        File::open(path).with_context(|| format!("failed to open GTF file: {}", path.display()))?;
    let mut reader = gtf::io::Reader::new(BufReader::new(file));
    build_annotation_index_raw(reader.record_bufs(), feature_type, name_attr)
}

/// Parse a GFF3 file and build a per-chromosome interval index.
///
/// `feature_type` filters by the type column (e.g. `"gene"`, `"mRNA"`);
/// `None` includes all types.  `name_attr` is the attribute tag used as the
/// feature name (e.g. `"ID"`, `"Name"`); falls back to `"chrom:start-end"`.
/// Coordinates are converted from 1-based inclusive to 0-based half-open.
pub fn parse_gff_file(
    path: &Path,
    feature_type: Option<&str>,
    name_attr: &str,
) -> Result<(RegionIndex, Vec<VarMeta>)> {
    let (intervals_by_chrom, var) = parse_gff_raw(path, feature_type, name_attr)?;
    Ok((build_region_index(intervals_by_chrom), var))
}

/// Parse a GFF3 file into per-chromosome **unsorted** interval lists plus var metadata.
fn parse_gff_raw(path: &Path, feature_type: Option<&str>, name_attr: &str) -> Result<RawIntervals> {
    let file =
        File::open(path).with_context(|| format!("failed to open GFF file: {}", path.display()))?;
    let mut reader = gff::io::Reader::new(BufReader::new(file));
    build_annotation_index_raw(reader.record_bufs(), feature_type, name_attr)
}

/// Parse one or more annotation files (BED, GTF, GFF3 — any mix) and merge
/// their intervals into a single [`RegionIndex`].
///
/// Format is detected from each file's extension:
///
/// | Extension(s)                | Parser | Default `feature_type`             | Default `name_attr`   |
/// |-----------------------------|--------|------------------------------------|-----------------------|
/// | `.gtf`, `.gtf.gz`           | GTF    | `"transcript"` / `"exon"`*         | `"transcript_id"` / `"gene_id"`* |
/// | `.gff`, `.gff3`, `.gff3.gz` | GFF3   | auto (`"transcript"`/`"mRNA"`)* / `"exon"`* | `"ID"` / `"Parent"`* |
/// | anything else               | BED    | —                                  | —                     |
///
/// \* When `metagene = true` the exon type and gene-grouping attribute are used
/// instead of the transcript-level defaults.
///
/// Feature names are concatenated in file order.  Each file's `val` indices are
/// shifted by the cumulative feature count so every feature in the combined
/// index maps to a unique slot in the returned name vector.
///
/// `feature_type` and `name_attr` only affect GTF / GFF3 files; BED files
/// always use the 4th column (or `"chrom:start-end"`) as the feature name.
/// `None` selects the per-format default shown above.
///
/// When `metagene = true`, exon intervals are grouped by their gene-id
/// attribute so that every exon belonging to the same gene shares one `val`
/// (= gene index).  This lets the counting loop count a read for a gene
/// regardless of how many of its exons the read overlaps.
pub fn parse_annotation_files<P: AsRef<Path>>(
    paths: impl IntoIterator<Item = P>,
    feature_type: Option<&str>,
    name_attr: Option<&str>,
    metagene: bool,
) -> Result<(RegionIndex, Vec<VarMeta>)> {
    let mut all_var: Vec<VarMeta> = Vec::new();
    let mut intervals_by_chrom: AHashMap<String, Vec<Interval>> = AHashMap::new();

    for path in paths {
        let path = path.as_ref();

        let filename = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let stem = filename.strip_suffix(".gz").unwrap_or(filename);
        let is_gtf = stem.ends_with(".gtf");
        let is_gff = stem.ends_with(".gff") || stem.ends_with(".gff3");

        // Metagene mode groups exon records by gene-id attribute so all exons
        // of the same gene share one val (= gene index in `all_var`).
        // The offset shift below still applies so multi-file merges work.
        if metagene && (is_gtf || is_gff) {
            let exon_type = feature_type.unwrap_or("exon");
            let gene_attr = if is_gtf {
                name_attr.unwrap_or("gene_id")
            } else {
                name_attr.unwrap_or("Parent")
            };
            let offset = all_var.len();
            let (file_ivs, var) = if is_gtf {
                let file = File::open(path)
                    .with_context(|| format!("failed to open GTF file: {}", path.display()))?;
                let mut reader = gtf::io::Reader::new(BufReader::new(file));
                build_metagene_index_raw(reader.record_bufs(), exon_type, gene_attr)
            } else {
                let file = File::open(path)
                    .with_context(|| format!("failed to open GFF file: {}", path.display()))?;
                let mut reader = gff::io::Reader::new(BufReader::new(file));
                build_metagene_index_raw(reader.record_bufs(), exon_type, gene_attr)
            }
            .with_context(|| format!("failed to parse annotation file: {}", path.display()))?;

            all_var.extend(var);
            for (chrom, ivs) in file_ivs {
                let bucket = intervals_by_chrom.entry(chrom).or_default();
                for mut iv in ivs {
                    iv.val += offset;
                    bucket.push(iv);
                }
            }
            continue;
        }

        let offset = all_var.len();

        // Per-format defaults.  GTF: region type "transcript", name attribute
        // "transcript_id".  GFF3: name attribute "ID"; the region type is
        // auto-detected (GENCODE-style uses "transcript", Ensembl/RefSeq use
        // "mRNA") when the caller didn't specify one.
        let (file_intervals, var) = if is_gtf {
            parse_gtf_raw(
                path,
                Some(feature_type.unwrap_or("transcript")),
                name_attr.unwrap_or("transcript_id"),
            )
        } else if is_gff {
            let transcript_type = match feature_type {
                Some(ft) => ft,
                None => detect_gff_transcript_type(path)?,
            };
            parse_gff_raw(path, Some(transcript_type), name_attr.unwrap_or("ID"))
        } else {
            parse_bed_raw(path)
        }
        .with_context(|| format!("failed to parse annotation file: {}", path.display()))?;

        all_var.extend(var);

        // Shift each file's feature indices by the running feature count so
        // every feature maps to a unique column in the combined var metadata.
        for (chrom, ivs) in file_intervals {
            let bucket = intervals_by_chrom.entry(chrom).or_default();
            for mut iv in ivs {
                iv.val += offset;
                bucket.push(iv);
            }
        }
    }

    Ok((build_region_index(intervals_by_chrom), all_var))
}

/// Group exon-type records by a gene-id attribute so all exons of the same
/// gene share one region index (`val`).
///
/// This is the metagene parsing path.  The resulting intervals have
/// `val` = gene index (into the returned `VarMeta` vec) rather than a
/// per-exon index, so `build_counting_index` (blacklist subtraction) and the
/// counting inner loop naturally treat all exons of one gene as one region.
fn build_metagene_index_raw<I>(iter: I, exon_type: &str, gene_id_attr: &str) -> Result<RawIntervals>
where
    I: Iterator<Item = io::Result<gff::feature::RecordBuf>>,
{
    let mut intervals_by_chrom: AHashMap<String, Vec<Interval>> = AHashMap::new();
    let mut var: Vec<VarMeta> = Vec::new();
    let mut gene_to_idx: AHashMap<String, usize> = AHashMap::new();

    for result in iter {
        let record = result.context("failed to read annotation record")?;

        if record.ty() != exon_type.as_bytes() {
            continue;
        }

        let chrom = record.reference_sequence_name().to_string();
        let start = record.start().get() - 1;
        let end = record.end().get();

        if end <= start {
            continue;
        }

        let gene_name = record
            .attributes()
            .get(gene_id_attr.as_bytes())
            .and_then(|v| {
                v.as_string()
                    .or_else(|| v.as_array().and_then(|a| a.first()).map(|s| s.as_ref()))
            })
            .map(|s| s.to_string())
            .unwrap_or_else(|| format!("{}:{}-{}", chrom, start, end));

        let gene_idx = match gene_to_idx.get(&gene_name).copied() {
            Some(idx) => {
                // Extend the gene's genomic span to cover this exon.
                var[idx].start = var[idx].start.min(start);
                var[idx].end = var[idx].end.max(end);
                idx
            }
            None => {
                let idx = var.len();
                gene_to_idx.insert(gene_name.clone(), idx);
                var.push(VarMeta {
                    chrom: chrom.clone(),
                    start,
                    end,
                    name: gene_name.clone(),
                    strand: gff_strand_to_char(record.strand()),
                });
                idx
            }
        };

        intervals_by_chrom.entry(chrom).or_default().push(Interval {
            start,
            end,
            val: gene_idx,
            name: gene_name,
        });
    }

    Ok((intervals_by_chrom, var))
}

/// Pick the region (transcript) feature type for a GFF3 file when the caller
/// didn't specify one.  GFF3 flavors disagree on the name: GENCODE uses
/// `transcript`, while Ensembl / RefSeq use `mRNA`.  We prefer `transcript`,
/// then `mRNA`, and fall back to `transcript` (yielding no features, so the
/// user knows to pass an explicit type) when neither is present.
///
/// Note: this scans the file's column-3 values once before parsing.
fn detect_gff_transcript_type(path: &Path) -> Result<&'static str> {
    let present = scan_feature_types(path)?;
    for candidate in ["transcript", "mRNA"] {
        if present.contains(candidate) {
            return Ok(candidate);
        }
    }
    Ok("transcript")
}

/// Collect the set of column-3 feature types present in a GFF/GTF file via a
/// fast line scan (no full record parsing).
fn scan_feature_types(path: &Path) -> Result<AHashSet<String>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open annotation file: {}", path.display()))?;
    let mut types: AHashSet<String> = AHashSet::new();
    for line in BufReader::new(file).lines() {
        let line = line.context("failed to read annotation line")?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(col3) = line.split('\t').nth(2)
            && !types.contains(col3)
        {
            types.insert(col3.to_string());
        }
    }
    Ok(types)
}

/// Build a "counting index" from a feature index by subtracting blacklisted
/// regions from every interval.
///
/// Each original interval `[start, end)` with feature-index `val` is replaced
/// by the sub-intervals that remain after removing any overlapping blacklist
/// segments.  All sub-intervals inherit the same `val`, so a single
/// `ChromIndex::find` query per read is sufficient to count into the correct
/// matrix column.
///
/// If `blacklist` is `None` the returned index is structurally identical to
/// `feature_index`.
pub fn build_counting_index(
    feature_index: &RegionIndex,
    blacklist: Option<&RegionIndex>,
) -> RegionIndex {
    feature_index
        .iter()
        .map(|(chrom, chrom_idx)| {
            let new_intervals: Vec<Interval> = chrom_idx
                .iter()
                .flat_map(|iv| {
                    valid_sub_intervals(iv.start, iv.end, chrom, blacklist)
                        .into_iter()
                        .map(move |(s, e)| Interval {
                            start: s,
                            end: e,
                            val: iv.val,
                            name: iv.name.clone(),
                        })
                })
                .collect();
            (chrom.clone(), ChromIndex::build(new_intervals))
        })
        .collect()
}

/// Return the portions of `[feat_start, feat_end)` not covered by the
/// blacklist on `chrom`.  Returns the original interval unchanged when there
/// is no blacklist or no overlap.
fn valid_sub_intervals(
    feat_start: usize,
    feat_end: usize,
    chrom: &str,
    blacklist: Option<&RegionIndex>,
) -> Vec<(usize, usize)> {
    let Some(bl_index) = blacklist else {
        return vec![(feat_start, feat_end)];
    };
    let Some(bl_chrom_idx) = bl_index.get(chrom) else {
        return vec![(feat_start, feat_end)];
    };

    // Collect blacklist segments overlapping the feature, clamped to its bounds.
    let mut bl_segs: Vec<(usize, usize)> = bl_chrom_idx
        .find(feat_start, feat_end)
        .map(|iv| (iv.start.max(feat_start), iv.end.min(feat_end)))
        .collect();

    if bl_segs.is_empty() {
        return vec![(feat_start, feat_end)];
    }

    // Sort then merge overlapping / adjacent blacklist segments.
    bl_segs.sort_unstable();
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (s, e) in bl_segs {
        match merged.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => merged.push((s, e)),
        }
    }

    // Valid sub-intervals are the gaps between merged blacklist segments.
    let mut valid: Vec<(usize, usize)> = Vec::new();
    let mut cursor = feat_start;
    for (bl_s, bl_e) in merged {
        if cursor < bl_s {
            valid.push((cursor, bl_s));
        }
        cursor = cursor.max(bl_e);
    }
    if cursor < feat_end {
        valid.push((cursor, feat_end));
    }
    valid
}
