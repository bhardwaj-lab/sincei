//! Open BAM files, and reuse readers across work chunks.
//!
//! Every reader here is BAI-indexed: the counting loops query one chromosome
//! chunk at a time rather than streaming whole files. `BamWorker` holds a
//! reader (and motif filter) per rayon thread, so a chunk does not pay to
//! reopen the file it was handed.
//!
//! Header reading falls back to the BAM binary reference dictionary when the
//! SAM header text does not satisfy noodles' strict parsing (like 10x Genomics
//! BAM files).

use std::ffi::CStr;
use std::fs::File;
use std::io::{self, Read};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::{Mutex, OnceLock};

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result};
use bstr::BString;
use noodles::bam;
use noodles::bgzf;
use noodles::sam::Header;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::header::record::value::{Map, map::ReadGroup, map::ReferenceSequence};
use noodles::sam::header::{ReadGroups, ReferenceSequences};
use twobit::TwoBitFile;

use super::filters::MotifFilter;

/// Alias for a BAI-indexed BAM reader opened from a file path.
pub(crate) type BamReader = bam::io::IndexedReader<bgzf::io::Reader<File>>;

/// Records inspected when checking that the BAM tags a run depends on are present.
///
/// Large enough that a tag used by the file will certainly appear, small enough
/// to be unnoticeable next to the run it guards.
const TAG_SAMPLE_READS: usize = 100_000;

/// Check the barcode tag, and the UMI tag when one was asked for, on every BAM.
pub(crate) fn ensure_barcode_tags_present(
    bam_paths: &[&Path],
    bc_tag: Tag,
    umi_tag: Option<Tag>,
) -> Result<()> {
    for path in bam_paths {
        ensure_tag_present(path, "--cellTag", bc_tag)?;
        if let Some(umi) = umi_tag {
            ensure_tag_present(path, "--umiTag", umi)?;
        }
    }
    Ok(())
}

/// Fail when `tag` is on none of the first [`TAG_SAMPLE_READS`] records.
fn ensure_tag_present(path: &Path, option: &str, tag: Tag) -> Result<()> {
    // Goes through `open_indexed_bam` rather than a plain reader so a BAM with a
    // non-compliant SAM header still gets checked; otherwise this would reject
    // every file the binary-dictionary fallback exists to support.
    let (mut reader, header) = open_indexed_bam(path)?;

    let mut seen = 0usize;
    for name in header.reference_sequences().keys() {
        let region: noodles::core::Region = String::from_utf8_lossy(name.as_ref())
            .parse()
            .with_context(|| format!("failed to build a query for reference {name:?}"))?;
        let Ok(query) = reader.query(&header, &region) else {
            continue;
        };

        for result in query.records() {
            let record = result.with_context(|| format!("failed to read {}", path.display()))?;
            for field in record.data().iter() {
                let (found, _) =
                    field.with_context(|| format!("failed to read {}", path.display()))?;
                if found == tag {
                    return Ok(());
                }
            }

            seen += 1;
            if seen >= TAG_SAMPLE_READS {
                break;
            }
        }
        if seen >= TAG_SAMPLE_READS {
            break;
        }
    }

    // Reaching here means the BAM held no such tag, because the records carry
    // other tags, carry none at all, or because there are no records.
    let name: &[u8; 2] = tag.as_ref();
    anyhow::bail!(
        "{} has no {} tag (from {}) in its first {} records.\n\
         Check which tag the BAM uses, e.g. with `samtools view {} | head -1`",
        path.display(),
        String::from_utf8_lossy(name),
        option,
        TAG_SAMPLE_READS,
        path.display()
    )
}

/// Fail unless the 2bit genome describes the same chromosomes as the BAMs.
///
/// The motif filter reads reference bases at each read's coordinates, so a
/// genome from another assembly returns the wrong sequence for every read. It
/// also puts reads past the end of a chromosome, where the reference
/// implementation dies inside py2bit and this one would quietly treat the
/// failed lookup as "no match" and drop the read.
///
/// Only the chromosomes named by both are compared.
pub(crate) fn ensure_genome_matches_bams(genome_path: &Path, bam_paths: &[&Path]) -> Result<()> {
    let genome = TwoBitFile::open(genome_path)
        .with_context(|| format!("failed to open 2bit genome: {}", genome_path.display()))?;
    let lengths: AHashMap<String, usize> = genome
        .chrom_names()
        .into_iter()
        .zip(genome.chrom_sizes())
        .collect();

    for path in bam_paths {
        let header = read_bam_header(path)?;
        let mut shared = 0usize;

        for (name, seq) in header.reference_sequences() {
            let chrom = String::from_utf8_lossy(name.as_ref()).into_owned();
            let Some(&genome_len) = lengths.get(&chrom) else {
                continue;
            };
            shared += 1;
            let bam_len = seq.length().get();
            anyhow::ensure!(
                bam_len == genome_len,
                "{} and {} disagree about {}: {} bases in the BAM, {} in the genome.\n\
                 They are different assemblies, so --motifFilter would read the wrong \
                 bases for every read.",
                path.display(),
                genome_path.display(),
                chrom,
                bam_len,
                genome_len
            );
        }

        anyhow::ensure!(
            shared > 0,
            "{} and {} name no chromosome in common, so --motifFilter has no \
             sequence to read.\nThe genome has {}; check whether one side uses a \
             `chr` prefix and the other does not.",
            path.display(),
            genome_path.display(),
            genome
                .chrom_names()
                .iter()
                .take(3)
                .cloned()
                .collect::<Vec<_>>()
                .join(", ")
        );
    }

    Ok(())
}

/// The `@RG` IDs a BAM declares, in header order.
///
/// This is the group axis for `--groupTag`: a merged BAM records which source
/// file each read came from, and the header lists the possible values. Reading
/// them here means the row space is known before a single record is parsed, so
/// the chunk-parallel counting loop needs only a read-only lookup.
///
/// Errors when the BAM declares none, because that is indistinguishable from an
/// un-merged file and counting one of those per-barcode would silently merge
/// cells that share a barcode across samples.
pub(crate) fn read_group_ids(header: &Header, path: &Path) -> Result<Vec<Vec<u8>>> {
    let ids: Vec<Vec<u8>> = header.read_groups().keys().map(|id| id.to_vec()).collect();

    anyhow::ensure!(
        !ids.is_empty(),
        "{} declares no @RG read groups, so --groupTag cannot tell its samples \
         apart.\nMerge the inputs with read groups, e.g. `samtools merge -r`",
        path.display()
    );
    Ok(ids)
}

/// Warn once per unrecognised group value, then stay quiet.
///
/// A read whose group tag is not among the BAM's `@RG` IDs is skipped. That is
/// worth saying, because it means data is being dropped -- but a merged file
/// with an unexpected group would otherwise emit the message on every record,
/// so each distinct value is reported exactly once.
pub(crate) fn warn_unknown_group(group: &[u8]) {
    static SEEN: OnceLock<Mutex<AHashSet<Vec<u8>>>> = OnceLock::new();
    let seen = SEEN.get_or_init(|| Mutex::new(AHashSet::new()));

    let mut seen = match seen.lock() {
        Ok(guard) => guard,
        // A poisoned lock only means some other thread panicked while warning;
        // losing the dedup is better than propagating that here.
        Err(poisoned) => poisoned.into_inner(),
    };
    if seen.insert(group.to_vec()) {
        eprintln!(
            "warning: read group {:?} is not declared in the BAM header; \
             those reads are being skipped",
            String::from_utf8_lossy(group)
        );
    }
}

/// Read just the header of a BAM file, tolerating non-compliant SAM header.
/// Builds an indexed reader so a missing `.bai` is reported as an error.
pub(crate) fn read_bam_header(path: &Path) -> Result<Header> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .with_context(|| {
            format!(
                "failed to open indexed BAM (does the .bai index file exist?): {}",
                path.display()
            )
        })?;
    match reader.read_header() {
        Ok(header) => Ok(header),
        Err(_) => read_header_from_binary_dict(path),
    }
}

/// Open a BAI-indexed BAM reader together with its header, tolerating
/// non-compliant SAM header text.
///
/// noodles validates the SAM header text strictly (per the hts-spec), which
/// makes it reject some real-world BAMs, such as 10x Genomics output, whose
/// `@HD` line lacks a conforming `VN` version field.
/// We fall back to reconstructing the header straight from the BAM **binary**
/// reference dictionary, which is a simple, well-defined block that 10x writes
/// correctly.
pub(crate) fn open_indexed_bam(path: &Path) -> Result<(BamReader, Header)> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(path)
        .with_context(|| {
            format!(
                "failed to open indexed BAM (does the .bai index file exist?): {}",
                path.display()
            )
        })?;
    let header = match reader.read_header() {
        Ok(header) => header,
        Err(_) => read_header_from_binary_dict(path)?,
    };
    Ok((reader, header))
}

/// Reconstruct a [`Header`] from the BAM binary reference dictionary, rather
/// than from the SAM header text.
///
/// BAM header layout (after BGZF decompression):
/// `magic[4] "BAM\1"`, `l_text: u32`, `text[l_text]`, `n_ref: u32`, then per
/// reference: `l_name: u32`, `name[l_name]` (NUL-terminated), `l_ref: u32`.
///
/// The header text is still read, but only to recover `@RG` records (see
/// [`read_groups_from_text`]); everything else comes from the binary dictionary.
fn read_header_from_binary_dict(path: &Path) -> Result<Header> {
    let file =
        File::open(path).with_context(|| format!("failed to open BAM: {}", path.display()))?;
    let mut reader = bgzf::io::Reader::new(file);

    let mut magic = [0u8; 4];
    reader
        .read_exact(&mut magic)
        .context("failed to read BAM magic number")?;
    anyhow::ensure!(&magic == b"BAM\x01", "not a BAM file: {}", path.display());

    // Read (rather than skip) the non-compliant SAM header text: the strict
    // parser generally rejects the @HD line, so the @RG records are generally
    // recoverable.
    let l_text = read_u32(&mut reader)? as usize;
    let mut text = vec![0u8; l_text];
    reader
        .read_exact(&mut text)
        .context("failed to read BAM header text")?;
    let read_groups = read_groups_from_text(&text)?;

    let n_ref = read_u32(&mut reader)?;
    let mut reference_sequences = ReferenceSequences::with_capacity(n_ref as usize);
    for _ in 0..n_ref {
        let l_name = read_u32(&mut reader)? as usize;
        let mut name_buf = vec![0u8; l_name];
        reader
            .read_exact(&mut name_buf)
            .context("failed to read reference name")?;
        let name: BString = CStr::from_bytes_with_nul(&name_buf)
            .map_err(|e| anyhow::anyhow!("invalid reference name in BAM header: {e}"))?
            .to_bytes()
            .into();

        let l_ref = read_u32(&mut reader)? as usize;
        let length = NonZeroUsize::new(l_ref)
            .with_context(|| format!("reference {name:?} has zero length"))?;
        reference_sequences.insert(name, Map::<ReferenceSequence>::new(length));
    }

    let mut header = Header::default();
    *header.reference_sequences_mut() = reference_sequences;
    *header.read_groups_mut() = read_groups;
    Ok(header)
}

/// Recover `@RG` records from raw SAM header text.
///
/// Order is first-seen, so the group index built from this is deterministic.
fn read_groups_from_text(text: &[u8]) -> Result<ReadGroups> {
    let mut read_groups = ReadGroups::new();

    for line in text.split(|&b| b == b'\n') {
        let Some(fields) = line.strip_prefix(b"@RG\t") else {
            continue;
        };
        let Some(id) = fields
            .split(|&b| b == b'\t')
            .find_map(|field| field.strip_prefix(b"ID:"))
        else {
            continue;
        };
        // Tolerate CRLF line endings.
        let id: BString = id.strip_suffix(b"\r").unwrap_or(id).into();

        anyhow::ensure!(
            !read_groups.contains_key(&id),
            "duplicate read group ID in BAM header: {id}"
        );
        read_groups.insert(id, Map::<ReadGroup>::default());
    }

    Ok(read_groups)
}

fn read_u32<R: Read>(reader: &mut R) -> io::Result<u32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

/// Per-worker (per-thread) reusable state for chunk-parallel counters.
///
/// Rayon's `map_init` hands one `BamWorker` to each thread, which then
/// processes many chunks. This amortizes two costs:
///   * opening the BAM file and parsing its header, and
///   * opening the 2bit genome for the motif filter.
///
/// The reader is keyed by path so a thread that alternates between input BAMs
/// only reopens when the path actually changes.
pub(crate) struct BamWorker<'a> {
    path: Option<&'a Path>,
    reader: Option<BamReader>,
    header: Option<Header>,
    motif: Option<MotifFilter>,
}

impl<'a> BamWorker<'a> {
    pub(crate) fn new() -> Self {
        Self {
            path: None,
            reader: None,
            header: None,
            motif: None,
        }
    }

    /// Ensure the reader is open for `path` (reopening only on a path change)
    /// and, when `motif_ingredients` is supplied, that the motif filter has
    /// been constructed once for this worker, then borrow the reader, header,
    /// and optional motif filter as disjoint parts so all three can be used
    /// simultaneously inside the per-chunk loop.
    pub(crate) fn prepare(
        &mut self,
        path: &'a Path,
        motif_ingredients: Option<(&Path, &[(String, String)])>,
    ) -> Result<(&mut BamReader, &Header, &mut Option<MotifFilter>)> {
        if self.path != Some(path) {
            let (reader, header) = open_indexed_bam(path)?;
            self.reader = Some(reader);
            self.header = Some(header);
            self.path = Some(path);
        }
        if self.motif.is_none()
            && let Some((genome, motifs)) = motif_ingredients
        {
            self.motif = Some(MotifFilter::new(genome, motifs.to_vec())?);
        }
        Ok((
            self.reader.as_mut().expect("just set above"),
            self.header.as_ref().expect("just set above"),
            &mut self.motif,
        ))
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn a_genome_whose_chromosome_lengths_match_the_bam_is_accepted() {
        // test_i1.bam is mm10 chr5 (151,834,684), so a genome accepted here has
        // to carry chromosome 5: the check rejects a pair naming no chromosome
        // in common.
        //
        // NOTE: `mm10_chr1.2bit` is the *extended* suite's genome
        // (benchmarks/pipeline/cases.py), whose BAMs are on chr1. It does not
        // match this BAM, so this test would fail rather than pass if that file
        // were ever placed in tests/testdata.
        let bam = testdata().join("test_i1.bam");
        let genome = testdata().join("mm10_chr1.2bit");
        if !genome.exists() {
            return; // the 2bit is not committed; the extended suite covers it
        }
        assert!(ensure_genome_matches_bams(&genome, &[bam.as_path()]).is_ok());
    }

    #[test]
    fn a_missing_genome_file_is_reported_with_its_path() {
        let bam = testdata().join("test_i1.bam");
        let missing = testdata().join("no_such_genome.2bit");
        let err = ensure_genome_matches_bams(&missing, &[bam.as_path()])
            .unwrap_err()
            .to_string();
        assert!(err.contains("no_such_genome.2bit"), "{err}");
    }

    use super::*;
    use std::path::PathBuf;
    use tempfile::TempDir;

    /// `test_i1.bam` is a coordinate-sorted, BAI-indexed slice of chromosome 5.
    /// Its reads carry the barcode in `BC` and the UMI in `RX`.
    fn testdata() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/testdata")
    }

    fn test_bam() -> PathBuf {
        testdata().join("test_i1.bam")
    }

    fn tag(a: u8, b: u8) -> Tag {
        Tag::new(a, b)
    }

    // Header reading

    #[test]
    fn the_header_carries_the_reference_sequences() {
        let header = read_bam_header(&test_bam()).unwrap();
        let refs = header.reference_sequences();

        assert!(!refs.is_empty(), "expected at least one reference sequence");
        let (name, seq) = refs.get_index(0).unwrap();
        assert_eq!(name.as_slice(), b"5");
        assert_eq!(usize::from(seq.length()), 151_834_684);
    }

    #[test]
    fn a_missing_bam_is_reported_with_its_path() {
        let err = read_bam_header(Path::new("/nonexistent/reads.bam"))
            .unwrap_err()
            .to_string();
        assert!(err.contains("/nonexistent/reads.bam"), "{err}");
    }

    #[test]
    fn a_bam_without_its_index_is_rejected_and_says_so() {
        // The fallback path must not hide a missing .bai: copy the BAM alone.
        let dir = TempDir::new().unwrap();
        let unindexed = dir.path().join("no_index.bam");
        std::fs::copy(test_bam(), &unindexed).unwrap();

        let err = read_bam_header(&unindexed).unwrap_err().to_string();
        assert!(err.contains(".bai"), "{err}");
    }

    // The binary-dictionary fallback

    #[test]
    fn the_binary_dictionary_reconstructs_the_same_references_as_the_header_text() {
        // 10x BAMs fail noodles' strict SAM header parsing, so the dictionary is
        // read instead. On a compliant BAM both routes must agree.
        let from_text = read_bam_header(&test_bam()).unwrap();
        let from_dict = read_header_from_binary_dict(&test_bam()).unwrap();

        let text_refs = from_text.reference_sequences();
        let dict_refs = from_dict.reference_sequences();
        assert_eq!(text_refs.len(), dict_refs.len());

        for (name, seq) in text_refs {
            let other = dict_refs
                .get(name)
                .unwrap_or_else(|| panic!("reference {name:?} missing from the binary dictionary"));
            assert_eq!(seq.length(), other.length());
        }
    }

    #[test]
    fn a_file_that_is_not_a_bam_is_rejected() {
        let dir = TempDir::new().unwrap();
        let not_bam = dir.path().join("plain.txt");
        std::fs::write(&not_bam, b"this is not a BAM file").unwrap();

        assert!(read_header_from_binary_dict(&not_bam).is_err());
    }

    // Tag checking

    #[test]
    fn the_tags_the_file_uses_are_accepted() {
        let bam = test_bam();
        let paths = [bam.as_path()];

        ensure_barcode_tags_present(&paths, tag(b'B', b'C'), None).unwrap();
        ensure_barcode_tags_present(&paths, tag(b'B', b'C'), Some(tag(b'R', b'X'))).unwrap();
    }

    #[test]
    fn a_barcode_tag_the_file_does_not_use_is_reported_with_advice() {
        let bam = test_bam();
        let err = ensure_barcode_tags_present(&[bam.as_path()], tag(b'Z', b'Z'), None)
            .unwrap_err()
            .to_string();

        assert!(err.contains("ZZ"), "{err}");
        assert!(err.contains("--cellTag"), "{err}");
        assert!(err.contains("samtools view"), "{err}");
    }

    #[test]
    fn a_missing_umi_tag_is_blamed_on_the_umi_option_not_the_cell_option() {
        let bam = test_bam();
        let err =
            ensure_barcode_tags_present(&[bam.as_path()], tag(b'B', b'C'), Some(tag(b'Z', b'Z')))
                .unwrap_err()
                .to_string();

        assert!(err.contains("--umiTag"), "{err}");
    }

    // Indexed reader and per-thread worker

    #[test]
    fn an_indexed_reader_comes_back_with_its_header() {
        let (_reader, header) = open_indexed_bam(&test_bam()).unwrap();
        assert!(!header.reference_sequences().is_empty());
    }

    #[test]
    fn a_worker_reuses_its_reader_across_calls_for_the_same_path() {
        let bam = test_bam();
        let mut worker = BamWorker::new();

        let first_len = {
            let (_reader, header, motif) = worker.prepare(&bam, None).unwrap();
            assert!(motif.is_none(), "no motif filter was asked for");
            header.reference_sequences().len()
        };

        // Second call for the same path must not reopen, and must agree.
        let (_reader, header, _motif) = worker.prepare(&bam, None).unwrap();
        assert_eq!(header.reference_sequences().len(), first_len);
    }

    #[test]
    fn a_fresh_worker_holds_nothing_until_it_is_prepared() {
        let worker = BamWorker::new();
        assert!(worker.path.is_none());
        assert!(worker.reader.is_none());
        assert!(worker.header.is_none());
        assert!(worker.motif.is_none());
    }

    #[test]
    fn preparing_a_worker_on_a_missing_bam_fails() {
        let mut worker = BamWorker::new();
        let missing = Path::new("/nonexistent/reads.bam");
        assert!(worker.prepare(missing, None).is_err());
    }

    // Read groups

    fn merged_bam() -> PathBuf {
        testdata().join("test_i1_i2.bam")
    }

    #[test]
    fn a_merged_bam_declares_its_source_samples() {
        let header = read_bam_header(&merged_bam()).unwrap();
        let ids = read_group_ids(&header, &merged_bam()).unwrap();
        assert_eq!(ids, vec![b"test_i1".to_vec(), b"test_i2".to_vec()]);
    }

    #[test]
    fn a_bam_without_read_groups_is_rejected_with_advice() {
        // test_i1.bam is a single sample, so it declares no @RG.
        let header = read_bam_header(&test_bam()).unwrap();
        let err = read_group_ids(&header, &test_bam())
            .unwrap_err()
            .to_string();

        assert!(err.contains("no @RG"), "{err}");
        assert!(err.contains("samtools merge -r"), "{err}");
    }

    #[test]
    fn read_groups_survive_a_header_the_strict_parser_rejects() {
        // test_i1_i2_badheader.bam has an @HD line with no VN field, so noodles
        // refuses the header text and the binary-dictionary path runs. The @RG
        // lines are intact, and must still come through.
        let bad = testdata().join("test_i1_i2_badheader.bam");
        let recovered = read_group_ids(&read_bam_header(&bad).unwrap(), &bad).unwrap();
        let clean =
            read_group_ids(&read_bam_header(&merged_bam()).unwrap(), &merged_bam()).unwrap();

        assert_eq!(recovered, clean);
    }

    // Parsing @RG out of raw header text

    #[test]
    fn ids_are_taken_in_first_seen_order() {
        let text = b"@HD\tVN:1.6\n@RG\tID:b\tSM:two\n@SQ\tSN:1\tLN:10\n@RG\tID:a\n";
        let groups = read_groups_from_text(text).unwrap();
        let ids: Vec<&[u8]> = groups.keys().map(|k| k.as_ref()).collect();
        assert_eq!(ids, vec![b"b".as_slice(), b"a".as_slice()]);
    }

    #[test]
    fn records_without_an_id_field_are_skipped() {
        let text = b"@RG\tSM:nameless\n@RG\tID:real\n";
        let groups = read_groups_from_text(text).unwrap();
        assert_eq!(groups.len(), 1);
        assert!(groups.contains_key(b"real".as_slice()));
    }

    #[test]
    fn a_header_with_no_read_groups_yields_none() {
        let text = b"@HD\tVN:1.6\n@SQ\tSN:1\tLN:10\n";
        assert!(read_groups_from_text(text).unwrap().is_empty());
    }

    #[test]
    fn carriage_returns_are_not_part_of_the_id() {
        let text = b"@RG\tID:windows\r\n";
        let groups = read_groups_from_text(text).unwrap();
        assert!(groups.contains_key(b"windows".as_slice()), "{groups:?}");
    }

    #[test]
    fn a_duplicate_read_group_id_is_rejected() {
        // Two groups sharing an ID would silently merge two samples into one row.
        let text = b"@RG\tID:same\n@RG\tID:same\n";
        let err = read_groups_from_text(text).unwrap_err().to_string();
        assert!(err.contains("duplicate read group"), "{err}");
    }
}
