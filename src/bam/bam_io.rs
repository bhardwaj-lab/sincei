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

use anyhow::{Context, Result};
use bstr::BString;
use noodles::bam;
use noodles::bgzf;
use noodles::sam::Header;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::header::ReferenceSequences;
use noodles::sam::header::record::value::{Map, map::ReferenceSequence};

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
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .with_context(|| format!("failed to open BAM: {}", path.display()))?;
    // The record iterator starts past the header, so it has to be consumed even
    // though nothing here depends on it.
    reader
        .read_header()
        .with_context(|| format!("failed to read BAM header: {}", path.display()))?;

    for result in reader.records().take(TAG_SAMPLE_READS) {
        let record = result.with_context(|| format!("failed to read {}", path.display()))?;
        for field in record.data().iter() {
            let (seen, _) = field.with_context(|| format!("failed to read {}", path.display()))?;
            if seen == tag {
                return Ok(());
            }
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

/// Reconstruct a [`Header`] from the BAM binary reference dictionary, skipping
/// the SAM header text entirely.
///
/// BAM header layout (after BGZF decompression):
/// `magic[4] "BAM\1"`, `l_text: u32`, `text[l_text]`, `n_ref: u32`, then per
/// reference: `l_name: u32`, `name[l_name]` (NUL-terminated), `l_ref: u32`.
fn read_header_from_binary_dict(path: &Path) -> Result<Header> {
    let file =
        File::open(path).with_context(|| format!("failed to open BAM: {}", path.display()))?;
    let mut reader = bgzf::io::Reader::new(file);

    let mut magic = [0u8; 4];
    reader
        .read_exact(&mut magic)
        .context("failed to read BAM magic number")?;
    anyhow::ensure!(&magic == b"BAM\x01", "not a BAM file: {}", path.display());

    // Skip the (non-compliant) SAM header text.
    let l_text = u64::from(read_u32(&mut reader)?);
    io::copy(&mut reader.by_ref().take(l_text), &mut io::sink())
        .context("failed to skip BAM header text")?;

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
    Ok(header)
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
}
