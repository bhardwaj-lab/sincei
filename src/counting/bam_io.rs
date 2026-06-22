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
use noodles::sam::header::ReferenceSequences;
use noodles::sam::header::record::value::{Map, map::ReferenceSequence};

use super::filters::MotifFilter;

/// Concrete type of a BAI-indexed BAM reader opened from a file path.
pub(crate) type BamReader = bam::io::IndexedReader<bgzf::io::Reader<File>>;

/// Read just the header of a BAM file, tolerating non-compliant SAM header text
/// (see module docs).  Builds an indexed reader so a missing `.bai` is still
/// reported here rather than later.
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
/// makes it reject some real-world BAMs, notably, 10x Genomics output, whose
/// `@HD` line lacks a conforming `VN` version field.
/// We fall back to reconstructing the header straight from the BAM **binary**
/// reference dictionary, which is a simple, well-defined block that 10x writes
/// correctly.
///
/// The fallback only triggers when the strict parse fails, so fully compliant
/// BAMs are read exactly as before (preserving `@RG`/`@PG`/etc.).
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

/// Per-worker (per-thread) reusable state for the chunk-parallel counters.
///
/// Rayon's `map_init` hands one `BamWorker` to each thread, which then
/// processes many chunks.  This amortizes two costs:
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
