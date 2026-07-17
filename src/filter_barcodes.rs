use std::path::{Path, PathBuf};

use ahash::{AHashMap, AHashSet};
use anyhow::{Context, Result};
use noodles::bam;
use noodles::sam::alignment::Record as AlignmentRecord;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rayon::prelude::*;
use triple_accel::hamming::hamming;

use crate::annotation::parse_annotation::parse_blacklist_bed;
use crate::annotation::region_index::{ChromIndex, GenomeIndex};
use crate::bam::bam_io::{open_indexed_bam, read_bam_header};

type BinsByBarcode = AHashMap<String, AHashSet<(String, usize)>>;
type BarcodeCounts = (Vec<(String, usize)>, Vec<String>);

fn run_filter_barcodes(
    bamfile: &Path,
    whitelist: Option<Vec<String>>,
    blacklist_file_name: Option<&Path>,
    cell_tag: &str,
    min_hamming_dist: usize,
    min_count: usize,
    min_mapping_quality: Option<u8>,
    bin_size: usize,
    chr_to_skip: &[String],
    num_threads: usize,
    chunk_size: usize,
) -> Result<BarcodeCounts> {
    if bin_size == 0 {
        anyhow::bail!("bin_size must be greater than zero");
    }
    anyhow::ensure!(chunk_size > 0, "chunk_size must be greater than zero");

    let whitelist = whitelist.unwrap_or_default();
    let whitelist_is_active = !whitelist.is_empty();
    let whitelist_set: AHashSet<String> = whitelist.iter().cloned().collect();
    let blacklist_index = if let Some(p) = blacklist_file_name {
        parse_blacklist_bed(p)?
    } else {
        GenomeIndex::new()
    };

    let tag = parse_tag(cell_tag)?;

    let header = read_bam_header(bamfile)?;

    let chrom_sizes: Vec<(String, usize)> = header
        .reference_sequences()
        .iter()
        .filter(|(name, _)| !chr_to_skip.contains(&name.to_string()))
        .map(|(name, seq)| (name.to_string(), seq.length().get()))
        .collect();

    // Build chunk work list sorted by descending size.
    let mut chunks: Vec<(String, usize, usize)> = chrom_sizes
        .iter()
        .flat_map(|(chrom, chrom_len)| {
            (0..*chrom_len)
                .step_by(chunk_size)
                .map(move |start| (chrom.clone(), start, (start + chunk_size).min(*chrom_len)))
        })
        .collect();
    chunks.sort_unstable_by_key(|b| std::cmp::Reverse(b.2 - b.1));

    let n_threads = if num_threads == 0 {
        rayon::current_num_threads()
    } else {
        num_threads
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .context("failed to build thread pool")?;

    let partial_maps: Vec<BinsByBarcode> = pool.install(|| {
        chunks
            .par_iter()
            .map(|(chrom, chunk_start, chunk_end)| -> Result<BinsByBarcode> {
                let (mut local_reader, local_header) = open_indexed_bam(bamfile)?;

                let region_str = format!("{}:{}-{}", chrom, chunk_start + 1, chunk_end);
                let region: noodles::core::Region = region_str
                    .parse()
                    .with_context(|| format!("failed to parse region: {}", region_str))?;

                let query = match local_reader.query(&local_header, &region) {
                    Ok(q) => q,
                    Err(_) => return Ok(AHashMap::new()),
                };

                let mut local_bins: BinsByBarcode = AHashMap::new();

                for result in query.records() {
                    let record = result.context("failed to read BAM record")?;

                    let flags = record.flags();
                    if flags.is_unmapped() {
                        continue;
                    }

                    if let Some(min_mq) = min_mapping_quality {
                        match record.mapping_quality() {
                            Some(mq) if mq.get() >= min_mq => {}
                            _ => continue,
                        }
                    }

                    let Some(aln_start) = record
                        .alignment_start()
                        .transpose()
                        .context("failed to read alignment start")?
                    else {
                        continue;
                    };
                    let Some(aln_end) = record
                        .alignment_end()
                        .transpose()
                        .context("failed to read alignment end")?
                    else {
                        continue;
                    };

                    let start = aln_start.get().saturating_sub(1);
                    // Ownership: a read belongs to the chunk that contains its
                    // alignment_start. The BAI query returns overlapping reads,
                    // so skip anything that started before this chunk.
                    if start < *chunk_start {
                        continue;
                    }

                    let end: usize = aln_end.get();
                    if end <= start {
                        continue;
                    }

                    if is_blacklisted(&blacklist_index, chrom, start, end) {
                        continue;
                    }

                    let Some(barcode) = read_cell_barcode(&record, &tag)? else {
                        continue;
                    };

                    if whitelist_is_active
                        && !barcode_is_whitelisted(
                            &barcode,
                            &whitelist_set,
                            &whitelist,
                            min_hamming_dist,
                        )
                    {
                        continue;
                    }

                    // A read counts once, in the bin holding its start.
                    let bin_idx = start / bin_size;
                    local_bins
                        .entry(barcode)
                        .or_default()
                        .insert((chrom.clone(), bin_idx));
                }

                Ok(local_bins)
            })
            .collect::<Result<Vec<_>>>()
    })?;

    let mut bins_by_barcode: BinsByBarcode = AHashMap::new();
    for partial in partial_maps {
        for (barcode, bins) in partial {
            bins_by_barcode.entry(barcode).or_default().extend(bins);
        }
    }

    // Return *all* detected barcodes with their non-zero-bin counts; `min_count`
    // only determines which are "selected" (it does not drop rows), matching the
    // Python `scFilterBarcodes` output where every detected barcode is listed
    // with a `selected` flag.
    let mut barcode_counts: Vec<(String, usize)> = bins_by_barcode
        .into_iter()
        .map(|(barcode, bins)| (barcode, bins.len()))
        .collect();

    barcode_counts.sort_by(|l, r| r.1.cmp(&l.1).then_with(|| l.0.cmp(&r.0)));

    let selected_barcodes = barcode_counts
        .iter()
        .filter(|(_, count)| *count >= min_count)
        .map(|(bc, _)| bc.clone())
        .collect();

    Ok((barcode_counts, selected_barcodes))
}

fn parse_tag(cell_tag: &str) -> Result<Tag> {
    let bytes = cell_tag.as_bytes();
    if bytes.len() != 2 {
        anyhow::bail!("barcode tag must be exactly two characters");
    }
    Ok(Tag::new(bytes[0], bytes[1]))
}

fn read_cell_barcode(record: &bam::Record, tag: &Tag) -> Result<Option<String>> {
    match record.data().get(tag) {
        Some(Ok(Value::String(value))) => {
            Ok(Some(String::from_utf8_lossy(value.as_ref()).into_owned()))
        }
        Some(Ok(_)) => Ok(None),
        Some(Err(error)) => Err(error.into()),
        None => Ok(None),
    }
}

fn barcode_is_whitelisted(
    barcode: &str,
    whitelist_set: &AHashSet<String>,
    whitelist: &[String],
    min_hamming_dist: usize,
) -> bool {
    if min_hamming_dist == 0 {
        return whitelist_set.contains(barcode);
    }
    whitelist.iter().any(|wl_bc| {
        wl_bc.len() == barcode.len()
            && hamming(barcode.as_bytes(), wl_bc.as_bytes()) as usize <= min_hamming_dist
    })
}

fn is_blacklisted(
    blacklist_index: &GenomeIndex,
    chromosome: &str,
    start: usize,
    end: usize,
) -> bool {
    blacklist_chrom_index(blacklist_index, chromosome)
        .map(|idx| idx.find(start, end).next().is_some())
        .unwrap_or(false)
}

fn blacklist_chrom_index<'a>(
    blacklist_index: &'a GenomeIndex,
    chromosome: &str,
) -> Option<&'a ChromIndex> {
    blacklist_index
        .get(chromosome)
        .or_else(|| {
            chromosome
                .strip_prefix("chr")
                .and_then(|c| blacklist_index.get(c))
        })
        .or_else(|| blacklist_index.get(&format!("chr{chromosome}")))
}

#[pyfunction(signature = (
    bamfile,
    whitelist = None,
    blacklist_file_name = None,
    cell_tag = "CB",
    min_hamming_dist = 0,
    min_count = 0,
    min_mapping_quality = None,
    bin_size = 100_000,
    chr_to_skip = vec![],
    num_threads = 0,
    chunk_size = 1_000_000,
))]
pub fn filter_barcodes(
    bamfile: PathBuf,
    whitelist: Option<Vec<String>>,
    blacklist_file_name: Option<PathBuf>,
    cell_tag: &str,
    min_hamming_dist: usize,
    min_count: usize,
    min_mapping_quality: Option<u8>,
    bin_size: usize,
    chr_to_skip: Vec<String>,
    num_threads: usize,
    chunk_size: usize,
) -> PyResult<BarcodeCounts> {
    run_filter_barcodes(
        bamfile.as_path(),
        whitelist,
        blacklist_file_name.as_deref(),
        cell_tag,
        min_hamming_dist,
        min_count,
        min_mapping_quality,
        bin_size,
        &chr_to_skip,
        num_threads,
        chunk_size,
    )
    .map_err(|e| PyRuntimeError::new_err(format!("{e:#}")))
}
