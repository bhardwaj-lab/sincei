use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use bed_utils::intervaltree::{Interval, Lapper};
use noodles::bam;
use noodles::sam::alignment::Record as AlignmentRecord;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use triple_accel::hamming::hamming;

type BinsByBarcode = HashMap<String, HashSet<(String, usize)>>;
type BlacklistIndex = HashMap<String, Lapper<usize, ()>>;

fn run_filter_barcodes(
    bamfile: &Path,
    whitelist: Option<Vec<String>>,
    blacklist_file_name: Option<&Path>,
    cell_tag: &str,
    min_hamming_dist: usize,
    min_count: usize,
    min_mapping_quality: Option<u8>,
    bin_size: usize,
) -> Result<(Vec<(String, usize)>, Vec<String>)> {
    if bin_size == 0 {
        anyhow::bail!("bin_size must be greater than zero");
    }

    let whitelist = whitelist.unwrap_or_default();
    let whitelist_is_active = !whitelist.is_empty();
    let whitelist_exact: HashSet<&str> = whitelist.iter().map(String::as_str).collect();
    let blacklist_index = if let Some(blacklist_file_name) = blacklist_file_name {
        load_blacklist_index(blacklist_file_name)?
    } else {
        HashMap::new()
    };

    let mut reader = bam::io::Reader::new(
        File::open(bamfile)
            .with_context(|| format!("failed to open BAM file {}", bamfile.display()))?,
    );
    let header = reader.read_header().context("failed to read BAM header")?;

    let tag = parse_tag(cell_tag)?;
    let mut bins_by_barcode: BinsByBarcode = HashMap::new();

    for result in reader.records() {
        let record = result.context("failed to read BAM record")?;

        if record.flags().is_unmapped() {
            continue;
        }

        if let Some(min_mapping_quality) = min_mapping_quality {
            match record.mapping_quality() {
                Some(mapping_quality) if mapping_quality.get() >= min_mapping_quality => {}
                _ => continue,
            }
        }

        let Some(chromosome) = record
            .reference_sequence(&header)
            .transpose()
            .context("failed to read record reference sequence")?
            .map(|(name, _)| name.to_string())
        else {
            continue;
        };

        let Some(start) = record
            .alignment_start()
            .transpose()
            .context("failed to read record alignment start")?
        else {
            continue;
        };
        let Some(end) = record
            .alignment_end()
            .transpose()
            .context("failed to read record alignment end")?
        else {
            continue;
        };

        let start = start.get().saturating_sub(1);
        let end: usize = end.get();

        if end <= start {
            continue;
        }

        if is_blacklisted(&blacklist_index, &chromosome, start, end) {
            continue;
        }

        let Some(barcode) = read_cell_barcode(&record, &tag)? else {
            continue;
        };

        if whitelist_is_active
            && !barcode_is_whitelisted(&barcode, &whitelist_exact, &whitelist, min_hamming_dist)
        {
            continue;
        }

        let first_bin = start / bin_size;
        let last_bin = (end - 1) / bin_size;
        let bin_set = bins_by_barcode.entry(barcode).or_default();

        for bin_index in first_bin..=last_bin {
            bin_set.insert((chromosome.clone(), bin_index));
        }
    }

    let mut barcode_counts: Vec<(String, usize)> = bins_by_barcode
        .into_iter()
        .filter_map(|(barcode, bins)| {
            let count = bins.len();
            (count >= min_count).then_some((barcode, count))
        })
        .collect();

    barcode_counts.sort_by(|left, right| right.1.cmp(&left.1).then_with(|| left.0.cmp(&right.0)));

    let selected_barcodes = barcode_counts
        .iter()
        .map(|(barcode, _)| barcode.clone())
        .collect();

    Ok((barcode_counts, selected_barcodes))
}

fn load_blacklist_index(path: &Path) -> Result<BlacklistIndex> {
    let file = File::open(path)
        .with_context(|| format!("failed to open blacklist file {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut intervals_by_chromosome: HashMap<String, Vec<Interval<usize, ()>>> = HashMap::new();

    for (line_number, line) in reader.lines().enumerate() {
        let line =
            line.with_context(|| format!("failed to read blacklist line {}", line_number + 1))?;
        let line = line.trim();

        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let mut fields = line.split('\t');
        let Some(chromosome) = fields.next() else {
            continue;
        };
        let Some(start) = fields.next() else {
            continue;
        };
        let Some(end) = fields.next() else {
            continue;
        };

        let start = start.parse::<usize>().with_context(|| {
            format!(
                "invalid blacklist start position on line {}",
                line_number + 1
            )
        })?;
        let end = end.parse::<usize>().with_context(|| {
            format!("invalid blacklist end position on line {}", line_number + 1)
        })?;

        if end <= start {
            continue;
        }

        intervals_by_chromosome
            .entry(chromosome.to_string())
            .or_default()
            .push(Interval {
                start,
                stop: end,
                val: (),
            });
    }

    Ok(intervals_by_chromosome
        .into_iter()
        .map(|(chromosome, intervals)| (chromosome, Lapper::new(intervals)))
        .collect())
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
    whitelist_exact: &HashSet<&str>,
    whitelist: &[String],
    min_hamming_dist: usize,
) -> bool {
    if min_hamming_dist == 0 {
        return whitelist_exact.contains(barcode);
    }

    whitelist.iter().any(|whitelist_barcode| {
        whitelist_barcode.len() == barcode.len()
            && hamming(barcode.as_bytes(), whitelist_barcode.as_bytes()) as usize
                <= min_hamming_dist
    })
}

fn is_blacklisted(
    blacklist_index: &BlacklistIndex,
    chromosome: &str,
    start: usize,
    end: usize,
) -> bool {
    blacklist_lapper(blacklist_index, chromosome)
        .map(|lapper| lapper.find(start, end).next().is_some())
        .unwrap_or(false)
}

fn blacklist_lapper<'a>(
    blacklist_index: &'a BlacklistIndex,
    chromosome: &str,
) -> Option<&'a Lapper<usize, ()>> {
    blacklist_index
        .get(chromosome)
        .or_else(|| {
            chromosome
                .strip_prefix("chr")
                .and_then(|chrom| blacklist_index.get(chrom))
        })
        .or_else(|| blacklist_index.get(&format!("chr{chromosome}")))
}

#[pyfunction(signature = (
	bamfile,
	whitelist=None,
	blacklist_file_name=None,
	cell_tag="CB",
	min_hamming_dist=0,
	min_count=0,
	min_mapping_quality=None,
	bin_size=100000,
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
) -> PyResult<(Vec<(String, usize)>, Vec<String>)> {
    run_filter_barcodes(
        bamfile.as_path(),
        whitelist,
        blacklist_file_name.as_deref(),
        cell_tag,
        min_hamming_dist,
        min_count,
        min_mapping_quality,
        bin_size,
    )
    .map_err(|error| PyRuntimeError::new_err(error.to_string()))
}
