//! Infer fragment size to extend reads to, from the first reads in a BAM file.
//!
//! Covers the behavior of `--extendReads` when given as a flag without an
//! integer value: each BAM's leading records are sampled to tell whether the
//! library is paired-end and, if so, what its median fragment length is.

use std::path::Path;

use anyhow::{Context, Result};
use noodles::bam;

/// `extend_reads` value meaning "estimate from the data".
///
/// A bare `--extendReads` (no value given) is parsed to this, which triggers
/// [`estimate_fragment_length`].
pub(crate) const ESTIMATE_FRAGMENT_LENGTH: usize = 0;

/// Records read from the start of each BAM when inferring its library layout.
///
/// A prefix rather than the whole file: fragment lengths are tightly
/// distributed, so a large sample buys nothing over a bounded one.
pub(crate) const SAMPLE_READS: usize = 100_000;

/// What the first [`SAMPLE_READS`] records of a BAM say about its library.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum LibraryLayout {
    /// No sampled read carries the paired flag.
    SingleEnd,
    /// Sampled reads are paired.  The median holds the TLEN of the properly
    /// paired fragments seen, or `None` when the sample contained none to
    /// measure (e.g. every mate unmapped) — paired, but nothing to estimate from.
    PairedEnd {
        median_fragment_length: Option<usize>,
    },
}

/// Resolve an `--extendReads` request to a concrete fragment length.
///
/// `None` leaves reads unextended; `Some(`[`ESTIMATE_FRAGMENT_LENGTH`]`)` triggers
/// [`estimate_fragment_length`] over `bam_paths`; any other `Some(n)` is used
/// verbatim.
pub(crate) fn resolve_extend_reads(
    extend_reads: Option<usize>,
    bam_paths: &[(&Path, &str)],
) -> Result<Option<usize>> {
    match extend_reads {
        Some(ESTIMATE_FRAGMENT_LENGTH) => Ok(Some(estimate_fragment_length(bam_paths)?)),
        other => Ok(other),
    }
}

/// Estimate a fragment length from the reads at the start of each BAM.
///
/// Every BAM must be paired-end since fragment length cannot be inferred from
/// single-end reads, so this errors and names the offending files rather than
/// guessing.  Otherwise each BAM contributes the median fragment length of its
/// sample, and their mean is the value used to extend reads — one value covers
/// all the BAMs, which are weighted equally regardless of depth.
fn estimate_fragment_length(bam_paths: &[(&Path, &str)]) -> Result<usize> {
    anyhow::ensure!(
        !bam_paths.is_empty(),
        "cannot estimate a fragment length without a BAM file"
    );

    let layouts = sample_library_layouts(bam_paths)?;

    let named = |wanted: fn(&LibraryLayout) -> bool| -> Vec<String> {
        layouts
            .iter()
            .zip(bam_paths)
            .filter(|(layout, _)| wanted(layout))
            .map(|(_, (path, _))| path.display().to_string())
            .collect()
    };

    let single_end = named(|l| matches!(l, LibraryLayout::SingleEnd));
    if !single_end.is_empty() {
        anyhow::bail!(
            "{} contain(s) only single-end reads, please provide a value to `--extendReads`",
            single_end.join(", ")
        );
    }

    let unmeasurable = named(|l| {
        matches!(
            l,
            LibraryLayout::PairedEnd {
                median_fragment_length: None
            }
        )
    });
    if !unmeasurable.is_empty() {
        anyhow::bail!(
            "{} contain(s) no properly paired fragments in the first {} reads, so a fragment \
             length cannot be estimated; please provide a value to `--extendReads`",
            unmeasurable.join(", "),
            SAMPLE_READS
        );
    }

    let medians: Vec<usize> = layouts
        .iter()
        .filter_map(|l| match l {
            LibraryLayout::PairedEnd {
                median_fragment_length,
            } => *median_fragment_length,
            LibraryLayout::SingleEnd => None,
        })
        .collect();

    // Mean of the per-BAM medians, rounded to the nearest base.
    let total: usize = medians.iter().sum();
    Ok((total + medians.len() / 2) / medians.len())
}

/// Read the first [`SAMPLE_READS`] records of each BAM and report its
/// library layout, index-aligned with `bam_paths`.
///
/// A BAM is single-end when no sampled read carries the paired flag.  Otherwise
/// the median is taken over the |TLEN| of properly paired primary reads; only
/// positive TLEN is counted, so each fragment contributes once.
pub(crate) fn sample_library_layouts(bam_paths: &[(&Path, &str)]) -> Result<Vec<LibraryLayout>> {
    bam_paths
        .iter()
        .map(|(path, _)| sample_library_layout(path))
        .collect()
}

fn sample_library_layout(path: &Path) -> Result<LibraryLayout> {
    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .with_context(|| format!("failed to open BAM: {}", path.display()))?;
    // The record iterator starts just past the header, so it must be consumed
    // first even though the layout does not depend on it.
    reader
        .read_header()
        .with_context(|| format!("failed to read BAM header: {}", path.display()))?;

    let mut any_paired = false;
    let mut fragment_lengths: Vec<usize> = Vec::new();

    for result in reader.records().take(SAMPLE_READS) {
        let record =
            result.with_context(|| format!("failed to read BAM record: {}", path.display()))?;
        let flags = record.flags();

        if flags.is_segmented() {
            any_paired = true;
        }
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        if !flags.is_properly_segmented() {
            continue;
        }
        let tlen = record.template_length();
        if tlen > 0 {
            fragment_lengths.push(tlen as usize);
        }
    }

    Ok(if any_paired {
        LibraryLayout::PairedEnd {
            median_fragment_length: median(fragment_lengths),
        }
    } else {
        LibraryLayout::SingleEnd
    })
}

/// Median of `values`, rounding a two-value midpoint up.  `None` when empty.
fn median(mut values: Vec<usize>) -> Option<usize> {
    if values.is_empty() {
        return None;
    }
    values.sort_unstable();
    let mid = values.len() / 2;
    Some(if values.len().is_multiple_of(2) {
        (values[mid - 1] + values[mid] + 1) / 2
    } else {
        values[mid]
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn data(name: &str) -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/testdata")
            .join(name)
    }

    #[test]
    fn median_of_odd_and_even_samples() {
        assert_eq!(median(vec![]), None);
        assert_eq!(median(vec![7]), Some(7));
        // Unsorted input, odd count: the middle value.
        assert_eq!(median(vec![9, 1, 5]), Some(5));
        // Even count: midpoint of the two middle values, rounded up.
        assert_eq!(median(vec![1, 2, 3, 4]), Some(3));
        assert_eq!(median(vec![10, 20]), Some(15));
        // A single outlier does not drag the median the way a mean would.
        assert_eq!(median(vec![1, 2, 3, 1000]), Some(3));
    }

    #[test]
    fn reports_the_median_fragment_length_of_a_paired_bam() {
        // test_i1.bam: TLENs 244 301 312 353 353 421 459 -> median 353.
        let p = data("test_i1.bam");
        assert_eq!(
            sample_library_layouts(&[(p.as_path(), "s1")]).unwrap(),
            vec![LibraryLayout::PairedEnd {
                median_fragment_length: Some(353)
            }]
        );
    }

    #[test]
    fn reports_a_layout_per_bam() {
        // test_i2.bam: TLENs 313 314 440 447 -> median (314+440)/2 = 377.
        let p1 = data("test_i1.bam");
        let p2 = data("test_i2.bam");
        let layouts =
            sample_library_layouts(&[(p1.as_path(), "s1"), (p2.as_path(), "s2")]).unwrap();
        assert_eq!(
            layouts,
            vec![
                LibraryLayout::PairedEnd {
                    median_fragment_length: Some(353)
                },
                LibraryLayout::PairedEnd {
                    median_fragment_length: Some(377)
                },
            ]
        );
    }

    #[test]
    fn estimate_averages_the_per_bam_medians() {
        // medians 353 and 377 -> mean 365.
        let p1 = data("test_i1.bam");
        let p2 = data("test_i2.bam");
        assert_eq!(
            estimate_fragment_length(&[(p1.as_path(), "s1"), (p2.as_path(), "s2")]).unwrap(),
            365
        );
        // A single BAM estimates to its own median.
        assert_eq!(
            estimate_fragment_length(&[(p1.as_path(), "s1")]).unwrap(),
            353
        );
    }

    #[test]
    fn estimating_without_any_bam_errors_rather_than_dividing_by_zero() {
        // The counting functions only reject an empty BAM list *after* resolving
        // the extension, so this path is reachable and must not panic.
        let err = estimate_fragment_length(&[])
            .err()
            .expect("estimating with no BAMs should fail");
        assert!(format!("{err:#}").contains("without a BAM file"));
    }

    #[test]
    fn resolve_passes_through_a_concrete_length_without_reading_bams() {
        // A real value is returned as-is; no path is touched, so a bogus path
        // is fine here.
        let bogus = Path::new("/does/not/exist.bam");
        let paths = [(bogus, "s1")];
        assert_eq!(resolve_extend_reads(Some(250), &paths).unwrap(), Some(250));
        assert_eq!(resolve_extend_reads(None, &paths).unwrap(), None);
    }
}
