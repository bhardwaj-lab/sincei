//! Genomic regions, and the indexes that find them by position.
//!
//! A [`Feature`] is one column of the count matrix; an [`Interval`] is one
//! stretch to search, pointing back at its feature through `var_idx`. The two
//! are deliberately not 1:1, the exons of a transcript, or the pieces a
//! blacklist cuts a feature into, are several intervals sharing one feature.
//!
//! Reads are then assigned by one of two strategies: a [`GenomeIndex`] (one
//! [`ChromIndex`] per chromosome) searches arbitrary annotation regions in
//! O(log N + k), while a [`BinIndex`] needs no search at all, because a uniform
//! tiling puts a read's bin one division away.

use ahash::AHashMap;

/// A genomic feature defined by a chromosome, half-open interval, name, and strand.
#[derive(Clone, Debug)]
pub struct Feature {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    /// Feature name or `chrom:start-end` for bins/unnamed features.
    pub name: String,
    /// Strand of the feature: `'+'`, `'-'`, or `'*'` (unstranded/unknown).
    pub strand: char,
}

/// A half-open genomic interval `[start, end)` and what feature it belongs to
/// (noted in `var_idx`).
///
/// Several intervals may share a `var_idx`. E.g., the exons of one transcript,
/// or the pieces a feature is split by blacklist subtraction.
#[derive(Clone, Debug)]
pub struct Interval {
    pub start: usize,
    pub end: usize,
    /// The index in `Vec<Feature>` to which this interval belongs.
    pub var_idx: usize,
}

/// Searchable index over all the intervals on a ** single chromosome**, as a
/// sorted flat array.
///
/// Intervals are sorted by `start` at build time. `max_end[i]` is the
/// prefix-maximum of `end` over `intervals[0..=i]`, enabling O(log N + k)
/// overlap queries with early exit (k = number of hits).
///
/// Correct for nested and overlapping intervals: the early-exit condition
/// `max_end[j] <= query_start` implies every interval at index ≤ j has
/// `end ≤ query_start`, so none can overlap the query regardless of nesting.
pub struct ChromIndex {
    intervals: Vec<Interval>,
    max_end: Vec<usize>,
}

impl ChromIndex {
    /// Build from an unsorted interval list. O(N log N).
    pub fn build(mut intervals: Vec<Interval>) -> Self {
        intervals.sort_unstable_by_key(|iv| iv.start);
        let max_end = intervals
            .iter()
            .scan(0usize, |m, iv| {
                *m = (*m).max(iv.end);
                Some(*m)
            })
            .collect();
        Self { intervals, max_end }
    }

    /// Return all intervals that overlap the half-open query `[start, end)`,
    /// in descending order of `start`.
    ///
    /// Since the intervals are sorted by `start`, `partition_point` finds the
    /// first one that begins at or after the query ends. It, and everything
    /// after it, starts too late to overlap, so only the intervals below it are
    /// candidates. Those are then walked back-to-front, which lets the scan
    /// stop as soon as it can:
    ///
    /// * `take_while` **ends** the scan at the first interval whose `max_end`
    ///   (the greatest `end` among it and everything before it) is at or before
    ///   the query's start. No earlier interval can reach the query either, so
    ///   there is nothing left to find. This keeps a lookup O(log N + k) rather
    ///   than a walk to index 0.
    /// * `filter` **skips** an individual interval that ends at or before the
    ///   query's start. Its neighbours may still overlap (they can be nested
    ///   inside a longer interval), so the scan carries on.
    ///
    /// Every interval not filtered out starts before the query ends and ends
    /// after it begins, a half-open overlap.
    pub fn find(&self, start: usize, end: usize) -> impl Iterator<Item = &Interval> {
        let candidates = self.intervals.partition_point(|iv| iv.start < end);

        self.intervals[..candidates]
            .iter()
            .zip(&self.max_end[..candidates])
            .rev()
            .take_while(move |(_, max_end)| **max_end > start)
            .map(|(interval, _)| interval)
            .filter(move |interval| interval.end > start)
    }

    /// Iterate over all intervals in start-sorted order.
    pub fn iter(&self) -> impl Iterator<Item = &Interval> {
        self.intervals.iter()
    }

    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }
}

/// Searchable index over the intervals of a **whole genome**, one
/// [`ChromIndex`] per chromosome, keyed by chromosome name.
///
/// Callers resolve `chrom -> ChromIndex` once per work chunk (each covers a
/// single chromosome) and then query that index per read.
///
/// `AHashMap` preserves the iteration order of chromosomes.
pub type GenomeIndex = AHashMap<String, ChromIndex>;

/// Uniform-bin tiling of a **whole genome**, covering every chromosome.
///
/// The bin equivalent of [`GenomeIndex`], for when the regions are a tiling
/// rather than arbitrary annotation features. Nothing is searched and no
/// [`Interval`] is stored: a bin's position is arithmetic, so the bin holding a
/// read is found in O(1) rather than O(log N + k).
///
/// Bin `i` on chromosome `c` covers the half-open interval
/// `[i * step_size,  min(i * step_size + bin_size, chrom_size))`.
///
/// `step_size == bin_size` -> non-overlapping tiles (most common).
/// `step_size < bin_size`  -> sliding windows.
pub struct BinIndex {
    pub bin_size: usize,
    pub step_size: usize,
    /// `chrom -> (first_feature_idx, n_bins_on_chrom)`, in BAM-header order.
    /// The first element turns a bin's index within its chromosome into its
    /// column in the count matrix.
    pub chrom_bins: AHashMap<String, (usize, usize)>,
    /// Total bins across all chromosomes, the number of matrix columns.
    pub n_bins: usize,
    /// `chrom -> first tiled base`. Zero for every chromosome unless a region
    /// restricted the tiling to a window, in which case bin `b` of that
    /// chromosome starts at `window_start + b * step_size` rather than at
    /// `b * step_size`.
    pub chrom_window_start: AHashMap<String, usize>,
}

impl BinIndex {
    /// First tiled base on `chrom`. Zero unless a region restricted the tiling,
    /// so callers can add it unconditionally.
    pub fn window_start(&self, chrom: &str) -> usize {
        self.chrom_window_start.get(chrom).copied().unwrap_or(0)
    }
}

/// Tile each chromosome into bins, without naming them.
///
/// A bin's position is arithmetic, so the tiling needs nothing per bin beyond
/// the count: chromosome `c` holds `ceil(chrom_size / step_size)` bins, since a
/// bin starts every `step_size` bp for as long as that start is inside the
/// chromosome.
///
/// Use this wherever the bins are only ever counted into and written out by
/// coordinate, coverage tracks, for one. [`build_bin_index`] additionally
/// materializes a [`Feature`] per bin, which only the count matrices need for
/// their `var` table, and which at a 100 bp bin size means tens of millions of
/// allocated names.
///
/// Chromosomes with size 0 are skipped.
pub fn build_bigwig_index(
    chrom_sizes: &[(String, usize)],
    bin_size: usize,
    step_size: usize,
) -> BinIndex {
    build_bigwig_index_in_window(&whole_chromosomes(chrom_sizes), bin_size, step_size)
}

/// Turn chromosome sizes into windows spanning each whole chromosome, so the
/// unrestricted tiling is the windowed tiling with nothing cut away.
fn whole_chromosomes(chrom_sizes: &[(String, usize)]) -> Vec<(String, usize, usize)> {
    chrom_sizes
        .iter()
        .map(|(chrom, size)| (chrom.clone(), 0, *size))
        .collect()
}

/// Tile a window of each chromosome into bins, without naming them.
///
/// `windows` holds `(chrom, window_start, window_end)` in 0-based half-open
/// coordinates, already clamped to the chromosome. Bins start at `window_start`
/// and step by `step_size` for as long as that start is inside the window, so
/// the window rather than the chromosome is the universe being tiled: a run
/// restricted to a region gets bins for that region only.
///
/// [`build_bigwig_index`] is this function with every window covering a whole
/// chromosome. Empty windows (`end <= start`) are skipped, as size-0
/// chromosomes are there.
pub fn build_bigwig_index_in_window(
    windows: &[(String, usize, usize)],
    bin_size: usize,
    step_size: usize,
) -> BinIndex {
    let mut chrom_bins: AHashMap<String, (usize, usize)> = AHashMap::new();
    let mut chrom_window_start: AHashMap<String, usize> = AHashMap::new();
    let mut n_bins = 0usize;

    for (chrom, window_start, window_end) in windows {
        let Some(span) = window_end.checked_sub(*window_start).filter(|s| *s > 0) else {
            continue;
        };
        let chrom_n_bins = span.div_ceil(step_size);
        chrom_bins.insert(chrom.clone(), (n_bins, chrom_n_bins));
        chrom_window_start.insert(chrom.clone(), *window_start);
        n_bins += chrom_n_bins;
    }

    BinIndex {
        bin_size,
        step_size,
        chrom_bins,
        n_bins,
        chrom_window_start,
    }
}

/// Tile each chromosome into bins and return a [`BinIndex`] together with one
/// [`Feature`] per bin, named `"chrom:start-end"`.
///
/// The last bin's `end` is clamped to the chromosome size.
/// Chromosomes with size 0 are skipped.
///
/// The tiling itself comes from [`build_bigwig_index`], so the two cannot
/// disagree on how many bins a chromosome has.
pub fn build_bin_index(
    chrom_sizes: &[(String, usize)],
    bin_size: usize,
    step_size: usize,
) -> (BinIndex, Vec<Feature>) {
    build_bin_index_in_window(&whole_chromosomes(chrom_sizes), bin_size, step_size)
}

/// Tile a window of each chromosome into bins and return a [`BinIndex`] with one
/// [`Feature`] per bin, named `"chrom:start-end"`.
///
/// The windowed counterpart of [`build_bin_index`]: `windows` is interpreted as
/// in [`build_bigwig_index_in_window`], and the last bin's `end` is clamped to
/// the window end rather than to the chromosome size.
pub fn build_bin_index_in_window(
    windows: &[(String, usize, usize)],
    bin_size: usize,
    step_size: usize,
) -> (BinIndex, Vec<Feature>) {
    let index = build_bigwig_index_in_window(windows, bin_size, step_size);
    let mut var: Vec<Feature> = Vec::with_capacity(index.n_bins);

    for (chrom, window_start, window_end) in windows {
        if window_end <= window_start {
            continue;
        }
        let mut bin_start = *window_start;
        while bin_start < *window_end {
            let bin_end = (bin_start + bin_size).min(*window_end);
            var.push(Feature {
                chrom: chrom.clone(),
                start: bin_start,
                end: bin_end,
                name: format!("{}:{}-{}", chrom, bin_start, bin_end),
                strand: '*',
            });
            bin_start += step_size;
        }
    }

    (index, var)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(start: usize, end: usize, var_idx: usize) -> Interval {
        Interval {
            start,
            end,
            var_idx,
        }
    }

    /// Feature indices overlapping `[start, end)`, sorted so assertions do not
    /// depend on the (descending-start) iteration order of `find`.
    fn hits(idx: &ChromIndex, start: usize, end: usize) -> Vec<usize> {
        let mut v: Vec<usize> = idx.find(start, end).map(|iv| iv.var_idx).collect();
        v.sort_unstable();
        v
    }

    #[test]
    fn find_returns_overlapping_intervals_only() {
        // Deliberately unsorted input: `build` must sort by start.
        let idx = ChromIndex::build(vec![iv(20, 30, 2), iv(0, 10, 0), iv(5, 15, 1)]);

        assert_eq!(hits(&idx, 12, 13), vec![1]);
        assert_eq!(hits(&idx, 0, 5), vec![0]);
        assert_eq!(hits(&idx, 0, 25), vec![0, 1, 2]);
        assert_eq!(hits(&idx, 30, 40), Vec::<usize>::new());
    }

    #[test]
    fn find_treats_intervals_as_half_open() {
        let idx = ChromIndex::build(vec![iv(10, 20, 0)]);

        // Touching at the ends is not an overlap.
        assert_eq!(hits(&idx, 0, 10), Vec::<usize>::new());
        assert_eq!(hits(&idx, 20, 30), Vec::<usize>::new());
        // One shared base is.
        assert_eq!(hits(&idx, 9, 11), vec![0]);
        assert_eq!(hits(&idx, 19, 21), vec![0]);
    }

    #[test]
    fn find_handles_nested_intervals() {
        // The max_end prefix-maximum exists so a long interval early in the
        // sorted order is not missed once a short later one ends before the query.
        let idx = ChromIndex::build(vec![iv(0, 1000, 0), iv(10, 20, 1), iv(500, 510, 2)]);

        assert_eq!(hits(&idx, 900, 950), vec![0]);
        assert_eq!(hits(&idx, 505, 506), vec![0, 2]);
        assert_eq!(hits(&idx, 15, 16), vec![0, 1]);
    }

    #[test]
    fn find_on_empty_index_yields_nothing() {
        let idx = ChromIndex::build(vec![]);
        assert!(idx.is_empty());
        assert_eq!(idx.len(), 0);
        assert_eq!(hits(&idx, 0, 100), Vec::<usize>::new());
    }

    #[test]
    fn find_returns_every_duplicate_at_the_same_start() {
        let idx = ChromIndex::build(vec![iv(100, 200, 0), iv(100, 200, 1), iv(100, 150, 2)]);
        assert_eq!(hits(&idx, 120, 130), vec![0, 1, 2]);
        assert_eq!(hits(&idx, 160, 170), vec![0, 1]);
    }

    #[test]
    fn the_name_free_index_tiles_exactly_like_the_named_one() {
        // `build_bigwig_index` derives the bin count arithmetically while
        // `build_bin_index` walks the tiling, so they have to be checked against
        // each other: a disagreement would silently shift every bin's column.
        let chrom_sizes = vec![
            ("chr1".to_string(), 250),
            ("chrEmpty".to_string(), 0),
            ("chr2".to_string(), 100),
            ("chrShort".to_string(), 30),
        ];

        for (bin_size, step_size) in [(100, 100), (50, 25), (100, 50), (1, 1), (1000, 1000)] {
            let bare = build_bigwig_index(&chrom_sizes, bin_size, step_size);
            let (named, var) = build_bin_index(&chrom_sizes, bin_size, step_size);

            assert_eq!(bare.n_bins, named.n_bins, "bin/step {bin_size}/{step_size}");
            assert_eq!(bare.n_bins, var.len(), "bin/step {bin_size}/{step_size}");
            assert_eq!(
                bare.chrom_bins, named.chrom_bins,
                "bin/step {bin_size}/{step_size}"
            );
        }
    }

    #[test]
    fn bins_tile_the_chromosome_and_clamp_the_last_bin() {
        let (index, var) = build_bin_index(&[("chr1".to_string(), 250)], 100, 100);

        assert_eq!(index.n_bins, 3);
        assert_eq!(index.chrom_bins.get("chr1"), Some(&(0, 3)));
        let spans: Vec<(usize, usize)> = var.iter().map(|v| (v.start, v.end)).collect();
        assert_eq!(spans, vec![(0, 100), (100, 200), (200, 250)]);
        assert_eq!(var[2].name, "chr1:200-250");
        assert!(var.iter().all(|v| v.strand == '*'));
    }

    #[test]
    fn step_size_below_bin_size_produces_sliding_windows() {
        let (index, var) = build_bin_index(&[("chr1".to_string(), 100)], 50, 25);

        assert_eq!(index.n_bins, 4);
        let spans: Vec<(usize, usize)> = var.iter().map(|v| (v.start, v.end)).collect();
        assert_eq!(spans, vec![(0, 50), (25, 75), (50, 100), (75, 100)]);
    }

    #[test]
    fn chromosome_offsets_are_cumulative_and_empty_chroms_are_skipped() {
        let (index, var) = build_bin_index(
            &[
                ("chr1".to_string(), 150),
                ("chrEmpty".to_string(), 0),
                ("chr2".to_string(), 100),
            ],
            100,
            100,
        );

        assert_eq!(index.chrom_bins.get("chr1"), Some(&(0, 2)));
        assert_eq!(index.chrom_bins.get("chr2"), Some(&(2, 1)));
        assert!(!index.chrom_bins.contains_key("chrEmpty"));
        assert_eq!(index.n_bins, 3);
        assert_eq!(var[2].chrom, "chr2");
    }

    #[test]
    fn chromosome_shorter_than_one_bin_gets_a_single_clamped_bin() {
        let (index, var) = build_bin_index(&[("chr1".to_string(), 30)], 100, 100);

        assert_eq!(index.n_bins, 1);
        assert_eq!((var[0].start, var[0].end), (0, 30));
        assert_eq!(var[0].name, "chr1:0-30");
    }
}
