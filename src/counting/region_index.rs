use ahash::AHashMap;

/// Per-feature metadata, ordered by feature (column) index.  Used to populate
/// the AnnData `var` DataFrame (`chrom`, `start`, `end`, `name`).
#[derive(Clone, Debug)]
pub struct VarMeta {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
    /// `chrom:start-end` for bins / unnamed BED features, or the feature name
    /// when the annotation provides one.
    pub name: String,
    /// Strand of the feature: `'+'`, `'-'`, or `'*'` when unstranded / unknown
    /// (always `'*'` for bins and BED features without a strand column).
    pub strand: char,
}

/// A half-open genomic interval `[start, end)` carrying an arbitrary value.
#[derive(Clone, Debug)]
pub struct Interval {
    pub start: usize,
    pub end: usize,
    /// Feature or bin index — used as the column key in the count matrix.
    pub val: usize,
    /// Feature name (e.g. gene ID, `"chrom:start-end"`). Empty for blacklist intervals.
    pub name: String,
}

/// Per-chromosome interval index backed by a sorted flat array.
///
/// Intervals are sorted by `start` at build time.  `max_end[i]` is the
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
    /// Build from an unsorted interval list.  O(N log N).
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

    /// Return all intervals that overlap the half-open query `[start, end)`.
    ///
    /// `partition_point` gives the first index where `iv.start >= end`; all
    /// indices below it are candidates.  Scanning backwards, we yield hits
    /// and stop early when `max_end[i] <= start`.
    pub fn find(&self, start: usize, end: usize) -> impl Iterator<Item = &Interval> {
        let p = self.intervals.partition_point(|iv| iv.start < end);
        FindIter {
            intervals: &self.intervals,
            max_end: &self.max_end,
            query_start: start,
            idx: p,
        }
    }

    /// Iterate all intervals in start-sorted order.
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

struct FindIter<'a> {
    intervals: &'a [Interval],
    max_end: &'a [usize],
    query_start: usize,
    idx: usize,
}

impl<'a> Iterator for FindIter<'a> {
    type Item = &'a Interval;

    fn next(&mut self) -> Option<&'a Interval> {
        loop {
            if self.idx == 0 {
                return None;
            }
            self.idx -= 1;
            // max_end[idx] is the max end of intervals[0..=idx].
            // If that max ≤ query_start, no remaining interval can overlap.
            if self.max_end[self.idx] <= self.query_start {
                return None;
            }
            let iv = &self.intervals[self.idx];
            if iv.end > self.query_start {
                return Some(iv);
            }
        }
    }
}

/// Per-chromosome interval index.
///
/// `AHashMap` preserves the insertion order of chromosome names, so iteration
/// order is deterministic (matches the order chromosomes were first seen in the
/// annotation or BAM header).
pub type RegionIndex = AHashMap<String, ChromIndex>;

/// Bin index
///
/// Precomputed geometry for a uniform-bin tiling of the genome.
///
/// Bin `i` on chromosome `c` covers the half-open interval
/// `[i * step_size,  min(i * step_size + bin_size, chrom_size))`.
///
/// `step_size == bin_size` → non-overlapping tiles (most common).
/// `step_size < bin_size`  → sliding windows.
///
/// The hot-loop lookup is O(1) arithmetic — no binary search needed.
pub struct BinIndex {
    pub bin_size: usize,
    pub step_size: usize,
    /// `chrom → (first_feature_idx, n_bins_on_chrom)`, in BAM-header order.
    pub chrom_bins: AHashMap<String, (usize, usize)>,
    pub n_bins: usize,
}

/// Tile each chromosome into bins and return a [`BinIndex`] together with
/// ordered feature names (`"chrom:start-end"`).
///
/// The last bin's `end` is clamped to the chromosome size.
/// Chromosomes with size 0 are skipped.
pub fn build_bin_index(
    chrom_sizes: &[(String, usize)],
    bin_size: usize,
    step_size: usize,
) -> (BinIndex, Vec<VarMeta>) {
    let mut var: Vec<VarMeta> = Vec::new();
    let mut chrom_bins: AHashMap<String, (usize, usize)> = AHashMap::new();

    for (chrom, chrom_size) in chrom_sizes {
        let chrom_size = *chrom_size;
        if chrom_size == 0 {
            continue;
        }
        let first_idx = var.len();
        let mut n_bins = 0usize;
        let mut bin_start = 0;
        while bin_start < chrom_size {
            let bin_end = (bin_start + bin_size).min(chrom_size);
            var.push(VarMeta {
                chrom: chrom.clone(),
                start: bin_start,
                end: bin_end,
                name: format!("{}:{}-{}", chrom, bin_start, bin_end),
                strand: '*',
            });
            n_bins += 1;
            bin_start += step_size;
        }
        if n_bins > 0 {
            chrom_bins.insert(chrom.clone(), (first_idx, n_bins));
        }
    }

    let n_bins = var.len();
    (
        BinIndex {
            bin_size,
            step_size,
            chrom_bins,
            n_bins,
        },
        var,
    )
}
