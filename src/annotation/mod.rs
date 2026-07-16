//! Parsers and abstractions for genomic region annotation files.
//!
//! [`region_index`] holds the region types ([`Feature`], [`Interval`]) and the
//! structures that make overlap queries fast ([`ChromIndex`], [`GenomeIndex`],
//! [`BinIndex`]).
//! [`parse_annotation`] reads BED / GTF / GFF3 files into those
//! structures.  The `counting` module consumes both to assign reads to features.

pub mod parse_annotation;
pub mod region_index;

pub use parse_annotation::{
    DEFAULT_EXON_TYPES, GENCODE_TRANSCRIPT_TYPE, parse_annotation_files, parse_blacklist_bed,
};
pub use region_index::{BinIndex, ChromIndex, Feature, GenomeIndex, Interval, build_bin_index};
