//! Read and filter BAM records into the structs the counting functions consume.
//!
//! [`bam_io`] opens indexed BAM files and drives per-chunk queries.
//! [`sc_record`] parses a raw record into an [`ScRecord`]. The barcode/UMI,
//! coordinates and optional QC fields needed to assign a read to a feature.
//! [`filters`] decides which records are kept (mapping quality, SAM flags,
//! QC thresholds, duplicates, blacklisted regions, and the 5′ motif check).
//! [`fragment_length`] samples a BAM's leading records to infer its library
//! layout and median fragment length, backing a valueless `--extendReads`.

pub mod bam_io;
pub mod filters;
pub mod fragment_length;
pub mod sc_record;

pub use filters::{
    BLACKLIST_MIN_OVERLAP_PERCENT, DupMethod, DuplicateFilter, MotifFilter, QcFilter,
    RawRecordFilter, is_blacklisted,
};
pub use sc_record::{AdjustRead, ScRecord};
