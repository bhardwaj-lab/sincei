//! Reading and filtering BAM records into the form the counting loops consume.
//!
//! [`bam_io`] opens indexed BAM files and drives per-chunk queries;
//! [`sc_record`] parses a raw record into an [`ScRecord`] — the barcode/UMI,
//! coordinates and optional QC fields needed to assign a read to a feature;
//! [`filters`] decides which records are kept (mapping quality, SAM flags,
//! QC thresholds, duplicates, and the 5′ motif check).

pub mod bam_io;
pub mod filters;
pub mod sc_record;

pub use filters::{DupMethod, DuplicateFilter, MotifFilter, QcFilter, RawRecordFilter};
pub use sc_record::{AdjustRead, ScRecord};
