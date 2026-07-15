//! Reading BAM records into the form the counting loops consume.
//!
//! [`bam_io`] opens indexed BAM files and drives per-chunk queries;
//! [`sc_record`] parses a raw record into an [`ScRecord`] — the barcode/UMI,
//! coordinates and optional QC fields needed to assign a read to a feature.

pub mod bam_io;
pub mod sc_record;

pub use sc_record::ScRecord;
