//! Make read count matrices from BAM files and annotations.
//!
//! [`count_bam_bins`] counts into a uniform tiling of the genome and
//! [`count_bam_features`] into the regions of annotation files; both write a
//! cell × region matrix as AnnData.
//! [`coverage`] instead pools cells into groups and writes pseudo-bulk signal tracks.
//!
//! The reads come from `crate::bam` and the regions from `crate::annotation`.

pub mod bin_matrix;
pub mod count_utils;
pub mod coverage;
pub mod feature_matrix;
pub mod params;

pub use bin_matrix::count_bam_bins;
pub use feature_matrix::count_bam_features;
pub use params::CountingParams;
