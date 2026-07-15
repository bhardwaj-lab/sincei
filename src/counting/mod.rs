pub mod bin_matrix;
pub mod count_utils;
pub mod coverage;
pub mod feature_matrix;
pub mod filters;
pub mod params;

pub use bin_matrix::count_bam_bins;
pub use feature_matrix::count_bam_features;
pub use filters::{DupMethod, DuplicateFilter, MotifFilter, QcFilter, RawRecordFilter};
pub use params::CountingParams;
