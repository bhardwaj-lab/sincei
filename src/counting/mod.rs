pub mod bin_matrix;
pub mod count_utils;
pub mod coverage;
pub mod feature_matrix;
pub mod params;

pub use bin_matrix::count_bam_bins;
pub use feature_matrix::count_bam_features;
pub use params::CountingParams;
