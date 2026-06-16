pub mod bam_io;
pub mod bin_matrix;
pub mod count_utils;
pub mod feature_matrix;
pub mod filters;
pub mod params;
pub mod parse_annotation;
pub mod region_index;
pub mod sc_record;

pub use bin_matrix::count_bam_bins;
pub use feature_matrix::count_bam_features;
pub use filters::{DupMethod, DuplicateFilter, MotifFilter, QcFilter, RawRecordFilter};
pub use params::CountingParams;
pub use parse_annotation::{
    parse_annotation_files, parse_bed_file, parse_gff_file, parse_gtf_file,
};
pub use region_index::{BinIndex, ChromIndex, Interval, RegionIndex, VarMeta, build_bin_index};
pub use sc_record::ScRecord;
