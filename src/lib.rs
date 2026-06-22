mod count_reads;
pub mod counting;
mod filter_barcodes;
mod filter_stats;
mod find_vcrs;
mod score_features;

use pyo3::prelude::*;
use pyo3::types::PyModule;

// Use jemalloc on non-MSVC targets
#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

// Version function
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[pymodule]
fn _sincei(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(version, m)?)?;

    // Preprocessing
    m.add_function(wrap_pyfunction!(filter_barcodes::filter_barcodes, m)?)?;
    m.add_function(wrap_pyfunction!(filter_stats::filter_stats, m)?)?;

    // Counting
    m.add_function(wrap_pyfunction!(count_reads::count_bins, m)?)?;
    m.add_function(wrap_pyfunction!(count_reads::count_features, m)?)?;

    // Export
    m.add_function(wrap_pyfunction!(counting::coverage::bulk_coverage, m)?)?;

    // Variable chromatin region detection
    m.add_function(wrap_pyfunction!(find_vcrs::find_vcrs, m)?)?;

    // Feature scoring
    m.add_function(wrap_pyfunction!(score_features::score_features, m)?)?;

    Ok(())
}
