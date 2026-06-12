mod filter_barcodes;
mod filter_stats;
mod count_reads;
mod utils;
pub mod counting;
pub mod export;
pub mod preprocessing;

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
    m.add_function(wrap_pyfunction!(export::coverage::bulk_coverage, m)?)?;

    Ok(())
}
