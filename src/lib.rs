mod filter_barcodes;
mod filter_stats;

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

    // Preprocessing functions
    m.add_function(wrap_pyfunction!(filter_barcodes::filter_barcodes, m)?)?;

    Ok(())
}
