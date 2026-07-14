//! Build script.
//!
//! The library ships as a Python extension module: it links no CPython symbols
//! and relies on the interpreter that imports it to supply them.  A test binary
//! has no such host, so `cargo test --no-default-features` links `libpython`
//! directly — and needs an rpath to find it at run time (conda / pyenv
//! interpreters live outside the loader's default search paths).
//!
//! The rpath is emitted for every linkable target because cargo only accepts
//! the test-only form (`rustc-link-arg-tests`) on nightly.  It is inert in the
//! shipped extension module, which links no `libpython` to search for.
fn main() {
    println!("cargo::rerun-if-changed=build.rs");

    if let Some(lib_dir) = pyo3_build_config::get().lib_dir.as_deref() {
        println!("cargo::rustc-link-arg=-Wl,-rpath,{lib_dir}");
    }
}
