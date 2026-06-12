//! Lightweight, opt-in phase timing for the counting pipeline.
//!
//! These are the only profiling *hooks* that must live inside the library
//! (phase boundaries can only be measured from the inside).  The actual
//! benchmark drivers live in the top-level `benchmarks/` directory.
//!
//! Timing is printed to stderr only when the `SINCEI_PROFILE` environment
//! variable is set, so it has no effect on normal runs.

use std::time::Instant;

/// Returns `true` when phase timing should be printed to stderr
/// (set the `SINCEI_PROFILE` environment variable to enable).
pub(crate) fn profile_enabled() -> bool {
    std::env::var_os("SINCEI_PROFILE").is_some()
}

/// Print `label` with the elapsed time since `since`, if profiling is enabled.
pub(crate) fn log_phase(enabled: bool, label: &str, since: Instant) {
    if enabled {
        eprintln!("[profile] {label:<22} {:>8.3} s", since.elapsed().as_secs_f64());
    }
}
