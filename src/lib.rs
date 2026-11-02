//! sylph: ultrafast genome ANI queries and taxonomic profiling.
//!
//! The library core (this crate compiled with `--no-default-features`) is
//! dependency-light and builds for `wasm32-unknown-unknown`: it provides the
//! sketch types ([`types`]), the seeding kernels ([`seeding`]) and push-style
//! sketchers ([`sketchers`]). The `cli` feature (on by default) adds the
//! command-line front end and all filesystem/threading/profiling code.
//!
//! The pure read-dedup/pair kernels are private; they are reachable through
//! their historical [`sketch`] paths (`sketch::pair_kmer`, ...).

pub mod constants;
mod kernels;
pub mod seeding;
pub mod sketchers;
pub mod sketch;
pub mod types;

#[cfg(target_arch = "x86_64")]
pub mod avx2_seeding;

// Command-line front end, I/O, parallelism and profiling. None of these are
// needed to reuse the sketching kernels.
#[cfg(feature = "cli")]
pub mod cmdline;
#[cfg(feature = "cli")]
pub mod contain;
#[cfg(feature = "cli")]
pub mod inference;
#[cfg(feature = "cli")]
pub mod inspect;
#[cfg(feature = "cli")]
pub mod parallel_sketch;
#[cfg(feature = "cli")]
pub mod twostage_db;
