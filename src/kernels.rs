//! Pure, dependency-light k-mer kernels shared by the CLI and the streaming
//! sketchers ([`crate::sketchers`]).
//!
//! Nothing in here touches the filesystem, threads or logging; every function
//! works on caller-provided slices and buffers. These were previously defined
//! in `sketch.rs` (which pulls in needletail/clap/rayon) and are re-exported
//! from [`crate::sketch`] so existing import paths keep working.

use crate::types::{Kmer, BYTE_TO_SEQ};
use fxhash::{FxHashMap, FxHashSet};
#[cfg(feature = "approx-dedup")]
use fxhash::FxHasher;

pub type Marker = u32;

/// Splits one read into the two "double k-mers" used as a PCR-duplicate
/// signature by the single-end sketcher.
#[inline]
pub fn pair_kmer_single(s1: &[u8]) -> Option<([Marker; 2], [Marker; 2])> {
    let k = std::mem::size_of::<Marker>() * 4;
    if s1.len() < 4 * k + 2 {
        return None;
    } else {
        let mut kmer_f = 0;
        let mut kmer_g = 0;
        let mut kmer_r = 0;
        let mut kmer_t = 0;
        let halfway = s1.len() / 2;
        // len(s1)/2 + (k-1)* 2 + 2 < len(s1)
        for i in 0..k {
            let nuc_1 = BYTE_TO_SEQ[s1[2 * i] as usize] as Marker;
            let nuc_2 = BYTE_TO_SEQ[s1[2 * i + halfway] as usize] as Marker;
            let nuc_3 = BYTE_TO_SEQ[s1[1 + 2 * i] as usize] as Marker;
            let nuc_4 = BYTE_TO_SEQ[s1[1 + 2 * i + halfway] as usize] as Marker;

            kmer_f <<= 2;
            kmer_f |= nuc_1;

            kmer_r <<= 2;
            kmer_r |= nuc_2;

            kmer_g <<= 2;
            kmer_g |= nuc_3;

            kmer_t <<= 2;
            kmer_t |= nuc_4;
        }
        return Some(([kmer_f, kmer_r], [kmer_g, kmer_t]));
    }
}

/// Builds a read pair's double k-mer signature from both mates.
#[inline]
pub fn pair_kmer(s1: &[u8], s2: &[u8]) -> Option<([Marker; 2], [Marker; 2])> {
    let k = std::mem::size_of::<Marker>() * 4;
    if s1.len() < 2 * k + 1 || s2.len() < 2 * k + 1 {
        return None;
    } else {
        let mut kmer_f = 0;
        let mut kmer_g = 0;
        let mut kmer_r = 0;
        let mut kmer_t = 0;
        for i in 0..k {
            let nuc_1 = BYTE_TO_SEQ[s1[2 * i] as usize] as Marker;
            let nuc_2 = BYTE_TO_SEQ[s2[2 * i] as usize] as Marker;
            let nuc_3 = BYTE_TO_SEQ[s1[1 + 2 * i] as usize] as Marker;
            let nuc_4 = BYTE_TO_SEQ[s2[1 + 2 * i] as usize] as Marker;

            kmer_f <<= 2;
            kmer_f |= nuc_1;

            kmer_r <<= 2;
            kmer_r |= nuc_2;

            kmer_g <<= 2;
            kmer_g |= nuc_3;

            kmer_t <<= 2;
            kmer_t |= nuc_4;
        }
        return Some(([kmer_f, kmer_r], [kmer_g, kmer_t]));
    }
}

/// Exact PCR-duplicate removal: a k-mer whose `(hash, signature)` was already
/// seen in the same read (pair) is dropped, unless `no_dedup` is set or its
/// count reached `threshold` (the single-end `MAX_DEDUP_COUNT` escape).
pub fn dup_removal_lsh_full_exact(
    kmer_counts: &mut FxHashMap<Kmer, u32>,
    kmer_to_pair_set: &mut FxHashSet<(u64, [Marker; 2])>,
    km: &u64,
    kmer_pair: Option<([Marker; 2], [Marker; 2])>,
    num_dup_removed: &mut usize,
    no_dedup: bool,
    threshold: Option<u32>,
) {
    let c = kmer_counts.entry(*km).or_insert(0);
    let mut c_threshold = u32::MAX;
    if let Some(t) = threshold {
        c_threshold = t;
    }
    if !no_dedup && *c < c_threshold {
        if let Some(doublepairs) = kmer_pair {
            let mut ret = false;
            if kmer_to_pair_set.contains(&(*km, doublepairs.0)) {
                //Need this when using approximate data structures
                if *c > 0 {
                    ret = true;
                }
            } else {
                kmer_to_pair_set.insert((*km, doublepairs.0));
            }
            if kmer_to_pair_set.contains(&(*km, doublepairs.1)) {
                if *c > 0 {
                    ret = true;
                }
            } else {
                kmer_to_pair_set.insert((*km, doublepairs.1));
            }
            if ret {
                *num_dup_removed += 1;
                return;
            }
        }
    }
    *c += 1;
}

/// Approximate PCR-duplicate removal backed by a `scalable_cuckoo_filter`.
///
/// Only available with the `approx-dedup` feature (enabled by `cli`): the
/// filter crate depends on `rand`'s default features, which pull in
/// `getrandom`, so a getrandom-free (`--no-default-features`) build uses
/// [`dup_removal_lsh_full_exact`] only.
#[cfg(feature = "approx-dedup")]
pub fn dup_removal_lsh_full<R: rand::Rng>(
    kmer_counts: &mut FxHashMap<Kmer, u32>,
    kmer_to_pair_set: &mut scalable_cuckoo_filter::ScalableCuckooFilter<
        (u64, [Marker; 2]),
        FxHasher,
        R,
    >,
    km: &u64,
    kmer_pair: Option<([Marker; 2], [Marker; 2])>,
    num_dup_removed: &mut usize,
    no_dedup: bool,
) {
    let c = kmer_counts.entry(*km).or_insert(0);
    if !no_dedup {
        if let Some(doublepairs) = kmer_pair {
            let mut ret = false;
            if kmer_to_pair_set.contains(&(*km, doublepairs.0)) {
                //Need this when using approximate data structures
                if *c > 0 {
                    ret = true;
                }
            } else {
                kmer_to_pair_set.insert(&(*km, doublepairs.0));
            }
            if kmer_to_pair_set.contains(&(*km, doublepairs.1)) {
                if *c > 0 {
                    ret = true;
                }
            } else {
                kmer_to_pair_set.insert(&(*km, doublepairs.1));
            }
            if ret {
                *num_dup_removed += 1;
                return;
            }
        }
    }
    *c += 1;
}
