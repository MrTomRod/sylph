use crate::types::*;

#[inline]
pub fn mm_hash64(kmer: u64) -> u64 {
    //TODO this is bugged. Fix after release
    let mut key = kmer;
    key = !key.wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    key
}

#[inline]
pub fn rev_hash_64(hashed_key: u64) -> u64 {
    let mut key = hashed_key;

    // Invert h_key = h_key.wrapping_add(h_key << 31)
    let mut tmp: u64 = key.wrapping_sub(key << 31);
    key = key.wrapping_sub(tmp << 31);

    // Invert h_key = h_key ^ h_key >> 28;
    tmp = key ^ key >> 28;
    key = key ^ tmp >> 28;

    // Invert h_key = h_key.wrapping_add(h_key << 2).wrapping_add(h_key << 4)
    key = key.wrapping_mul(14933078535860113213u64);

    // Invert h_key = h_key ^ h_key >> 14;
    tmp = key ^ key >> 14;
    tmp = key ^ tmp >> 14;
    tmp = key ^ tmp >> 14;
    key = key ^ tmp >> 14;

    // Invert h_key = h_key.wrapping_add(h_key << 3).wrapping_add(h_key << 8)
    key = key.wrapping_mul(15244667743933553977u64);

    // Invert h_key = h_key ^ h_key >> 24
    tmp = key ^ key >> 24;
    key = key ^ tmp >> 24;

    // Invert h_key = (!h_key).wrapping_add(h_key << 21)
    tmp = !key;
    tmp = !(key.wrapping_sub(tmp << 21));
    tmp = !(key.wrapping_sub(tmp << 21));
    key = !(key.wrapping_sub(tmp << 21));

    key
}

pub fn decode(byte: u64) -> u8 {
    if byte == 0 {
        return b'A';
    } else if byte == 1 {
        return b'C';
    } else if byte == 2 {
        return b'G';
    } else if byte == 3 {
        return b'T';
    } else {
        panic!("decoding failed")
    }
}
pub fn print_string(kmer: u64, k: usize) {
    let mut bytes = vec![];
    let mask = 3;
    for i in 0..k {
        let val = kmer >> 2 * i;
        let val = val & mask;
        bytes.push(decode(val));
    }
    dbg!(std::str::from_utf8(&bytes.into_iter().rev().collect::<Vec<u8>>()).unwrap());
}
#[inline]
fn _position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

pub fn fmh_seeds(
    string: &[u8],
    kmer_vec: &mut Vec<u64>,
    c: usize,
    k: usize
) {
    type MarkerBits = u64;
    if string.len() < k {
        return;
    }

    let marker_k = k;
    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    //    let threshold = i64::MIN + (u64::MAX / (c as u64)) as i64;
    //    let threshold_marker = i64::MIN + (u64::MAX / sketch_params.marker_c as u64) as i64;

    let threshold_marker = u64::MAX / (c as u64);
    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        //        let nuc_f = KmerEnc::encode(string[i]
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        //        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
    }
    for i in marker_k-1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker &= marker_rev_mask;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
        //        rolling_kmer_r &= max_mask;
        //        KmerEnc::print_string(rolling_kmer_f, k);
        //        KmerEnc::print_string(rolling_kmer_r, k);
        //

        let canonical_marker = rolling_kmer_f_marker < rolling_kmer_r_marker;
        let canonical_kmer_marker = if canonical_marker {
            rolling_kmer_f_marker
        } else {
            rolling_kmer_r_marker
        };
        let hash_marker = mm_hash64(canonical_kmer_marker);

        if hash_marker < threshold_marker {
            kmer_vec.push(hash_marker as u64);
        }
    }
}

pub fn fmh_seeds_positions(
    string: &[u8],
    kmer_vec: &mut Vec<(usize,usize,u64)>,
    c: usize,
    k: usize,
    contig_number: usize
) {
    type MarkerBits = u64;
    if string.len() < k {
        return;
    }

    let marker_k = k;
    let mut rolling_kmer_f_marker: MarkerBits = 0;
    let mut rolling_kmer_r_marker: MarkerBits = 0;

    let marker_reverse_shift_dist = 2 * (marker_k - 1);
    let marker_mask = MarkerBits::MAX >> (std::mem::size_of::<MarkerBits>() * 8 - 2 * marker_k);
    let marker_rev_mask = !(3 << (2 * marker_k - 2));
    let len = string.len();
    //    let threshold = i64::MIN + (u64::MAX / (c as u64)) as i64;
    //    let threshold_marker = i64::MIN + (u64::MAX / sketch_params.marker_c as u64) as i64;

    let threshold_marker = u64::MAX / (c as u64);
    for i in 0..marker_k - 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize] as u64;
        //        let nuc_f = KmerEnc::encode(string[i]
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        //        rolling_kmer_r = KmerEnc::rc(rolling_kmer_f, k);
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
    }
    for i in marker_k-1..len {
        let nuc_byte = string[i] as usize;
        let nuc_f = BYTE_TO_SEQ[nuc_byte] as u64;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_f_marker <<= 2;
        rolling_kmer_f_marker |= nuc_f;
        rolling_kmer_f_marker &= marker_mask;
        rolling_kmer_r_marker >>= 2;
        rolling_kmer_r_marker &= marker_rev_mask;
        rolling_kmer_r_marker |= nuc_r << marker_reverse_shift_dist;
        //        rolling_kmer_r &= max_mask;
        //        KmerEnc::print_string(rolling_kmer_f, k);
        //        KmerEnc::print_string(rolling_kmer_r, k);
        //

        let canonical_marker = rolling_kmer_f_marker < rolling_kmer_r_marker;
        let canonical_kmer_marker = if canonical_marker {
            rolling_kmer_f_marker
        } else {
            rolling_kmer_r_marker
        };
        let hash_marker = mm_hash64(canonical_kmer_marker);

        if hash_marker < threshold_marker {
            kmer_vec.push((contig_number, i, hash_marker as u64));
        }
    }
}

/// Which seeding kernel `extract_markers` / `extract_markers_positions`
/// (thin wrappers over `extract_markers_with`, re-exported from
/// [`crate::sketch`]) should use.
///
/// Both kernels require `c > 0` (the k-mer downsampling threshold is
/// `u64::MAX / c`) and only support `k` up to 32 on the scalar path. In
/// practice every sylph build so far uses `k = 21` or `k = 31`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum Seeding {
    /// Whatever this build does natively: sylph's AVX2 kernel on x86_64 CPUs
    /// that support AVX2, the scalar kernel everywhere else. This is the
    /// behaviour the `sylph` CLI has always had.
    ///
    /// On an AVX2-capable host, `k` must be 21 or 31: the AVX2 kernel derives
    /// its shift constants from `2 * (k - 1)` and panics for any other `k`
    /// once a sequence is long enough to enter the kernel body (at least
    /// `k + 1` bases for marker extraction, at least `2 * k` for marker
    /// positions). Shorter sequences return no markers before that check, so
    /// they do not panic. On hosts without AVX2 the scalar kernel is used,
    /// which works for `k <= 32` but yields different k-mers than the AVX2
    /// kernel.
    #[default]
    Native,
    /// Portable scalar emulation of the x86_64 AVX2 kernel.
    ///
    /// It produces exactly the k-mers the AVX2 kernel produces, in exactly the
    /// same order, on every architecture (so browser/wasm and aarch64 builds
    /// can reproduce x86_64 sketches). See `fmh_seeds_avx2_compat` and
    /// `fmh_seeds_avx2_compat_positions`.
    ///
    /// Only `k = 21` and `k = 31` are guaranteed to match the AVX2 kernel
    /// (the AVX2 kernel accepts no other value; the scalar kernel itself
    /// supports `k <= 32`).
    Avx2Compat,
}

/// As [`extract_markers`], with an explicit seeding strategy.
pub(crate) fn extract_markers_with(
    seeding: Seeding,
    string: &[u8],
    kmer_vec: &mut Vec<u64>,
    c: usize,
    k: usize,
) {
    match seeding {
        Seeding::Avx2Compat => fmh_seeds_avx2_compat(string, kmer_vec, c, k),
        Seeding::Native => {
            #[cfg(target_arch = "x86_64")]
            {
                if is_x86_feature_detected!("avx2") {
                    use crate::avx2_seeding::*;
                    unsafe {
                        extract_markers_avx2(string, kmer_vec, c, k);
                    }
                    return;
                }
            }
            fmh_seeds(string, kmer_vec, c, k);
        }
    }
}

/// Appends the hashed canonical k-mers of `string` to `kmer_vec` using this
/// build's native seeding kernel (AVX2 on x86_64 when available, scalar
/// otherwise).
pub(crate) fn extract_markers(string: &[u8], kmer_vec: &mut Vec<u64>, c: usize, k: usize) {
    extract_markers_with(Seeding::Native, string, kmer_vec, c, k);
}

/// As [`extract_markers_positions`], with an explicit seeding strategy.
pub(crate) fn extract_markers_positions_with(
    seeding: Seeding,
    string: &[u8],
    kmer_vec: &mut Vec<(usize, usize, u64)>,
    c: usize,
    k: usize,
    contig_number: usize,
) {
    match seeding {
        Seeding::Avx2Compat => {
            fmh_seeds_avx2_compat_positions(string, kmer_vec, c, k, contig_number)
        }
        Seeding::Native => {
            #[cfg(target_arch = "x86_64")]
            {
                if is_x86_feature_detected!("avx2") {
                    use crate::avx2_seeding::*;
                    unsafe {
                        extract_markers_avx2_positions(string, kmer_vec, c, k, contig_number);
                    }
                    return;
                }
            }
            fmh_seeds_positions(string, kmer_vec, c, k, contig_number);
        }
    }
}

/// Appends `(contig_number, k-mer end position, hash)` markers using this
/// build's native seeding kernel.
pub(crate) fn extract_markers_positions(
    string: &[u8],
    kmer_vec: &mut Vec<(usize, usize, u64)>,
    c: usize,
    k: usize,
    contig_number: usize,
) {
    extract_markers_positions_with(Seeding::Native, string, kmer_vec, c, k, contig_number);
}

/// Scalar emulation of sylph's x86_64 `avx2_seeding::extract_markers_avx2`.
///
/// The AVX2 kernel splits the `n - k + 1` k-mer start positions of a read into
/// four equal segments of `(n - k + 1) / 4` positions and processes them in
/// lock-step, one k-mer of each segment per step, in segment order. As a
/// consequence it
///
/// * never sees the last `(n - k + 1) % 4` k-mers of a read,
/// * yields nothing for a read of exactly `k` bases (or shorter), and
/// * appends k-mers in the interleaved segment order above (which matters for
///   sylph's order-dependent read deduplication).
///
/// This function reproduces that output exactly and portably for `k = 21` and
/// `k = 31`, the only two values the AVX2 kernel accepts: it picks its shift
/// constants from `2 * (k - 1)` and `panic!()`s for any other `k` once a
/// sequence is long enough to reach the kernel body (at least `k + 1` bases
/// here; at least `2 * k` for [`fmh_seeds_avx2_compat_positions`]). Shorter
/// inputs return early, before the AVX2 kernel's `k` check, so they are
/// silently empty rather than a panic. For other `k` values this function
/// still runs, but its output is not a compatibility guarantee for any real
/// sylph build.
fn fmh_seeds_avx2_compat(string: &[u8], kmer_vec: &mut Vec<u64>, c: usize, k: usize) {
    if string.len() < k + 1 {
        return;
    }
    let seg = (string.len() - k + 1) / 4;
    if seg == 0 {
        return;
    }
    let reverse_shift = 2 * (k - 1);
    let mask = u64::MAX >> (64 - 2 * k);
    let rev_mask = !(3u64 << (2 * k - 2));
    let threshold = u64::MAX / (c as u64);

    let mut f = [0u64; 4];
    let mut r = [0u64; 4];
    for i in 0..k - 1 {
        for s in 0..4 {
            let nuc_f = BYTE_TO_SEQ[string[s * seg + i] as usize] as u64;
            f[s] = (f[s] << 2) | nuc_f;
            r[s] = (r[s] >> 2) | ((3 - nuc_f) << reverse_shift);
        }
    }
    for i in k - 1..seg + k - 1 {
        for s in 0..4 {
            let nuc_f = BYTE_TO_SEQ[string[s * seg + i] as usize] as u64;
            f[s] = ((f[s] << 2) | nuc_f) & mask;
            r[s] = ((r[s] >> 2) & rev_mask) | ((3 - nuc_f) << reverse_shift);
            let hash = mm_hash64(f[s].min(r[s]));
            if hash < threshold {
                kmer_vec.push(hash);
            }
        }
    }
}

/// Scalar emulation of sylph's x86_64
/// `avx2_seeding::extract_markers_avx2_positions`.
///
/// Same segment/tail behaviour as [`fmh_seeds_avx2_compat`], plus the AVX2
/// positions kernel's additional rule: a read shorter than `2 * k` yields
/// nothing at all. Positions are the k-mer end offsets (`segment start + i`),
/// as in the AVX2 kernel.
fn fmh_seeds_avx2_compat_positions(
    string: &[u8],
    kmer_vec: &mut Vec<(usize, usize, u64)>,
    c: usize,
    k: usize,
    contig_number: usize,
) {
    if string.len() < k {
        return;
    }
    let seg = (string.len() - k + 1) / 4;
    if string.len() < 2 * k || seg == 0 {
        return;
    }
    let reverse_shift = 2 * (k - 1);
    let mask = u64::MAX >> (64 - 2 * k);
    let rev_mask = !(3u64 << (2 * k - 2));
    let threshold = u64::MAX / (c as u64);

    let mut f = [0u64; 4];
    let mut r = [0u64; 4];
    for i in 0..k - 1 {
        for s in 0..4 {
            let nuc_f = BYTE_TO_SEQ[string[s * seg + i] as usize] as u64;
            f[s] = (f[s] << 2) | nuc_f;
            r[s] = (r[s] >> 2) | ((3 - nuc_f) << reverse_shift);
        }
    }
    for i in k - 1..seg + k - 1 {
        for s in 0..4 {
            let nuc_f = BYTE_TO_SEQ[string[s * seg + i] as usize] as u64;
            f[s] = ((f[s] << 2) | nuc_f) & mask;
            r[s] = ((r[s] >> 2) & rev_mask) | ((3 - nuc_f) << reverse_shift);
            let hash = mm_hash64(f[s].min(r[s]));
            if hash < threshold {
                kmer_vec.push((contig_number, s * seg + i, hash));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random read (xorshift).
    fn random_read(seed: u64, len: usize) -> Vec<u8> {
        let mut state = seed | 1;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        (0..len).map(|_| b"ACGT"[(next() % 4) as usize]).collect()
    }

    /// Includes ambiguous `N`s, lower-case bases and IUPAC codes. `BYTE_TO_SEQ`
    /// maps `N` and the other non-ACGT bytes to `A`, while lower-case
    /// `a`/`c`/`g`/`t` map to the same codes as their upper-case forms
    /// (so a lower-case `c`/`g`/`t` is *not* turned into `A`).
    fn noisy_read(seed: u64, len: usize) -> Vec<u8> {
        let mut read = random_read(seed, len);
        for (i, b) in read.iter_mut().enumerate() {
            if i % 37 == 0 {
                *b = b'N';
            } else if i % 53 == 0 {
                *b = b.to_ascii_lowercase();
            } else if i % 89 == 0 {
                *b = b'R';
            }
        }
        read
    }

    #[test]
    #[cfg(target_arch = "x86_64")]
    fn avx2_compat_matches_avx2_seeding() {
        if !is_x86_feature_detected!("avx2") {
            eprintln!("AVX2 not available; skipping AVX2 parity test");
            return;
        }
        for k in [21usize, 31] {
            let lengths = [
                k - 1,
                k,
                k + 1,
                k + 2,
                2 * k - 1,
                2 * k,
                2 * k + 1,
                4 * k - 1,
                4 * k,
                100,
                103,
                105,
                157,
                1000,
                4099,
            ];
            for c in [1usize, 5, 100, 200] {
                for (n, &len) in lengths.iter().enumerate() {
                    for (kind, read) in [
                        ("acgt", random_read(len as u64 * 31 + k as u64 + n as u64, len)),
                        ("noisy", noisy_read(len as u64 * 17 + k as u64 + n as u64, len)),
                    ] {
                        let mut native = Vec::new();
                        let mut compat = Vec::new();
                        unsafe {
                            crate::avx2_seeding::extract_markers_avx2(
                                &read,
                                &mut native,
                                c,
                                k,
                            );
                        }
                        fmh_seeds_avx2_compat(&read, &mut compat, c, k);
                        assert_eq!(native, compat, "markers k={k} len={len} c={c} kind={kind}");

                        let mut native = Vec::new();
                        let mut compat = Vec::new();
                        unsafe {
                            crate::avx2_seeding::extract_markers_avx2_positions(
                                &read,
                                &mut native,
                                c,
                                k,
                                /*contig_number=*/ 3,
                            );
                        }
                        fmh_seeds_avx2_compat_positions(&read, &mut compat, c, k, 3);
                        assert_eq!(
                            native, compat,
                            "positions k={k} len={len} c={c} kind={kind}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn avx2_compat_drops_the_tail_and_short_reads() {
        for k in [21usize, 31] {
            let read = random_read(k as u64, 4 * k + 3);
            // shorter than k: nothing
            let mut out = Vec::new();
            fmh_seeds_avx2_compat(&read[..k - 1], &mut out, 1, k);
            assert!(out.is_empty());
            // exactly k: nothing
            let mut out = Vec::new();
            fmh_seeds_avx2_compat(&read[..k], &mut out, 1, k);
            assert!(out.is_empty());
            // positions: nothing below 2k
            let mut out = Vec::new();
            fmh_seeds_avx2_compat_positions(&read[..2 * k - 1], &mut out, 1, k, 0);
            assert!(out.is_empty());

            // n - k + 1 not divisible by 4: exactly 4 * ((n - k + 1) / 4)
            // k-mers of the scalar path remain (before c-subsampling).
            let mut scalar = Vec::new();
            fmh_seeds(&read, &mut scalar, 1, k);
            let mut compat = Vec::new();
            fmh_seeds_avx2_compat(&read, &mut compat, 1, k);
            assert_eq!(compat.len(), 4 * ((read.len() - k + 1) / 4));
            // every compat k-mer is also a scalar k-mer (as a multiset)
            let mut scalar_sorted = scalar.clone();
            scalar_sorted.sort_unstable();
            let mut compat_sorted = compat.clone();
            compat_sorted.sort_unstable();
            assert!(compat_sorted.iter().all(|h| scalar_sorted.contains(h)));
        }
    }
}
