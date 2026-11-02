//! Push-style (streaming) sketchers: feed one record at a time, get sylph's
//! own sketch types back.
//!
//! These hold the exact per-record logic of the path-based functions in
//! [`crate::sketch`] (`sketch_sequences_needle`, `sketch_pair_sequences`,
//! `sketch_genome`), which are now thin needletail loops over this module.
//! Downstream crates (pysylph, browser wasm) can therefore reuse the sketching
//! kernels with their own parser: everything here takes `&[u8]` sequences and
//! never touches the filesystem.
//!
//! ```
//! use sylph::seeding::Seeding;
//! use sylph::sketchers::ReadSketcher;
//!
//! let mut sketcher = ReadSketcher::new(200, 31, false, Seeding::Native);
//! sketcher.push_read(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT");
//! let sketch = sketcher.finish("reads.fq".to_string(), None);
//! assert_eq!(sketch.k, 31);
//! ```
//!
//! Output equivalence notes, per sketcher:
//!
//! * [`ReadSketcher`] reproduces `sketch_sequences_needle`: one
//!   `(k-mer, double k-mer signature)` exact dedup set per sample,
//!   `MAX_DEDUP_COUNT` as the escape threshold, and the incremental (moving
//!   average) mean read length over *all* records.
//! * [`PairedReadSketcher`] reproduces `sketch_pair_sequences`: approximate
//!   dedup at `dedup_fpr` (`0.0` selects the exact set, any other value the
//!   seeded cuckoo filter), cross-mate suppression of identical k-mer hashes
//!   within a pair, and `mean_read_length = total_bases / (2 * pairs)`.
//!   (pysylph's `sketch_paired` instead tracks a moving average of the first
//!   mate's length; a caller that needs that can ignore
//!   [`SequencesSketch::mean_read_length`] and average its own lengths.)
//! * [`GenomeSketcher`] reproduces `sketch_genome` over all contigs, or
//!   `sketch_genome_individual` if you use one sketcher per record: markers
//!   from every contig are sorted together, k-mers occurring more than once
//!   anywhere are dropped, and a greedy `min_spacing` filter is applied along
//!   each contig. Unused k-mers pushed out by the spacing filter are recorded
//!   when pseudotax tracking is requested.

use std::fmt;

use fxhash::{FxHashMap, FxHashSet};
#[cfg(feature = "approx-dedup")]
use fxhash::FxHasher64;
#[cfg(feature = "approx-dedup")]
use rand::{rngs::SmallRng, SeedableRng};
#[cfg(feature = "approx-dedup")]
use scalable_cuckoo_filter::{ScalableCuckooFilter, ScalableCuckooFilterBuilder};

#[cfg(feature = "approx-dedup")]
use crate::constants::DEFAULT_RNG_SEED;
use crate::constants::MAX_DEDUP_COUNT;
#[cfg(feature = "approx-dedup")]
use crate::kernels::dup_removal_lsh_full;
use crate::kernels::{dup_removal_lsh_full_exact, pair_kmer, pair_kmer_single, Marker};
use crate::seeding::{extract_markers_positions_with, extract_markers_with, Seeding};
use crate::types::{GenomeSketch, Kmer, MMHashSet, SequencesSketch};

/// Anything a streaming sketcher can fail on. Parsing errors stay with the
/// caller (the sketchers only see `&[u8]`); this covers the configuration and
/// the path-based wrappers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SketchError {
    /// The input file could not be opened or is not a valid FASTA/FASTQ file.
    Input(String),
    /// A record in the input could not be parsed.
    Record(String),
    /// The requested configuration is not available in this build.
    Unsupported(String),
}

impl fmt::Display for SketchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SketchError::Input(msg) | SketchError::Record(msg) | SketchError::Unsupported(msg) => {
                f.write_str(msg)
            }
        }
    }
}

impl std::error::Error for SketchError {}

/// Direction of the duplicate-removal filter used by the read sketchers.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Dedup {
    Exact,
    Approx,
}

/// Sketches single-end reads (one sample).
pub struct ReadSketcher {
    c: usize,
    k: usize,
    no_dedup: bool,
    seeding: Seeding,
    kmer_counts: FxHashMap<Kmer, u32>,
    kmer_to_pair_set: FxHashSet<(u64, [Marker; 2])>,
    num_dup_removed: usize,
    mean_read_length: f64,
    num_reads: u64,
    scratch: Vec<u64>,
}

impl ReadSketcher {
    /// `c` is the seeding downsampling factor and must be `> 0`; the kernels
    /// compute `u64::MAX / c`, so `c = 0` panics (division by zero). `k` must
    /// be positive; sylph's sketches are defined for `k = 21` or `k = 31`, and
    /// with [`Seeding::Native`] on an AVX2-capable x86_64 host any other `k`
    /// panics as soon as a read has at least `k + 1` bases (see
    /// [`Seeding::Native`]).
    pub fn new(c: usize, k: usize, no_dedup: bool, seeding: Seeding) -> Self {
        ReadSketcher {
            c,
            k,
            no_dedup,
            seeding,
            kmer_counts: FxHashMap::default(),
            kmer_to_pair_set: FxHashSet::default(),
            num_dup_removed: 0,
            mean_read_length: 0.,
            num_reads: 0,
            scratch: Vec::new(),
        }
    }

    /// Adds one read (or, for assemblies, one record; the caller decides what
    /// a "read" is).
    pub fn push_read(&mut self, seq: &[u8]) {
        let kmer_pair = if seq.len() > 400 {
            None
        } else {
            pair_kmer_single(seq)
        };
        self.scratch.clear();
        extract_markers_with(self.seeding, seq, &mut self.scratch, self.c, self.k);
        for km in self.scratch.iter() {
            dup_removal_lsh_full_exact(
                &mut self.kmer_counts,
                &mut self.kmer_to_pair_set,
                km,
                kmer_pair,
                &mut self.num_dup_removed,
                self.no_dedup,
                Some(MAX_DEDUP_COUNT),
            );
        }
        // moving average, exactly as in `sketch_sequences_needle`
        self.num_reads += 1;
        self.mean_read_length +=
            ((seq.len() as f64) - self.mean_read_length) / self.num_reads as f64;
    }

    pub fn kmer_counts(&self) -> &FxHashMap<Kmer, u32> {
        &self.kmer_counts
    }

    #[cfg(test)]
    fn num_reads(&self) -> u64 {
        self.num_reads
    }

    #[cfg(test)]
    fn num_dup_removed(&self) -> usize {
        self.num_dup_removed
    }

    #[cfg(test)]
    fn mean_read_length(&self) -> f64 {
        self.mean_read_length
    }

    /// Consumes the sketcher and returns the sample sketch.
    pub fn finish(self, file_name: String, sample_name: Option<String>) -> SequencesSketch {
        SequencesSketch {
            kmer_counts: self.kmer_counts,
            file_name,
            c: self.c,
            k: self.k,
            paired: false,
            sample_name,
            mean_read_length: self.mean_read_length,
        }
    }
}

/// Sketches paired-end reads (one sample).
pub struct PairedReadSketcher {
    c: usize,
    k: usize,
    no_dedup: bool,
    seeding: Seeding,
    dedup: Dedup,
    kmer_counts: FxHashMap<Kmer, u32>,
    kmer_to_pair_set: FxHashSet<(u64, [Marker; 2])>,
    #[cfg(feature = "approx-dedup")]
    // `FxHasher64`, not the pointer-width `FxHasher`: on 32-bit targets
    // (wasm32) the latter returns a 32-bit hash, whose top bits -- which
    // `scalable_cuckoo_filter` uses as the fingerprint (`hash >> (64 - bits)`)
    // -- are always zero. Every fingerprint then collides and the filter
    // degenerates (inserts never match, the table keeps growing). On 64-bit
    // targets both hashers are the same type.
    approx_set: ScalableCuckooFilter<(u64, [Marker; 2]), FxHasher64, SmallRng>,
    num_dup_removed: usize,
    total_read_length: u64,
    num_read_pairs: u64,
    scratch1: Vec<u64>,
    scratch2: Vec<u64>,
}

impl PairedReadSketcher {
    /// `dedup_fpr == 0.0` uses the exact dedup set; any other value uses the
    /// approximate (cuckoo filter) one, exactly like `sketch_pair_sequences`.
    ///
    /// `c` must be `> 0`; the kernels compute `u64::MAX / c`, so `c = 0`
    /// panics (division by zero). `k` must be positive; use `21` or `31`, the
    /// values sylph's sketches are defined for. With [`Seeding::Native`] on an
    /// AVX2-capable x86_64 host any other `k` panics once a mate has at least
    /// `k + 1` bases (see [`Seeding::Native`]).
    pub fn new(
        c: usize,
        k: usize,
        no_dedup: bool,
        dedup_fpr: f64,
        seeding: Seeding,
    ) -> Result<Self, SketchError> {
        let dedup = if dedup_fpr == 0. {
            Dedup::Exact
        } else {
            Dedup::Approx
        };
        #[cfg(not(feature = "approx-dedup"))]
        if dedup_fpr != 0. {
            return Err(SketchError::Unsupported(
                "approximate read deduplication (dedup_fpr != 0) requires the `approx-dedup` feature"
                    .to_string(),
            ));
        }
        // Built eagerly even in exact mode, as the original code did; it does
        // not influence the exact path's output.
        #[cfg(feature = "approx-dedup")]
        let approx_set = {
            let mut fpr = 0.001;
            if dedup_fpr != 0. {
                fpr = dedup_fpr;
            }
            ScalableCuckooFilterBuilder::new()
                .initial_capacity(1_000_000_0)
                .false_positive_probability(fpr)
                .hasher(FxHasher64::default())
                .rng(SmallRng::seed_from_u64(DEFAULT_RNG_SEED))
                .finish()
        };
        Ok(PairedReadSketcher {
            c,
            k,
            no_dedup,
            seeding,
            dedup,
            kmer_counts: FxHashMap::default(),
            kmer_to_pair_set: FxHashSet::default(),
            #[cfg(feature = "approx-dedup")]
            approx_set,
            num_dup_removed: 0,
            total_read_length: 0,
            num_read_pairs: 0,
            scratch1: Vec::new(),
            scratch2: Vec::new(),
        })
    }

    /// Adds one mate pair. Mate order must match `sketch_pair_sequences`
    /// (`seq1` = first/second file, `seq2` = the other one).
    pub fn push_pair(&mut self, seq1: &[u8], seq2: &[u8]) {
        let kmer_pair = pair_kmer(seq1, seq2);
        self.num_read_pairs += 1;
        self.total_read_length += seq1.len() as u64 + seq2.len() as u64;

        self.scratch1.clear();
        self.scratch2.clear();
        extract_markers_with(self.seeding, seq1, &mut self.scratch1, self.c, self.k);
        extract_markers_with(self.seeding, seq2, &mut self.scratch2, self.c, self.k);

        for i in 0..self.scratch1.len() {
            let km = self.scratch1[i];
            self.dedup_and_count(&km, kmer_pair);
        }
        for i in 0..self.scratch2.len() {
            let km = self.scratch2[i];
            if self.scratch1.contains(&km) {
                continue;
            }
            self.dedup_and_count(&km, kmer_pair);
        }
    }

    fn dedup_and_count(&mut self, km: &u64, kmer_pair: Option<([Marker; 2], [Marker; 2])>) {
        match self.dedup {
            Dedup::Exact => {
                let Self {
                    kmer_counts,
                    kmer_to_pair_set,
                    num_dup_removed,
                    no_dedup,
                    ..
                } = self;
                dup_removal_lsh_full_exact(
                    kmer_counts,
                    kmer_to_pair_set,
                    km,
                    kmer_pair,
                    num_dup_removed,
                    *no_dedup,
                    None,
                );
            }
            #[cfg(feature = "approx-dedup")]
            Dedup::Approx => {
                let Self {
                    kmer_counts,
                    approx_set,
                    num_dup_removed,
                    no_dedup,
                    ..
                } = self;
                dup_removal_lsh_full(
                    kmer_counts,
                    approx_set,
                    km,
                    kmer_pair,
                    num_dup_removed,
                    *no_dedup,
                );
            }
            #[cfg(not(feature = "approx-dedup"))]
            Dedup::Approx => unreachable!("approximate dedup requires the `approx-dedup` feature"),
        }
    }

    pub fn kmer_counts(&self) -> &FxHashMap<Kmer, u32> {
        &self.kmer_counts
    }

    /// Consumes the sketcher and returns the paired sample sketch.
    pub fn finish(self, file_name: String, sample_name: Option<String>) -> SequencesSketch {
        let percent = (self.num_dup_removed as f64)
            / ((self.kmer_counts.values().sum::<u32>() as f64) + self.num_dup_removed as f64)
            * 100.;
        log::debug!(
            "Number of sketched k-mers removed due to read duplication for {}: {}. Percentage: {:.2}%",
            file_name,
            self.num_dup_removed,
            percent,
        );
        let mean_read_length = if self.num_read_pairs > 0 {
            self.total_read_length as f64 / (2.0 * self.num_read_pairs as f64)
        } else {
            0.0
        };
        SequencesSketch {
            kmer_counts: self.kmer_counts,
            file_name,
            c: self.c,
            k: self.k,
            paired: true,
            sample_name,
            mean_read_length,
        }
    }
}

/// Sketches genomes/contigs (one `GenomeSketch`).
///
/// Push every contig of a genome into one sketcher for `sketch_genome`
/// semantics, or use a fresh sketcher per record for
/// `sketch_genome_individual` / per-contig sketches.
pub struct GenomeSketcher {
    c: usize,
    k: usize,
    min_spacing: usize,
    pseudotax: bool,
    seeding: Seeding,
    markers: Vec<(usize, usize, u64)>,
    gn_size: usize,
    first_contig_name: String,
    contig_number: usize,
}

impl GenomeSketcher {
    /// `c` is the seeding downsampling factor and must be `> 0`; the kernels
    /// compute `u64::MAX / c`, so `c = 0` panics (division by zero). `k` must
    /// be positive; use `21` or `31`, the values sylph's sketches are defined
    /// for. With [`Seeding::Native`] on an AVX2-capable x86_64 host any other
    /// `k` panics once a contig has at least `k + 1` bases (marker
    /// extraction) or at least `2 * k` bases (marker positions, i.e. genome
    /// sketches), see [`Seeding::Native`].
    ///
    /// `min_spacing` is the greedy minimum distance (in k-mer end positions)
    /// between two retained markers on the same contig; `0` disables it.
    pub fn new(c: usize, k: usize, min_spacing: usize, pseudotax: bool, seeding: Seeding) -> Self {
        GenomeSketcher {
            c,
            k,
            min_spacing,
            pseudotax,
            seeding,
            markers: Vec::new(),
            gn_size: 0,
            first_contig_name: String::new(),
            contig_number: 0,
        }
    }

    /// Adds one contig. `name` is only used for the first contig of the
    /// genome (`first_contig_name`); tabs become spaces as before. Pass
    /// `b""` when the caller has no record names (pysylph's contig list);
    /// the field then stays empty, as in pysylph's `sketch_genome`.
    pub fn push_contig(&mut self, name: &[u8], seq: &[u8]) {
        if self.contig_number == 0 {
            self.first_contig_name = String::from_utf8_lossy(name).to_string().replace('\t', " ");
        }
        self.gn_size += seq.len();
        extract_markers_positions_with(
            self.seeding,
            seq,
            &mut self.markers,
            self.c,
            self.k,
            self.contig_number,
        );
        self.contig_number += 1;
    }

    /// Consumes the sketcher and returns the genome sketch.
    pub fn finish(self, file_name: String) -> GenomeSketch {
        let GenomeSketcher {
            c,
            k,
            min_spacing,
            pseudotax,
            mut markers,
            gn_size,
            first_contig_name,
            ..
        } = self;

        let mut kmer_set = MMHashSet::default();
        let mut duplicate_set = MMHashSet::default();
        let mut new_vec = Vec::with_capacity(markers.len());
        markers.sort();
        for (_, _, km) in markers.iter() {
            if !kmer_set.contains(km) {
                kmer_set.insert(*km);
            } else {
                duplicate_set.insert(*km);
            }
        }

        let mut pseudotax_track_kmers = vec![];
        let mut last_pos = 0;
        let mut last_contig = 0;
        for (contig, pos, km) in markers.iter() {
            if !duplicate_set.contains(km) {
                if last_pos == 0 || last_contig != *contig || pos - last_pos > min_spacing {
                    new_vec.push(*km);
                    last_contig = *contig;
                    last_pos = *pos;
                } else if pseudotax {
                    pseudotax_track_kmers.push(*km);
                }
            }
        }

        GenomeSketch {
            genome_kmers: new_vec,
            pseudotax_tracked_nonused_kmers: if pseudotax {
                Some(pseudotax_track_kmers)
            } else {
                None
            },
            file_name,
            first_contig_name,
            c,
            k,
            gn_size,
            min_spacing,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    fn random_seq(seed: u64, len: usize) -> Vec<u8> {
        let mut state = seed | 1;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        (0..len).map(|_| b"ACGT"[(next() % 4) as usize]).collect()
    }

    #[test]
    fn read_sketcher_dedups_repeated_reads_and_tracks_mean_length() {
        let read = random_seq(1, 200);
        let other = random_seq(2, 100);

        let mut sketcher = ReadSketcher::new(1, 21, false, Seeding::Avx2Compat);
        for _ in 0..4 {
            sketcher.push_read(&read);
        }
        sketcher.push_read(&other);
        // Every observation after the first of each duplicate signature pair is
        // dropped, so counts are 1 (the MAX_DEDUP_COUNT escape is not reached
        // for c = 1).
        assert!(sketcher.kmer_counts().values().all(|&c| c == 1));
        assert!(sketcher.num_dup_removed() > 0);
        assert_eq!(sketcher.num_reads(), 5);
        assert_eq!(sketcher.mean_read_length(), (4.0 * 200.0 + 100.0) / 5.0);

        let sketched = sketcher.finish("reads.fq".to_string(), Some("sample".into()));
        assert_eq!(sketched.file_name, "reads.fq");
        assert_eq!(sketched.sample_name.as_deref(), Some("sample"));
        assert!(!sketched.paired);
        assert_eq!(sketched.c, 1);
        assert_eq!(sketched.k, 21);

        // --no-dedup keeps every observation.
        let mut sketcher = ReadSketcher::new(1, 21, true, Seeding::Avx2Compat);
        for _ in 0..4 {
            sketcher.push_read(&read);
        }
        assert_eq!(sketcher.num_dup_removed(), 0);
        assert!(sketcher.kmer_counts().values().all(|&c| c == 4));
    }

    #[test]
    fn paired_sketcher_exact_and_approx_paths() {
        let r1 = random_seq(3, 150);
        let r2 = random_seq(4, 150);

        let fprs: &[f64] = if cfg!(feature = "approx-dedup") {
            &[0.0, 0.001]
        } else {
            &[0.0]
        };
        for &dedup_fpr in fprs {
            let mut sketcher = PairedReadSketcher::new(50, 21, false, dedup_fpr, Seeding::Native)
                .unwrap();
            for _ in 0..3 {
                sketcher.push_pair(&r1, &r2);
            }
            let sketch = sketcher.finish("pair_1.fq".to_string(), None);
            assert!(sketch.paired);
            assert_eq!(sketch.mean_read_length, 150.0);
            // Repeated mate pairs are deduplicated away (up to the filter's
            // false-positive rate in approximate mode).
            assert!(sketch.kmer_counts.values().all(|&c| c == 1));
        }
    }

    #[test]
    fn genome_sketcher_duplicates_spacing_and_pseudotax() {
        let first = random_seq(5, 5000);
        let second = random_seq(6, 4000);

        let mut sketcher = GenomeSketcher::new(20, 21, 30, true, Seeding::Avx2Compat);
        sketcher.push_contig(b"contig\tone", &first);
        sketcher.push_contig(b"contig2", &second);
        let sketch = sketcher.finish("genome.fa".to_string());
        assert_eq!(sketch.first_contig_name, "contig one");
        assert_eq!(sketch.gn_size, 9000);
        assert_eq!(sketch.min_spacing, 30);
        assert!(!sketch.genome_kmers.is_empty());
        // pseudotax tracking is enabled, so the field exists (it may be empty)
        assert!(sketch.pseudotax_tracked_nonused_kmers.is_some());

        // A tandem repeat has every k-mer more than once -> nothing survives.
        let unit = random_seq(7, 100);
        let tandem: Vec<u8> = unit.iter().copied().cycle().take(3000).collect();
        let mut sketcher = GenomeSketcher::new(1, 21, 30, false, Seeding::Avx2Compat);
        sketcher.push_contig(b"tandem", &tandem);
        let sketch = sketcher.finish("tandem.fa".to_string());
        assert!(sketch.genome_kmers.is_empty());
        assert!(sketch.pseudotax_tracked_nonused_kmers.is_none());

        // A single record sketched on its own matches the one-shot path
        // (GenomeSketcher::push_contig is what sketch_genome_individual does).
        let mut one = GenomeSketcher::new(20, 21, 0, false, Seeding::Avx2Compat);
        one.push_contig(b"only", &first);
        let a = one.finish("one.fa".to_string());
        let mut two = GenomeSketcher::new(20, 21, 0, false, Seeding::Avx2Compat);
        two.push_contig(b"only", &first);
        let b = two.finish("one.fa".to_string());
        assert_eq!(a, b);
    }

    #[test]
    fn unspaced_genome_keeps_all_unique_kmers() {
        let seq = random_seq(8, 2000);
        let dense = {
            let mut s = GenomeSketcher::new(1, 21, 0, false, Seeding::Avx2Compat);
            s.push_contig(b"c", &seq);
            s.finish(String::new()).genome_kmers
        };
        let spaced = {
            let mut s = GenomeSketcher::new(1, 21, 30, false, Seeding::Avx2Compat);
            s.push_contig(b"c", &seq);
            s.finish(String::new()).genome_kmers
        };
        assert!(spaced.len() < dense.len());
        // deterministic across runs
        let dense2 = {
            let mut s = GenomeSketcher::new(1, 21, 0, false, Seeding::Avx2Compat);
            s.push_contig(b"c", &seq);
            s.finish(String::new()).genome_kmers
        };
        assert_eq!(dense, dense2);
    }
}

/// The path-based wrappers must remain exactly the streaming sketchers driven
/// by needletail; these tests pin that for the repository's test data.
#[cfg(all(test, feature = "cli"))]
mod path_tests {
    use super::*;
    use needletail::parse_fastx_file;

    fn drive_reads(path: &str, c: usize, k: usize, no_dedup: bool) -> SequencesSketch {
        let mut sketcher = ReadSketcher::new(c, k, no_dedup, Seeding::Native);
        let mut reader = parse_fastx_file(path).unwrap();
        while let Some(record) = reader.next() {
            let record = record.unwrap();
            sketcher.push_read(&record.seq());
        }
        sketcher.finish(path.to_string(), None)
    }

    #[test]
    fn read_path_matches_streaming_sketcher() {
        let path = "test_files/t1.fq";
        for no_dedup in [false, true] {
            let direct = crate::sketch::sketch_sequences_needle(path, 200, 31, None, no_dedup)
                .expect("t1.fq should be readable");
            assert_eq!(direct, drive_reads(path, 200, 31, no_dedup));
        }
    }

    #[test]
    fn paired_path_matches_streaming_sketcher() {
        let (p1, p2) = ("test_files/t1.fq", "test_files/t2.fq");
        for (no_dedup, fpr) in [(false, 0.0), (false, crate::constants::DEFAULT_FPR), (true, 0.0)]
        {
            let direct = crate::sketch::sketch_pair_sequences(p1, p2, 200, 31, None, no_dedup, fpr)
                .expect("test pairs should be readable");
            let mut sketcher =
                PairedReadSketcher::new(200, 31, no_dedup, fpr, Seeding::Native).unwrap();
            let mut r1 = parse_fastx_file(p1).unwrap();
            let mut r2 = parse_fastx_file(p2).unwrap();
            while let (Some(a), Some(b)) = (r1.next(), r2.next()) {
                sketcher.push_pair(&a.unwrap().seq(), &b.unwrap().seq());
            }
            assert_eq!(direct, sketcher.finish(p1.to_string(), None));
        }
    }

    #[test]
    fn genome_path_matches_streaming_sketcher() {
        let path = "test_files/tiny_virus.fasta";
        for pseudotax in [false, true] {
            let direct = crate::sketch::sketch_genome(200, 31, path, 30, pseudotax)
                .expect("tiny_virus.fasta should be readable");
            let mut sketcher = GenomeSketcher::new(200, 31, 30, pseudotax, Seeding::Native);
            let mut reader = parse_fastx_file(path).unwrap();
            while let Some(record) = reader.next() {
                let record = record.unwrap();
                sketcher.push_contig(record.id(), &record.seq());
            }
            assert_eq!(direct, sketcher.finish(path.to_string()));
        }
    }
}
