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
//! * [`GenomeSketcher`] reproduces `sketch_genome` over all contigs, or
//!   `sketch_genome_individual` if you use one sketcher per record: markers
//!   from every contig are sorted together, k-mers occurring more than once
//!   anywhere are dropped, and a greedy `min_spacing` filter is applied along
//!   each contig. Unused k-mers pushed out by the spacing filter are recorded
//!   when pseudotax tracking is requested.

use std::fmt;

use fxhash::{FxHashMap, FxHashSet, FxHasher};
use rand::{rngs::SmallRng, SeedableRng};
use scalable_cuckoo_filter::{ScalableCuckooFilter, ScalableCuckooFilterBuilder};

use crate::constants::{DEFAULT_RNG_SEED, MAX_DEDUP_COUNT};
use crate::kernels::{dup_removal_lsh_full, dup_removal_lsh_full_exact, pair_kmer, pair_kmer_single, Marker};
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

    pub fn num_reads(&self) -> u64 {
        self.num_reads
    }

    pub fn num_dup_removed(&self) -> usize {
        self.num_dup_removed
    }

    pub fn mean_read_length(&self) -> f64 {
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
    approx_set: ScalableCuckooFilter<(u64, [Marker; 2]), FxHasher, SmallRng>,
    num_dup_removed: usize,
    total_read_length: u64,
    num_read_pairs: u64,
    scratch1: Vec<u64>,
    scratch2: Vec<u64>,
}

impl PairedReadSketcher {
    /// `dedup_fpr == 0.0` uses the exact dedup set; any other value uses the
    /// approximate (cuckoo filter) one, exactly like `sketch_pair_sequences`.
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
        let mut fpr = 0.001;
        if dedup_fpr != 0. {
            fpr = dedup_fpr;
        }
        // Built eagerly even in exact mode, as the original code did; it does
        // not influence the exact path's output.
        let approx_set = ScalableCuckooFilterBuilder::new()
            .initial_capacity(1_000_000_0)
            .false_positive_probability(fpr)
            .hasher(FxHasher::default())
            .rng(SmallRng::seed_from_u64(DEFAULT_RNG_SEED))
            .finish();
        Ok(PairedReadSketcher {
            c,
            k,
            no_dedup,
            seeding,
            dedup,
            kmer_counts: FxHashMap::default(),
            kmer_to_pair_set: FxHashSet::default(),
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
        }
    }

    pub fn kmer_counts(&self) -> &FxHashMap<Kmer, u32> {
        &self.kmer_counts
    }

    pub fn num_read_pairs(&self) -> u64 {
        self.num_read_pairs
    }

    pub fn num_dup_removed(&self) -> usize {
        self.num_dup_removed
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
    /// genome (`first_contig_name`); tabs become spaces as before.
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

    pub fn num_contigs(&self) -> usize {
        self.contig_number
    }

    pub fn genome_size(&self) -> usize {
        self.gn_size
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
