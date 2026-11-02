//! Golden parity tests for the sketching pipeline.
//!
//! `tests/golden/sketches.golden.txt` was produced by the **pre-refactor**
//! path-based sketch functions (the `expose-library-api` branch, commit
//! `4b5c2be`, i.e. before push-style sketchers existed) on the small fixtures
//! committed in `tests/golden/`. This test recomputes every case with
//! (a) the current path-based wrappers and (b) the push-style sketchers fed
//! directly, and requires byte-identical digests.
//!
//! Why this is not circular: the wrappers are now thin loops over the
//! sketchers, so comparing the two alone can only catch a divergence between
//! the loop and the sketcher, never a regression that affects both. The
//! goldens are an independent, frozen record of what the old implementation
//! produced.
//!
//! Host dependence: the old code (and `Seeding::Native`) uses the x86_64 AVX2
//! seeding kernel on AVX2-capable CPUs, which yields a different k-mer set
//! than the scalar kernel. The goldens were generated on such a host, so the
//! `Seeding::Native` checks are skipped (with a message) elsewhere. The
//! `Seeding::Avx2Compat` check is portable: the compat kernel emulates the
//! AVX2 kernel bit-for-bit, so it must reproduce the same goldens on **any**
//! architecture.
//!
//! Fixtures are regenerated deterministically by `tests/golden/gen_fixtures.py`.

#![cfg(feature = "cli")]

use sylph::constants::DEFAULT_FPR;
use sylph::seeding::Seeding;
use sylph::sketchers::{GenomeSketcher, PairedReadSketcher, ReadSketcher};
use sylph::types::{GenomeSketch, Kmer, SequencesSketch};

const GOLDEN: &str = include_str!("golden/sketches.golden.txt");

// ---------------------------------------------------------------------------
// Canonical digest. This code is byte-for-byte the same as the generator's
// (scratch worktree, `examples/gen_goldens.rs`); keep the two in sync.
// ---------------------------------------------------------------------------

fn seq_block(id: &str, sk: &SequencesSketch) -> String {
    let mut kmers: Vec<(Kmer, u32)> = sk.kmer_counts.iter().map(|(k, v)| (*k, *v)).collect();
    kmers.sort_unstable();
    let list = kmers
        .iter()
        .map(|(k, c)| format!("{}:{}", k, c))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "case={}\nkind=reads\nfile={}\nsample={}\nc={}\nk={}\npaired={}\nmean_bits={}\nn_kmers={}\nkmers={}\n",
        id,
        sk.file_name,
        sk.sample_name.as_deref().unwrap_or("-"),
        sk.c,
        sk.k,
        sk.paired as u8,
        sk.mean_read_length.to_bits(),
        kmers.len(),
        list
    )
}

fn genome_block(id: &str, g: &GenomeSketch) -> String {
    let list = g
        .genome_kmers
        .iter()
        .map(|k| k.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let pseudotax = match &g.pseudotax_tracked_nonused_kmers {
        None => "none".to_string(),
        Some(v) => {
            let mut v = v.clone();
            v.sort_unstable();
            if v.is_empty() {
                "empty".to_string()
            } else {
                v.iter().map(|k| k.to_string()).collect::<Vec<_>>().join(",")
            }
        }
    };
    format!(
        "case={}\nkind=genome\nfile={}\nfirst_contig={}\nc={}\nk={}\ngn_size={}\nmin_spacing={}\nn_kmers={}\nkmers={}\npseudotax={}\n",
        id,
        g.file_name,
        g.first_contig_name,
        g.c,
        g.k,
        g.gn_size,
        g.min_spacing,
        g.genome_kmers.len(),
        list,
        pseudotax
    )
}

// ---------------------------------------------------------------------------
// Case table. MUST stay identical to the generator's.
// ---------------------------------------------------------------------------

enum Case {
    Reads {
        id: &'static str,
        fq: &'static str,
        c: usize,
        k: usize,
        no_dedup: bool,
    },
    Pairs {
        id: &'static str,
        r1: &'static str,
        r2: &'static str,
        c: usize,
        k: usize,
        no_dedup: bool,
        fpr: f64,
    },
    Genome {
        id: &'static str,
        fa: &'static str,
        c: usize,
        k: usize,
        min_spacing: usize,
        pseudotax: bool,
    },
    Individual {
        id: &'static str,
        fa: &'static str,
        c: usize,
        k: usize,
        min_spacing: usize,
        pseudotax: bool,
    },
}

fn cases() -> Vec<Case> {
    use Case::*;
    vec![
        Reads { id: "reads_mixed_c1_k21_dedup", fq: "tests/golden/reads_mixed.fq", c: 1, k: 21, no_dedup: false },
        Reads { id: "reads_mixed_c5_k31_dedup", fq: "tests/golden/reads_mixed.fq", c: 5, k: 31, no_dedup: false },
        Reads { id: "reads_mixed_c100_k21_dedup", fq: "tests/golden/reads_mixed.fq", c: 100, k: 21, no_dedup: false },
        Reads { id: "reads_mixed_c100_k31_nodedup", fq: "tests/golden/reads_mixed.fq", c: 100, k: 31, no_dedup: true },
        Reads { id: "reads_mixed_c5_k21_nodedup", fq: "tests/golden/reads_mixed.fq", c: 5, k: 21, no_dedup: true },
        Reads { id: "reads_dups_c2_k21_dedup", fq: "tests/golden/reads_dups.fq", c: 2, k: 21, no_dedup: false },
        Reads { id: "reads_dups_c2_k31_dedup", fq: "tests/golden/reads_dups.fq", c: 2, k: 31, no_dedup: false },
        Reads { id: "reads_dups_c2_k21_nodedup", fq: "tests/golden/reads_dups.fq", c: 2, k: 21, no_dedup: true },
        Reads { id: "reads_dups_c5_k31_dedup", fq: "tests/golden/reads_dups.fq", c: 5, k: 31, no_dedup: false },
        Reads { id: "reads_repeats_c5_k21_dedup", fq: "tests/golden/reads_repeats.fq", c: 5, k: 21, no_dedup: false },
        Reads { id: "reads_repeats_c5_k31_dedup", fq: "tests/golden/reads_repeats.fq", c: 5, k: 31, no_dedup: false },
        Reads { id: "reads_repeats_c10_k21_nodedup", fq: "tests/golden/reads_repeats.fq", c: 10, k: 21, no_dedup: true },
        Pairs { id: "pair_c5_k21_dedup", r1: "tests/golden/pair_r1.fq", r2: "tests/golden/pair_r2.fq", c: 5, k: 21, no_dedup: false, fpr: 0.0 },
        Pairs { id: "pair_c20_k31_dedup", r1: "tests/golden/pair_r1.fq", r2: "tests/golden/pair_r2.fq", c: 20, k: 31, no_dedup: false, fpr: 0.0 },
        Pairs { id: "pair_c20_k21_nodedup", r1: "tests/golden/pair_r1.fq", r2: "tests/golden/pair_r2.fq", c: 20, k: 21, no_dedup: true, fpr: 0.0 },
        Pairs { id: "pair_c5_k31_approx", r1: "tests/golden/pair_r1.fq", r2: "tests/golden/pair_r2.fq", c: 5, k: 31, no_dedup: false, fpr: DEFAULT_FPR },
        Pairs { id: "pair_trunc_c5_k21_dedup", r1: "tests/golden/pair_trunc_r1.fq", r2: "tests/golden/pair_trunc_r2.fq", c: 5, k: 21, no_dedup: false, fpr: 0.0 },
        Pairs { id: "pair_trunc_c10_k31_dedup", r1: "tests/golden/pair_trunc_r1.fq", r2: "tests/golden/pair_trunc_r2.fq", c: 10, k: 31, no_dedup: false, fpr: 0.0 },
        Genome { id: "genome_multi_c5_k21_ms30", fa: "tests/golden/genome_multi.fasta", c: 5, k: 21, min_spacing: 30, pseudotax: false },
        Genome { id: "genome_multi_c20_k21_ms30_pt", fa: "tests/golden/genome_multi.fasta", c: 20, k: 21, min_spacing: 30, pseudotax: true },
        Genome { id: "genome_multi_c5_k31_ms0", fa: "tests/golden/genome_multi.fasta", c: 5, k: 31, min_spacing: 0, pseudotax: false },
        Genome { id: "genome_multi_c20_k31_ms1_pt", fa: "tests/golden/genome_multi.fasta", c: 20, k: 31, min_spacing: 1, pseudotax: true },
        Genome { id: "genome_spacing_c2_k21_ms0", fa: "tests/golden/genome_spacing.fasta", c: 2, k: 21, min_spacing: 0, pseudotax: false },
        Genome { id: "genome_spacing_c2_k21_ms1", fa: "tests/golden/genome_spacing.fasta", c: 2, k: 21, min_spacing: 1, pseudotax: false },
        Genome { id: "genome_spacing_c2_k21_ms2_pt", fa: "tests/golden/genome_spacing.fasta", c: 2, k: 21, min_spacing: 2, pseudotax: true },
        Genome { id: "genome_spacing_c5_k21_ms15_pt", fa: "tests/golden/genome_spacing.fasta", c: 5, k: 21, min_spacing: 15, pseudotax: true },
        Genome { id: "genome_spacing_c5_k31_ms30_pt", fa: "tests/golden/genome_spacing.fasta", c: 5, k: 31, min_spacing: 30, pseudotax: true },
        Genome { id: "genome_spacing_c20_k31_ms61", fa: "tests/golden/genome_spacing.fasta", c: 20, k: 31, min_spacing: 61, pseudotax: false },
        Genome { id: "genome_short_c1_k21_ms0", fa: "tests/golden/genome_short.fasta", c: 1, k: 21, min_spacing: 0, pseudotax: false },
        Genome { id: "genome_short_c1_k31_ms0", fa: "tests/golden/genome_short.fasta", c: 1, k: 31, min_spacing: 0, pseudotax: false },
        Individual { id: "individual_multi_c5_k21_ms30", fa: "tests/golden/genome_multi.fasta", c: 5, k: 21, min_spacing: 30, pseudotax: false },
        Individual { id: "individual_multi_c20_k31_ms0_pt", fa: "tests/golden/genome_multi.fasta", c: 20, k: 31, min_spacing: 0, pseudotax: true },
        Individual { id: "individual_short_c5_k21_ms2_pt", fa: "tests/golden/genome_short.fasta", c: 5, k: 21, min_spacing: 2, pseudotax: true },
    ]
}

// ---------------------------------------------------------------------------
// Minimal FASTA/FASTQ reader for the fixtures.
//
// The path-based wrappers are exercised through needletail above; this reader
// exists so the push-style sketchers can be driven directly without adding
// needletail as a dev-dependency. It mirrors needletail for these fixtures:
// headers are kept whole (so a tab survives into `push_contig`), and a FASTQ
// record whose quality length disagrees with its sequence is reported as
// malformed -- which the path-based pair loop skips, continuing with the next
// record.
// ---------------------------------------------------------------------------

struct Record {
    name: Vec<u8>,
    seq: Vec<u8>,
}

/// `None` stands for a record needletail would report as an error.
fn parse_records(path: &str) -> Vec<Option<Record>> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut out = Vec::new();
    if text.starts_with('>') {
        let mut cur: Option<Record> = None;
        for line in text.lines() {
            if let Some(header) = line.strip_prefix('>') {
                if let Some(r) = cur.take() {
                    out.push(Some(r));
                }
                cur = Some(Record {
                    name: header.as_bytes().to_vec(),
                    seq: Vec::new(),
                });
            } else if let Some(r) = cur.as_mut() {
                r.seq.extend_from_slice(line.trim_end().as_bytes());
            }
        }
        if let Some(r) = cur.take() {
            out.push(Some(r));
        }
    } else {
        let lines: Vec<&str> = text.lines().collect();
        let mut i = 0;
        while i + 1 < lines.len() {
            let header = lines[i];
            let seq = lines[i + 1];
            let qual = lines.get(i + 3).copied().unwrap_or("");
            assert!(
                header.starts_with('@'),
                "{path}: expected a FASTQ header at line {}",
                i + 1
            );
            assert_eq!(lines.get(i + 2), Some(&"+"), "{path}: missing '+' line");
            if qual.len() == seq.len() {
                out.push(Some(Record {
                    name: header[1..].as_bytes().to_vec(),
                    seq: seq.as_bytes().to_vec(),
                }));
            } else {
                // needletail reports this record as an error and its reader
                // yields nothing afterwards, so the pair loop stops here.
                out.push(None);
                break;
            }
            i += 4;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Drivers: path-based wrappers vs. push-style sketchers.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
enum Source {
    Path,
    Sketcher(Seeding),
}

fn drive_reads(
    fq: &str,
    c: usize,
    k: usize,
    no_dedup: bool,
    seeding: Seeding,
) -> SequencesSketch {
    let mut sketcher = ReadSketcher::new(c, k, no_dedup, seeding);
    for record in parse_records(fq) {
        if let Some(record) = record {
            sketcher.push_read(&record.seq);
        }
    }
    sketcher.finish(fq.to_string(), None)
}

fn drive_pair(
    r1: &str,
    r2: &str,
    c: usize,
    k: usize,
    no_dedup: bool,
    fpr: f64,
    seeding: Seeding,
) -> SequencesSketch {
    let mut sketcher = PairedReadSketcher::new(c, k, no_dedup, fpr, seeding).unwrap();
    let recs1 = parse_records(r1);
    let recs2 = parse_records(r2);
    for i in 0..recs1.len() {
        let (a, b) = match (recs1[i].as_ref(), recs2.get(i).and_then(|r| r.as_ref())) {
            (Some(a), Some(b)) => (a, b),
            _ => continue,
        };
        sketcher.push_pair(&a.seq, &b.seq);
    }
    sketcher.finish(r1.to_string(), None)
}

fn drive_genome(
    fa: &str,
    c: usize,
    k: usize,
    min_spacing: usize,
    pseudotax: bool,
    seeding: Seeding,
) -> GenomeSketch {
    let mut sketcher = GenomeSketcher::new(c, k, min_spacing, pseudotax, seeding);
    for record in parse_records(fa) {
        if let Some(record) = record {
            sketcher.push_contig(&record.name, &record.seq);
        }
    }
    sketcher.finish(fa.to_string())
}

fn drive_individual(
    fa: &str,
    c: usize,
    k: usize,
    min_spacing: usize,
    pseudotax: bool,
    seeding: Seeding,
) -> Vec<GenomeSketch> {
    parse_records(fa)
        .into_iter()
        .map(|record| {
            let record = record.expect("fixtures have no malformed FASTA records");
            let mut sketcher = GenomeSketcher::new(c, k, min_spacing, pseudotax, seeding);
            sketcher.push_contig(&record.name, &record.seq);
            sketcher.finish(fa.to_string())
        })
        .collect()
}

fn blocks_for(case: &Case, source: Source) -> Vec<String> {
    use sylph::sketch::{
        sketch_genome, sketch_genome_individual, sketch_pair_sequences, sketch_sequences_needle,
    };
    match case {
        Case::Reads { id, fq, c, k, no_dedup } => {
            let sk = match source {
                Source::Path => {
                    sketch_sequences_needle(fq, *c, *k, None, *no_dedup).expect("readable fixture")
                }
                Source::Sketcher(seeding) => drive_reads(fq, *c, *k, *no_dedup, seeding),
            };
            vec![seq_block(id, &sk)]
        }
        Case::Pairs { id, r1, r2, c, k, no_dedup, fpr } => {
            let sk = match source {
                Source::Path => {
                    sketch_pair_sequences(r1, r2, *c, *k, None, *no_dedup, *fpr)
                        .expect("readable fixture")
                }
                Source::Sketcher(seeding) => {
                    drive_pair(r1, r2, *c, *k, *no_dedup, *fpr, seeding)
                }
            };
            vec![seq_block(id, &sk)]
        }
        Case::Genome { id, fa, c, k, min_spacing, pseudotax } => {
            let sk = match source {
                Source::Path => sketch_genome(*c, *k, fa, *min_spacing, *pseudotax)
                    .expect("readable fixture"),
                Source::Sketcher(seeding) => {
                    drive_genome(fa, *c, *k, *min_spacing, *pseudotax, seeding)
                }
            };
            vec![genome_block(id, &sk)]
        }
        Case::Individual { id, fa, c, k, min_spacing, pseudotax } => {
            let sketches = match source {
                Source::Path => {
                    sketch_genome_individual(*c, *k, fa, *min_spacing, *pseudotax)
                }
                Source::Sketcher(seeding) => {
                    drive_individual(fa, *c, *k, *min_spacing, *pseudotax, seeding)
                }
            };
            assert!(!sketches.is_empty(), "{id}: no records in {fa}");
            sketches
                .iter()
                .enumerate()
                .map(|(i, sk)| genome_block(&format!("{}#{}", id, i), sk))
                .collect()
        }
    }
}

fn render(source: Source) -> Vec<String> {
    cases().iter().flat_map(|c| blocks_for(c, source)).collect()
}

fn compare_blocks(blocks: &[String]) {
    let got = blocks.join("\n");
    if got == GOLDEN {
        return;
    }
    let expected_blocks: Vec<&str> = GOLDEN.split("\n\n").collect();
    let got_blocks: Vec<&str> = got.split("\n\n").collect();
    for (i, (want, have)) in expected_blocks.iter().zip(got_blocks.iter()).enumerate() {
        if want != have {
            panic!(
                "golden block #{i} differs\n--- golden ---\n{want}\n--- current code ---\n{have}"
            );
        }
    }
    panic!(
        "golden has {} blocks / {} bytes, current code produced {} blocks / {} bytes",
        expected_blocks.len(),
        GOLDEN.len(),
        got_blocks.len(),
        got.len()
    );
}

#[test]
fn golden_file_has_the_expected_shape() {
    assert_eq!(
        GOLDEN.split("\n\n").count(),
        39,
        "sketches.golden.txt was regenerated or truncated"
    );
    // Cheap sanity check that the digest is not vacuous.
    assert!(GOLDEN.lines().any(|l| l.starts_with("n_kmers=")));
    let duplicates = GOLDEN
        .lines()
        .filter(|l| l.starts_with("kmers="))
        .any(|l| l.contains(":2") || l.contains(":3") || l.contains(":5"));
    assert!(duplicates, "goldens should cover k-mer counts > 1");
}

#[cfg(target_arch = "x86_64")]
fn native_is_avx2() -> bool {
    is_x86_feature_detected!("avx2")
}

#[cfg(not(target_arch = "x86_64"))]
fn native_is_avx2() -> bool {
    false
}

#[test]
fn goldens_match_path_based_functions_native() {
    if !native_is_avx2() {
        eprintln!(
            "SKIPPED: goldens were generated with the x86_64 AVX2 seeding kernel; \
             this host does not run it (target_arch = {}, avx2 = false). \
             The portable check is `goldens_match_streaming_sketchers_avx2_compat`.",
            std::env::consts::ARCH
        );
        return;
    }
    compare_blocks(&render(Source::Path));
}

#[test]
fn goldens_match_streaming_sketchers_native() {
    if !native_is_avx2() {
        eprintln!("SKIPPED: Seeding::Native is not the AVX2 kernel on this host");
        return;
    }
    compare_blocks(&render(Source::Sketcher(Seeding::Native)));
}

/// Portable: `Seeding::Avx2Compat` emulates the AVX2 kernel on every
/// architecture, so it must reproduce the AVX2-generated goldens everywhere.
#[test]
fn goldens_match_streaming_sketchers_avx2_compat() {
    compare_blocks(&render(Source::Sketcher(Seeding::Avx2Compat)));
}
