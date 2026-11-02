# Golden sketching fixtures

`sketches.golden.txt` is a frozen record of what the **pre-refactor**,
path-based sketching functions produce on the fixtures in this directory.
`tests/golden_parity.rs` requires the current code (path wrappers *and* the
push-style sketchers) to reproduce it byte for byte.

## Why

The equivalence tests in `src/sketchers.rs` compare the path-based wrappers
with the push-style sketchers they are now implemented with, so they cannot
see a bug that affects both. These goldens come from `expose-library-api`
(commit `4b5c2be`, before the sketchers existed) and are independent of the
current implementation.

## Fixtures

Generated deterministically (fixed-seed xorshift64, Python standard library
only) by `gen_fixtures.py`:

    python3 tests/golden/gen_fixtures.py

| file | covers |
| --- | --- |
| `reads_mixed.fq` | reads shorter than / exactly / slightly longer than k (for k = 21 and 31), N bases, lower-case bases, an IUPAC code, an all-N read, reads over 400 bp (no double-k-mer signature) |
| `reads_dups.fq` | 8 identical copies of one read, 3 of another; reads sharing a common core behind different prefixes (same k-mers, different signatures) including a repeat of the signature seen exactly when the count reaches `MAX_DEDUP_COUNT`, and a unique read |
| `reads_repeats.fq` | tandem unit repeat, internal duplicated segment, duplicated 30 bp block, low-complexity `ACGT` repeat |
| `pair_r1.fq`, `pair_r2.fq` | unequal mate lengths, mates below the 33 bp signature minimum, a duplicated mate pair |
| `pair_trunc_r1.fq`, `pair_trunc_r2.fq` | malformed (too-short quality) second mate and fewer records in mate 2 |
| `genome_multi.fasta` | three contigs, tab in the first name, N/lower-case contig, tandem-repeat contig |
| `genome_spacing.fasta` | two contigs, used for `min_spacing` 0/1/2/15/30/61 and cross-contig boundaries |
| `genome_short.fasta` | contigs shorter than / equal to / just over `2 * k` |

## Goldens

`sketches.golden.txt` was produced on an x86_64 host with AVX2, with a small
generator binary (`examples/gen_goldens.rs`) that was run **from the old
checkout** in a scratch worktree and is deliberately not part of this branch.
It called
`sketch_sequences_needle`, `sketch_pair_sequences`, `sketch_genome` and
`sketch_genome_individual` over the case table in `tests/golden_parity.rs` and
wrote the canonical digest below. To regenerate (only needed if the case table
or a fixture changes):

1. `git worktree add --detach <scratch>/old 4b5c2be`
2. symlink `<scratch>/old/tests/golden` to this directory,
3. add `examples/gen_goldens.rs` to `<scratch>/old` with the same case table
   and the `seq_block` / `genome_block` functions from
   `tests/golden_parity.rs` (they are the canonicalisation),
4. `cargo run --release --example gen_goldens` from `<scratch>/old`, and
5. move the new `sketches.golden.txt` back here.

### Digest format

One block per sketch, blocks separated by a blank line, `key=value` lines:

* reads: `file_name`, `sample_name` (`-` for none), `c`, `k`, `paired`,
  `mean_bits` (`f64::to_bits` of `mean_read_length`, so it is exact and
  host-independent), `n_kmers`, `kmers` = `hash:count` pairs sorted by hash.
* genome: `file_name`, `first_contig`, `c`, `k`, `gn_size`, `min_spacing`,
  `n_kmers`, `kmers` = the retained k-mers in output order (the spacing filter
  is deterministic, so the order is part of the contract), `pseudotax` =
  `none`, `empty`, or the tracked k-mers sorted by hash.
* individual mode emits one genome block per record, with the case id suffixed
  `#<record index>`.

The `pair_c5_k31_approx` case pins the approximate (cuckoo-filter) dedup path.
Its output depends on `scalable_cuckoo_filter` / `rand` internals, so a
dependency update that changes the filter's false-positive behaviour will
(legitimately) fail that block; the exact-dedup cases are stable.

## Host dependence

The goldens were generated with the x86_64 AVX2 seeding kernel, whose k-mer set
differs from the scalar kernel's. The `Seeding::Native` tests therefore skip
(with a message) unless the host is x86_64 with AVX2. The
`Seeding::Avx2Compat` test is portable: that kernel emulates the AVX2 kernel
exactly, so it must match the same goldens on any architecture.
