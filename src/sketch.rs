//! CLI-side sketching: argument parsing, the file walker and the path-based
//! wrappers around needletail.
//!
//! This module is always available so that the historical
//! `sylph::sketch::{extract_markers, pair_kmer, ...}` paths keep working, but
//! everything that needs needletail/clap/rayon/filesystem access is behind the
//! `cli` feature. The actual sketching logic lives in [`crate::sketchers`].

// The pure kernels live in the dependency-light `crate::kernels` module;
// re-export them here so that `sylph::sketch::{extract_markers, pair_kmer,
// dup_removal_lsh_full, ...}` keeps working for existing users.
pub use crate::kernels::{dup_removal_lsh_full_exact, pair_kmer, pair_kmer_single, Marker};
#[cfg(feature = "approx-dedup")]
pub use crate::kernels::dup_removal_lsh_full;
pub use crate::seeding::{extract_markers, extract_markers_positions, Seeding};
pub use crate::sketchers::{GenomeSketcher, PairedReadSketcher, ReadSketcher, SketchError};

#[cfg(feature = "cli")]
mod cli {
    use super::{GenomeSketcher, PairedReadSketcher, ReadSketcher, Seeding, SketchError};
    use crate::cmdline::*;
    use crate::constants::*;
    use crate::types::*;
    use log::*;
    use memory_stats::memory_stats;
    use needletail::parse_fastx_file;
    use rayon::prelude::*;
    use std::fs;
    use std::fs::File;
    use std::io::BufWriter;
    use std::io::{prelude::*, BufReader};
    use std::path::Path;
    use std::sync::Mutex;
    use std::thread;
    use std::time::Duration;

    pub fn check_vram_and_block(max_ram: usize, file: &str) {
        if let Some(usage) = memory_stats() {
            let mut gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
            if (max_ram as f64) < gb_usage_curr {
                log::debug!(
                    "Max memory reached. Blocking sketch for {}. Curr memory {}, max mem {}",
                    file,
                    gb_usage_curr,
                    max_ram
                );
            }
            while (max_ram as f64) < gb_usage_curr {
                let five_second = Duration::from_secs(1);
                thread::sleep(five_second);
                if let Some(usage) = memory_stats() {
                    gb_usage_curr = usage.virtual_mem as f64 / 1_000_000_000 as f64;
                    if (max_ram as f64) >= gb_usage_curr {
                        log::debug!("Sketching for {} freed", file);
                    }
                } else {
                    break;
                }
            }
        }
    }

    pub fn is_fastq(file: &str) -> bool {
        if file.ends_with(".fq")
            || file.ends_with(".fnq")
            || file.ends_with(".fastq")
            || file.ends_with(".fq.gz")
            || file.ends_with(".fnq.gz")
            || file.ends_with(".fastq.gz")
        {
            return true;
        } else {
            return false;
        }
    }

    pub fn is_fasta(file: &str) -> bool {
        if file.ends_with(".fa")
            || file.ends_with(".fna")
            || file.ends_with(".fasta")
            || file.ends_with(".fa.gz")
            || file.ends_with(".fna.gz")
            || file.ends_with(".fasta.gz")
        {
            return true;
        } else {
            return false;
        }
    }

    fn check_args_valid(args: &SketchArgs) {
        let level;
        if args.trace {
            level = log::LevelFilter::Trace;
        } else if args.debug {
            level = log::LevelFilter::Debug;
        } else {
            level = log::LevelFilter::Info;
        }

        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();

        simple_logger::SimpleLogger::new()
            .with_level(level)
            .init()
            .unwrap();

        if args.files.is_empty()
            && args.list_sequence.is_none()
            && args.first_pair.is_empty()
            && args.second_pair.is_empty()
            && args.genomes.is_none()
            && args.reads.is_none()
            && args.list_genomes.is_none()
            && args.list_reads.is_none()
            && args.list_first_pair.is_none()
            && args.list_second_pair.is_none()
        {
            error!("No input sequences found; see sylph sketch -h for help. Exiting.");
            std::process::exit(1);
        }

        if args.fpr < 0. || args.fpr >= 1. {
            error!("Invalid FPR for sketching. Must be in [0,1).");
            std::process::exit(1);
        }
    }

    fn parse_ambiguous_files(
        args: &SketchArgs,
        read_inputs: &mut Vec<String>,
        genome_inputs: &mut Vec<String>,
    ) {
        let mut all_files = vec![];
        if args.list_sequence.is_some() {
            let file_list = args.list_sequence.as_ref().unwrap();
            parse_line_file(file_list, &mut all_files);
        }

        all_files.extend(args.files.clone());

        for file in all_files {
            if is_fastq(&file) {
                read_inputs.push(file);
            } else if is_fasta(&file) {
                genome_inputs.push(file);
            } else {
                warn!(
                    "{} does not have a fasta/fastq/gzip type extension; skipping",
                    file
                );
            }
        }
    }

    fn parse_reads_and_genomes(
        args: &SketchArgs,
        read_inputs: &mut Vec<String>,
        genome_inputs: &mut Vec<String>,
    ) {
        if let Some(genomes_syl_in) = args.genomes.clone() {
            for gn_file in genomes_syl_in {
                genome_inputs.push(gn_file);
            }
        }
        if let Some(reads_syl_in) = args.reads.clone() {
            for rd_file in reads_syl_in {
                read_inputs.push(rd_file);
            }
        }

        if args.list_reads.is_some() {
            let file_reads = args.list_reads.as_ref().unwrap();
            parse_line_file(file_reads, read_inputs);
        }

        if args.list_genomes.is_some() {
            let file_genomes = args.list_genomes.as_ref().unwrap();
            parse_line_file(file_genomes, genome_inputs);
        }
    }

    fn parse_paired_end_reads(
        args: &SketchArgs,
        first_pairs: &mut Vec<String>,
        second_pairs: &mut Vec<String>,
    ) {
        if args.first_pair.len() != args.second_pair.len() {
            error!("Different number of paired sequences. Exiting.");
            std::process::exit(1);
        }

        for f in args.first_pair.iter() {
            first_pairs.push(f.clone());
        }

        for f in args.second_pair.iter() {
            second_pairs.push(f.clone());
        }

        if args.list_first_pair.is_some() {
            let file_first_pair = args.list_first_pair.as_ref().unwrap();
            parse_line_file(file_first_pair, first_pairs);
        }

        if args.list_second_pair.is_some() {
            let file_second_pair = args.list_second_pair.as_ref().unwrap();
            parse_line_file(file_second_pair, second_pairs)
        }

        if first_pairs.len() != second_pairs.len() {
            error!("Different number of paired sequences. Exiting.");
            std::process::exit(1);
        }
    }

    fn parse_line_file(file_name: &str, vec: &mut Vec<String>) {
        let file = File::open(file_name).unwrap();
        let reader = BufReader::new(file);
        for line in reader.lines() {
            vec.push(line.unwrap());
        }
    }

    fn parse_sample_names(args: &SketchArgs) -> Option<Vec<String>> {
        if args.list_sample_names.is_none() && args.sample_names.is_none() {
            return None;
        } else {
            let mut sample_names = vec![];
            if let Some(file) = &args.list_sample_names {
                parse_line_file(file, &mut sample_names);
                return Some(sample_names);
            }
            if let Some(vec) = &args.sample_names {
                sample_names.extend(vec.clone());
            }
            return Some(sample_names);
        }
    }

    pub fn sketch(args: SketchArgs) {
        let mut read_inputs = vec![];
        let mut genome_inputs = vec![];
        let mut first_pairs = vec![];
        let mut second_pairs = vec![];

        check_args_valid(&args);
        parse_ambiguous_files(&args, &mut read_inputs, &mut genome_inputs);
        parse_reads_and_genomes(&args, &mut read_inputs, &mut genome_inputs);
        parse_paired_end_reads(&args, &mut first_pairs, &mut second_pairs);

        let sample_names = parse_sample_names(&args);
        if let Some(names) = &sample_names {
            if names.len() != first_pairs.len() + read_inputs.len() {
                log::error!("Sample name length is not equal to the number of reads. Exiting");
                std::process::exit(1);
            }
        }

        let mut max_ram = usize::MAX;
        if args.max_ram.is_some() {
            max_ram = args.max_ram.unwrap();
            if max_ram < 7 {
                log::error!("Max ram must be >= 7. Exiting.");
                std::process::exit(1);
            }
        }

        if genome_inputs.is_empty() && args.db_out_name != "database" {
            log::warn!(
                "-o is set but no genomes are present. -o only applies to genomes; see -d for reads"
            );
        }

        if !first_pairs.is_empty() && !second_pairs.is_empty() {
            info!("Sketching paired sequences...");
            let pipeline_params = crate::parallel_sketch::PipelineParams {
                batch_records: args.sketch_batch_size,
                batch_max_bytes: args.sketch_batch_max_bytes,
                channel_depth: args.sketch_channel_depth,
                num_shards: args.sketch_shards,
            };
            crate::parallel_sketch::run_file_slots(first_pairs.len(), args.threads, args.sample_threads, |i, slot_threads| {
                let read_file1 = &first_pairs[i];
                let read_file2 = &second_pairs[i];

                let mut sample_name = None;
                if let Some(name) = &sample_names {
                    sample_name = Some(name[i].clone());
                }
                let read_sketch_opt = if crate::parallel_sketch::should_use_pipeline(
                    args.no_sketch_pipeline,
                    slot_threads,
                    read_file1,
                ) {
                    crate::parallel_sketch::sketch_pair_sequences_parallel(
                        read_file1,
                        read_file2,
                        args.c,
                        args.k,
                        sample_name.clone(),
                        args.no_dedup,
                        args.fpr,
                        slot_threads,
                        &pipeline_params,
                    )
                } else {
                    match sketch_pair_sequences(
                        read_file1,
                        read_file2,
                        args.c,
                        args.k,
                        sample_name.clone(),
                        args.no_dedup,
                        args.fpr,
                    ) {
                        Ok(sketch) => Some(sketch),
                        // The library already logged the problem; keep the CLI's
                        // historical hard failure for unreadable files.
                        Err(SketchError::Input(_)) => std::process::exit(1),
                        Err(_) => None,
                    }
                };
                if read_sketch_opt.is_some() {
                    let res = fs::create_dir_all(&args.sample_output_dir);
                    if res.is_err() {
                        error!("Could not create directory at {}", args.sample_output_dir);
                        std::process::exit(1);
                    }
                    let pref = Path::new(&args.sample_output_dir);
                    let read_sketch = read_sketch_opt.unwrap();

                    let sketch_name;
                    if sample_name.is_some() {
                        sketch_name = read_sketch.sample_name.as_ref().unwrap();
                    } else {
                        sketch_name = &read_sketch.file_name;
                    }

                    let read_file_path = Path::new(&sketch_name).file_name().unwrap();
                    let file_path = pref.join(&read_file_path);

                    let file_path_str = format!(
                        "{}.paired{}",
                        file_path.to_str().unwrap(),
                        SAMPLE_FILE_SUFFIX
                    );

                    let mut read_sk_file = BufWriter::new(
                        File::create(&file_path_str)
                            .expect(&format!("{} path not valid; exiting ", file_path_str)),
                    );

                    bincode::serialize_into(&mut read_sk_file, &read_sketch).unwrap();
                    info!("Sketching {} complete.", file_path_str);
                }
            });
        }

        if !read_inputs.is_empty() {
            info!("Sketching non-paired sequences...");
        }

        let pipeline_params = crate::parallel_sketch::PipelineParams {
            batch_records: args.sketch_batch_size,
            batch_max_bytes: args.sketch_batch_max_bytes,
            channel_depth: args.sketch_channel_depth,
            num_shards: args.sketch_shards,
        };
        crate::parallel_sketch::run_file_slots(read_inputs.len(), args.threads, args.sample_threads, |i, slot_threads| {
            let pref = Path::new(&args.sample_output_dir);
            std::fs::create_dir_all(pref)
                .expect("Could not create directory for output sample files (-d). Exiting...");

            let read_file = &read_inputs[i];

            check_vram_and_block(max_ram, read_file);
            let mut sample_name = None;
            if let Some(name) = &sample_names {
                sample_name = Some(name[i + first_pairs.len()].clone());
            }

            let read_sketch_opt = if crate::parallel_sketch::should_use_pipeline(
                args.no_sketch_pipeline,
                slot_threads,
                read_file,
            ) {
                crate::parallel_sketch::sketch_sequences_needle_parallel(
                    read_file,
                    args.c,
                    args.k,
                    sample_name.clone(),
                    args.no_dedup,
                    slot_threads,
                    &pipeline_params,
                )
            } else {
                sketch_sequences_needle(read_file, args.c, args.k, sample_name.clone(), args.no_dedup)
            };

            if read_sketch_opt.is_some() {
                let read_sketch = read_sketch_opt.unwrap();
                let sketch_name;
                if sample_name.is_some() {
                    sketch_name = read_sketch.sample_name.as_ref().unwrap();
                } else {
                    sketch_name = &read_sketch.file_name;
                }
                let read_file_path = Path::new(&sketch_name).file_name().unwrap();
                let file_path = pref.join(&read_file_path);

                let file_path_str = format!("{}{}", file_path.to_str().unwrap(), SAMPLE_FILE_SUFFIX);

                let mut read_sk_file = BufWriter::new(
                    File::create(&file_path_str)
                        .expect(&format!("{} path not valid; exiting.", file_path_str)),
                );

                bincode::serialize_into(&mut read_sk_file, &read_sketch).unwrap();
                info!("Sketching {} complete.", file_path_str);
            }
        });

        if !genome_inputs.is_empty() {
            info!("Sketching genomes...");
            let iter_vec: Vec<usize> = (0..genome_inputs.len()).into_iter().collect();
            let counter: Mutex<usize> = Mutex::new(0);
            let pref = Path::new(&args.db_out_name);
            let file_path_str = format!("{}{}", pref.to_str().unwrap(), QUERY_FILE_SUFFIX);
            let path = std::path::Path::new(&file_path_str);
            let prefix = path.parent().unwrap();
            std::fs::create_dir_all(prefix)
                .expect("Could not create directory for output database file (-o). Exiting...");
            let all_genome_sketches = Mutex::new(vec![]);

            iter_vec.into_par_iter().for_each(|i| {
                let genome_file = &genome_inputs[i];
                if args.individual {
                    let indiv_gn_sketches = sketch_genome_individual(
                        args.c,
                        args.k,
                        genome_file,
                        args.min_spacing_kmer,
                        !args.no_pseudotax,
                    );
                    all_genome_sketches
                        .lock()
                        .unwrap()
                        .extend(indiv_gn_sketches);
                } else {
                    let genome_sketch = sketch_genome(
                        args.c,
                        args.k,
                        genome_file,
                        args.min_spacing_kmer,
                        !args.no_pseudotax,
                    );
                    if genome_sketch.is_some() {
                        all_genome_sketches
                            .lock()
                            .unwrap()
                            .push(genome_sketch.unwrap());
                    }
                }
                let mut c = counter.lock().unwrap();
                *c += 1;
                if *c % 100 == 0 && *c != 0 {
                    info!("{} genomes processed.", *c);
                }
            });

            if all_genome_sketches.lock().unwrap().is_empty() {
                warn!(
                    "No valid genomes to sketch; {} is not output",
                    file_path_str
                );
            } else {
                let mut genome_sk_file = BufWriter::new(
                    File::create(&file_path_str).expect(&format!("{} not valid ", file_path_str)),
                );
                info!("Wrote all genome sketches to {}", file_path_str);
                bincode::serialize_into(&mut genome_sk_file, &all_genome_sketches).unwrap();
            }
        }

        info!("Finished.");
    }

    pub fn sketch_genome_individual(
        c: usize,
        k: usize,
        ref_file: &str,
        min_spacing: usize,
        pseudotax: bool,
    ) -> Vec<GenomeSketch> {
        let reader = parse_fastx_file(&ref_file);
        if !reader.is_ok() {
            warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
            return vec![];
        }
        let mut reader = reader.unwrap();
        let mut return_vec = vec![];
        while let Some(record) = reader.next() {
            if record.is_ok() {
                let record = record.expect(&format!("Invalid record for file {} ", ref_file));
                let mut sketcher =
                    GenomeSketcher::new(c, k, min_spacing, pseudotax, Seeding::Native);
                sketcher.push_contig(record.id(), &record.seq());
                return_vec.push(sketcher.finish(ref_file.to_string()));
            } else {
                warn!("File {} is not a valid fasta/fastq file", ref_file);
                return vec![];
            }
        }
        return return_vec;
    }

    pub fn sketch_genome(
        c: usize,
        k: usize,
        ref_file: &str,
        min_spacing: usize,
        pseudotax: bool,
    ) -> Option<GenomeSketch> {
        let reader = parse_fastx_file(&ref_file);
        if !reader.is_ok() {
            warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
            return None;
        }
        let mut reader = reader.unwrap();
        let mut sketcher = GenomeSketcher::new(c, k, min_spacing, pseudotax, Seeding::Native);
        while let Some(record) = reader.next() {
            if record.is_ok() {
                let record = record.expect(&format!("Invalid record for file {} ", ref_file));
                sketcher.push_contig(record.id(), &record.seq());
            } else {
                warn!("File {} is not a valid fasta/fastq file", ref_file);
                return None;
            }
        }
        Some(sketcher.finish(ref_file.to_string()))
    }

    pub fn sketch_pair_sequences(
        read_file1: &str,
        read_file2: &str,
        c: usize,
        k: usize,
        sample_name: Option<String>,
        no_dedup: bool,
        dedup_fpr: f64,
    ) -> Result<SequencesSketch, SketchError> {
        let r1o = parse_fastx_file(&read_file1);
        let r2o = parse_fastx_file(&read_file2);
        if r1o.is_err() || r2o.is_err() {
            log::error!("Paired end reading failed for '{}' and '{}'. Make sure the files are present or the sequences are valid.", read_file1, read_file2);
            return Err(SketchError::Input(format!(
                "Paired end reading failed for '{}' and '{}'. Make sure the files are present or the sequences are valid.",
                read_file1, read_file2
            )));
        }

        let mut reader1 = r1o.unwrap();
        let mut reader2 = r2o.unwrap();
        let mut sketcher =
            PairedReadSketcher::new(c, k, no_dedup, dedup_fpr, Seeding::Native)?;

        loop {
            let n1 = reader1.next();
            let n2 = reader2.next();
            if let Some(rec1_o) = n1 {
                if let Some(rec2_o) = n2 {
                    if let Ok(rec1) = rec1_o {
                        if let Ok(rec2) = rec2_o {
                            sketcher.push_pair(&rec1.seq(), &rec2.seq());
                        }
                    } else {
                        return Err(SketchError::Record(format!(
                            "Invalid record for file {} ",
                            read_file1
                        )));
                    }
                }
            } else {
                break;
            }
        }
        Ok(sketcher.finish(read_file1.to_string(), sample_name))
    }

    pub fn sketch_sequences_needle(
        read_file: &str,
        c: usize,
        k: usize,
        sample_name: Option<String>,
        no_dedup: bool,
    ) -> Option<SequencesSketch> {
        let ref_file = &read_file;
        let reader = parse_fastx_file(&ref_file);
        if !reader.is_ok() {
            warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
            return None;
        }
        let mut reader = reader.unwrap();
        let mut sketcher = ReadSketcher::new(c, k, no_dedup, Seeding::Native);
        while let Some(record) = reader.next() {
            if record.is_ok() {
                let record = record.expect(&format!("Invalid record for file {} ", ref_file));
                sketcher.push_read(&record.seq());
            } else {
                warn!("File {} is not a valid fasta/fastq file", ref_file);
            }
        }
        Some(sketcher.finish(read_file.to_string(), sample_name))
    }
}

#[cfg(feature = "cli")]
pub use cli::*;
