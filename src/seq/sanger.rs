use bio::alphabets::dna::revcomp;
/// Simulation of Sanger/capillary sequencing reads
use clap::Args as ClapArgs;
use rand::distributions::Uniform;
use rand::prelude::Distribution;
use rand_xoshiro::Xoshiro256Plus;
use statrs::distribution::Normal;

use crate::seq::Args as SeqArgs;
use crate::seq::CigarElement;
use crate::seq::ReadFromFragment;
use crate::seq::ReadLengthModel;

use super::{CigarOp, Direction, Strand};

#[derive(ClapArgs, Debug, Clone)]
pub struct Args {
    /// The model to use for read length distribution
    #[clap(id =  "sanger-read-length-model", long, value_enum, default_value_t = ReadLengthModel::Normal)]
    length_model: ReadLengthModel,
    /// The minimal read length if uniformly distributed
    #[clap(id = "sanger-read-length-min", long, default_value_t = 400)]
    read_length_min: usize,
    /// The maximal read length if uniformly distributed
    #[clap(id = "sanger-read-length-max", long, default_value_t = 600)]
    read_length_max: usize,
    /// The read length mean if normally distributed
    #[clap(id = "sanger-read-length-mean", long, default_value_t = 400.0)]
    read_length_mean: f64,
    /// The read length standard deviation if normally distributed
    #[clap(id = "sanger-read-length-stddev", long, default_value_t = 40.0)]
    read_length_stddev: f64,
    /// Mismatch probability ramp start
    #[clap(id = "sanger-prob-mismatch-begin", long, default_value_t = 0.005)]
    prob_mismatch_begin: f64,
    /// Mismatch probability ramp end
    #[clap(id = "sanger-prob-mismatch-end", long, default_value_t = 0.01)]
    prob_mismatch_end: f64,
    /// Insert probability ramp start
    #[clap(id = "sanger-prob-insertion-begin", long, default_value_t = 0.001)]
    prob_ins_begin: f64,
    /// Insert probability ramp end
    #[clap(id = "sanger-prob-insertion-end", long, default_value_t = 0.0025)]
    prob_ins_end: f64,
    /// Insert probability ramp start
    #[clap(id = "sanger-prob-deletion-begin", long, default_value_t = 0.001)]
    prob_del_begin: f64,
    /// Insert probability ramp end
    #[clap(id = "sanger-prob-deletion-end", long, default_value_t = 0.0025)]
    prob_del_end: f64,
    /// Quality mean at start for matches
    #[clap(id = "sanger-quality-match-start-mean", long, default_value_t = 40.0)]
    qual_match_start_mean: f64,
    /// Quality mean at end for matches
    #[clap(id = "sanger-quality-match-end-mean", long, default_value_t = 39.5)]
    qual_match_end_mean: f64,
    /// Quality standard deviation at start for matches
    #[clap(id = "sanger-quality-match-start-stddev", long, default_value_t = 0.1)]
    qual_match_start_stddev: f64,
    /// Quality standard deviation at end for matches
    #[clap(id = "sanger-quality-match-end-stddev", long, default_value_t = 2.0)]
    qual_match_end_stddev: f64,
    /// Quality mean at start for errors
    #[clap(id = "sanger-quality-error-start-mean", long, default_value_t = 30.0)]
    qual_err_start_mean: f64,
    /// Quality mean at end for errors
    #[clap(id = "sanger-quality-error-end-mean", long, default_value_t = 20.0)]
    qual_err_end_mean: f64,
    /// Quality standard deviation at start for errors
    #[clap(id = "sanger-quality-error-start-stddev", long, default_value_t = 2.0)]
    qual_err_start_stddev: f64,
    /// Quality standard deviation at end for errors
    #[clap(id = "sanger-quality-error-end-stddev", long, default_value_t = 5.0)]
    qual_err_end_stddev: f64,
}

impl Args {
    #[allow(dead_code)] // only used in tests
    pub fn new() -> Self {
        Self {
            length_model: ReadLengthModel::Normal,
            read_length_min: 400,
            read_length_max: 600,
            read_length_mean: 400.0,
            read_length_stddev: 40.0,
            prob_mismatch_begin: 0.005,
            prob_mismatch_end: 0.01,
            prob_ins_begin: 0.001,
            prob_ins_end: 0.0025,
            prob_del_begin: 0.001,
            prob_del_end: 0.0025,
            qual_match_start_mean: 40.0,
            qual_match_end_mean: 39.5,
            qual_match_start_stddev: 0.1,
            qual_match_end_stddev: 2.0,
            qual_err_start_mean: 30.0,
            qual_err_end_mean: 20.0,
            qual_err_start_stddev: 2.0,
            qual_err_end_stddev: 5.0,
        }
    }
}

/// The `ReadFromFragment` implementation for Sanger reads
pub struct SangerFromFragment {
    /// Generic sequencing simulation options
    seq_args: SeqArgs,
    /// Sanger-specific simulation options
    args: Args,
}

impl SangerFromFragment {
    pub fn new(seq_args: &SeqArgs, args: &Args) -> Self {
        Self {
            seq_args: seq_args.clone(),
            args: args.clone(),
        }
    }
}

impl ReadFromFragment for SangerFromFragment {
    fn simulate_read(
        &self,
        seq: &mut Vec<u8>,
        quals: &mut Vec<u8>,
        info: &mut super::SeqInfo,
        frag: &[u8],
        dir: super::Direction,
        strand: super::Strand,
        rng: &mut Xoshiro256Plus,
    ) {
        info.clear();

        // pick read length
        let sample_length = self.pick_read_length(rng);
        if sample_length > frag.len() {
            panic!(
                "Sanger read is too long ({}) for fragment ({})",
                sample_length,
                frag.len()
            );
        }

        // simulate CIGAR string
        self.simulate_cigar(&mut info.cigar, rng, sample_length);

        // simulate sequence (materialize mismatches and insertions)
        self.materialize_cigar(seq, rng, frag, sample_length, &info.cigar, dir, strand);

        // simulate qualities
        self.simulate_qualities(quals, &info.cigar, rng, sample_length);

        // reverse qualities if necessary
        if strand == Strand::Reverse {
            quals.reverse();
        }

        // write out sequencing information if instructed to do so
        if self.seq_args.embed_read_info {
            let len_in_ref = info.cigar.len_in_ref();
            info.orig_seq = match dir {
                Direction::Left => {
                    if strand == Strand::Reverse {
                        revcomp(&frag[0..len_in_ref])
                    } else {
                        frag[0..len_in_ref].into()
                    }
                }
                Direction::Right => {
                    if strand == Strand::Reverse {
                        revcomp(&frag[frag.len() - len_in_ref..frag.len()])
                    } else {
                        frag[frag.len() - len_in_ref..frag.len()].into()
                    }
                }
            };
            info.strand = strand;
        }
    }
}

// Supporting functions for `impl ReadFromFragment for SangerFromFragment`
impl SangerFromFragment {
    fn pick_read_length(&self, rng: &mut Xoshiro256Plus) -> usize {
        match self.args.length_model {
            ReadLengthModel::Uniform => {
                let dist = Uniform::new(
                    self.args.read_length_min as f64,
                    self.args.read_length_max as f64,
                );
                dist.sample(rng).round() as usize
            }
            ReadLengthModel::Normal => {
                let dist = Normal::new(
                    self.args.read_length_mean as f64,
                    self.args.read_length_stddev as f64,
                )
                .expect("invalid normal distribution");
                dist.sample(rng).round() as usize
            }
        }
    }

    /// Simulate CIGAR string.  We can do this with position specific parameters only andn thus
    /// independent of the context.
    fn simulate_cigar(
        &self,
        cigar: &mut super::CigarString,
        rng: &mut Xoshiro256Plus,
        sample_length: usize,
    ) {
        let mut i = 0;
        let dist = Uniform::new(0.0, 1.0);
        while i < sample_length {
            let x = dist.sample(rng);
            let pos = (i as f64) / ((sample_length as f64) - 1.0);
            let p_mismatch = self.args.prob_mismatch_begin
                + pos * (self.args.prob_mismatch_end - self.args.prob_mismatch_begin);
            let p_ins =
                self.args.prob_ins_begin + pos * (self.args.prob_ins_end - self.args.prob_ins_begin);
            let p_del =
                self.args.prob_del_begin + pos * (self.args.prob_del_end - self.args.prob_del_begin);
            let p_match = 1.0 - p_mismatch - p_ins - p_del;

            if x < p_match {
                // match
                let j = cigar.append(CigarElement::new(CigarOp::Match, 1)).1;
                let k = j as usize;
                i += k;
            } else if x < p_match + p_mismatch {
                // single nucleotide mismatch
                i += cigar.append(CigarElement::new(CigarOp::Diff, 1)).1 as usize;
            } else if x < p_match + p_mismatch + p_ins {
                // insertion
                if i > 0 && (i + 1) < sample_length {
                    // no indel at beginning / end
                    i += cigar.append(CigarElement::new(CigarOp::Ins, 1)).1 as usize;
                }
            } else {
                // deletion
                if i > 0 && (i + 1) < sample_length {
                    // no indel at beginning / end
                    i += cigar.append(CigarElement::new(CigarOp::Del, 1)).1 as usize;
                }
            }
        }
    }

    /// Materialize the CIGAR string in `cigar` and write results to `seq`.
    fn materialize_cigar(
        &self,
        seq: &mut Vec<u8>,
        rng: &mut Xoshiro256Plus,
        frag: &[u8],
        sample_length: usize,
        cigar: &super::CigarString,
        dir: super::Direction,
        strand: super::Strand,
    ) {
        // Note: can be re-used for Illumina
        let mut frag_segment = Vec::with_capacity(sample_length);
        if dir == Direction::Left {
            // Take from the left
            frag_segment.extend_from_slice(&frag[0..sample_length]);
        } else {
            // Take from the right
            let end = frag.len();
            let begin = end - sample_length;
            frag_segment.extend_from_slice(&frag[begin..end]);
        }
        if strand == Strand::Reverse {
            // Reverse-complement because we are on the reverse strand
            frag_segment = bio::alphabets::dna::revcomp(&frag_segment);
        }

        self.simulate_sequence(seq, rng, &frag_segment, cigar);
    }

    /// Materialize CIGAR into `seq` from `frag_segment` (after reverse-complementing when necessary)
    fn simulate_sequence(
        &self,
        seq: &mut Vec<u8>,
        rng: &mut Xoshiro256Plus,
        frag_segment: &[u8],
        cigar: &super::CigarString,
    ) {
        // Note: this function can be re-used for Illumina

        // Position in `frag_segment`
        let mut pos = 0;

        let chars = b"ACGT";
        let uniform = Uniform::from(0..chars.len());

        for CigarElement { op, count } in &cigar.elements {
            match op {
                CigarOp::Match => {
                    (0..*count).for_each(|_j| {
                        seq.push(frag_segment[pos]);
                        pos += 1;
                    });
                }
                CigarOp::Del => pos += count,
                CigarOp::Ins => {
                    (0..*count).for_each(|_j| {
                        seq.push(chars[uniform.sample(rng)]);
                    });
                }
                CigarOp::Diff => {
                    (0..*count).for_each(|_j| {
                        let x = uniform.sample(rng);
                        let xu = b"ACGT"[x];
                        let xl = b"acgt"[x];
                        let c = frag_segment[pos];
                        let offset = (c == xu) || (c == xl);
                        seq.push(b"ACGTN"[x + offset as usize]);
                        pos += 1;
                    });
                }
                _ => panic!("CIGAR operation not generated by mason"),
            }
        }
    }

    /// Simulate the qualities for the read
    fn simulate_qualities(
        &self,
        quals: &mut Vec<u8>,
        cigar: &super::CigarString,
        rng: &mut Xoshiro256Plus,
        sample_length: usize,
    ) {
        let mut pos = 0; // in fragment

        for CigarElement { op, count } in &cigar.elements {
            (0..*count).for_each(|_j| {
                if *op != CigarOp::Del {
                    pos += 1;
                }

                let rel_pos = 1.0 * (pos as f64) / (sample_length as f64);

                let mean_stddev = match *op {
                    CigarOp::Del => {
                        // no quality to give out
                        None
                    }
                    CigarOp::Ins | CigarOp::Diff => {
                        // quality for insertion/mismatch
                        Some((
                            self.args.qual_err_start_mean
                                + rel_pos
                                    * (self.args.qual_err_end_mean - self.args.qual_err_start_mean),
                            self.args.qual_err_start_stddev
                                + rel_pos
                                    * (self.args.qual_err_end_stddev
                                        - self.args.qual_err_start_stddev),
                        ))
                    }
                    CigarOp::Match => {
                        // quality for match
                        Some((
                            self.args.qual_match_start_mean
                                + rel_pos
                                    * (self.args.qual_match_end_mean
                                        - self.args.qual_match_start_mean),
                            self.args.qual_match_start_stddev
                                + rel_pos
                                    * (self.args.qual_match_end_stddev
                                        - self.args.qual_match_start_stddev),
                        ))
                    }
                    _ => panic!("CIGAR operation not generated by mason"),
                };
                if let Some((mean, std_dev)) = mean_stddev {
                    let dist = Normal::new(mean, std_dev).expect("invalid normal distribution");
                    let val = dist.sample(rng);
                    let phred = (-10.0 * val.log10()).round() as i64;
                    let phred = std::cmp::max(0, std::cmp::min(phred, 40)) as u8;
                    quals.push(b'!' + phred);
                }
            })
        }
    }
}
