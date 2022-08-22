/// Simulation of Illumina short reads
use clap::Args as ClapArgs;
use rand_xoshiro::Xoshiro256Plus;

use crate::seq::Args as SeqArgs;
use crate::seq::ReadFromFragment;

#[derive(ClapArgs, Debug, Clone)]
pub struct Args {
    /// Length of the reads to simulate
    #[clap(id = "illumina-read-length", long, default_value_t = 100)]
    read_length: usize,
    /// Optional path to file with positional error probabilities (floating point numbers
    /// separated by space, each giving a positional error rate).
    #[clap(id = "illumina-error-profile-file")]
    probability_mismatch_file: Option<String>,
    /// Probability for an insertion (indepenent of position)
    #[clap(id = "illumina-prob-insert", long, default_value_t = 0.00005)]
    err_prob_ins: f64,
    /// Probability for a deleteion (independent of position)
    #[clap(id = "illumina-prob-deletion", long, default_value_t = 0.00005)]
    err_prob_del: f64,
    /// Scale factor toa pply to mismatch probabilities
    #[clap(id = "illumina-prob-mismatch-scale", long, default_value_t = 1.0)]
    err_prob_mismatch_scale: f64,
    /// Probability of a mismatch (single-nucleotide exchange)
    #[clap(id = "illumina-prob-mismatch", long, default_value_t = 1.0)]
    err_prob_mismatch: f64,
    /// Probability of a mismatch in the first base
    #[clap(id = "illumina-prob-mismatch-begin", long, default_value_t = 0.002)]
    err_prob_mismatch_begin: f64,
    /// Probability of a mismatch in the last base
    #[clap(id = "illumina-prob-mismatch-end", long, default_value_t = 0.012)]
    err_prob_mismatch_end: f64,
    /// Relative position in the read betwen 0 and 1 where the steeper curve begins
    #[clap(id = "illumina-position-raise", long, default_value_t = 0.66)]
    position_raise: f64,
    /// Optional path to template FASTQ files to use for positional qualities for the first
    /// mate.  Patterns of Ns will be applied to the simulated reads.  If set, this will be
    /// used instead of the built-in model.
    #[clap(id = "illumina-left-template-fastq")]
    template_fastq_left: Option<String>,
    /// Optional path to template FASTQ files for the left mate, see documentation above.
    #[clap(id = "illumina-right-template-fastq")]
    template_fastq_right: Option<String>,
    /// Mean quality for non-mismatches at the first base
    #[clap(id = "illumina-quality-mean-begin", long, default_value_t = 40.0)]
    quality_mean_begin: f64,
    /// Mean quality for non-mismatches at the last base
    #[clap(id = "illumina-quality-mean-end", long, default_value_t = 39.5)]
    quality_mean_end: f64,
    /// Standard deviation quality for non-mismatches at the first base
    #[clap(id = "illumina-quality-stddev-begin", long, default_value_t = 0.005)]
    quality_stddev_begin: f64,
    /// Standard deviation quality for non-mismatches at the last base
    #[clap(id = "illumina-quality-stddev-end", long, default_value_t = 10.0)]
    quality_stddev_end: f64,
    /// Mean quality for mismatches at the first base
    #[clap(
        id = "illumina-mismatch-quality-mean-begin",
        long,
        default_value_t = 40.0
    )]
    quality_mismatch_mean_begin: f64,
    /// Mean quality for mismatches at the last base
    #[clap(
        id = "illumina-mismatch-quality-mean-end",
        long,
        default_value_t = 30.0
    )]
    quality_mismatch_mean_end: f64,
    /// Standard deviation quality for mismatches at the first base
    #[clap(
        id = "illumina-mismatch-quality-stddev-begin",
        long,
        default_value_t = 3.0
    )]
    quality_mismatch_stddev_begin: f64,
    /// Standard deviation quality for mismatches at the lastbase
    #[clap(
        id = "illumina-mismatch-quality-stddev-end",
        long,
        default_value_t = 15.0
    )]
    quality_mismatch_stddev_end: f64,
}

impl Args {
    #[allow(dead_code)] // only used in tests
    pub fn new() -> Self {
        Self {
            read_length: 100,
            probability_mismatch_file: None,
            err_prob_ins: 0.00005,
            err_prob_del: 0.00005,
            err_prob_mismatch_scale: 1.0,
            err_prob_mismatch: 1.0,
            err_prob_mismatch_begin: 0.002,
            err_prob_mismatch_end: 0.012,
            position_raise: 0.66,
            template_fastq_left: None,
            template_fastq_right: None,
            quality_mean_begin: 40.0,
            quality_mean_end: 39.5,
            quality_stddev_begin: 0.005,
            quality_stddev_end: 10.0,
            quality_mismatch_mean_begin: 40.0,
            quality_mismatch_mean_end: 30.0,
            quality_mismatch_stddev_begin: 3.0,
            quality_mismatch_stddev_end: 15.0,
        }
    }
}

/// The `ReadFromFragment` implementation for Illumina reads
pub struct IlluminaFromFragment {
    /// Generic sequencing simulation options
    seq_args: SeqArgs,
    /// Illumina-specific simulation options
    args: Args,
}

impl IlluminaFromFragment {
    pub fn new(seq_args: &SeqArgs, args: &Args) -> Self {
        Self {
            seq_args: seq_args.clone(),
            args: args.clone(),
        }
    }
}

impl ReadFromFragment for IlluminaFromFragment {
    fn simulate_read(
        &self,
        _seq: &mut Vec<u8>,
        _quals: &mut Vec<u8>,
        _info: &mut super::SeqInfo,
        _frag: &[u8],
        _dir: super::Direction,
        _strand: super::Strand,
        _rng: &mut Xoshiro256Plus,
    ) {
        todo!()
    }
}
