/// Simulation of Sanger/capillary sequencing reads
use clap::Args as ClapArgs;

use crate::seq::Args as SeqArgs;
use crate::seq::ReadFromFragment;
use crate::seq::ReadLengthModel;

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
        _seq: &mut Vec<u8>,
        _quals: &mut Vec<u8>,
        _info: &mut super::SeqInfo,
        _frag: &[u8],
        _dir: super::Direction,
        _strand: super::Strand,
    ) {
        todo!()
    }
}
