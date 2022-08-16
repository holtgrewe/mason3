/// Simulation of Roche 454 reads
use clap::builder::ArgAction as ClapArgAction;
use clap::Args as ClapArgs;

use crate::seq::ReadLengthModel;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// The model to use for read length distribution
    #[clap(id = "454-read-length-model", long, value_enum, default_value_t = ReadLengthModel::Uniform)]
    length_model: ReadLengthModel,
    /// The minimal read length if uniformly distributed
    #[clap(id = "454-read-length-min", long, default_value_t = 10)]
    read_length_min: usize,
    /// The maximal read length if uniformly distributed
    #[clap(id = "454-read-length-max", long, default_value_t = 600)]
    read_length_max: usize,
    /// The read length mean if normally distributed
    #[clap(id = "454-read-length-mean", long, default_value_t = 400.0)]
    read_length_mean: f64,
    /// The read length standard deviation if normally distributed
    #[clap(id = "454-read-length-stddev", long, default_value_t = 40.0)]
    read_length_stddev: f64,
    /// If true, then `sigma=k*sqrt(r)`, otherwise `sigma=k*r` is used
    #[clap(id = "454-no-sqrt-in-std-dev", long, action = ClapArgAction::SetFalse, default_value_t = true)]
    model_sqrt_in_stddev: bool,
    /// Proprtionality factor for calcuating standard deviation proporitional to
    /// `sqrt(homopolymer length)`
    #[clap(id = "454-proportionality-factor", long, default_value_t = 0.15)]
    model_k: f64,
    /// The mean of the lognormal distribution for the noise
    #[clap(id = "454-background-noise-mean", long, default_value_t = 0.23)]
    model_bg_noise_mean: f64,
    /// The standard deviation of the lognromal distribution for the noise
    #[clap(id = "454-background-noise-stddev", long, default_value_t = 0.15)]
    model_bg_noise_stddev: f64,
}
