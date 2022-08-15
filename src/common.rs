/// Commonly shared data structures
use clap::Args as ClapArgs;
use clap_verbosity_flag::Verbosity;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Verbosity of the program
    #[clap(flatten)]
    pub verbose: Verbosity,
    /// The seed to use for the random number generator
    #[clap(short = 's', long = "seed", default_value_t = 42)]
    pub seed: u64,
}
