/// Sequencing simulation
pub mod illumina;
pub mod roche454;
pub mod sanger;

use clap::Args as ClapArgs;

#[derive(ClapArgs, Debug)]
pub struct Args {
}
