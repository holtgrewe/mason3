/// `mason frag_sequencing` -- simulate reads from fragments
use fastrand::Rng;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use clap::Args as ClapArgs;
use console::{Emoji, Term};
use indicatif::{ProgressBar, ProgressStyle};

use crate::common::Args as CommonArgs;
use crate::common::prefix_lines;
use crate::seq::Args as SequencingArgs;
use crate::seq::illumina::Args as IlluminaArgs;
use crate::seq::roche454::Args as Roche454Args;
use crate::seq::sanger::Args as SangerArgs;

/// Configuration for `methylation` sub command
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Input file name
    #[clap(short = 'i', long = "in")]
    input_filename: String,
    /// Output file name for left reads
    #[clap(short = 'o', long = "out")]
    output_filename_left: String,
    /// Output file name for right reads
    #[clap(short = 'o', long = "out-right")]
    output_filename_right: String,

    /// Common sequencing simulation options
    #[clap(flatten)]
    sequencing: SequencingArgs,
    /// Illumina sequencing simulation options
    #[clap(flatten)]
    illumina: IlluminaArgs,
    /// Roche 454 sequencing simulation options
    #[clap(flatten)]
    roche454: Roche454Args,
    /// Sanger sequencing simulation options
    #[clap(flatten)]
    sanger: SangerArgs,
}

pub fn run(term: &Term, common_args: &CommonArgs, args: &Args) -> Result<(), anyhow::Error> {
    term.write_line(&format!(
        "{}mason frag_sequencing -- simulate sequencing of fragments",
        Emoji("🧬 ", "")
    ))?;
    term.write_line(&format!("{}configuration:", Emoji("🔧 ", "")))?;
    term.write_line(&prefix_lines(
        "   ",
        &format!("common = {:#?}", &common_args),
    ))?;
    term.write_line(&prefix_lines("   ", &format!("sequencing = {:#?}", &args)))?;
    // simulate_levels(term, common_args, args)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
}