/// Main entry point for mason
mod cmd_genome;
mod cmd_methylation;
mod common;
mod err;
mod methylation;

use clap::{Parser, Subcommand};
use console::{Emoji, Term};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    /// Commonly used arguments
    #[clap(flatten)]
    common: common::Args,

    /// The sub command to run
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Create contigs with synthetic sequence
    Genome(cmd_genome::Args),
    /// Simulate methylation levels
    Methylation(cmd_methylation::Args),
}

fn main() -> Result<(), anyhow::Error> {
    let cli = Cli::parse();

    let term = Term::stderr();
    match &cli.command {
        Commands::Genome(args) => {
            cmd_genome::run(&term, &cli.common, args)?;
        }
        Commands::Methylation(args) => {
            cmd_methylation::run(&term, &cli.common, args)?;
        }
    }
    term.write_line(&format!("All done. Have a nice day!{}", Emoji(" 😃", "")))?;

    Ok(())
}
