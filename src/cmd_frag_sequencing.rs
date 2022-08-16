/// `mason frag-sequencing` -- simulate reads from fragments
use bio::io::fasta;
use clap::Args as ClapArgs;
use console::{Emoji, Term};
use fastrand::Rng;

use crate::common::prefix_lines;
use crate::common::Args as CommonArgs;
use crate::seq::illumina::Args as IlluminaArgs;
use crate::seq::illumina::IlluminaFromFragment;
use crate::seq::roche454::Args as Roche454Args;
use crate::seq::roche454::Roche454FromFragment;
use crate::seq::sanger::Args as SangerArgs;
use crate::seq::sanger::SangerFromFragment;
use crate::seq::Args as SequencingArgs;
use crate::seq::ReadFromFragment;
use crate::seq::ReadSimulator;

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
    #[clap(long = "out-right")]
    output_filename_right: Option<String>,
    /// Force single-end sequencing although `--out-right` is given
    #[clap(long = "force-single-end", action, default_value_t = false)]
    force_single_end: bool,

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

/// Run the read simulation from the fragments
fn run_simulation(term: &Term, common_args: &CommonArgs, args: &Args) -> Result<(), anyhow::Error> {
    let reader = fasta::Reader::from_file(&args.input_filename)?;
    let writer_left = fasta::Writer::to_file(&args.output_filename_left);
    let writer_right = if let Some(path) = &args.output_filename_right {
        if !args.force_single_end {
            Some(fasta::Writer::to_file(&path)?)
        } else {
            None
        }
    } else {
        None
    };
    let mut rng = Rng::with_seed(common_args.seed);
    let mut rng_meth = Rng::with_seed(common_args.seed);

    let simulator = ReadSimulator::new(&args.sequencing, &mut rng, &mut rng_meth);
    let read_gen: Box<dyn ReadFromFragment> = match args.sequencing.technology {
        crate::seq::SequencingTechnology::Illumina => {
            Box::new(IlluminaFromFragment::new(&args.sequencing, &args.illumina))
        }
        crate::seq::SequencingTechnology::Roche454 => {
            Box::new(Roche454FromFragment::new(&args.sequencing, &args.roche454))
        }
        crate::seq::SequencingTechnology::Sanger => {
            Box::new(SangerFromFragment::new(&args.sequencing, &args.sanger))
        }
    };

    for (i, frag_record) in reader.records().enumerate() {
        let frag_record = frag_record?;
        let frag_id = frag_record
            .id()
            .split_whitespace()
            .next()
            .expect("fragment must have an ID");
        if args.sequencing.embed_read_info {
        } else {
        };
    }

    Ok(())
}

/// Main entry point for `mason frag-sequencing`
pub fn run(term: &Term, common_args: &CommonArgs, args: &Args) -> Result<(), anyhow::Error> {
    term.write_line(&format!(
        "{}mason frag-sequencing -- simulate sequencing of fragments",
        Emoji("🧬 ", "")
    ))?;
    term.write_line(&format!("{}configuration:", Emoji("🔧 ", "")))?;
    term.write_line(&prefix_lines(
        "   ",
        &format!("common = {:#?}", &common_args),
    ))?;
    term.write_line(&prefix_lines("   ", &format!("sequencing = {:#?}", &args)))?;
    run_simulation(term, common_args, args)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
}
