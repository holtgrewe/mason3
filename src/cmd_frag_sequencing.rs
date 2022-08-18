/// `mason frag-sequencing` -- simulate reads from fragments
use bio::io::fasta;
use clap::Args as ClapArgs;
use console::{Emoji, Term};
use fastrand::Rng;

use crate::common::prefix_lines;
use crate::common::Args as CommonArgs;
use crate::seq::bs::Args as BSArgs;
use crate::seq::illumina::Args as IlluminaArgs;
use crate::seq::illumina::IlluminaFromFragment;
use crate::seq::roche454::Args as Roche454Args;
use crate::seq::roche454::Roche454FromFragment;
use crate::seq::sanger::Args as SangerArgs;
use crate::seq::sanger::SangerFromFragment;
use crate::seq::Args as SequencingArgs;
use crate::seq::ReadFromFragment;
use crate::seq::ReadSimulator;
use crate::seq::SeqInfo;

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
fn run_simulation(
    _term: &Term,
    common_args: &CommonArgs,
    args: &Args,
) -> Result<(), anyhow::Error> {
    // Construct readers and writers for file I/o
    let reader = fasta::Reader::from_file(&args.input_filename)?;
    let mut writer_left = fasta::Writer::to_file(&args.output_filename_left)?;
    let mut writer_right = if let Some(path) = &args.output_filename_right {
        if !args.force_single_end {
            Some(fasta::Writer::to_file(&path)?)
        } else {
            None
        }
    } else {
        None
    };
    // Initialize random number generators
    let mut rng = Rng::with_seed(common_args.seed);
    let mut rng_meth = Rng::with_seed(common_args.seed);
    // Create dummy for BS treatment simulation (deactivated)
    let dummy_bs_args = BSArgs::new();

    // Create appropriate helpers for the actual read simulation
    let mut simulator =
        ReadSimulator::new(&args.sequencing, &dummy_bs_args, &mut rng, &mut rng_meth);
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

    // Allocate buffers
    let mut seq_l = Vec::new();
    let mut seq_r = Vec::new();
    let mut quals_l = Vec::new();
    let mut quals_r = Vec::new();
    let mut info_l = SeqInfo::new();
    let mut info_r = SeqInfo::new();
    let mut meth_frag_buffer_dummy = Vec::new();

    // Perform the actual read simulation
    for (i, frag_record) in reader.records().enumerate() {
        let frag_record = frag_record?;
        // Build read ID, optionally with fragment ID as additional read info
        let read_id = format!("{}{}", args.sequencing.read_name_prefix, i + 1);

        // Simulate SE/PE read and write to file
        if let Some(ref mut writer_right) = writer_right.as_mut() {
            simulator.simulate_paired_end(
                &mut seq_l,
                &mut quals_l,
                &mut info_l,
                &mut seq_r,
                &mut quals_r,
                &mut info_r,
                &frag_record.seq(),
                None,
                &mut meth_frag_buffer_dummy,
                &read_gen,
            );

            let (desc_l, desc_r) = if args.sequencing.embed_read_info {
                (
                    Some(format!(
                        "{} FRAG_ID={}",
                        info_l.comment_string(),
                        &frag_record.id()
                    )),
                    Some(format!(
                        "{} FRAG_ID={}",
                        info_r.comment_string(),
                        &frag_record.id()
                    )),
                )
            } else {
                (None, None)
            };

            writer_left.write(&read_id, desc_l.as_deref(), &seq_l)?;
            writer_right.write(&read_id, desc_r.as_deref(), &seq_r)?;
        } else {
            simulator.simulate_single_end(
                &mut seq_l,
                &mut quals_l,
                &mut info_l,
                &frag_record.seq(),
                None,
                &mut meth_frag_buffer_dummy,
                &read_gen,
            );

            let desc_l = Some(format!(
                "{} FRAG_ID={}",
                info_l.comment_string(),
                &frag_record.id()
            ));

            writer_left.write(&read_id, desc_l.as_deref(), &seq_l)?;
        }
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
