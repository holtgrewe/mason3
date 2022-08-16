/// `mason methylation` - simulate methylation levels
use bio::io::fasta;
use clap::Args as ClapArgs;
use console::{Emoji, Term};
use fastrand::Rng;

use crate::common::prefix_lines;
use crate::common::Args as CommonArgs;
use crate::methylation::simulate_methylation_levels;
use crate::methylation::Args as MethylationArgs;

/// Configuration for `methylation` sub command
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Input file name
    #[clap(short = 'i', long = "in")]
    input_filename: String,
    /// Output file name
    #[clap(short = 'o', long = "out")]
    output_filename: String,

    /// Methylation level simulation
    #[clap(flatten)]
    methylation: MethylationArgs,
}

/// Simulate methylation levels according to configuration
fn simulate_levels(
    term: &Term,
    common_args: &CommonArgs,
    args: &Args,
) -> Result<(), anyhow::Error> {
    term.write_line(&format!(
        "{} Starting methylation level simulation...",
        Emoji("🛫 ", "")
    ))?;
    let mut rng = Rng::with_seed(common_args.seed);
    let reader = fasta::Reader::from_file(&args.input_filename)?;
    let mut writer = fasta::Writer::to_file(&args.output_filename)?;

    let mut records = reader.records();
    while let Some(Ok(record)) = records.next() {
        let mut levels_top: Vec<u8> = Vec::with_capacity(record.seq().len());
        let mut levels_bot: Vec<u8> = Vec::with_capacity(record.seq().len());

        simulate_methylation_levels(
            record.id(),
            record.seq(),
            &args.methylation,
            &mut levels_top,
            &mut levels_bot,
            &mut rng,
        );

        writer.write(&format!("{}/TOP", &record.id()), None, &levels_top)?;
        writer.write(&format!("{}/BOT", &record.id()), None, &levels_bot)?;
    }

    term.write_line(&format!(
        "{}done with methylation level simulation",
        Emoji("🛬 ", "")
    ))?;
    Ok(())
}

pub fn run(term: &Term, common_args: &CommonArgs, args: &Args) -> Result<(), anyhow::Error> {
    term.write_line(&format!(
        "{}mason methylation -- simulate methylation levels",
        Emoji("📈 ", "")
    ))?;
    term.write_line(&format!("{}configuration:", Emoji("🔧 ", "")))?;
    term.write_line(&prefix_lines(
        "   ",
        &format!("common = {:#?}", &common_args),
    ))?;
    term.write_line(&prefix_lines("   ", &format!("methylation = {:#?}", &args)))?;
    simulate_levels(term, common_args, args)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use clap_verbosity_flag;
    use console::Term;
    use file_diff::diff;
    use tempdir::TempDir;

    #[test]
    fn test_run_smoke_test() {
        let tmp_dir = TempDir::new("test").unwrap();
        let common_args = CommonArgs {
            verbose: clap_verbosity_flag::Verbosity::new(0, 0),
            seed: 0,
        };
        let args = Args {
            input_filename: "./tests/methylation/input.fa".to_string(),
            output_filename: tmp_dir
                .path()
                .join("out.fa")
                .into_os_string()
                .into_string()
                .unwrap(),
            methylation: MethylationArgs {
                enabled: true,
                cg_mu: 0.6,
                cg_sigma: 0.03,
                chg_mu: 0.08,
                chg_sigma: 0.008,
                chh_mu: 0.05,
                chh_sigma: 0.005,
            },
        };
        let term = Term::stderr();

        run(&term, &common_args, &args).unwrap();

        assert!(!diff(
            "./tests/methylation/expected.fa",
            &args.output_filename
        ));
    }
}
