/// `mason genome` - simulate random sequence
use fastrand::Rng;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use clap::Args as ClapArgs;
use console::{Emoji, Term};
use indicatif::{ProgressBar, ProgressStyle};

use crate::common::Args as CommonArgs;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Output file name
    #[clap(short = 'o', long = "out-file")]
    output_filename: String,
    /// Lengths of the contigs
    #[clap(short = 'l', long = "contig-length", required = true)]
    contig_lengths: Vec<u64>,
    /// The length of lines in FASTA file to write out to.  Use 0 to write one
    /// line per contig.
    #[clap(long = "line-length", default_value_t = 70)]
    line_length: u64,
}

/// Simulate genome to output file.
fn simulate_genome_to_file(
    term: &Term,
    common_args: &CommonArgs,
    args: &Args,
    file: &mut File,
) -> Result<(), anyhow::Error> {
    term.write_line(&format!(
        "{} Starting genome simulation...",
        Emoji("🛫 ", "")
    ))?;
    let rng = Rng::with_seed(common_args.seed);
    let mut writer = BufWriter::new(file);

    for (i, contig_len) in args.contig_lengths.iter().enumerate() {
        let bar = ProgressBar::new(*contig_len);
        bar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{msg} | {wide_bar:.cyan/blue} {pos:>7}/{len:7} [{elapsed_precise}/{eta_precise}]",
                )
                .unwrap()
                .progress_chars("##-"),
        );
        bar.set_message(format!("contig {}/{}", i + 1, args.contig_lengths.len()));
        bar.set_position(0);

        let chars = b"ACGT";

        writeln!(&mut writer, ">{}", i + 1)?;
        for i in 0..*contig_len {
            // write line breaks
            if i > 0 && args.line_length > 0 && i % args.line_length == 0 {
                writeln!(&mut writer)?;
            }
            // sample and write character
            let c = rng.usize(0..4);
            writer.write_all(&[chars[c]])?;
            // advance progress bar
            if i % 10_000 == 0 {
                bar.set_position(i);
            }
        }
        writeln!(&mut writer)?;

        bar.finish_and_clear();
        bar.println(&format!(
            "{}contig {}/{} done",
            Emoji("✓ ", ""),
            i + 1,
            args.contig_lengths.len()
        ));
    }

    term.write_line(&format!("{}done with genome simulation", Emoji("🛬 ", "")))?;
    Ok(())
}

pub fn run(term: &Term, common_args: &CommonArgs, args: &Args) -> Result<(), anyhow::Error> {
    term.write_line(&format!(
        "{}mason genome -- genome simulator",
        Emoji("🧬 ", "")
    ))?;
    term.write_line(&format!(
        "{}configuration:\ncommon = {:#?}\ngenome = {:#?}",
        Emoji("🔧 ", ""),
        &common_args,
        &args
    ))?;
    let mut output_file = File::create(&args.output_filename)?;
    simulate_genome_to_file(term, common_args, args, &mut output_file)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{run, Args};
    use crate::common::Args as CommonArgs;

    use clap_verbosity_flag;
    use console::Term;
    use file_diff::diff;
    use tempdir::TempDir;

    #[test]
    fn test_no_line_length() {
        let tmp_dir = TempDir::new("test").unwrap();
        let common_args = CommonArgs {
            verbose: clap_verbosity_flag::Verbosity::new(0, 0),
            seed: 0,
        };
        let args = Args {
            output_filename: tmp_dir
                .path()
                .join("out.fa")
                .into_os_string()
                .into_string()
                .unwrap(),
            contig_lengths: vec![100, 100],
            line_length: 0,
        };
        let term = Term::stderr();

        run(&term, &common_args, &args).unwrap();

        assert!(!diff("./tests/genome/expected.fa", &args.output_filename));
    }

    #[test]
    fn test_line_length_80() {
        let tmp_dir = TempDir::new("test").unwrap();
        let common_args = CommonArgs {
            verbose: clap_verbosity_flag::Verbosity::new(0, 0),
            seed: 0,
        };
        let args = Args {
            output_filename: tmp_dir
                .path()
                .join("out.fa")
                .into_os_string()
                .into_string()
                .unwrap(),
            contig_lengths: vec![100, 100],
            line_length: 80,
        };
        let term = Term::stderr();

        run(&term, &common_args, &args).unwrap();

        assert!(!diff(
            "./tests/genome/expected-ll-80.fa",
            &args.output_filename
        ));
    }
}
