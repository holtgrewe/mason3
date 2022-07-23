/// mason_genome - simulate random sequence
mod err;

use fastrand::Rng;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

use clap::Parser;
use clap_verbosity_flag::Verbosity;
use console::{Emoji, Term};
use indicatif::{ProgressBar, ProgressStyle};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Verbosity of the program.
    #[clap(flatten)]
    verbose: Verbosity,
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
    /// The seed to use for the random number generator.
    #[clap(short = 's', long = "seed", default_value_t = 42)]
    seed: u64,
}

/// Simulate genome to output file.
fn simulate_genome_to_file(term: &Term, args: &Args, file: &mut File) -> Result<(), anyhow::Error> {
    term.write_line(&format!(
        "{} Starting genome simulation...",
        Emoji("🛫 ", "")
    ))?;
    let rng = Rng::with_seed(args.seed);
    let mut writer = BufWriter::new(file);

    for (i, contig_len) in args.contig_lengths.iter().enumerate() {
        let bar = ProgressBar::new(*contig_len);
        bar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{msg} | {wide_bar:.cyan/blue} {pos:>7}/{len:7} [{elapsed_precise}/{eta_precise}]",
                )
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

fn run(term: &Term, args: &Args) -> Result<(), anyhow::Error> {
    term.write_line(&format!("{}configuration:\n{:#?}", Emoji("🔧 ", ""), &args))?;
    let mut output_file = File::create(&args.output_filename)?;
    simulate_genome_to_file(term, args, &mut output_file)?;
    Ok(())
}

fn main() -> Result<(), anyhow::Error> {
    let args = Args::parse();

    let term = Term::stderr();
    term.write_line(&format!(
        "{}mason_genome -- genome simulator",
        Emoji("🧬 ", "")
    ))?;
    run(&term, &args)?;
    term.write_line(&format!("All done. Have a nice day!{}", Emoji(" 😃", "")))?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{run, Args};

    use console::Term;
    use file_diff::diff;
    use tempdir::TempDir;

    #[test]
    fn test_no_line_length() {
        let tmp_dir = TempDir::new("test").unwrap();
        let args = Args {
            verbose: clap_verbosity_flag::Verbosity::new(0, 0),
            output_filename: tmp_dir
                .path()
                .join("out.fa")
                .into_os_string()
                .into_string()
                .unwrap(),
            contig_lengths: vec![100, 100],
            line_length: 0,
            seed: 0,
        };
        let term = Term::stderr();

        run(&term, &args).unwrap();

        assert!(!diff("./tests/genome/expected.fa", &args.output_filename));
    }

    #[test]
    fn test_line_length_80() {
        let tmp_dir = TempDir::new("test").unwrap();
        let args = Args {
            verbose: clap_verbosity_flag::Verbosity::new(0, 0),
            output_filename: tmp_dir
                .path()
                .join("out.fa")
                .into_os_string()
                .into_string()
                .unwrap(),
            contig_lengths: vec![100, 100],
            line_length: 80,
            seed: 0,
        };
        let term = Term::stderr();

        run(&term, &args).unwrap();

        assert!(!diff(
            "./tests/genome/expected-ll-80.fa",
            &args.output_filename
        ));
    }
}
