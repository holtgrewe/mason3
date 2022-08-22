/// Implementation of methylation level simulation
///
/// Methylation levels are represented by characters in resolution of `0.0125` steps. ASCII
/// characters are used starting with `'!'` for `0.0`, `'>'` is skipped so there cannot
/// be any conflict in FASTA representation.
///
/// The methylation levels for the forward strand are written out as `${name}/TOP` and
/// for the reverse strand, that are written out as `${name}/BOT` where `${name}` is the
/// name of the original contig.
use bio::alphabets;
use bio::utils::TextSlice;
use clap::Args as ClapArgs;
use console::Emoji;
use indicatif::{ProgressBar, ProgressStyle};
use rand::prelude::Distribution;
use rand_xoshiro::Xoshiro256Plus;
use statrs::distribution::{Beta, ContinuousCDF};

const M_CG: usize = 10;

const M_CCG: usize = 74;
const M_CAG: usize = 66;
const M_CTG: usize = 98;
const M_CAA: usize = 64;
const M_CAC: usize = 65;
const M_CAT: usize = 68;
const M_CCA: usize = 72;
const M_CCC: usize = 73;
const M_CCT: usize = 76;
const M_CTA: usize = 96;
const M_CTC: usize = 97;
const M_CTT: usize = 100;
const M_AAG: usize = 2;
const M_AGG: usize = 18;
const M_ATG: usize = 34;
const M_GAG: usize = 130;
const M_GGG: usize = 146;
const M_GTG: usize = 162;
const M_TAG: usize = 258;
const M_TGG: usize = 274;
const M_TTG: usize = 290;

/// Configuration for simulating methylation levels
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Whether to enable methylation level simulation
    #[clap(default_value_t = true)]
    pub enabled: bool,
    /// Mean of beta distribution for methylation at CpG loci
    #[clap(long = "meth-cg-mu", default_value_t = 0.6)]
    pub cg_mu: f64,
    /// Standard deviation of beta distribution for methylation at CpG loci
    #[clap(long = "meth-cg-sigma", default_value_t = 0.03)]
    pub cg_sigma: f64,
    /// Mean of beta distribution for methylation at CHG loci
    #[clap(long = "meth-chg-mu", default_value_t = 0.08)]
    pub chg_mu: f64,
    /// Standard deviation of beta distribution for methylation at CHG loci
    #[clap(long = "meth-chg-sigma", default_value_t = 0.008)]
    pub chg_sigma: f64,
    /// Mean of beta distribution for methylation at CHH loci
    #[clap(long = "meth-chh-mu", default_value_t = 0.05)]
    pub chh_mu: f64,
    /// Standard deviation of beta distribution for methylation at CHH loci
    #[clap(long = "meth-chh-sigma", default_value_t = 0.005)]
    pub chh_sigma: f64,
}

/// Collect the PDFs that we need for the simulation
struct MethPDFs {
    /// PDF for CpG loci
    cg: Beta,
    /// PDF for CHG loci
    chg: Beta,
    /// PDF for CHH loci
    chh: Beta,
}

impl MethPDFs {
    /// Create with mu/sigmas for beta distributions
    fn new(cg_a: f64, cg_b: f64, chg_a: f64, chg_b: f64, chh_a: f64, chh_b: f64) -> Self {
        Self {
            cg: Beta::new(cg_a, cg_b).unwrap(),
            chg: Beta::new(chg_a, chg_b).unwrap(),
            chh: Beta::new(chh_a, chh_b).unwrap(),
        }
    }
}

/// Compute alpha for beta distribution with mu and sigma
fn mu_sigma_to_alpha(mu: f64, sigma: f64) -> f64 {
    ((1.0 - mu) / sigma / sigma - 1.0 / mu) * mu * mu
}

/// Compute beta for beta distribution with mu and sigma
fn mu_sigma_to_beta(mu: f64, sigma: f64) -> f64 {
    mu_sigma_to_alpha(mu, sigma) * (1.0 / mu - 1.0)
}

/// Convert methylation level (in `[0.0, 1.0]`) to character
fn level_to_char(level: f64) -> u8 {
    let c = (level / 0.0125).round() as u8;
    if level > 1.0 {
        println!("level = {}, c = {}", level, c);
    }
    if c + b'!' < b'>' {
        c + 33
    } else {
        c + 34
    }
}

// /// Convert character to level
#[allow(dead_code)]
pub fn char_to_level(char: u8) -> f64 {
    if char < b'>' {
        ((char - 33) as f64) * 0.0125
    } else {
        ((char - 34) as f64) * 0.0125
    }
}

/// Update slice `vals` at position `i` with `max(vals[i], val)`
fn update_slice(pos: usize, val: u8, vals: &mut [u8]) {
    vals[pos] = std::cmp::max(vals[pos], val);
}

/// Compute random level character
fn random_level_char(beta: &Beta, rng: &mut Xoshiro256Plus) -> u8 {
    let x = beta.sample(rng);
    let level = beta.cdf(x);
    level_to_char(level)
}

/// Handle 2-mer at position.
fn handle_two_mer(
    levels_top: &mut [u8],
    levels_bot: &mut [u8],
    rng: &mut Xoshiro256Plus,
    pos: usize,
    hash_value: usize,
    pdfs: &MethPDFs,
) {
    if hash_value == M_CG {
        // CpG is symmetric
        update_slice(pos, random_level_char(&pdfs.cg, rng), levels_top);
        update_slice(pos + 1, random_level_char(&pdfs.cg, rng), levels_bot);
    }
}

/// Handle 3-mer at position.
fn handle_three_mer(
    levels_top: &mut [u8],
    levels_bot: &mut [u8],
    rng: &mut Xoshiro256Plus,
    pos: usize,
    hash_value: usize,
    pdfs: &MethPDFs,
) {
    match hash_value {
        M_CAG | M_CTG => {
            update_slice(pos, random_level_char(&pdfs.chg, rng), levels_top);
            update_slice(pos + 2, random_level_char(&pdfs.chg, rng), levels_bot);
        }
        M_CCG => {
            update_slice(pos, random_level_char(&pdfs.chg, rng), levels_top);
            update_slice(pos + 2, random_level_char(&pdfs.chg, rng), levels_bot);
        }
        M_CAA | M_CAC | M_CAT | M_CCA | M_CCC | M_CCT | M_CTA | M_CTC | M_CTT => {
            update_slice(pos, random_level_char(&pdfs.chh, rng), levels_top);
        }
        M_AAG | M_AGG | M_ATG | M_GAG | M_GGG | M_GTG | M_TAG | M_TGG | M_TTG => {
            update_slice(pos + 2, random_level_char(&pdfs.chh, rng), levels_bot);
        }
        _ => (),
    }
}

/// Runs the methylation level simulation
pub fn simulate_methylation_levels(
    id: &str,
    seq: TextSlice<'_>,
    args: &Args,
    levels_top: &mut Vec<u8>,
    levels_bot: &mut Vec<u8>,
    rng: &mut Xoshiro256Plus,
) {
    levels_top.resize(seq.len(), b'!');
    levels_bot.resize(seq.len(), b'!');

    let contig_len = seq.len();
    let bar = ProgressBar::new(contig_len.try_into().unwrap());
    bar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{msg} | {wide_bar:.cyan/blue} {pos:>7}/{len:7} [{elapsed_precise}/{eta_precise}]",
            )
            .unwrap()
            .progress_chars("##-"),
    );
    bar.set_message(id.to_owned());
    bar.set_position(0);

    let pdfs = MethPDFs::new(
        mu_sigma_to_alpha(args.cg_mu, args.cg_sigma),
        mu_sigma_to_beta(args.cg_mu, args.cg_sigma),
        mu_sigma_to_alpha(args.chg_mu, args.chg_sigma),
        mu_sigma_to_beta(args.chg_mu, args.chg_sigma),
        mu_sigma_to_alpha(args.chh_mu, args.chh_sigma),
        mu_sigma_to_beta(args.chh_mu, args.chh_sigma),
    );

    // Create local DNA5 alphabet and rank transform
    let dna_alphabet = alphabets::Alphabet::new(b"ACGTN");
    let dna_ranks = alphabets::RankTransform::new(&dna_alphabet);

    let mut qgrams2 = dna_ranks.qgrams(2, seq);
    let qgrams3 = dna_ranks.qgrams(3, seq);

    let mut pos: usize = 0;
    for hash_value_3 in qgrams3 {
        let hash_value_2 = qgrams2.next().unwrap();

        handle_two_mer(levels_top, levels_bot, rng, pos, hash_value_2, &pdfs);
        handle_three_mer(levels_top, levels_bot, rng, pos, hash_value_3, &pdfs);

        // advance progress bar
        if pos % 10_000 == 0 {
            bar.set_position(pos.try_into().unwrap());
        }
        // advance position
        pos += 1;
    }
    // Handle trailing qgram
    if let Some(hash_value_2) = qgrams2.next() {
        handle_two_mer(levels_top, levels_bot, rng, pos, hash_value_2, &pdfs);
    }

    bar.finish_and_clear();
    bar.println(&format!("{}{} done", Emoji("✓ ", ""), id,));
}

#[cfg(test)]
mod tests {
    use super::*;

    use statrs::distribution::Beta;
    use statrs::prec::almost_eq;
    use statrs::statistics::Distribution;

    #[test]
    fn test_level_to_char() {
        assert_eq!(level_to_char(0.0), b'!');
        assert_eq!(level_to_char(1.0), b'r');

        assert_eq!(level_to_char(0.0062), b'!');
        assert_eq!(level_to_char(0.0063), b'"');
    }

    #[test]
    fn test_char_to_level() {
        assert_eq!(char_to_level(b'!'), 0.0);
        assert_eq!(char_to_level(b'r'), 1.0);
    }

    #[test]
    fn test_mu_sigma_beta() {
        {
            let beta1 =
                Beta::new(mu_sigma_to_alpha(0.5, 0.05), mu_sigma_to_beta(0.5, 0.05)).unwrap();
            assert!(almost_eq(beta1.mean().unwrap(), 0.5, 0.001));
            assert!(almost_eq(beta1.variance().unwrap(), 0.05 * 0.05, 0.001));
        }

        {
            let beta2 = Beta::new(
                mu_sigma_to_alpha(0.08, 0.008),
                mu_sigma_to_beta(0.08, 0.008),
            )
            .unwrap();
            assert!(almost_eq(beta2.mean().unwrap(), 0.08, 0.001));
            assert!(almost_eq(beta2.variance().unwrap(), 0.008 * 0.008, 0.001));
        }

        {
            let beta3 = Beta::new(
                mu_sigma_to_alpha(0.06, 0.005),
                mu_sigma_to_beta(0.06, 0.005),
            )
            .unwrap();
            assert!(almost_eq(beta3.mean().unwrap(), 0.06, 0.001));
            assert!(almost_eq(beta3.variance().unwrap(), 0.005 * 0.005, 0.001));
        }
    }

    #[test]
    fn test_hash_two_mers() {
        let dna_alphabet = alphabets::Alphabet::new(b"ACGTN");
        let dna_ranks = alphabets::RankTransform::new(&dna_alphabet);

        let seqs = vec![(b"CG", M_CG)];
        for (seq, expected) in seqs {
            let mut qgrams = dna_ranks.qgrams(2, seq);
            let hash_value = qgrams.next().unwrap();
            assert_eq!(hash_value, expected);
        }
    }

    #[test]
    fn test_hash_three_mers() {
        let dna_alphabet = alphabets::Alphabet::new(b"ACGTN");
        let dna_ranks = alphabets::RankTransform::new(&dna_alphabet);

        let seqs = vec![
            (b"CCG", M_CCG),
            (b"CAG", M_CAG),
            (b"CTG", M_CTG),
            (b"CAA", M_CAA),
            (b"CAC", M_CAC),
            (b"CAT", M_CAT),
            (b"CCA", M_CCA),
            (b"CCC", M_CCC),
            (b"CCT", M_CCT),
            (b"CTA", M_CTA),
            (b"CTC", M_CTC),
            (b"CTT", M_CTT),
            (b"AAG", M_AAG),
            (b"AGG", M_AGG),
            (b"ATG", M_ATG),
            (b"GAG", M_GAG),
            (b"GGG", M_GGG),
            (b"GTG", M_GTG),
            (b"TAG", M_TAG),
            (b"TGG", M_TGG),
            (b"TTG", M_TTG),
        ];
        for (seq, expected) in seqs {
            let mut qgrams = dna_ranks.qgrams(3, seq);
            let hash_value = qgrams.next().unwrap();
            assert_eq!(hash_value, expected);
        }
    }
}
