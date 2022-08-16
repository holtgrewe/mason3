/// Sequencing simulation
pub mod illumina;
pub mod roche454;
pub mod sanger;

use clap::Args as ClapArgs;

/// Enum for selecting sequencing technology
#[derive(PartialEq, Eq, Debug, Clone, clap::ValueEnum, enum_utils::FromStr)]
pub enum SequencingTechnology {
    /// Illumina sequencing
    Illumina,
    /// Roche 454 sequencing
    Roche454,
    /// Sanger/capillary sequencing
    Sanger,
}

/// Enum for the mate orientation in case of paired sequencing
#[derive(PartialEq, Eq, Debug, Clone, clap::ValueEnum, enum_utils::FromStr)]
pub enum MateOrientation {
    /// R1 --> <-- R2
    ForwardReverse,
    /// R1 <-- --> R1
    ReverseForward,
    /// R1 --> --> R2
    ForwardForward,
    /// R2 --> --> R1
    ForwardForward2,
}

/// Select strands to simulate sequencing from
#[derive(PartialEq, Eq, Debug, Clone, clap::ValueEnum, enum_utils::FromStr)]
pub enum Strands {
    /// Both strands
    Both,
    /// Forward strand
    Forward,
    /// Reverse strand
    Reverse,
}

#[derive(PartialEq, Eq, Debug, Clone, clap::ValueEnum, enum_utils::FromStr)]
pub enum ReadLengthModel {
    /// Uniform distribution of read lengths
    Uniform,
    /// Normal distribution of read lengths
    Normal,
}

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// The sequencing technology used
    #[clap(long = "seq-technology", value_enum, default_value_t = SequencingTechnology::Illumina)]
    technology: SequencingTechnology,
    /// The default orientation for Illumina paired-end reads
    #[clap(long = "seq-mate-orientation", value_enum, default_value_t = MateOrientation::ForwardReverse)]
    default_orientation: MateOrientation,
    /// The strands to simulate sequences from
    #[clap(long = "seq-strands", value_enum, default_value_t = Strands::Reverse)]
    strands: Strands,
    /// Whether or not to embed read simulation information
    #[clap(long = "embed-read-info", action, default_value_t = false)]
    embed_read_info: bool,
    /// The read name prefix to use
    #[clap(long = "read-name-prefix", default_value = "simulated.")]
    read_name_prefix: String,
}
