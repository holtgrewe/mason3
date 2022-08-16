/// Sequencing simulation
pub mod illumina;
pub mod roche454;
pub mod sanger;

use clap::Args as ClapArgs;
use fastrand::Rng;

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

/// Select a single strand to simulate from
pub enum Strand {
    /// Forward strand
    Forward,
    /// Reverse strand
    Reverse,
}

/// Select sequencing direction (from left or right side)
pub enum Direction {
    /// Left to right
    Left,
    /// Right to left
    Right,
}

#[derive(PartialEq, Eq, Debug, Clone, clap::ValueEnum, enum_utils::FromStr)]
pub enum ReadLengthModel {
    /// Uniform distribution of read lengths
    Uniform,
    /// Normal distribution of read lengths
    Normal,
}

#[derive(ClapArgs, Debug, Clone)]
pub struct Args {
    /// The sequencing technology used
    #[clap(long = "seq-technology", value_enum, default_value_t = SequencingTechnology::Illumina)]
    pub technology: SequencingTechnology,
    /// The default orientation for Illumina paired-end reads
    #[clap(long = "seq-mate-orientation", value_enum, default_value_t = MateOrientation::ForwardReverse)]
    pub default_orientation: MateOrientation,
    /// The strands to simulate sequences from
    #[clap(long = "seq-strands", value_enum, default_value_t = Strands::Reverse)]
    pub strands: Strands,
    /// Whether or not to embed read simulation information
    #[clap(long = "embed-read-info", action, default_value_t = false)]
    pub embed_read_info: bool,
    /// The read name prefix to use
    #[clap(long = "read-name-prefix", default_value = "simulated.")]
    pub read_name_prefix: String,
}

/// Information from sequencing
pub struct SeqInfo {
    /// Original sample sequence, that together with the errors introduced by sequencing gives the
    /// read sequence
    pub orig_seq: Vec<u8>,
    /// The CIGAR string, `MXID` for matches, mismatches, insertions, deletions with respect to
    /// the reference
    pub cigar: String,
    /// Whether or not this comes from the forward strand
    pub is_forward: bool,
    /// The contig ID
    pub contig_id: usize,
    /// The haplotype ID
    pub haplotype_id: usize,
    /// The begin position of the sequence
    pub pos: usize,
    /// Number of bases covering SNVs
    pub snv_count: usize,
    /// NUmber of bases covering indels
    pub indel_count: usize,
}

/// Trait for simulating a read from a fragment
pub trait ReadFromFragment {
    /// Simulate a read
    fn simulate_read(
        &self,
        seq: &mut Vec<u8>,
        quals: &mut Vec<u8>,
        info: &mut SeqInfo,
        frag: &[u8],
        dir: Direction,
        strand: Strand,
    );
}

/// Generic read simulation algorithm that provides a rich interface for simulating reads given a
/// `ReadFromFragment`
///
/// Construction of these structs is fairly cheap as it only stores only references and only
/// initialize a buffer once required.
pub struct ReadSimulator<'a> {
    /// Random number generation for read simulation
    rng: &'a mut Rng,
    /// Random number generator for methylation level simulation
    meth_rng: &'a mut Rng,
    /// Overall sequencing configuration
    seq_args: &'a Args,
    /// Buffer for the materialization of BS-seq treated fragments
    meth_frag: Option<Vec<u8>>,
}

/// Rich interface that `ReadSimulator` implements on top of `ReadFromFragment`
impl<'a> ReadSimulator<'a> {
    /// Construct read simulator with methylation simulation
    pub fn new(seq_args: &'a Args, rng: &'a mut Rng, meth_rng: &'a mut Rng) -> Self {
        Self {
            rng: rng,
            meth_rng: meth_rng,
            seq_args: seq_args,
            meth_frag: None,
        }
    }

    /// Simulate single-end reads
    pub fn simulate_single_end(
        &mut self,
        _seq: &mut Vec<u8>,
        _quals: &mut Vec<u8>,
        _info: &mut SeqInfo,
        _frag: &[u8],
        _levels: Option<&[u8]>,
        _sim: &dyn ReadFromFragment,
    ) {
        todo!()
    }

    /// Simulate paired-end reads
    pub fn simulate_paired_end(
        &mut self,
        _seq_l: &mut Vec<u8>,
        _quals_l: &mut Vec<u8>,
        _info_l: &mut SeqInfo,
        _seq_r: &mut Vec<u8>,
        _quals_r: &mut Vec<u8>,
        _info_r: &mut SeqInfo,
        _frag: &[u8],
        _levels: Option<&[u8]>,
        _sim: &dyn ReadFromFragment,
    ) {
        todo!()
    }
}
