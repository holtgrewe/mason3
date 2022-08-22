pub mod bs;
/// Sequencing simulation
pub mod illumina;
pub mod roche454;
pub mod sanger;

use clap::Args as ClapArgs;
use rand::{distributions::Uniform, prelude::Distribution};
use rand_xoshiro::Xoshiro256Plus;

use bs::Args as BSArgs;

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
#[derive(PartialEq, Eq, Debug, Clone, Copy, clap::ValueEnum, enum_utils::FromStr)]
pub enum Strands {
    /// Both strands
    Both,
    /// Forward strand
    Forward,
    /// Reverse strand
    Reverse,
}

/// Select a single strand to simulate from
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub enum Strand {
    /// Forward strand
    Forward,
    /// Reverse strand
    Reverse,
}

/// Constant for selecting strand with integer
static STRANDS: &[Strand] = &[Strand::Forward, Strand::Reverse];

/// Select sequencing direction (from left or right side)
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
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
    pub cigar: CigarString,
    /// The strand that the sequence comes from
    pub strand: Strand,
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

impl SeqInfo {
    /// Create new, empty `SeqInfo` record
    pub fn new() -> Self {
        Self {
            orig_seq: Vec::new(),
            cigar: CigarString::new(),
            strand: Strand::Forward,
            contig_id: 0,
            haplotype_id: 0,
            pos: 0,
            snv_count: 0,
            indel_count: 0,
        }
    }

    /// Clear `SeqInfo` record
    pub fn clear(&mut self) {
        self.orig_seq.clear();
        self.cigar.clear();
        self.strand = Strand::Forward;
        self.contig_id = 0;
        self.haplotype_id = 0;
        self.pos = 0;
        self.snv_count = 0;
        self.indel_count = 0;
    }

    /// Convert to string, suitable for embedding in FASTA header
    pub fn comment_string(&self) -> String {
        format!(
            "SEQUENCE={} HAPLOTYPE={} BEGIN_POS={} SAMPLE_SEQUENCE={} CIGAR={} STRAND={} \
            NUM_SNPS={} NUM_INDELS={}",
            self.contig_id,
            self.haplotype_id,
            self.pos,
            String::from_utf8_lossy(&self.orig_seq),
            &self.cigar.to_str(),
            match self.strand {
                Strand::Forward => 'F',
                Strand::Reverse => 'R',
            },
            self.snv_count,
            self.indel_count,
        )
    }
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
        rng: &mut Xoshiro256Plus,
    );
}

/// Represents a CIGAR operation
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub enum CigarOp {
    Match,
    Ins,
    Del,
    // RefSkip,
    SoftClip,
    HardClip,
    Pad,
    // Equal,
    Diff,
}

/// Represents a CIGAR element
#[derive(PartialEq, Eq, Debug, Clone, Copy)]
pub struct CigarElement {
    /// The operations
    pub op: CigarOp,
    /// Number of operations
    pub count: usize,
}

impl CigarElement {
    /// Create a new CIGAR element
    pub fn new(op: CigarOp, count: usize) -> Self {
        Self { op, count }
    }

    /// Convert a CIGAR element to a `String`
    pub fn to_str(&self) -> String {
        match self.op {
            CigarOp::Match => format!("{}M", self.count),
            CigarOp::Ins => format!("{}I", self.count),
            CigarOp::Del => format!("{}D", self.count),
            // CigarOp::RefSkip => format!("{}N", self.count),
            CigarOp::SoftClip => format!("{}S", self.count),
            CigarOp::HardClip => format!("{}H", self.count),
            CigarOp::Pad => format!("{}P", self.count),
            // CigarOp::Equal => format!("{}=", self.count),
            CigarOp::Diff => format!("{}X", self.count),
        }
    }
}

/// Wrapper for a "CIGAR string"
pub struct CigarString {
    pub elements: Vec<CigarElement>,
}

impl CigarString {
    /// Construct a new, empty `CigarString`
    pub fn new() -> Self {
        Self {
            elements: Vec::new(),
        }
    }

    /// Convert a CIGAR string to its `String` representation
    pub fn to_str(&self) -> String {
        self.elements
            .iter()
            .map(|x| x.to_str())
            .collect::<Vec<String>>()
            .join("")
    }

    /// Compute length in reference
    pub fn len_in_ref(&self) -> usize {
        self.elements
            .iter()
            .filter(|elem| {
                elem.op != CigarOp::Ins
                    && elem.op != CigarOp::SoftClip
                    && elem.op != CigarOp::HardClip
                    && elem.op != CigarOp::Pad
            })
            .map(|elem| elem.count)
            .sum()
    }

    /// Clear the `CigarString`
    pub fn clear(&mut self) {
        self.elements.clear();
    }

    /// Append CIGAR operation to the `CigarString`.
    ///
    /// Returns `(a, b)` where `a` is the difference in resulting read length and `b` is the
    /// difference in used input/reference sequence length.
    ///
    /// Note: the count in `e` must be `1`!
    pub fn append(&mut self, e: CigarElement) -> (isize, isize) {
        if e.count != 1 {
            panic!("invalid CIGAR element {:?}", &e);
        }
        if let Some(last) = self.elements.last_mut() {
            // Trailing operation and new operation cancel each other out
            if (last.op == CigarOp::Ins && e.op == CigarOp::Del)
                || (last.op == CigarOp::Del && e.op == CigarOp::Ins)
            {
                if last.count > 1 {
                    last.count -= 1;
                } else {
                    self.elements.pop();
                }
                if e.op == CigarOp::Del {
                    return (-1, 0);
                } else {
                    return (0, -1);
                }
            }
        }

        // No canceling out of events.  The read length increases by one if the operation is
        // no deletion and one base of input sequence is used up if the operation is not an
        // insertion.
        if !self.elements.is_empty() && self.elements.last().unwrap().op == e.op {
            self.elements.last_mut().unwrap().count += e.count;
        } else {
            self.elements.push(e);
        }
        (
            (e.op != CigarOp::Del) as isize,
            (e.op != CigarOp::Ins) as isize,
        )
    }
}

/// Generic read simulation algorithm that provides a rich interface for simulating reads given a
/// `ReadFromFragment`
///
/// Construction of these structs is fairly cheap as it only stores only references and only
/// initialize a buffer once required.
pub struct ReadSimulator<'a> {
    /// Random number generation for read simulation
    rng: &'a mut Xoshiro256Plus,
    /// Random number generator for methylation level simulation
    meth_rng: &'a mut Xoshiro256Plus,
    /// Overall sequencing configuration
    args: &'a Args,
    /// Configuration for BS treatment
    bs_args: &'a BSArgs,
}

/// Rich interface that `ReadSimulator` implements on top of `ReadFromFragment`
impl<'a> ReadSimulator<'a> {
    /// Construct read simulator with methylation simulation
    pub fn new(
        args: &'a Args,
        bs_args: &'a BSArgs,
        rng: &'a mut Xoshiro256Plus,
        meth_rng: &'a mut Xoshiro256Plus,
    ) -> Self {
        Self {
            rng,
            meth_rng,
            args,
            bs_args,
        }
    }

    /// Simulate single-end reads
    pub fn simulate_single_end(
        &mut self,
        seq: &mut Vec<u8>,
        quals: &mut Vec<u8>,
        info: &mut SeqInfo,
        frag: &[u8],
        levels: Option<&[u8]>,
        meth_frag_buffer: &mut Vec<u8>,
        read_gen: &Box<dyn ReadFromFragment>,
    ) {
        let dist = Uniform::from(0..2);
        let strand = match self.args.strands {
            Strands::Both => STRANDS[dist.sample(self.rng)],
            Strands::Forward => Strand::Forward,
            Strands::Reverse => Strand::Reverse,
        };

        if !self.bs_args.bs_sim_enabled {
            // Forward to PE read simulation
            self.simulate_single_end_impl(seq, quals, info, frag, read_gen, strand);
        } else {
            levels.expect("missing methylation levels for BS read simulation");
            // Pick strandedness of the BS-treated fragment
            let bs_strand = if self.bs_args.protocol == bs::Protocol::Directional {
                strand
            } else {
                STRANDS[dist.sample(self.meth_rng)]
            };
            // Simulate the BS treatment of the fragment
            self.simulate_bs_treatment(
                frag,
                levels.unwrap(),
                bs_strand == Strand::Reverse,
                meth_frag_buffer,
            );
            // Simulate the actual SE read
            self.simulate_single_end_impl(seq, quals, info, meth_frag_buffer, read_gen, strand);
        }
    }

    /// Simulate paired-end reads
    pub fn simulate_paired_end(
        &mut self,
        seq_l: &mut Vec<u8>,
        quals_l: &mut Vec<u8>,
        info_l: &mut SeqInfo,
        seq_r: &mut Vec<u8>,
        quals_r: &mut Vec<u8>,
        info_r: &mut SeqInfo,
        frag: &[u8],
        levels: Option<&[u8]>,
        meth_frag_buffer: &mut Vec<u8>,
        read_gen: &Box<dyn ReadFromFragment>,
    ) {
        let dist = Uniform::from(0..2);
        let strand = match self.args.strands {
            Strands::Both => STRANDS[dist.sample(self.rng)],
            Strands::Forward => Strand::Forward,
            Strands::Reverse => Strand::Reverse,
        };

        if !self.bs_args.bs_sim_enabled {
            // Forward to PE read simulation
            self.simulate_paired_end_impl(
                seq_l, quals_l, info_l, seq_r, quals_r, info_r, frag, read_gen, strand,
            );
        } else {
            levels.expect("missing methylation levels for BS read simulation");
            // Pick strandedness of the BS-treated fragment
            let bs_strand = if self.bs_args.protocol == bs::Protocol::Directional {
                strand
            } else {
                STRANDS[dist.sample(self.meth_rng)]
            };
            // Simulate the BS treatment of the fragment
            self.simulate_bs_treatment(
                frag,
                levels.unwrap(),
                bs_strand == Strand::Reverse,
                meth_frag_buffer,
            );
            // Simulate the actual PE read
            self.simulate_paired_end_impl(
                seq_l,
                quals_l,
                info_l,
                seq_r,
                quals_r,
                info_r,
                meth_frag_buffer,
                read_gen,
                strand,
            );
        }
    }

    /// Implementation of actual single-end read simulation without BS treatment
    fn simulate_single_end_impl(
        &mut self,
        seq: &mut Vec<u8>,
        quals: &mut Vec<u8>,
        info: &mut SeqInfo,
        frag: &[u8],
        read_gen: &Box<dyn ReadFromFragment>,
        strand: Strand,
    ) {
        read_gen.simulate_read(
            seq,
            quals,
            info,
            frag,
            Direction::Left,
            strand,
            &mut self.rng,
        );
    }

    /// Implementation of actual paired-end read simulation without BS treatment
    fn simulate_paired_end_impl(
        &mut self,
        seq_l: &mut Vec<u8>,
        quals_l: &mut Vec<u8>,
        info_l: &mut SeqInfo,
        seq_r: &mut Vec<u8>,
        quals_r: &mut Vec<u8>,
        info_r: &mut SeqInfo,
        frag: &[u8],
        read_gen: &Box<dyn ReadFromFragment>,
        strand: Strand,
    ) {
        use Direction::*;
        use MateOrientation::*;
        use Strand::*;

        let (direction_l, strand_l, direction_r, strand_r) =
            match (&self.args.default_orientation, strand) {
                (ForwardReverse, Forward) => (Left, Forward, Right, Reverse),
                (ForwardReverse, Reverse) => (Right, Reverse, Left, Forward),
                (ReverseForward, Forward) => (Left, Reverse, Right, Forward),
                (ReverseForward, Reverse) => (Right, Forward, Left, Reverse),
                (ForwardForward, Forward) => (Left, Forward, Right, Forward),
                (ForwardForward, Reverse) => (Right, Reverse, Left, Reverse),
                (ForwardForward2, Forward) => (Right, Forward, Left, Forward),
                (ForwardForward2, Reverse) => (Left, Reverse, Right, Reverse),
            };

        read_gen.simulate_read(
            seq_l,
            quals_l,
            info_l,
            frag,
            direction_l,
            strand_l,
            &mut self.rng,
        );
        read_gen.simulate_read(
            seq_r,
            quals_r,
            info_r,
            frag,
            direction_r,
            strand_r,
            &mut self.rng,
        );
    }

    /// Simulate BS treatment to fragment
    fn simulate_bs_treatment(
        &self,
        _frag: &[u8],
        _levels: &[u8],
        _reverse: bool,
        _meth_frag: &mut Vec<u8>,
    ) {
        todo!()
    }
}
