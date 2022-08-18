pub mod bs;
/// Sequencing simulation
pub mod illumina;
pub mod roche454;
pub mod sanger;

use clap::Args as ClapArgs;
use fastrand::Rng;

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
static STRANDS: &'static [Strand] = &[Strand::Forward, Strand::Reverse];

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

impl SeqInfo {
    /// Create new, empty `SeqInfo` record
    pub fn new() -> Self {
        Self {
            orig_seq: Vec::new(),
            cigar: String::new(),
            is_forward: false,
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
        self.is_forward = false;
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
            &self.cigar,
            if self.is_forward { 'F' } else { 'R' },
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
    /// Configuration for BS treatment
    bs_args: &'a BSArgs,
}

/// Rich interface that `ReadSimulator` implements on top of `ReadFromFragment`
impl<'a> ReadSimulator<'a> {
    /// Construct read simulator with methylation simulation
    pub fn new(
        seq_args: &'a Args,
        bs_args: &'a BSArgs,
        rng: &'a mut Rng,
        meth_rng: &'a mut Rng,
    ) -> Self {
        Self {
            rng,
            meth_rng,
            seq_args,
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
        let strand = match self.seq_args.strands {
            Strands::Both => STRANDS[self.rng.usize(0..2)],
            Strands::Forward => Strand::Forward,
            Strands::Reverse => Strand::Reverse,
        };

        if self.bs_args.bs_sim_enabled {
            // Forward to PE read simulation
            self.simulate_single_end_impl(seq, quals, info, frag, read_gen, strand);
        } else {
            levels.expect("missing methylation levels for BS read simulation");
            // Pick strandedness of the BS-treated fragment
            let bs_strand = if self.bs_args.protocol == bs::Protocol::Directional {
                strand
            } else {
                STRANDS[self.meth_rng.usize(0..2)]
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
        let strand = match self.seq_args.strands {
            Strands::Both => STRANDS[self.rng.usize(0..2)],
            Strands::Forward => Strand::Forward,
            Strands::Reverse => Strand::Reverse,
        };

        if self.bs_args.bs_sim_enabled {
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
                STRANDS[self.meth_rng.usize(0..2)]
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
        &self,
        seq: &mut Vec<u8>,
        quals: &mut Vec<u8>,
        info: &mut SeqInfo,
        _frag: &[u8],
        _read_gen: &Box<dyn ReadFromFragment>,
        _strand: Strand,
    ) {
        seq.clear();
        quals.clear();
        info.clear();
        todo!()
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
        _frag: &[u8],
        _read_gen: &Box<dyn ReadFromFragment>,
        _strand: Strand,
    ) {
        seq_l.clear();
        quals_l.clear();
        info_l.clear();
        seq_r.clear();
        quals_r.clear();
        info_r.clear();
        todo!();
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
