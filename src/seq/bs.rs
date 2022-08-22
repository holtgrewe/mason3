/// Simulation of bisulfite treatment on reads
use clap::Args as ClapArgs;

/// Enum for the BS treatment protocol
#[derive(PartialEq, Eq, Debug, Clone, Copy, clap::ValueEnum, enum_utils::FromStr)]
pub enum Protocol {
    /// Directional protocol variant
    Directional,
    /// Undirectional protocol variant
    Undirectional,
}

const DEFAULT_ENABLE_BS_SIM: bool = false;
const DEFAULT_CONVERSION_RATE: f64 = 0.99;
const DEFAULT_BS_SEQ_PROTOCOL: Protocol = Protocol::Directional;

/// Command line arguments for BS treatment
#[derive(ClapArgs, Debug, Clone)]
pub struct Args {
    /// Whether or not BS-treatement simulation is enabled
    #[clap(long = "enable-bs-seq", action, default_value_t = DEFAULT_ENABLE_BS_SIM)]
    pub bs_sim_enabled: bool,
    /// The rate that unmethylated Cs become Ts
    #[clap(long = "bs-seq-conversion-rate", default_value_t = DEFAULT_CONVERSION_RATE)]
    pub bs_conversion_rate: f64,
    /// The protocol to use for the simulation
    #[clap(long = "bs-seq-protocol", value_enum, default_value_t = DEFAULT_BS_SEQ_PROTOCOL)]
    pub protocol: Protocol,
}

impl Args {
    /// Create new args with same defaults as those from clap
    pub fn new() -> Self {
        Self {
            bs_sim_enabled: DEFAULT_ENABLE_BS_SIM,
            bs_conversion_rate: DEFAULT_CONVERSION_RATE,
            protocol: DEFAULT_BS_SEQ_PROTOCOL,
        }
    }
}
