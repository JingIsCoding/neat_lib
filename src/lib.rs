mod loader;
mod network;

pub use network::population::Population;
pub use network::genome::Genome;
pub use network::gene::Gene;
pub use network::connection_gene::ConnectionGene;
pub use network::config::*;
pub use network::errors::*;
pub use loader::save_load;
