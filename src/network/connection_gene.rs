use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ConnectionGene {
    pub from: usize,
    pub to: usize,
    pub innovation: usize,
    pub weight: f64,
    pub enabled: bool,
}

impl ConnectionGene {
    pub fn new(from: usize, to: usize, innovation: usize, weight: f64, enabled: bool) -> Self {
        ConnectionGene { from, to, innovation, weight, enabled }
    }

    pub fn new_from(connected_gene: &ConnectionGene) -> Self {
        ConnectionGene { 
            from: connected_gene.from, 
            to: connected_gene.to, 
            innovation: connected_gene.innovation, 
            weight: connected_gene.weight, 
            enabled: connected_gene.enabled,
        }
    }
}
