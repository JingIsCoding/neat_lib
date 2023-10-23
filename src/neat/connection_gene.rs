use serde::{Deserialize, Serialize};
use crate::{neat::attribute::FloatAttribute, Config};

#[derive(Clone, Serialize, Deserialize)]
pub struct ConnectionGene {
    pub from: usize,
    pub to: usize,
    pub innovation: usize,
    pub weight: FloatAttribute,
    pub enabled: bool,
    config: Config
}

impl ConnectionGene {
    pub fn new(config: &Config, from: usize, to: usize, innovation: usize, weight: f64, enabled: bool) -> Self {
        ConnectionGene { config: config.clone(), from, to, innovation, weight: FloatAttribute::new(weight), enabled }
    }

    pub fn mutate(&mut self) {
        self.weight.mutate(self.config.weight_replace_chance, self.config.weight_mutation_chance, self.config.weight_mutation_power, self.config.weight_min_value, self.config.weight_max_value, self.config.weight_init_mean, self.config.weight_init_stdev);
    }
}

impl std::fmt::Debug for ConnectionGene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ConnectionGene[ from:{:?} \tto:{:?} \tinnoovation:{:?} \tweight:{:?} \tenabled:{:?} ]", self.from, self.to, self.innovation, self.weight, self.enabled)
    }
}

