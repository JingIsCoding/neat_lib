use serde::{Deserialize, Serialize};
use crate::{network::attribute::FloatAttribute, Config};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ConnectionGene {
    pub from: usize,
    pub to: usize,
    pub innovation: usize,
    pub weight: FloatAttribute,
    pub enabled: bool,
    config: Config
}

impl ConnectionGene {
    pub fn new(config: Config, from: usize, to: usize, innovation: usize, weight: f64, enabled: bool) -> Self {
        ConnectionGene { config, from, to, innovation, weight: FloatAttribute::new(weight), enabled }
    }

    pub fn mutate(&mut self) {
        let config = self.config;
        self.weight.mutate(config.weight_replace_chance, config.weight_mutation_chance, config.weight_mutation_power, config.weight_min_value, config.weight_max_value, config.weight_init_mean, config.weight_init_stdev);
    }

    pub fn debug(&self) {
        println!("ConnectionGene[ from:{:?} \tto:{:?} \tinnoovation:{:?} \tweight:{:?} \tenabled:{:?} ]", self.from, self.to, self.innovation, self.weight, self.enabled);
    }
}
