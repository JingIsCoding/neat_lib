use serde::{Deserialize, Serialize};
use crate::loader::save_load::{Loader, FileSaverLoader};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Config {
    pub inputs: u32,
    pub outputs: u32,
    pub population: u32,
    pub hidden_nodes: u32,
    pub weight_mutation_step: f64,
    pub perturbed_chance: f64,
    pub weight_chance: f64,
    pub weight_mutation_chance: f64,
    pub node_mutation_chance: f64,
    pub connection_mutation_chance: f64,
    pub bias_mutation_chance: f64,
    pub disable_mutation_chance: f64,
    pub enable_mutation_chance: f64,
    pub remove_connection_mutataion_chance: f64,
    pub crossover_chance: f64,

    pub compatibility_threshold: f64,
    pub excess_coefficent: f64,
    pub disjoint_coefficent: f64,
    pub weight_coefficent: f64,

    pub stale_species_threshold: u32,
    pub stale_population_threshold: u32,
}

impl Config {
    pub fn new_from_path(path: &str) -> Self {
        let loader = FileSaverLoader::new(path);
        loader.load().expect("can not get config")
    }

    pub fn new(inputs: u32, outputs: u32, population: u32) -> Self {
        Config { 
            inputs,
            outputs,
            population,
            hidden_nodes: 1000,
            weight_mutation_step: 0.1,
            perturbed_chance: 0.9,
            weight_chance: 0.3,
            weight_mutation_chance: 0.9,
            node_mutation_chance: 0.03,
            connection_mutation_chance: 0.05,
            bias_mutation_chance: 0.15,
            disable_mutation_chance: 0.1,
            enable_mutation_chance: 0.2,
            remove_connection_mutataion_chance: 0.05,
            crossover_chance: 0.75,

            compatibility_threshold: 0.5,
            excess_coefficent: 2.0,
            disjoint_coefficent: 2.0,
            weight_coefficent: 0.4,

            stale_species_threshold: 15,
            stale_population_threshold: 20,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_config_new() {
        Config::new_from_path("./config.json");
    }
}
