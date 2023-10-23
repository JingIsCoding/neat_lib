use serde::{Deserialize, Serialize};
use crate::loader::save_load::{Loader, FileSaverLoader};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct Config {
    pub inputs: u32,
    pub outputs: u32,
    pub population: u32,
    pub conn_add_chance: f64,
    pub conn_delete_chance: f64,

    pub bias_init_mean: f64,
    pub bias_init_stdev: f64,
    pub bias_max_value: f64,
    pub bias_min_value: f64,
    pub bias_mutate_power: f64,
    pub bias_mutate_rate: f64,
    pub bias_replace_rate: f64,

    pub response_init_mean: f64,
    pub response_init_stdev: f64,
    pub response_max_value: f64,
    pub response_min_value: f64,
    pub response_mutate_power: f64,
    pub response_mutate_rate : f64,
    pub response_replace_rate: f64,

    pub weight_init_mean: f64,
    pub weight_init_stdev: f64,
    pub weight_min_value: f64,
    pub weight_max_value: f64,
    pub weight_mutation_step: f64,
    pub weight_perturbed_chance: f64,
    pub weight_mutation_chance: f64,
    pub weight_mutation_power: f64,
    pub weight_replace_chance: f64,

    pub node_add_chance: f64,
    pub node_delete_chance: f64,

    pub disable_mutation_chance: f64,
    pub enabled_mutation_chance: f64,
    pub crossover_chance: f64,

    pub compatibility_threshold: f64,
    pub excess_coefficent: f64,
    pub disjoint_coefficent: f64,
    pub weight_coefficent: f64,

    pub stale_species_threshold: u32,
    pub stale_population_threshold: u32,

    pub elitism: usize,
    pub survival_threshold: f64,
    pub min_species_size: usize,

    pub activation_func: String,
    pub aggregation_func: String,
}

impl Default for Config {
    fn default() -> Self {
        return Config { 
            inputs: 2, 
            outputs: 1, 
            population: 150, 
            conn_add_chance: 0.5, 
            conn_delete_chance: 0.5, 

            bias_init_mean: 0.0, 
            bias_init_stdev: 1.0, 
            bias_max_value: 30.0, 
            bias_min_value: -30.0, 
            bias_mutate_power: 0.5, 
            bias_mutate_rate: 0.7, 
            bias_replace_rate: 0.1, 

            response_init_mean: 1.0, 
            response_init_stdev: 0.0, 
            response_max_value: 30.0, 
            response_min_value: -30.0, 
            response_mutate_power: 0.0, 
            response_mutate_rate: 0.0, 
            response_replace_rate: 0.0, 

            weight_init_mean: 0.0, 
            weight_init_stdev: 1.0, 
            weight_max_value: 30.0, 
            weight_min_value: -30.0, 
            weight_mutation_step: 0.0, 
            weight_perturbed_chance: 0.0, 
            weight_mutation_chance: 0.5, 
            weight_mutation_power: 0.8, 
            weight_replace_chance: 0.1, 

            node_add_chance: 0.3, 
            node_delete_chance: 0.3, 

            disable_mutation_chance: 0.05, 
            enabled_mutation_chance: 0.05, 
            crossover_chance: 0.75, 
            compatibility_threshold: 1.0, 
            excess_coefficent: 1.0, 
            disjoint_coefficent: 1.0, 
            weight_coefficent: 0.5, 
            stale_species_threshold: 15, 
            stale_population_threshold: 20,

            elitism: 2,
            survival_threshold: 0.2,
            min_species_size: 2,

            activation_func: "sigmod".to_owned(),
            aggregation_func: "sum".to_owned(),
        }
    }
}

impl Config {
    pub fn new_from_path(path: &str) -> Self {
        let loader = FileSaverLoader::new(path);
        loader.load().expect("can not get config")
    }
}

#[cfg(test)]
mod tests {
    use super::Config;
    #[test]
    fn test_default() {
        let config = Config::default();
        assert_eq!(config.inputs, 2);
        assert_eq!(config.node_add_chance, 0.3);
    }

    #[test]
    fn test_new_from_path() {
        let config = Config::new_from_path("src/neat/test_config_1.json");
        assert_eq!(config.inputs, 2);
        assert_eq!(config.node_add_chance, 0.5);
        assert_eq!(config.activation_func, "sigmod");
    }
}
