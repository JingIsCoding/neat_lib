use std::collections::{HashMap, HashSet};
use std::cmp::{max, min};
use rand::Rng;
use serde::{Deserialize, Serialize};

use uuid;
use super::gene::GeneType;
use super::innovation_number;
use super::attribute::FloatAttribute;
use super::{config::Config, gene::Gene, connection_gene::ConnectionGene};
use crate::network::errors::Errors;
use crate::network::activation::sigmod;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Genome {
    pub ref_id: uuid::Uuid,
    pub config: Config,
    pub fitness: f64,
    pub adjusted_fitness: f64,
    // index to gene
    pub genes: HashMap<usize, Gene>,
    // innovation number to connection
    pub connections: HashMap<usize, ConnectionGene>,
}

impl Genome {
    pub fn new(config: Config) -> Self {
        Genome { 
            config,
            ref_id: uuid::Uuid::new_v4(),
            fitness: 0.0,
            adjusted_fitness: 0.0,
            genes: HashMap::new(),
            connections: HashMap::new(),
        }
    }

    pub fn copy_from(genome: &Genome, keep_id: bool) -> Self {
        let id = if keep_id {
            genome.ref_id
        } else {
            uuid::Uuid::new_v4()
        };
        Genome { 
            config: genome.config,
            ref_id: id,
            fitness: genome.fitness,
            adjusted_fitness: genome.adjusted_fitness,
            genes: genome.genes.clone(),
            connections: genome.connections.clone(),
        }
    }

    pub fn is_same_species(parent1: &Genome, parent2: &Genome) -> bool {
        let mut matching = 0.0;
        let mut disjoint = 0.0;
        let mut excess = 0.0;
        let mut weight = 0.0;
        let mut low_max_innov = 0;
        let mut delta = 0.0;
        let max_innov_1 = parent1.connections.keys().max().unwrap_or(&(0 as usize));
        let max_innov_2 = parent2.connections.keys().max().unwrap_or(&(0 as usize));
        low_max_innov = min(*max_innov_1, *max_innov_2);
        let mut innovation_keys: HashSet<&usize> = HashSet::new();
        innovation_keys.extend(parent1.connections.keys().collect::<Vec<&usize>>());
        innovation_keys.extend(parent2.connections.keys().collect::<Vec<&usize>>());
        for innov in innovation_keys {
            if parent1.connections.contains_key(innov) && parent2.connections.contains_key(innov) {
                matching += 1.0;
                weight += (parent1.connections.get(innov).unwrap().weight - parent2.connections.get(innov).unwrap().weight).abs();
            } else {
                if *innov <= low_max_innov {
                    disjoint += 1.0;
                } else {
                    excess += 1.0;
                }
            }
        }
        let config = parent1.config;
        let n = matching + disjoint + excess;
        if n > 0.0 {
            delta = (config.excess_coefficent * excess + config.disjoint_coefficent * disjoint) / n + (config.weight_coefficent * weight) / matching
        }
        delta < config.compatibility_threshold
    }

    pub fn cross_over(parent1: &Genome, parent2: &Genome) -> Self {
        let mut rng = rand::thread_rng();
        let (best, second_best) = if parent1.fitness > parent2.fitness {
            (parent1, parent2)
        } else {
            (parent2, parent1)
        };
        let mut genome = Self::new(parent1.config);
        let mut innovation_keys: HashSet<&usize> = HashSet::new();
        innovation_keys.extend(parent1.connections.keys().collect::<Vec<&usize>>());
        innovation_keys.extend(parent2.connections.keys().collect::<Vec<&usize>>());
        for innov in innovation_keys {
            let conn = if best.connections.contains_key(innov) && second_best.connections.contains_key(innov) {
                let (best_node, second_node) = (best.connections.get(innov).unwrap(), second_best.connections.get(innov).unwrap());
                let mut conn = if rng.gen() {
                    best_node.clone()
                } else {
                    second_node.clone()
                };
                if best_node.enabled != second_node.enabled {
                    conn.enabled = rng.gen::<f32>() < 0.75;
                }
                conn
            } else if best.connections.contains_key(innov) {
                best.connections.get(innov).unwrap().clone()
            } else {
                second_best.connections.get(innov).unwrap().clone()
            };
            genome.connections.insert(conn.innovation, conn);
        }
        genome.generate_network();
        genome
    }

    pub fn generate_network(&mut self) {
        let mut rng = rand::thread_rng();
        let (inputs, outputs) = (self.config.inputs as usize, self.config.outputs as usize);
        self.genes.clear();
        // input
        for i in 0..inputs {
            self.genes.insert(i,Gene::new(self.config,i, GeneType::Input));
        }
        // Bias
        // self.genes.insert(self.config.inputs as usize, Gene::new(GeneType::Input));
        //
        // output
        for i in inputs..(inputs + outputs) {
            self.genes.insert(i,Gene::new(self.config,i, GeneType::Output));
        }

        //minimal network
        if self.connections.is_empty() {
            for from in 0..inputs {
                for to in inputs..(inputs + outputs) {
                    let ran_weight = rng.gen_range(-0.1..0.1);
                    self.add_connection(from, to, ran_weight)
                }
            }
        }
    }

    fn add_connection(&mut self, from: usize, to: usize, weight: f64) {
        let conn = ConnectionGene::new(self.config, from as usize, to as usize, innovation_number::next(), weight, true);
        self.connections.insert(conn.innovation, conn);
        self.genes.get_mut(&to).unwrap().add_incomming_conn(conn.innovation);
    }

    fn remove_connection(&mut self, innov: usize) {
        self.connections.remove(&innov);
        self.genes.iter_mut().for_each(|(_, gene)| {
            gene.remove_incoming_conn(&innov);
        })
    }

    pub fn evalute_network(&mut self, input_values: &Vec<f64>) -> Result<Vec<f64>, Errors> {
        if input_values.len() != self.config.inputs as usize {
            return Err(Errors::InputSizeNotMatch("evalute_network inputs number does not match".to_owned()));
        }

        for (_ ,gene) in &mut self.genes {
            gene.value.set_value(0.0)
        }

        let (inputs, outputs) = (self.config.inputs as usize, self.config.outputs as usize);
        let mut output_values = vec![0.0; self.config.outputs as usize];
        // put inputs into input nodes
        for index in 0..inputs {
            if let Some(gene) = self.genes.get_mut(&index) {
                gene.value.set_value(input_values[index]);
            }
        }

        // evaluate hidden layer first
        for i in (inputs + outputs)..self.genes.len() {
            if let Err(err) = Self::evalute_gene_at(&i, &mut self.genes, &self.connections) {
                return Err(err)
            }
        }

        // evaluate outputs
        for i in inputs..(inputs + outputs) {
            if let Err(err) = Self::evalute_gene_at(&i, &mut self.genes, &self.connections) {
                return Err(err)
            }
        }

        // fill in output values
        for i in inputs..(inputs + outputs) {
            output_values[i - inputs] = self.genes.get(&i).unwrap().value.value();
        }

        Ok(output_values)
    }

    fn evalute_gene_at(index: &usize, genes: &mut HashMap<usize, Gene>, connections: &HashMap<usize, ConnectionGene>) -> Result<(), Errors> {
        let mut aggregation = FloatAttribute::new(0.0);
        if let Some(gene) = genes.get(index) {
            for conn_index in &gene.incoming_conns {
                if let Some(conn) = connections.get(conn_index) {
                    if conn.enabled {
                        aggregation += genes[&conn.from].value * conn.weight;
                    }
                } else {
                    return Err(Errors::ConnectionGeneNotExists(*conn_index));
                }
            }
            if let Some(gene) = genes.get_mut(index){
                gene.value.set_value(sigmod((gene.bias + (gene.response * aggregation)).value()));
            }
        } else {
            return Err(Errors::GeneNotExists(*index));
        }
        Ok(())
    }

    pub fn mutate(&mut self) {
        let mut rng = rand::thread_rng();
        'single_structural_mutate: {
            if rng.gen_range(0.0..1.0) < self.config.conn_add_chance {
                self.mutate_add_conn();
                break 'single_structural_mutate
            }

            if rng.gen_range(0.0..1.0) < self.config.conn_delete_chance {
                self.mutate_remove_conn();
                break 'single_structural_mutate
            }

            if rng.gen_range(0.0..1.0) < self.config.node_add_chance {
                self.mutate_add_node();
                break 'single_structural_mutate
            }

            if rng.gen_range(0.0..1.0) < self.config.node_delete_chance {
                self.mutate_remove_node();
                break 'single_structural_mutate
            }
        }

        if rng.gen_range(0.0..1.0) < self.config.weight_mutation_chance{
            self.mutate_weight();
        }

        if rng.gen_range(0.0..1.0) < self.config.enabled_mutation_chance {
            self.mutate_enable();
        }
    }

    fn mutate_weight(&mut self) {
        self.genes.iter_mut().for_each(| (_, gene) | gene.mutate());
        self.connections.iter_mut().for_each(| (_, conn) | conn.mutate());
    }

    fn mutate_add_conn(&mut self) {
        let mut rng = rand::thread_rng();
        let (inputs, _, size) = (self.config.inputs as usize, self.config.outputs as usize, self.genes.len());
        let ran_node1 = rng.gen_range(0..size);
        let ran_node2 = rng.gen_range(inputs..size);
        if ran_node1 == ran_node2 {
            return
        }
        let (from, to) = (min(ran_node1, ran_node2), max(ran_node1, ran_node2));
        if self.try_connect(from, to) {
            return
        } 
    }

    fn mutate_remove_conn(&mut self) {
        let mut rng = rand::thread_rng();
        let innovs: Vec<&usize> = self.connections.keys().collect();
        if innovs.is_empty() {
            return
        }
        self.remove_connection(innovs[rng.gen_range(0..innovs.len())].clone());
    }

    fn mutate_add_node(&mut self) {
        let next_node = self.genes.len().clone();
        let mut from: usize = 0;
        let mut to: usize = 0;
        let mut weight: FloatAttribute = FloatAttribute::new(0.0);
        if let Some(conn) = self.find_random_connection() {
            conn.enabled = false;
            self.add_node_on_conn(conn.from.clone(), conn.to.clone(), conn.weight.value())
        } else {
            return
        }
    }

    fn add_node_on_conn(&mut self, from: usize, to: usize, weight: f64) {
        let next_node = self.genes.len().clone();
        self.genes.insert(next_node, Gene::new(self.config, next_node, GeneType::Hidden));
        self.add_connection(from, next_node, 1.0);
        self.add_connection(next_node, to, weight);
    }

    fn mutate_remove_node(&mut self) {
        let mut rng = rand::thread_rng();
        if self.genes.len() == (self.config.inputs + self.config.outputs) as usize {
            return;
        }
        let to_delete = rng.gen_range((self.config.inputs as usize + self.config.outputs as usize)..self.genes.len());
        self.remove_node(to_delete)
    }

    fn remove_node(&mut self, index: usize) {
        let mut to_delete_conns = vec![];
        for (innov, conn) in &self.connections {
            if conn.to == index || conn.from == index {
                to_delete_conns.push(innov.clone());
            }
        }
        for innov in to_delete_conns {
            self.remove_connection(innov);
        }
        self.genes.remove(&index);
    }

    fn try_connect(&mut self, from: usize, to: usize) -> bool {
        if from == to {
            return false;
        }
        if !self.genes.contains_key(&from) || !self.genes.contains_key(&to) {
            return false
        }
        let mut rng = rand::thread_rng();
        let gene = self.genes.get(&to).unwrap();
        for innov in &gene.incoming_conns {
            let conn = self.connections.get(innov).unwrap();
            if conn.from == from {
                return false;
            }
        }
        self.add_connection(from, to, rng.gen_range(-1.0..1.0));
        true
    }

    fn find_random_connection(&mut self) -> Option<&mut ConnectionGene> {
        let mut rng = rand::thread_rng();
        let keys = self.connections.keys().collect::<Vec<&usize>>();
        if keys.is_empty() {
            return None;
        }
        let randome_key = keys[rng.gen_range(0..keys.len())];
        if let Some(connection) = self.connections.get_mut(&randome_key.clone()) {
            return Some(connection);
        }
        None
    }

    fn mutate_enable(&mut self) {
        let mut rng = rand::thread_rng();
        let keys = self.connections.keys().collect::<Vec<&usize>>();
        if keys.is_empty() {
            return;
        }
        let randome_key = keys[rng.gen_range(0..keys.len())];
        if let Some(connection) = self.connections.get_mut(&randome_key.clone()) {
            connection.enabled = rng.gen_bool(0.5)
        }
    }

    pub fn debug(&self) {
        for (_ ,gene) in &self.genes {
            gene.debug();
        }
        for (_ ,conn) in &self.connections {
            conn.debug();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::config::Config;
    #[test]
    fn test_genrates_networks() {
        let config = Config::default();
        let mut genome = Genome::new(config);
        genome.generate_network();
        assert!(genome.genes.len() == 2);
        assert!(genome.connections.len() == 1);
    }

    #[test]
    fn test_evalute_minimal_network() {
        let config = Config::default();
        let mut genome = Genome::new(config);
        genome.generate_network();
        let inputs = vec![2.0, 1.0];
        let outputs = genome.evalute_network(&inputs);
        assert_eq!(outputs.unwrap().len(), 2);
    }

    #[test]
    fn test_evalute_network() {
        let config = Config::default();
        let mut genome = Genome::new(config);
        genome.connections = HashMap::new();
        genome.connections.insert(1, ConnectionGene::new(config, 0, 3, 1, -0.41132283, true));
    }

    #[test]
    fn test_minimal_network_crossover() {
        let mut config = Config::default();
        config.inputs = 2;
        config.outputs = 2;
        let parent1 = Genome::new(config);
        let parent2 = Genome::copy_from(&parent1, false);
        let child = Genome::cross_over(&parent1, &parent2);
        assert_eq!(child.genes.len(), 4);
        assert_eq!(child.connections.len(), 4);
    }

    #[test]
    fn test_network_crossover() {
        let mut config = Config::default();
        config.inputs = 2;
        config.outputs = 2;
        let mut parent1 = Genome::new(config);
        let mut parent2 = Genome::copy_from(&parent1, false);
        parent1.connections.insert(1, ConnectionGene::new(config,0, 1, 1, 0.1, true));
        parent1.connections.insert(2, ConnectionGene::new(config,0, 2, 2, 0.1, true));

        parent2.connections.insert(1, ConnectionGene::new(config,0, 1, 1, 0.2, true));
        parent2.connections.insert(2, ConnectionGene::new(config,0, 2, 2, 0.2, true));
        parent2.connections.insert(3, ConnectionGene::new(config,0, 3, 3, 0.2, true));

        let child = Genome::cross_over(&parent1, &parent2);
        child.debug();
        assert_eq!(child.genes.len(), 4);
        assert_eq!(child.connections.len(), 3);
    }

    #[test]
    fn test_try_connect() {
        let mut config = Config::default();
        config.inputs = 2;
        config.outputs = 2;
        let mut genome = Genome::new(config);
        genome.generate_network();
        assert_eq!(genome.genes.len(), 4);
        assert_eq!(genome.connections.len(), 4);

        genome.remove_connection(0);
        genome.debug();
        assert_eq!(genome.connections.len(), 3);
        let expected_incomming_conn: HashSet<usize> = vec![2 as usize].into_iter().clone().collect();
        assert_eq!(genome.genes.get(&2).unwrap().incoming_conns, expected_incomming_conn, "");

        assert!(!genome.try_connect(0, 0));
        assert!(!genome.try_connect(0, 3));
        assert!(genome.try_connect(0, 2));
    }
}
