use std::collections::{HashMap, HashSet};
use std::cmp::min;
use rand::Rng;
use serde::{Deserialize, Serialize};

use uuid;
use super::gene::GeneType;
use super::innovation_number;
use super::attribute::FloatAttribute;
use super::{config::Config, gene::Gene, connection_gene::ConnectionGene};
use crate::neat::errors::Errors;
use crate::neat::activation::sigmod;

#[derive(Clone, Serialize, Deserialize)]
pub struct Genome {
    pub ref_id: uuid::Uuid,
    pub config: Config,
    pub fitness: f64,
    // index to gene
    pub genes: HashMap<usize, Gene>,
    // innovation number to connection
    pub connections: HashMap<usize, ConnectionGene>,
}

impl Genome {
    pub fn new(config: &Config) -> Self {
        let mut genes = HashMap::new();
        for i in 0..config.inputs {
            let gene = Gene::new(config, i as usize, GeneType::Input);
            genes.insert(gene.key, gene);
        }
        for i in config.inputs..(config.inputs + config.outputs) {
            let gene = Gene::new(config, i as usize, GeneType::Output);
            genes.insert(gene.key, gene);
        }
        Genome { 
            config: config.clone(),
            ref_id: uuid::Uuid::new_v4(),
            fitness: 0.0,
            genes,
            connections: HashMap::new(),
        }
    }

    pub fn new_from(genome: &Genome, keep_id: bool) -> Self {
        let id = if keep_id {
            genome.ref_id
        } else {
            uuid::Uuid::new_v4()
        };
        Genome { 
            config: genome.config.clone(),
            ref_id: id,
            fitness: 0.0,
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
        let config = &parent1.config;
        let n = matching + disjoint + excess;
        if n > 0.0 {
            delta = (config.excess_coefficent * excess + config.disjoint_coefficent * disjoint) / n + (config.weight_coefficent * weight) / matching
        }
        delta < config.compatibility_threshold
    }

    pub fn cross_over(parent1: &Genome, parent2: &Genome) -> Self {
        let (best, second_best) = if parent1.fitness > parent2.fitness {
            (parent1, parent2)
        } else {
            (parent2, parent1)
        };
        let mut child = Self::new(&parent1.config);
        Self::cross_over_genes(&mut child, best, second_best);
        Self::cross_over_connections(&mut child, best, second_best);
        child.minimal_network();
        child.check_connections();
        child
    }

    fn cross_over_genes(child: &mut Genome, best: &Genome, second_best: &Genome) {
        let mut rng = rand::thread_rng();
        let mut gene_keys: HashSet<&usize> = HashSet::new();
        gene_keys.extend(best.genes.keys().collect::<Vec<&usize>>());
        gene_keys.extend(second_best.genes.keys().collect::<Vec<&usize>>());
        for key in gene_keys {
            let gene = if best.genes.contains_key(key) && second_best.genes.contains_key(key) {
                if rng.gen() {
                    best.genes.get(key).unwrap().clone()
                } else {
                    second_best.genes.get(key).unwrap().clone()
                }
            } else if best.genes.contains_key(key) {
                best.genes.get(key).unwrap().clone()
            } else {
                second_best.genes.get(key).unwrap().clone()
            };
            child.genes.insert(*key, gene);
        }
    }

    fn cross_over_connections(child: &mut Genome, best: &Genome, second_best: &Genome) {
        let mut rng = rand::thread_rng();
        let mut innovation_keys: HashSet<&usize> = HashSet::new();
        innovation_keys.extend(best.connections.keys().collect::<Vec<&usize>>());
        innovation_keys.extend(second_best.connections.keys().collect::<Vec<&usize>>());
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
            child.add_connection(conn);
        }
    }

    pub fn minimal_network(&mut self) {
        //minimal network
        if self.connections.is_empty() {
            let mut rng = rand::thread_rng();
            let inputs: Vec<usize> = self.genes.iter().filter(|(_, gene)| { gene.gene_type == GeneType::Input }).map(| (key, _) | { key.clone() }).collect();
            let outputs: Vec<usize> = self.genes.iter().filter(|(_, gene)| { gene.gene_type == GeneType::Output }).map(| (key, _) | { key.clone() }).collect();
            for from in &inputs {
                for to in outputs.iter() {
                    let ran_weight = rng.gen_range(-0.1..0.1);
                    self.add_new_connection(*from, *to, ran_weight)
                }
            }
        }
    }

    pub(crate) fn add_connection(&mut self, conn: ConnectionGene) {
        self.genes.get_mut(&conn.to).unwrap().add_incomming_conn(conn.innovation);
        self.connections.insert(conn.innovation, conn);
    }

    fn add_new_connection(&mut self, from: usize, to: usize, weight: f64) {
        let conn = ConnectionGene::new(&self.config, from as usize, to as usize, innovation_number::next(), weight, true);
        self.add_connection(conn);
    }

    fn remove_connection(&mut self, innov: usize) {
        self.connections.remove(&innov);
        self.genes.iter_mut().for_each(|(_, gene)| {
            gene.remove_incoming_conn(&innov);
        })
    }

    pub fn mutate(&mut self) {
        let mut rng = rand::thread_rng();
        'single_structural_mutate: {
            if rng.gen_range(0.0..1.0) < self.config.conn_add_chance {
                self.mutate_add_conn();
                self.check_connections();
                break 'single_structural_mutate
            }

            if rng.gen_range(0.0..1.0) < self.config.conn_delete_chance {
                self.mutate_remove_conn();
                self.check_connections();
                break 'single_structural_mutate
            }

            if rng.gen_range(0.0..1.0) < self.config.node_add_chance {
                self.mutate_add_node();
                self.check_connections();
                break 'single_structural_mutate
            }

            if rng.gen_range(0.0..1.0) < self.config.node_delete_chance {
                self.mutate_remove_node();
                self.check_connections();
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
        let inputs: Vec<&usize> = self.genes.keys().collect();
        let outputs: Vec<(&usize, &Gene)> = self.genes.iter().filter(| (_, gene) | {
            gene.gene_type != GeneType::Input
        }).collect();
        let ran_index_1 = inputs[rng.gen_range(0..inputs.len())];
        let ran_index_2 = outputs[rng.gen_range(0..outputs.len())].0;
        if ran_index_1 == ran_index_2 {
            return
        }
        let (from, to) = (ran_index_1.min(ran_index_2), ran_index_1.max(ran_index_2));
        if self.connections.iter_mut().any(| (_ , conn) | {
            if conn.from == *from && conn.to == *to {
                conn.enabled = true;
                return true;
            }
            return false;
        }) {
            return
        }
        if self.try_connect(*from, *to) {
            return
        } 
    }

    fn mutate_remove_conn(&mut self) {
        let innovs: Vec<&usize> = self.connections.keys().collect();
        if innovs.is_empty() {
            return
        }
        let mut rng = rand::thread_rng();
        self.remove_connection(innovs[rng.gen_range(0..innovs.len())].clone());
    }

    fn mutate_add_node(&mut self) {
        if let Some(innov) = self.find_random_connection() {
            let conn = self.connections.get_mut(&innov).unwrap();
            conn.enabled = false;
            let from = conn.from.clone();
            let to = conn.to.clone();
            let weight = conn.weight.value().clone();
            self.add_node_on_conn(from, to, weight);
        }     
    }

    fn add_node_on_conn(&mut self, from: usize, to: usize, weight: f64) {
        let next_node = *self.genes.keys().max().unwrap() + 1;
        self.genes.insert(next_node, Gene::new(&self.config, next_node, GeneType::Hidden));
        self.add_new_connection(from, next_node, 1.0);
        self.add_new_connection(next_node, to, weight);
    }

    fn mutate_remove_node(&mut self) {
        let possible_nodes: Vec<(&usize, &Gene)> = self.genes.iter().filter(|(_, gene)| { gene.gene_type == GeneType::Hidden }).collect();
        if possible_nodes.is_empty() {
            return
        }
        let mut rng = rand::thread_rng();
        let to_delete = possible_nodes[rng.gen_range(0..possible_nodes.len())].0;
        self.remove_node(*to_delete)
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
        self.add_new_connection(from, to, rng.gen_range(-1.0..1.0));
        true
    }

    fn find_random_connection(&mut self) -> Option<usize> {
        let mut rng = rand::thread_rng();
        let keys = self.connections.keys().collect::<Vec<&usize>>();
        if keys.is_empty() {
            return None;
        }
        let randome_key = keys[rng.gen_range(0..keys.len())];
        if self.connections.contains_key(randome_key) {
            return Some(randome_key.clone());
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

    pub fn check_connections(&self) {
        for (_, conn) in &self.connections {
            if !self.genes.contains_key(&conn.from) {
                panic!("from key not exists >> {:?}\n{:?}", conn.innovation, self);
            } 
            if !self.genes.contains_key(&conn.to) {
                panic!("to key not exists >> {:?}\n{:?}", conn.innovation, self);
            } 
            if !self.genes.get(&conn.to).unwrap().incoming_conns.contains(&conn.innovation) {
                panic!("gene does not have incoming conn >> {:?}\n{:?}", conn.innovation, self);
            }
        }
    }
}

impl std::fmt::Debug for Genome {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut debug_content = String::new();

        debug_content.push_str(format!("id:\t{:?}\nfitness:\t{:?}\n", self.ref_id, self.fitness).as_str());

        let mut genes = self.genes.values().collect::<Vec<&Gene>>();
        genes.sort_by_key(| gene | gene.key);
        genes.iter().for_each(|gene| {
            debug_content.push_str(format!("{:?}\n", gene).as_str());
        });
        let mut connections = self.connections.values().collect::<Vec<&ConnectionGene>>();
        connections.sort_by_key(| gene | gene.innovation);
        connections.iter().for_each(|conn| {
            debug_content.push_str(format!("{:?}\n", conn).as_str());
        });
        write!(f, "{}", debug_content)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::config::Config;
    #[test]
    fn test_genrates_networks() {
        let mut config = Config::default();
        config.inputs = 1;
        config.outputs = 1;
        let mut genome = Genome::new(&config);
        genome.minimal_network();
        assert!(genome.genes.len() == 2);
        assert!(genome.connections.len() == 1);
    }

    #[test]
    fn test_evalute_network() {
        let config = Config::default();
        let mut genome = Genome::new(&config);
        genome.connections = HashMap::new();
        genome.connections.insert(1, ConnectionGene::new(&config, 0, 3, 1, -0.41132283, true));
    }

    #[test]
    fn test_minimal_network_crossover() {
        let mut config = Config::default();
        config.inputs = 2;
        config.outputs = 2;
        let parent1 = Genome::new(&config);
        let parent2 = Genome::new_from(&parent1, false);
        let child = Genome::cross_over(&parent1, &parent2);
        assert_eq!(child.genes.len(), 4);
        assert_eq!(child.connections.len(), 4);
    }

    #[test]
    fn test_network_crossover() {
        innovation_number::reset();
        let mut config = Config::default();
        config.inputs = 2;
        config.outputs = 2;
        let mut parent1 = Genome::new(&config);
        let mut parent2 = Genome::new_from(&parent1, false);
        parent1.connections.insert(1, ConnectionGene::new(&config,0, 1, 1, 0.1, true));
        parent1.connections.insert(2, ConnectionGene::new(&config,0, 2, 2, 0.1, true));

        parent2.connections.insert(1, ConnectionGene::new(&config,0, 1, 1, 0.2, true));
        parent2.connections.insert(2, ConnectionGene::new(&config,0, 2, 2, 0.2, true));
        parent2.connections.insert(3, ConnectionGene::new(&config,0, 3, 3, 0.2, true));

        let child = Genome::cross_over(&parent1, &parent2);
        assert_eq!(child.genes.len(), 4);
        assert_eq!(child.connections.len(), 3);
    }

    #[test]
    fn test_try_connect() {
        innovation_number::reset();
        let mut config = Config::default();
        config.inputs = 2;
        config.outputs = 2;
        let mut genome = Genome::new(&config);
        assert!(!genome.try_connect(0, 0));
        assert!(genome.try_connect(0, 2));
        assert_eq!(genome.connections.len(), 1);
        assert_eq!(genome.connections.get(&1).unwrap().from, 0);
        assert_eq!(genome.connections.get(&1).unwrap().to, 2);
        assert_eq!(genome.genes.get(&(2 as usize)).unwrap().incoming_conns.len(), 1);
        assert!(genome.genes.get(&(2 as usize)).unwrap().incoming_conns.contains(&1));
    }
}
