use std::collections::{HashMap, HashSet};
use std::cmp::{max, min};
use rand::Rng;
use serde::{Deserialize, Serialize};

use uuid;
use super::gene::GeneType;
use super::innovation_number;
use super::{config::Config, gene::Gene, connection_gene::ConnectionGene};

const MAX_RETRIES: u8 = 3;

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy, Serialize, Deserialize)]
pub enum Mutations {
    Steps,
    PerturbedChance,
    WeightChance,
    WeightMutationChance,
    NodeMutationChange,
    ConnectionMutationChance,
    BiasMutationChance,
    DisableMutationChance,
    EnableMutationChance
}

type MutationMap = HashMap<Mutations, f64>;

#[derive(Debug, Serialize, Deserialize)]
pub struct Genome {
    pub ref_id: uuid::Uuid,
    pub config: Config,
    pub mutations: MutationMap,
    pub fitness: f64,
    pub adjusted_fitness: f64,
    pub genes: HashMap<usize, Gene>,
    pub connections: Vec<ConnectionGene>,
}

impl Genome {
    pub fn new(config: Config) -> Self {
        let mut mutations = MutationMap::new();
        mutations.insert(Mutations::Steps, config.weight_mutation_step);
        mutations.insert(Mutations::PerturbedChance, config.perturbed_chance);
        mutations.insert(Mutations::WeightChance, config.weight_chance);
        mutations.insert(Mutations::WeightMutationChance, config.weight_mutation_chance);
        mutations.insert(Mutations::NodeMutationChange, config.node_mutation_chance);
        mutations.insert(Mutations::ConnectionMutationChance, config.connection_mutation_chance);
        mutations.insert(Mutations::BiasMutationChance, config.bias_mutation_chance);
        mutations.insert(Mutations::DisableMutationChance, config.disable_mutation_chance);
        mutations.insert(Mutations::EnableMutationChance, config.enable_mutation_chance);

        Genome { 
            config,
            mutations,
            ref_id: uuid::Uuid::new_v4(),
            fitness: 0.0,
            adjusted_fitness: 0.0,
            genes: HashMap::new(),
            connections: vec![],
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
            mutations: genome.mutations.clone(),
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
        let mut max_innov_1 = 0;
        let mut max_innov_2 = 0;
        let mut low_max_innov = 0;
        let mut delta = 0.0;
        let mut parent1_map = HashMap::new();
        let mut parent2_map = HashMap::new();

        for connection in parent1.connections.iter() {
            let innov = connection.innovation;
            parent1_map.insert(innov, connection);
            max_innov_1 = max(max_innov_1, innov);
        }

        for connection in parent2.connections.iter() {
            let innov = connection.innovation;
            parent2_map.insert(innov, connection);
            max_innov_2 = max(max_innov_2, innov);
        }

        if parent1_map.is_empty() || parent2_map.is_empty() {
            low_max_innov = 0;
        } else {
            low_max_innov = min(max_innov_1, max_innov_2);
        }

        let innovation_keys1= parent1_map.keys().fold(HashSet::new(), Self::agg);
        let innovation_keys2= parent2_map.keys().fold(HashSet::new(), Self::agg);
        let innovation_keys = innovation_keys1.union(&innovation_keys2).clone().collect::<HashSet<&usize>>();

        for innov in innovation_keys {
            if parent1_map.contains_key(innov) && parent2_map.contains_key(innov) {
                matching += 1.0;
                weight += (parent1_map.get(innov).unwrap().weight - parent2_map.get(innov).unwrap().weight).abs();
            } else {
                if *innov <= low_max_innov {
                    disjoint += 1.0;
                } else {
                    excess += 1.0;
                }
            }
        }

        let config = parent1.config;
        let n = matching+disjoint+excess;
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
        let mut max_best = 0;
        let mut max_second_best = 0;
        let mut best_map = HashMap::new();
        let mut second_best_map = HashMap::new();

        for connection in best.connections.iter() {
            let innov = connection.innovation;
            best_map.insert(innov, connection);
            max_best = max(max_best, innov);
        }

        for connection in second_best.connections.iter() {
            let innov = connection.innovation;
            second_best_map.insert(innov, connection);
            max_second_best = max(max_second_best, innov);
        }

        let innovation_keys1= best_map.keys().fold(HashSet::new(), Self::agg);
        let innovation_keys2= second_best_map.keys().fold(HashSet::new(), Self::agg);
        let innovation_keys = innovation_keys1.union(&innovation_keys2).collect::<HashSet<&usize>>();
        for innov in innovation_keys {
            let conn = if best_map.contains_key(innov) && second_best_map.contains_key(innov) {
                let (best_node, second_node) = (best_map.get(innov).unwrap(), second_best_map.get(innov).unwrap());
                let mut conn = if rng.gen() {
                    ConnectionGene::new_from(&best_node)
                } else {
                    ConnectionGene::new_from(&second_node)
                };
                if best_node.enabled != second_node.enabled {
                    conn.enabled = rng.gen::<f32>() < 0.75;
                }
                conn
            } else  {
                if best_map.contains_key(innov) {
                    ConnectionGene::new_from(best_map.get(innov).unwrap())
                } else {
                    ConnectionGene::new_from(second_best_map.get(innov).unwrap())
                }
            }; 
            genome.connections.push(conn);
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
            self.genes.insert(i,Gene::new(GeneType::Input));
        }
        // Bias
        // self.genes.insert(self.config.inputs as usize, Gene::new(GeneType::Input));
        //
        // output
        for i in inputs..(inputs + outputs) {
            self.genes.insert(i,Gene::new(GeneType::Output));
        }

        //minimal network
        if self.connections.is_empty() {
            for from in 0..inputs {
                for to in inputs..(inputs + outputs) {
                    let ran_weight = rng.gen_range(-0.1..0.1);
                    let conn = ConnectionGene::new(from as usize, to as usize, innovation_number::next(), ran_weight, true);
                    self.connections.push(conn);
                }
            }
        }

        for (index, conn) in self.connections.iter().enumerate() {
            if !self.genes.contains_key(&conn.from) {
                self.genes.insert(conn.from, Gene::new(GeneType::Hidden));
            }
            if !self.genes.contains_key(&conn.to) {
                self.genes.insert(conn.to, Gene::new(GeneType::Hidden));
            }
            self.genes.get_mut(&conn.to).unwrap().add_incomming_conn(index);
        }

    }

    pub fn evalute_network(&mut self, input_values: &Vec<f64>) -> Vec<f64> {
        assert!(input_values.len() == self.config.inputs as usize);
        for (_ ,gene) in &mut self.genes {
            gene.value = 0.0;
        }
        let (inputs, outputs) = (self.config.inputs as usize, self.config.outputs as usize);
        let mut output_values = vec![0.0; self.config.outputs as usize];
        // put inputs into input nodes
        for index in 0..inputs {
            if let Some(gene) = self.genes.get_mut(&index) {
                gene.value = input_values[index];
            }
        }
        let cloned_genes = self.genes.clone();
        // evaluate hidden layer first
        for i in (inputs + outputs)..cloned_genes.len() {
            let mut value = 0.0;
            if let Some(gene) = cloned_genes.get(&i) {
                for conn_index in &gene.incoming_conns {
                    let conn = self.connections[*conn_index];
                    if conn.enabled {
                        value += self.genes[&conn.from].value * conn.weight;
                    }
                }
                if let Some(gene) = self.genes.get_mut(&i){
                    gene.value = Self::sigmod(value);
                }
            }
        }

        // evaluate outputs
        for i in inputs..(inputs + outputs) {
            let mut value = 0.0;
            if let Some(gene) = cloned_genes.get(&i) {
                for conn_index in &gene.incoming_conns {
                    let conn = self.connections[*conn_index];
                    if conn.enabled {
                        value += self.genes[&conn.from].value * conn.weight;
                    }
                }
                if let Some(gene) = self.genes.get_mut(&i){
                    gene.value = Self::sigmod(value);
                }
            }
        }

        for i in inputs..(inputs + outputs) {
            output_values[i - inputs] = self.genes.get(&i).unwrap().value;
        }
        output_values
    }

    pub fn mutate(&mut self) {
        let mut rng = rand::thread_rng();
        for (_, mutation_value) in self.mutations.iter_mut() {
            if rng.gen_bool(0.5) {
                *mutation_value = 0.95 * *mutation_value;
            } else {
                *mutation_value = 1.05263 * *mutation_value;
            }
        }
        if rng.gen_range(0.0..1.0) < *self.mutations.get(&Mutations::WeightMutationChance).unwrap() {
            self.mutate_weight();
        }
        if rng.gen_range(0.0..1.0) < *self.mutations.get(&Mutations::EnableMutationChance).unwrap() {
            self.mutate_enable();
        }

        if rng.gen_range(0.0..1.0) < *self.mutations.get(&Mutations::DisableMutationChance).unwrap() {
            self.mutate_disable();
        }

        if rng.gen_range(0.0..1.0) < *self.mutations.get(&Mutations::ConnectionMutationChance).unwrap() {
            self.mutate_add_conn();
        }

        if rng.gen_range(0.0..1.0) < *self.mutations.get(&Mutations::NodeMutationChange).unwrap() {
            self.mutate_add_node();
        }
    }

    fn mutate_weight(&mut self) {
        let mut rng = rand::thread_rng();
        for conn in self.connections.iter_mut() {
            if rng.gen_range(0.0..1.0) < *self.mutations.get(&Mutations::WeightChance).unwrap() {
                if rng.gen_range(0.0..1.0) < *self.mutations.get(&Mutations::PerturbedChance).unwrap() {
                    conn.weight += rng.gen_range(-2.0..2.0) * self.mutations.get(&Mutations::Steps).unwrap();
                } else {
                    conn.weight = rng.gen_range(-4.0..4.0);
                }
            }
        }
    }

    fn mutate_add_conn(&mut self) {
        let mut rng = rand::thread_rng();
        let (inputs, _, size) = (self.config.inputs as usize, self.config.outputs as usize, self.genes.len());
        for _ in 0..MAX_RETRIES {
            let ran_node1 = rng.gen_range(0..size);
            let ran_node2 = rng.gen_range(inputs..size);
            if ran_node1 == ran_node2 {
                continue
            }
            let (from, to) = (min(ran_node1, ran_node2), max(ran_node1, ran_node2));
            if self.try_connect(from, to) {
                return
            } 
        }
    }

    fn try_connect(&mut self, from: usize, to: usize) -> bool {
        let mut rng = rand::thread_rng();
        if let Some(gene) = self.genes.get(&to) {
            for index in &gene.incoming_conns {
                let conn = self.connections[*index];
                if conn.from == from {
                    return false
                }
            }
            self.connections.push(ConnectionGene::new(from, to, innovation_number::next(), rng.gen_range(-1.0..1.0), true));
            return true;
        }
        false
    }

    fn mutate_add_node(&mut self) {
        if let Some(index) = self.find_random_enabled_node_index() {
            let next_node = self.genes.len();
            self.connections[index].enabled = false;
            let conn = self.connections[index];
            self.connections.push(ConnectionGene::new(conn.from, next_node, innovation_number::next(), 1.0, true));
            self.connections.push(ConnectionGene::new(next_node, conn.to, innovation_number::next(), conn.weight, true));
            self.generate_network();
        }
    }

    fn find_random_enabled_node_index(&mut self) -> Option<usize> {
        let mut rng = rand::thread_rng();
        let size = self.connections.len();
        if size == 0 {
            return None;
        }
        for _ in 0..self.config.hidden_nodes {
            let index = rng.gen_range(0..size);
            if self.connections[index].enabled {
                return Some(index);
            }
        }
        None
    }

    fn mutate_enable(&mut self) {
        let mut rng = rand::thread_rng();
        let mut retries = 0;
        if self.connections.is_empty() {
            return
        }
        let size = self.connections.len();
        while retries < MAX_RETRIES {
            let index = rng.gen_range(0..size);
            let conn = self.connections[index];
            if !conn.enabled {
                self.connections[index].enabled = true;
                return
            }
            retries += 1;
        }
    }

    fn mutate_disable(&mut self) {
        let mut rng = rand::thread_rng();
        let mut retries = 0;
        if self.connections.is_empty() {
            return
        }
        let size = self.connections.len();
        while retries < MAX_RETRIES {
            let index = rng.gen_range(0..size);
            let conn = self.connections[index];
            if conn.enabled {
                self.connections[index].enabled = false;
                return
            }
            retries += 1;
        }
    }

    fn sigmod(value: f64) -> f64 {
        1.0 / (1.0 + (-4.9 * value).exp())
    }

    fn agg(mut acc: HashSet<usize>, each: &usize) -> HashSet<usize> {
        acc.insert(*each);
        acc
    }

    pub fn debug(&self) {
        println!("gene {:?} connections {:?}", self.genes.len(), self.connections.len());
        let mut connections = self.connections.clone();
        connections.sort_by(|conn1, conn2| {
            if conn1.from < conn2.from {
                return std::cmp::Ordering::Less
            }
            return std::cmp::Ordering::Greater
        });
        for conn in connections {
            println!("{:?}", conn);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::config::Config;
    #[test]
    fn test_genrates_networks() {
        let config = Config::new(1, 1,10);
        let mut genome = Genome::new(config);
        genome.generate_network();
        assert!(genome.genes.len() == 2);
        assert!(genome.connections.len() == 1);
    }

    #[test]
    fn test_evalute_minimal_network() {
        let config = Config::new(2, 2,10);
        let mut genome = Genome::new(config);
        genome.generate_network();
        let inputs = vec![2.0, 1.0];
        let outputs = genome.evalute_network(&inputs);

        assert_eq!(outputs.len(), 2);
    }

    #[test]
    fn test_evalute_network() {
        let config = Config::new(2, 1,10);
        let mut genome = Genome::new(config);
        genome.connections = vec![
            ConnectionGene::new(0, 3, 1, -0.41132283, true),
            ConnectionGene::new(0, 4, 2, 0.7442026, true),
            ConnectionGene::new(0, 5, 2, 3.9696338, true),
            ConnectionGene::new(0, 6, 2, -3.487875, true),

            ConnectionGene::new(1, 3, 3, -0.48198652, true),
            ConnectionGene::new(1, 4, 4, -0.253062, true),
            ConnectionGene::new(1, 5, 4, -3.969754, true),
            ConnectionGene::new(1, 6, 4, 3.4878979, true),

            ConnectionGene::new(3, 2, 5, 0.4171474, true),
            ConnectionGene::new(4, 2, 6, -0.28525555, true),
            ConnectionGene::new(5, 2, 6, 3.68314, true),
            ConnectionGene::new(6, 2, 6, 4.206077, true),
        ];
        genome.generate_network();
        let inputs = vec![1.0, 1.0];
        let outputs = genome.evalute_network(&inputs);
        println!("{:?} {:?} > {:?}", 1, 1, outputs[0]);

        let inputs = vec![0.0, 0.0];
        let outputs = genome.evalute_network(&inputs);
        println!("{:?} {:?} > {:?}", 0, 0, outputs[0]);

        let inputs = vec![1.0, 0.0];
        let outputs = genome.evalute_network(&inputs);
        println!("{:?} {:?} > {:?}", 1, 0, outputs[0]);

        let inputs = vec![0.0, 0.0];
        let outputs = genome.evalute_network(&inputs);
        println!("{:?} {:?} > {:?}", 1, 0, outputs[0]);

        assert_eq!(outputs[0], 0.0);
    }

    #[test]
    fn test_minimal_network_crossover() {
        let config = Config::new(2, 2,10);
        let parent1 = Genome::new(config);
        let parent2 = Genome::copy_from(&parent1, false);
        let child = Genome::cross_over(&parent1, &parent2);
        assert_eq!(child.genes.len(), 4);
        assert_eq!(child.connections.len(), 4);
    }

    #[test]
    fn test_network_crossover() {
        let config = Config::new(2, 2,10);
        let mut parent1 = Genome::new(config);
        let mut parent2 = Genome::copy_from(&parent1, false);
        parent1.connections = vec![ConnectionGene::new(0, 2, 1, 0.5, true)];
        parent2.connections = vec![ConnectionGene::new(1, 2, 2, 0.5, true)];
        let child = Genome::cross_over(&parent1, &parent2);
        assert_eq!(child.genes.len(), 4);
        assert_eq!(child.connections.len(), 2);
    }

    #[test]
    fn test_try_connect() {
        let config = Config::new(2, 2,10);
        let mut genome = Genome::new(config);
        genome.generate_network();
        genome.genes.insert(4, Gene::new(GeneType::Hidden));
        assert_eq!(genome.try_connect(0, 3), false);
        assert_eq!(genome.try_connect(0, 4), true);
        assert_eq!(genome.connections.len(), 5);
        assert_eq!(genome.try_connect(0, 5), false);
    }
}
