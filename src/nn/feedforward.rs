use std::collections::HashMap;
use std::collections::HashSet;

use crate::ConnectionGene;
use super::super::neat::genome::*;
use super::super::neat::gene::*;
use super::super::neat::aggregation::*;
use super::super::neat::activation::*;
use super::super::neat::config::*;
use super::layer::*;
use super::Network;

#[derive(Debug)]
struct NodeEval {
    key: usize,
    bias: f64,
    response: f64,
    incomings: Vec<(usize, f64)> // node key to weight
}

pub struct FeedForwardNetwork {
    node_evals: Vec<NodeEval>,
    layers: Vec<Layer>,
    values: HashMap<usize, f64>,
    inputs: Vec<usize>,
    outputs: Vec<usize>,
    aggregation_func: AggregationFunc,
    activation_func: ActivationFunc,
    genome: Genome
}

impl std::fmt::Debug for FeedForwardNetwork {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "layers: {:?} \n values are: {:?} \n evals are: {:?}\n {:?}", self.layers, self.values, self.node_evals, self.genome)
    }
}

impl FeedForwardNetwork {
    pub fn new(config: &Config, genome: &Genome) -> Self {
        let conns: Vec<&ConnectionGene> = genome.connections.values().filter(| conn | conn.enabled ).collect();
        let mut inputs: Vec<usize> = genome.genes.values()
            .filter(| gene | gene.gene_type == GeneType::Input)
            .map(| gene | gene.key).collect();
        let mut outputs: Vec<usize> = genome.genes.values()
            .filter(| gene | gene.gene_type == GeneType::Output)
            .map(| gene | gene.key).collect();
        inputs.sort();
        outputs.sort();

        let layers = Self::create_layers(&inputs, &outputs, conns);
        let mut node_evals = vec![];
        for layer in &layers {
            for node in &layer.keys {
                let mut incomings = vec![];
                let node = genome.genes.get(node).unwrap();
                for incoming in &node.incoming_conns {
                    let conn = genome.connections.get(incoming).unwrap();
                    if conn.enabled {
                        incomings.push((conn.from, conn.weight.value()));
                    }
                }
                node_evals.push(NodeEval{ key: node.key,  bias: node.bias.value(), response: node.bias.value(), incomings });
            }
        }
        FeedForwardNetwork { 
            genome: genome.clone(),
            inputs, 
            outputs, 
            node_evals, 
            values: HashMap::new(), 
            layers,
            aggregation_func: get_aggregation(&config.aggregation_func), 
            activation_func: get_activation(&config.activation_func)
        }
    }

    fn create_layers(inputs: &Vec<usize>, outputs: &Vec<usize>, conns: Vec<&ConnectionGene>) -> Vec<Layer> {
        let required = Self::required_for_outputs(inputs, outputs, conns.clone());
        let mut layers = vec![];
        let mut s: HashSet<usize> = inputs.iter().map(| node | *node).collect();
        loop {
            let c: HashSet<usize> = conns.iter()
                .filter(| conn | s.contains(&conn.from) && !s.contains(&conn.to))
                .map(| conn | conn.to).collect();
            let mut t = HashSet::new();
            for n in c {
                if required.contains(&n) && conns.iter().filter(| conn | conn.to == n).all(| conn | s.contains(&conn.from)) {
                    t.insert(n);
                }
            }
            if t.is_empty() {
                break;
            }
            layers.push(Layer::new(t.clone()));
            s.extend(t);
        }
        layers
    }

    fn required_for_outputs(inputs: &Vec<usize>, outputs: &Vec<usize>, conns: Vec<&ConnectionGene>) -> HashSet<usize> {
        let inputs = inputs.iter().map(| node| *node ).collect::<HashSet<usize>>();
        let mut required = outputs.iter().map(| node| *node ).collect::<HashSet<usize>>();
        let mut s = required.clone();
        loop {
            let t: HashSet<usize> = conns.iter().filter(| conn | s.contains(&conn.to) && !s.contains(&conn.from)).map(| conn | conn.from).collect();
            if t.is_empty() {
                break;
            }
            let layer_nodes: HashSet<&usize> = t.iter().filter(| x | !inputs.contains(x)).collect();
            if layer_nodes.is_empty() {
                break;
            }
            required.extend(layer_nodes);
            s.extend(t);
        }
        required
    }

}

impl Network for FeedForwardNetwork {
    fn evalute(&mut self, inputs: Vec<f64>) -> Vec<f64> {
        assert!(inputs.len() == self.inputs.len());
        let zipped: Vec<(&f64, &usize)> = inputs.iter().zip(self.inputs.iter()).collect();
        // initilize inputs
        for (value, key) in zipped {
            self.values.insert(*key, *value);
        }

        for eval in &self.node_evals {
            let mut input_values = vec![];
            for (key, weight) in &eval.incomings {
                if let Some(value) = self.values.get(key) {
                    input_values.push(value * weight)
                } else {
                    panic!("can not find value at {:?} \n {:?}", key, self)
                }
            }
            let agg_value = (self.aggregation_func)(input_values);
            let value = (self.activation_func)(eval.response * agg_value + eval.bias);
            self.values.insert(eval.key, value);
        }

        // if output node does not have connections to it default to 0.0
        self.outputs.iter().map(| node | *self.values.get(node).unwrap_or(&0.0)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ConnectionGene;
    use crate::Config;

    #[test]
    fn test_required_for_outputs(){
        let config = Config::default();
        let inputs = vec![1, 2, 3];
        // assume we have 5 ,6,7 as hidden layer
        let outputs = vec![8, 9, 10];
        let conns = vec![
            ConnectionGene::new(&config, 1, 5, 1, 0.1, true), 
            ConnectionGene::new(&config, 2, 6, 1, 0.1, true),
            ConnectionGene::new(&config, 3, 7, 1, 0.1, true),
            ConnectionGene::new(&config, 5, 8, 1, 0.1, true),
            ConnectionGene::new(&config, 6, 9, 1, 0.1, true),
        ];
        let conns = conns.iter().map(|conn| conn).collect();
        let required = FeedForwardNetwork::required_for_outputs(&inputs, &outputs, conns);
        assert_eq!(required.len(), 5);
    }

    #[test]
    fn test_create_layers(){
        let config = Config::default();
        let inputs = vec![1, 2, 3];
        // assume we have 5 ,6,7 as hidden layer
        let outputs = vec![8, 9, 10];
        let conns = vec![
            ConnectionGene::new(&config, 1, 5, 1, 0.1, true), 
            ConnectionGene::new(&config, 2, 6, 1, 0.1, true),
            ConnectionGene::new(&config, 3, 7, 1, 0.1, true),
            ConnectionGene::new(&config, 5, 8, 1, 0.1, true),
            ConnectionGene::new(&config, 6, 9, 1, 0.1, true),
        ];
        let conns = conns.iter().map(|conn| conn).collect();
        let layers = FeedForwardNetwork::create_layers(&inputs, &outputs, conns);
        assert_eq!(layers.len(), 2);
    }

    #[test]
    fn test_new_from_genome(){
        let config = Config::default();
        let mut genome = Genome::new(&config);
        genome.minimal_network();
        let forward_nn = FeedForwardNetwork::new(&config, &genome);
        assert_eq!(forward_nn.layers.len(), 1);
        
        let mut genome = Genome::new(&config);
        genome.genes.insert(3, Gene::new(&config, 3, GeneType::Hidden));
        genome.genes.insert(4, Gene::new(&config, 4, GeneType::Hidden));
        genome.add_connection(ConnectionGene::new(&config, 0, 3, 1, 0.5, true));
        genome.add_connection(ConnectionGene::new(&config, 1, 4, 2, 0.5, true));
        genome.add_connection(ConnectionGene::new(&config, 3, 2, 3, 0.5, true));
        genome.add_connection(ConnectionGene::new(&config, 4, 2, 4, 0.5, true));
        let forward_nn = FeedForwardNetwork::new(&config, &genome);
        assert_eq!(forward_nn.layers.len(), 2);
    }

    #[test]
    fn test_evaluate(){
        let config = Config::default();
        let genome = Genome::new(&config);
        let mut forward_nn = FeedForwardNetwork::new(&config, &genome);
        let outpus = forward_nn.evalute(vec![0.1, 0.2]);
        assert_eq!(outpus.len(), 1);
        assert_eq!(outpus[0], 0.0);
    }
}
