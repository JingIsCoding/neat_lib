use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use crate::neat::attribute::FloatAttribute;
use crate::neat::config::Config;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum GeneType {
    Input,
    Output,
    Hidden,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Gene {
    config: Config,
    pub key: usize,
    pub value: FloatAttribute,
    pub bias: FloatAttribute,
    pub response: FloatAttribute,
    pub gene_type: GeneType,
    pub incoming_conns: HashSet<usize>
}

impl Gene {
    pub fn new(config: &Config, key: usize, gene_type: GeneType) -> Self {
        Self::new_with_value(config, key,0.0, gene_type)
    }

    pub fn new_with_value(config: &Config, key: usize, value: f64, gene_type: GeneType) -> Self {
        Gene { 
            config: config.clone(),
            key,
            value: FloatAttribute::new(value), 
            bias: FloatAttribute::new_norm(0.0, 1.0), 
            response: FloatAttribute::new_norm(0.0, 1.0), 
            incoming_conns: HashSet::new(), 
            gene_type 
        }
    }

    pub fn mutate(&mut self) {
        self.bias.mutate(self.config.bias_replace_rate, self.config.bias_mutate_rate, self.config.bias_mutate_power, self.config.bias_min_value, self.config.bias_max_value, self.config.bias_init_mean, self.config.bias_init_stdev);
        self.response.mutate(self.config.response_replace_rate, self.config.response_mutate_rate, self.config.response_mutate_power, self.config.response_min_value, self.config.response_max_value, self.config.response_init_mean, self.config.response_init_stdev);
    }

    pub fn add_incomming_conn(&mut self, incoming_conn: usize) -> bool {
        self.incoming_conns.insert(incoming_conn)
    }

    pub fn remove_incoming_conn(&mut self, incoming_conn: &usize) -> bool {
        self.incoming_conns.remove(incoming_conn)
    }
}

impl std::fmt::Debug for Gene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Gene[ key:{:?} type:{:?} \tvalue:{:?} \tbias:{:?} \tresponse:{:?} \tincoming_conns:{:?} ]", self.key, self.gene_type, self.value, self.bias, self.response, self.incoming_conns)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_with_value() {
        let config = Config::default();
        let gene = Gene::new_with_value(&config, 1, 0.0, GeneType::Input);
    }
}
