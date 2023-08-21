use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum GeneType {
    Input,
    Output,
    Hidden,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gene {
    pub value: f64,
    pub gene_type: GeneType,
    pub incoming_conns: Vec<usize>,
}

impl Gene {
    pub fn new(gene_type: GeneType) -> Self {
        Gene { value: 0.0, incoming_conns: vec![], gene_type }
    }

    pub fn new_with_value(value: f64, gene_type: GeneType) -> Self {
        Gene { value, incoming_conns: vec![], gene_type }
    }

    pub fn add_incomming_conn(&mut self, incoming_conn: usize) {
        self.incoming_conns.push(incoming_conn)
    }

    pub fn set_incomming_conns(&mut self, incoming_conns: Vec<usize>) {
        self.incoming_conns = incoming_conns;
    }
}
