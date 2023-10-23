use std::collections::HashSet;

#[derive(Debug)]
pub struct Layer {
    pub keys: HashSet<usize>
}

impl Layer {
    pub fn new(keys: HashSet<usize>) -> Self {
        Layer { keys }
    }
}
