pub(crate) mod feedforward;
pub(crate) mod layer;

pub trait Network {
    fn evalute(&mut self, inputs: Vec<f64>) -> Vec<f64>;
}

pub use feedforward::*;

