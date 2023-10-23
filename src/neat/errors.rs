use std::error::Error;
use std::fmt;
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub enum Errors {
    InputSizeNotMatch(String),
    GeneNotExists(usize),
    ConnectionGeneNotExists(usize),
    PopulationExtinction(),
    CanNotFindSolution()
}

impl Error for Errors {}

// Implement the Display trait for your custom error type
impl fmt::Display for Errors {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Errors::InputSizeNotMatch(details) => write!(f, "Input error: {}", details),
            Errors::GeneNotExists(index) => write!(f, "Gene #{} Does Not Exist", index),
            Errors::ConnectionGeneNotExists(innov) => write!(f, "Connection Gene #{} Does Not Exist", innov),
            Errors::PopulationExtinction() => write!(f, "Population has extincted"),
            Errors::CanNotFindSolution() => write!(f, "Can not find solution"),
        }
    }
}
