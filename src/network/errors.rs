use std::error::Error;
use std::fmt;
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub enum Errors {
    InputSizeNotMatch(String)
}

impl Error for Errors {}

// Implement the Display trait for your custom error type
impl fmt::Display for Errors {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Errors::InputSizeNotMatch(details) => write!(f, "Input error: {}", details),
        }
    }
}
