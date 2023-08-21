use crate::network::*;
use serde_json;
use std::fs::{File, OpenOptions, read_to_string};
use std::io::prelude::*;
use std::io::{BufRead, BufReader};

pub trait Saver<T> {
    fn save(&self, data: &T) -> std::io::Result<()>;
}

pub trait Loader<T> {
    fn load(&self) -> std::io::Result<T>;
}

pub struct FileSaverLoader {
    pub path: String
}

impl FileSaverLoader {
    pub fn new(path: &str) -> Self {
        FileSaverLoader { path: path.to_owned() }
    }
}

impl Saver<population::Population> for FileSaverLoader {
    fn save(&self, population: &population::Population) -> std::io::Result<()> {
        let mut file = OpenOptions::new().create(true).append(true).open(&self.path)?;
        let json = serde_json::to_string(population)?;
        file.write_all(format!("{}\n", json).as_bytes())
            
    }
}

impl Loader<population::Population> for FileSaverLoader {
    fn load(&self) -> std::io::Result<population::Population> {
        let file = File::open(&self.path)?;
        let reader = BufReader::new(file);
        let last_line = reader.lines().fold(None, |acc, line| Some(line.unwrap_or(acc.unwrap_or_default())));
        match last_line {
            Some(line) => Ok(serde_json::from_str(&line)?),
            None => Err(std::io::Error::new(std::io::ErrorKind::Other, "empty file"))
        }
    }
}

impl Loader<config::Config> for FileSaverLoader {
    fn load(&self) -> std::io::Result<config::Config> {
        let data = read_to_string(&self.path)?;
        let config: config::Config = serde_json::from_str(&data)?;
        Ok(config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_saver() {
        let config = config::Config::new(1, 2, 10);
        let population = population::Population::new(config, "test");
        let saver = FileSaverLoader::new("./tmp/test");
        saver.save(&population).expect("should work");
        saver.save(&population).expect("should work");
        let _: population::Population = saver.load().expect("should load");
    }

}
