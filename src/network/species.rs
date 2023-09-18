use std::vec;

use super::{genome::{Genome, self}, config::Config};
use serde::{Deserialize, Serialize};
use rand::Rng;

#[derive(Debug, Serialize, Deserialize)]
pub struct Species {
    config: Config,
    pub genomes: Vec<Genome>,
    pub top_fitness: f64,
    pub adjusted_fitness: f64,
    pub staleness: u32,
}

impl Species {
    pub fn new(config: Config) -> Self {
        Species { 
            config,
            genomes: vec![],
            top_fitness: 0.0,
            adjusted_fitness: 0.0,
            staleness: 0,
        }
    }

    pub fn new_from(species: &Species) -> Self {
        Self { 
            config: species.config, 
            genomes: species.genomes.clone(), 
            top_fitness: species.top_fitness, 
            adjusted_fitness: species.adjusted_fitness, 
            staleness: species.staleness 
        }
    }

    pub fn is_empty(&self) -> bool {
        self.genomes.is_empty()
    }

    pub fn get_genome_by_index(&mut self, index: usize) -> Option<&mut Genome> {
        if index >= self.genomes.len() {
            return None;
        }
        Some(&mut self.genomes[index])
    }

    pub fn get_top_genome_mut(&mut self) -> Option<&mut Genome> {
        return self.genomes.iter_mut().reduce(| a, b | if a.fitness > b.fitness { a } else { b });
    }

    pub fn get_top_genome(&self) -> Option<&Genome> {
        return self.genomes.iter().reduce(| a, b | if a.fitness > b.fitness { a } else { b });
    }

    pub fn add_genome(&mut self, genome: Genome)  {
        self.genomes.push(genome)
    }

    pub(crate) fn sort_genomes(&mut self) {
        self.genomes.sort_by(|genome1, genome2| {
            if genome2.fitness > genome1.fitness {
                return std::cmp::Ordering::Greater
            } else if genome2.fitness == genome1.fitness{
                return std::cmp::Ordering::Equal
            }
            return std::cmp::Ordering::Less
        })
    }

    pub fn remove_weak_genomes(&mut self) {
        self.sort_genomes();
        let survive_count = (self.genomes.len() / 2) + 1;
        for i in 0..self.genomes.len() {
            self.genomes.retain(| _ | i < survive_count);
        }
    }

    pub fn check_progress(&mut self) {
        let top_fitness = self.top_fitness;
        if let Some(genomo) = self.get_top_genome_mut() {
            if genomo.fitness > top_fitness {
                self.top_fitness = genomo.fitness;
                self.staleness = 0;
            } else {
                self.staleness += 1;
            }
        }
    }

    // resize species to the target size and return new genomes
    pub fn reproduce_to_size(&mut self, target_size: usize) -> Vec<Genome> {
        let mut children = vec![];
        let cutoff = (self.genomes.len() as f64 * self.config.survival_threshold).ceil().max(self.config.min_species_size as f64) as usize;
        let cutoff = cutoff.min(target_size);
        if cutoff < self.genomes.len() {
            self.sort_genomes();
            self.genomes = self.genomes[..cutoff].to_vec();
        }
        let diff = target_size - self.genomes.len();
        for _ in 0..diff {
            children.push(self.breed_new_child());
        } 
        children
    }

    // assume the genomes are sorted
    pub fn breed_new_child(&mut self) -> Genome {
        let mut rng = rand::thread_rng();
        let length = self.genomes.len();
        self.sort_genomes();
        let mut child = if rng.gen_range(0.0..1.0) < self.config.crossover_chance {
            let parent1 = &self.genomes[rng.gen_range(0..length)];
            let parent2 = &self.genomes[rng.gen_range(0..length)];
            Genome::cross_over(parent1, parent2)
        } else {
            Genome::new_from(&self.genomes[0], true)
        };
        child.mutate();
        child
    }

    pub fn calculate_adjusted_fitness(&mut self) -> f64 {
        self.adjusted_fitness = (self.genomes.iter().map(|genome| genome.fitness).sum::<f64>() / self.genomes.len() as f64).max(0.0);
        self.adjusted_fitness
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use super::super::config::Config;
    #[test]
    fn test_sort_genomes(){
        let config = Config::default();
        let mut species = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 3.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species.genomes = vec![genome1, genome2, genome3];

        species.sort_genomes();

        assert_eq!(species.genomes[0].fitness, 3.0);
        assert_eq!(species.genomes[1].fitness, 2.0);
        assert_eq!(species.genomes[2].fitness, 1.0);
    }

    #[test]
    fn test_check_progress_making_progress() {
        let config = Config::default();
        // higher fitness
        let mut species = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 3.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species.genomes = vec![genome1, genome2, genome3];
        species.check_progress();
        assert_eq!(species.staleness, 0);
        assert_eq!(species.top_fitness, 3.0);
    }

    #[test]
    fn test_check_progress_not_making_progress() {
        let config = Config::default();
        // higher fitness
        let mut species = Species::new(config);
        species.top_fitness = 4.0;
        let mut genome1 = Genome::new(config);
        genome1.fitness = 3.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species.genomes = vec![genome1, genome2, genome3];
        species.check_progress();
        assert_eq!(species.staleness, 1);
        assert_eq!(species.top_fitness, 4.0);
    }

    #[test]
    fn test_check_calculate_adjusted_fitness() {
        let config = Config::default();
        // higher fitness
        let mut species = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 3.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species.genomes = vec![genome1, genome2, genome3];
        species.calculate_adjusted_fitness();

        assert_eq!(species.calculate_adjusted_fitness(), 2.0);
    }
}
