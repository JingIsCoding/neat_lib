use std::vec;
use serde::{Deserialize, Serialize};

use crate::loader::save_load::{FileSaverLoader, Saver};

use super::species::Species;
use super::config::Config;
use super::genome::Genome;
use super::innovation_number::reset;
use crate::network::errors::*;

#[derive(Debug, Serialize, Deserialize)]
pub struct Population {
    pub name: String,
    pub generation: u32,
    pub top_fitness: f64,
    pub species: Vec<Species>,
    pub staleness: u32,
    config: Config,
}

impl Population {
    pub fn new(config: Config, name: &str) -> Self {
        reset();
        let mut pool = Population { 
            name: name.to_owned(),
            config,
            generation: 1,
            species: vec![],
            top_fitness: 0.0,
            staleness: 0,
        };
        pool.init_species();
        pool
    }

    fn init_species(&mut self) {
        let mut seed = Genome::new(self.config);
        seed.generate_network();
        for _ in 0..self.config.population {
            let copy_gene = Genome::copy_from(&seed, false);
            self.add_to_species(copy_gene);
        }
    }

    fn add_to_species(&mut self, genome: Genome) {
        for i in 0..self.species.len() {
            let species = &mut self.species[i];
            if let Some(existing_genome) = species.get_genome_by_index(0) {
                if Genome::is_same_species(&genome, existing_genome) {
                    species.add_genome(genome);
                    return
                };
            }
        }
        let mut new_species = Species::new(self.config);
        new_species.add_genome(genome);
        self.species.push(new_species);
    }

    pub fn evaluate(&mut self, inputs: Vec<Vec<f64>>) -> Result<Vec<Vec<f64>>, Errors> {
        let mut outputs = vec![vec![0.0; self.config.outputs as usize]; inputs.len()];
        let mut genomes = self.all_genomes();
        if inputs.len() != genomes.len() {
            return Err(Errors::InputSizeNotMatch("Evalute inputs size not match genomes size".to_owned()))
        }
        for (index, inputs) in inputs.iter().enumerate() {
            match genomes[index].evalute_network(inputs){
                Ok(outputs_at_index) => {
                    outputs[index] = outputs_at_index;
                },
                Err(err) => {
                    return Err(err)
                }
            }
        }
        Ok(outputs)
    }

    pub fn set_fitness(&mut self, fitnesses: Vec<f64>) -> Result<(), Errors> {
        let mut genomes = self.all_genomes();
        if fitnesses.len() != genomes.len() {
            return Err(Errors::InputSizeNotMatch("Set fitnesses size not match genomes size".to_owned()))
        }
        for (index, genome) in genomes.iter_mut().enumerate() {
            genome.fitness = fitnesses[index];
        }
        Ok(())
    }

    pub fn all_genomes(&mut self) -> Vec<&mut Genome> {
        self.species.iter_mut().fold(vec![], |mut acc, species| {
            for genomo in &mut species.genomes {
                acc.push(genomo);
            }
            acc
        })
    }

    fn calculate_genome_adjusted_fitness(&mut self) {
        self.species.iter_mut()
            .for_each(| species | species.calculate_adjusted_fitness());
    }
    
    fn remove_weak_genome_from_species(&mut self) {
        self.species.iter_mut().for_each(| species | species.remove_weak_genomes());
    }

    fn remove_stale_species(&mut self) {
        self.species.iter_mut().for_each(| species | species.check_progress());
        self.check_progress();

        self.species.retain(| species | species.staleness < self.config.stale_species_threshold || species.top_fitness >= self.top_fitness );
        if self.staleness >= self.config.stale_population_threshold {
            self.species.retain(| species | species.top_fitness >= self.top_fitness )
        }
    }

    fn check_progress(&mut self) {
        let (top_fitness, _) = self.get_top_fitness();
        if top_fitness > self.top_fitness {
            self.top_fitness = top_fitness;
            self.staleness = 0;
        } else {
            self.staleness += 1;
        }
    }

    pub fn breed_new_generation(&mut self) {
        let population = self.config.population;
        self.calculate_genome_adjusted_fitness();
        self.remove_stale_species();
        self.remove_weak_genome_from_species();
        let global_adjusted_fitness = self.species.iter().fold(0.0, | acc, species | acc + species.total_adjusted_fitness());
        let mut children = vec![];
        let mut new_species = vec![];
        let species_size = self.species.len();
        for species in self.species.iter_mut() {
            // keep the top genome in the species
            let mut cloned_species = Species::new_from(species);
            cloned_species.add_genome(Genome::copy_from(species.get_top_genome().unwrap(), true));
            new_species.push(cloned_species);
            let ratio: f64 = if global_adjusted_fitness == 0.0 {
                1.0 / species_size as f64
            } else {
                species.total_adjusted_fitness() / global_adjusted_fitness
            };
            println!("raito is {:?}", ratio);
            let child_num = (population as f64 * ratio) as i32;
            for _ in 0..(child_num - 1) {
                let child = species.breed_new_child();
                children.push(child);
            }
        }
        let diff = population  - children.len() as u32;
        for i in 0..diff {
            let size = self.species.len();
            children.push(self.species[i as usize % size].breed_new_child());
        }
        self.species = new_species;
        for child in children {
            self.add_to_species(child);
        }
        assert!(self.all_genomes().len() == population as usize);
        self.generation += 1;
    }

    fn calculate_pop_size_for_each_species(pop_size: usize, adjusted_fitnesses: Vec<f64>, species_pop_size: Vec<usize>, min_species_size: usize) -> Vec<usize> {
        let sum = adjusted_fitnesses.iter().fold(0.0, |acc, fitness| {
             acc + fitness
        });
        let mut new_species_pos_size = vec![];
        for (adjusted_fitness, pre_size) in adjusted_fitnesses.iter().zip(species_pop_size.iter()) {
            let mut species_size: i64 = 0;
            if sum > 0.0 {
                species_size = min_species_size.max((adjusted_fitness / sum * pop_size as f64) as usize) as i64;
            }         
            let d = (species_size - *pre_size as i64) as f64 * 0.5;
            let c = d.round() as i64;
            species_size = *pre_size as i64;
            if c.abs() > 0 {
                species_size += c;
            } else if d > 0.0 {
                species_size += 1;
            } else if d < 0.0 {
                species_size -= 1;
            }
            new_species_pos_size.push(species_size as usize);
        }
        let diff = pop_size as i64 - new_species_pos_size.iter().fold(0, | acc, size | acc + size ) as i64;
        for i in 0..diff.abs() {
            let index = i as usize % new_species_pos_size.len();
            if diff > 0 {
                new_species_pos_size[index] += 1;
            } else {
                new_species_pos_size[index] -= 1;
            }
        }
        new_species_pos_size
    }

    pub fn get_top_fitness(&mut self) -> (f64, Option<&mut Genome>) {
        let mut top_fitness = 0.0;
        let mut top_genomo = None;
        for species in &mut self.species {
            if let Some(genome) = species.get_top_genome() {
                if genome.fitness >= top_fitness {
                    top_fitness = genome.fitness;
                    top_genomo = Some(genome);
                }
            }
        }
        return (top_fitness, top_genomo);
    }

    pub fn save_to(&self, file: &str) -> std::io::Result<()> {
        let saver = FileSaverLoader::new(file);
        saver.save(self)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use super::super::config::Config;
    #[test]
    fn test_init_species() {
        let config = Config::default();
        let population = Population::new(config, "test");
        assert_eq!(population.species.len(), 1);
    }

    #[test]
    fn test_evaluate() {
        let config = Config::default();
        let mut population = Population::new(config, "test");
        assert_eq!(population.species.len(), 1);
        let outputs = population.evaluate(vec![vec![1.0], vec![0.4], vec![0.3], vec![0.5], vec![0.2]]);
        assert_eq!(outputs.unwrap().len(), 5);
    }

    #[test]
    fn test_set_fitnesses() {
        let config = Config::default();
        let mut population = Population::new(config, "test");
        let fitnesses = vec![0.1, 0.2, 0.3, 0.4, 0.5];
        population.set_fitness(vec![0.1, 0.2, 0.3, 0.4, 0.5]);
        let genomes = population.all_genomes();
        for (i, genome) in genomes.iter().enumerate() {
            assert_eq!(genome.fitness, fitnesses[i]);
        }
    }

    #[test]
    fn test_check_progress_making_progress() {
        let config = Config::default();
        // higher fitness
        let mut species1 = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 4.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species1.genomes = vec![genome1, genome2, genome3];

        let mut species2 = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 5.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species2.genomes = vec![genome1, genome2, genome3];

        let mut population = Population{
            name: "test".to_owned(),
            config,
            generation: 0,
            top_fitness: 0.0,
            staleness: 0,
            species: vec![species1, species2]
        };

        population.check_progress();
        assert_eq!(population.staleness, 0);
        assert_eq!(population.top_fitness, 5.0);
    }


    #[test]
    fn test_check_progress_not_making_progress() {
        let config = Config::default();
        // higher fitness
        let mut species1 = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 4.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species1.genomes = vec![genome1, genome2, genome3];

        let mut species2 = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 5.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species2.genomes = vec![genome1, genome2, genome3];

        let mut population = Population{
            config,
            name: "test".to_owned(),
            generation: 0,
            top_fitness: 6.0,
            staleness: 0,
            species: vec![species1, species2]
        };

        population.check_progress();
        assert_eq!(population.staleness, 1);
        assert_eq!(population.top_fitness, 6.0);
    }

    #[test]
    fn test_remove_stale_species() {
        let mut config = Config::default();
        config.stale_species_threshold = 1;
        config.stale_population_threshold = 1;

        // should become stale
        let mut species1 = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 4.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species1.genomes = vec![genome1, genome2, genome3];
        species1.top_fitness = 5.0;

        // should not become stale
        let mut species2 = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 5.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species2.genomes = vec![genome1, genome2, genome3];
        species2.top_fitness = 4.0;

        // should not become stale
        let mut species3 = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 6.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species3.genomes = vec![genome1, genome2, genome3];
        species3.top_fitness = 4.0;

        let mut population = Population{
            config,
            name: "test".to_owned(),
            generation: 0,
            top_fitness: 5.9,
            staleness: 0,
            species: vec![species1, species2, species3]
        };

        population.remove_stale_species();
        assert_eq!(population.species.len(), 2);
        assert_eq!(population.top_fitness, 6.0);
        assert_eq!(population.species[0].top_fitness, 5.0);
        assert_eq!(population.species[1].top_fitness, 6.0);
    }

    #[test]
    fn test_remove_stale_species_when_population_is_stale() {
        let mut config = Config::default();
        config.stale_species_threshold = 1;
        config.stale_population_threshold = 1;

        // should become stale
        let mut species1 = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 4.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species1.genomes = vec![genome1, genome2, genome3];
        species1.top_fitness = 5.0;

        // should not become stale
        let mut species2 = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 5.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species2.genomes = vec![genome1, genome2, genome3];
        species2.top_fitness = 4.0;

        // should not become stale
        let mut species3 = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 6.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species3.genomes = vec![genome1, genome2, genome3];
        species3.top_fitness = 4.0;

        let mut population = Population{
            config,
            name: "test".to_owned(),
            generation: 0,
            top_fitness: 6.0,
            staleness: 1,
            species: vec![species1, species2, species3]
        };

        population.remove_stale_species();
        assert_eq!(population.species.len(), 1);
        assert_eq!(population.species[0].top_fitness, 6.0);
    }

    #[test]
    fn test_breed_new_generation() {
        let config = Config::default();
        let mut population = Population::new(config, "test");
        population.breed_new_generation();
        assert_eq!(population.generation, 2);
    }
}

