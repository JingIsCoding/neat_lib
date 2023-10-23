use std::vec;
use serde::{Deserialize, Serialize};

use crate::loader::save_load::{FileSaverLoader, Saver};

use super::species::Species;
use super::config::Config;
use super::genome::Genome;
use super::innovation_number::reset;
use crate::neat::errors::*;

#[derive(Serialize, Deserialize)]
pub struct Population {
    pub name: String,
    pub generation: u32,
    pub top_fitness: f64,
    pub top_genome: Option<Genome>,
    pub species: Vec<Species>,
    pub staleness: u32,
    config: Config,
}

pub struct EvaluteOption  {
    pub epochs: u64,
    pub fitness_target: Option<f64>
}

impl EvaluteOption {
    pub fn new(epochs: u64) -> Self {
        EvaluteOption { epochs, fitness_target: None }
    }
}

pub type EvaluteFunc = fn(config: &Config, genome: &Genome) -> f64;

impl Population {
    pub fn new(config: &Config, name: &str) -> Self {
        reset();
        let mut pool = Population { 
            name: name.to_owned(),
            config: config.clone(),
            generation: 1,
            species: vec![],
            top_fitness: 0.0,
            top_genome: None,
            staleness: 0,
        };
        pool.init_species();
        pool
    }

    fn init_species(&mut self) {
        let mut seed = Genome::new(&self.config);
        seed.minimal_network();
        for _ in 0..self.config.population {
            let copy_gene = Genome::new_from(&seed, false);
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
        let mut new_species = Species::new(&self.config);
        new_species.add_genome(genome);
        self.species.push(new_species);
    }

    pub fn run(&mut self, evalute_func: EvaluteFunc, option: EvaluteOption) -> Result<Genome, Errors> {
        let mut k = 0;
        let config = self.config.clone();
        while k < option.epochs {
            k += 1;
            let genomes = self.all_genomes();
            for genome in genomes {
                let fitness = evalute_func(&config, genome);
                genome.fitness = fitness;
                if let Some(fitness_target) = option.fitness_target {
                    if fitness >= fitness_target {
                        return Ok(genome.clone());
                    }
                }
            }
            self.breed_new_generation()?;
            println!("round {:?} {:?}", k, self.get_top_genome())
        }
        let (_, genome) = self.get_top_genome();
        genome.ok_or(Errors::CanNotFindSolution())
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

    fn calculate_species_adjusted_fitness(&mut self) {
        self.species.iter_mut()
            .for_each(| species | {
                species.calculate_adjusted_fitness();
            });
    }
    
    fn remove_stale_species(&mut self) {
        self.species.iter_mut().for_each(| species | species.check_progress());
        self.check_progress();
        self.species.retain(| species | species.staleness < self.config.stale_species_threshold || species.top_fitness >= self.top_fitness );
        if self.staleness >= self.config.stale_population_threshold {
            self.species = vec![];
            self.staleness = 0;
            //self.species.retain(| species | species.top_fitness >= self.top_fitness )
        }
    }

    fn check_progress(&mut self) {
        let (top_fitness, genome) = self.get_top_genome();
        if top_fitness > self.top_fitness {
            self.top_fitness = top_fitness;
            self.staleness = 0;
            self.top_genome = genome;
        } else {
            self.staleness += 1;
        }
    }

    pub fn breed_new_generation(&mut self) -> Result<(), Errors> {
        let mut children = vec![];
        self.remove_stale_species();
        if self.species.is_empty() {
            if let Some(genome) = self.top_genome.clone() {
                for _ in 0..self.config.population {
                    children.push(genome.clone());
                }
            } else {
                return Err(Errors::PopulationExtinction());
            }
        } else {
            self.calculate_species_adjusted_fitness();
            let species_sizes = Self::calculate_pop_size_for_each_species(&self.species, self.config.population as usize, self.config.min_species_size);
            for (target_size, species) in species_sizes.iter().zip(&mut self.species) {
                children.extend(species.reproduce_to_size(*target_size));
            }
        }
        for child in children {
            self.add_to_species(child);
        }
        self.generation += 1;
        Ok(())
    }

    fn calculate_pop_size_for_each_species(species: &Vec<Species>, pop_size: usize, min_species_size: usize) -> Vec<usize> {
        let adjusted_fitness_sum = species.iter().fold(0.0, |acc, species| {
             acc + species.adjusted_fitness
        });
        let mut new_species_pos_size = vec![];
        for s in species {
            let adjusted_fitness = s.adjusted_fitness;
            let mut species_size: i64 = min_species_size as i64;
            if adjusted_fitness_sum > 0.0 {
                species_size = min_species_size.max((adjusted_fitness / adjusted_fitness_sum * pop_size as f64) as usize) as i64;
            };
            new_species_pos_size.push(species_size as usize);
        }
        let mut diff = pop_size as i64 - new_species_pos_size.iter().fold(0, | acc, size | acc + size ) as i64;
        let mut counter = 0;
        while diff != 0 {
            let index = counter as usize % new_species_pos_size.len();
            if diff > 0 {
                new_species_pos_size[index] += 1;
                diff -= 1;
            } else {
                if new_species_pos_size[index] > min_species_size {
                    new_species_pos_size[index] -= 1;
                    diff += 1;
                } else {
                    // every species is at min species size then give up.
                    if !new_species_pos_size.iter().any(| size | *size > min_species_size) {
                        break;
                    }
                }
            }
            counter += 1;
        }
        new_species_pos_size
    }

    pub fn get_top_genome(&self) -> (f64, Option<Genome>) {
        let mut top_fitness = 0.0;
        let mut top_genomo = None;
        for species in &self.species {
            if let Some(genome) = species.get_top_genome() {
                if genome.fitness >= top_fitness {
                    top_fitness = genome.fitness;
                    top_genomo = Some(genome.clone());
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

impl std::fmt::Debug for Population {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (_, genome) = self.get_top_genome();
        write!(f, "Generation {:?} \nTop fitness {:?}\nNum of species {:?}\nTop Genome {:?}\n", self.generation, self.top_fitness, self.species.len(), genome.unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::config::Config;
    #[test]
    fn test_init_species() {
        let config = Config::default();
        let population = Population::new(&config, "test");
        assert_eq!(population.species.len(), 1);
    }

    #[test]
    fn test_set_fitnesses() {
        let mut config = Config::default();
        config.population = 5;
        let mut population = Population::new(&config, "test");
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
        let mut species1 = Species::new(&config);
        let mut genome1 = Genome::new(&config);
        genome1.fitness = 4.0;
        let mut genome2 = Genome::new(&config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(&config);
        genome3.fitness = 2.0;
        species1.genomes = vec![genome1, genome2, genome3];

        let mut species2 = Species::new(&config);
        let mut genome1 = Genome::new(&config);
        genome1.fitness = 5.0;
        let mut genome2 = Genome::new(&config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(&config);
        genome3.fitness = 2.0;
        species2.genomes = vec![genome1, genome2, genome3];

        let mut population = Population{
            top_genome: None,
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
        let mut species1 = Species::new(&config);
        let mut genome1 = Genome::new(&config);
        genome1.fitness = 4.0;
        let mut genome2 = Genome::new(&config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(&config);
        genome3.fitness = 2.0;
        species1.genomes = vec![genome1, genome2, genome3];

        let mut species2 = Species::new(&config);
        let mut genome1 = Genome::new(&config);
        genome1.fitness = 5.0;
        let mut genome2 = Genome::new(&config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(&config);
        genome3.fitness = 2.0;
        species2.genomes = vec![genome1, genome2, genome3];

        let mut population = Population{
            top_genome: None,
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
        let mut species1 = Species::new(&config);
        let mut genome1 = Genome::new(&config);
        genome1.fitness = 4.0;
        let mut genome2 = Genome::new(&config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(&config);
        genome3.fitness = 2.0;
        species1.genomes = vec![genome1, genome2, genome3];
        species1.top_fitness = 5.0;

        // should not become stale
        let mut species2 = Species::new(&config);
        let mut genome1 = Genome::new(&config);
        genome1.fitness = 5.0;
        let mut genome2 = Genome::new(&config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(&config);
        genome3.fitness = 2.0;
        species2.genomes = vec![genome1, genome2, genome3];
        species2.top_fitness = 4.0;

        // should not become stale
        let mut species3 = Species::new(&config);
        let mut genome1 = Genome::new(&config);
        genome1.fitness = 6.0;
        let mut genome2 = Genome::new(&config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(&config);
        genome3.fitness = 2.0;
        species3.genomes = vec![genome1, genome2, genome3];
        species3.top_fitness = 4.0;

        let mut population = Population{
            config,
            name: "test".to_owned(),
            generation: 0,
            top_fitness: 5.9,
            staleness: 0,
            species: vec![species1, species2, species3],
            top_genome: None
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
        let mut species1 = Species::new(&config);
        let mut genome1 = Genome::new(&config);
        genome1.fitness = 4.0;
        let mut genome2 = Genome::new(&config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(&config);
        genome3.fitness = 2.0;
        species1.genomes = vec![genome1, genome2, genome3];
        species1.top_fitness = 5.0;

        // should not become stale
        let mut species2 = Species::new(&config);
        let mut genome1 = Genome::new(&config);
        genome1.fitness = 5.0;
        let mut genome2 = Genome::new(&config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(&config);
        genome3.fitness = 2.0;
        species2.genomes = vec![genome1, genome2, genome3];
        species2.top_fitness = 4.0;

        // should not become stale
        let mut species3 = Species::new(&config);
        let mut genome1 = Genome::new(&config);
        genome1.fitness = 6.0;
        let mut genome2 = Genome::new(&config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(&config);
        genome3.fitness = 2.0;
        species3.genomes = vec![genome1, genome2, genome3];
        species3.top_fitness = 4.0;

        let mut population = Population::new(&config, "test");
        population.species = vec![species1, species2, species3];
        population.remove_stale_species();
        assert_eq!(population.species.len(), 2);
        assert_eq!(population.species[0].top_fitness, 6.0);
    }

    #[test]
    fn test_breed_new_generation() {
        let config = Config::default();
        let mut population = Population::new(&config, "test");
        population.breed_new_generation();
        assert_eq!(population.generation, 2);
    }

    #[test]
    fn test_calculate_pop_size_for_each_species() {
        let mut config = Config::default();
        config.population = 100;

        let mut population = Population::new(&config, "test");

        // only on species
        population.species = vec![Species::new(&config)];
        population.species[0].adjusted_fitness = 0.0;
        population.species[0].genomes = vec![Genome::new(&config), Genome::new(&config)];
        let species_sizes = Population::calculate_pop_size_for_each_species(&population.species, 100, 2);
        assert_eq!(species_sizes[0], 100);

        // two species
        population.species = vec![Species::new(&config), Species::new(&config)];
        population.species[0].adjusted_fitness = 0.1;
        population.species[0].genomes = vec![Genome::new(&config), Genome::new(&config)];

        population.species[1].adjusted_fitness = 0.1;
        population.species[1].genomes = vec![Genome::new(&config), Genome::new(&config)];
        let species_sizes = Population::calculate_pop_size_for_each_species(&population.species, 100, 2);
        assert_eq!(species_sizes[0], 50);
        assert_eq!(species_sizes[1], 50);

        // two species
        population.species = vec![Species::new(&config), Species::new(&config)];
        population.species[0].adjusted_fitness = 0.0;
        population.species[0].genomes = vec![Genome::new(&config), Genome::new(&config)];

        population.species[1].adjusted_fitness = 0.1;
        population.species[1].genomes = vec![Genome::new(&config), Genome::new(&config)];
        let species_sizes = Population::calculate_pop_size_for_each_species(&population.species, 100, 2);
        assert_eq!(species_sizes[0], 2);
        assert_eq!(species_sizes[1], 98);

        // two species
        population.species = vec![Species::new(&config), Species::new(&config)];
        population.species[0].adjusted_fitness = 0.0;
        population.species[0].genomes = vec![Genome::new(&config), Genome::new(&config)];

        population.species[1].adjusted_fitness = 0.1;
        population.species[1].genomes = vec![Genome::new(&config), Genome::new(&config)];
        let species_sizes = Population::calculate_pop_size_for_each_species(&population.species, 4, 2);
        assert_eq!(species_sizes[0], 2);
        assert_eq!(species_sizes[1], 2);

        // two species
        population.species = vec![Species::new(&config), Species::new(&config)];
        population.species[0].adjusted_fitness = 0.0;
        population.species[0].genomes = vec![Genome::new(&config), Genome::new(&config)];

        population.species[1].adjusted_fitness = 0.1;
        population.species[1].genomes = vec![Genome::new(&config), Genome::new(&config)];
        let species_sizes = Population::calculate_pop_size_for_each_species(&population.species, 3, 2);
        assert_eq!(species_sizes[0], 2);
        assert_eq!(species_sizes[1], 2);

        population.species = vec![Species::new(&config), Species::new(&config)];
        population.species[0].adjusted_fitness = 0.0;
        population.species[0].genomes = vec![Genome::new(&config), Genome::new(&config)];

        population.species[1].adjusted_fitness = 0.1;
        population.species[1].genomes = vec![Genome::new(&config), Genome::new(&config)];
        let species_sizes = Population::calculate_pop_size_for_each_species(&population.species, 3, 1);
        assert_eq!(species_sizes[0], 1);
        assert_eq!(species_sizes[1], 2);

        population.species = vec![Species::new(&config), Species::new(&config)];
        population.species[0].adjusted_fitness = 0.0;
        population.species[0].genomes = vec![Genome::new(&config), Genome::new(&config)];

        population.species[1].adjusted_fitness = 0.1;
        population.species[1].genomes = vec![Genome::new(&config), Genome::new(&config)];
        let species_sizes = Population::calculate_pop_size_for_each_species(&population.species, 99, 2);
        assert_eq!(species_sizes[0], 2);
        assert_eq!(species_sizes[1], 97);
    }
}
