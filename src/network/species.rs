use super::{genome::Genome, config::Config};
use serde::{Deserialize, Serialize};
use rand::Rng;

#[derive(Debug, Serialize, Deserialize)]
pub struct Species {
    config: Config,
    pub genomes: Vec<Genome>,
    pub top_fitness: f64,
    pub staleness: u32,
}

impl Species {
    pub fn new(config: Config) -> Self {
        Species { 
            config,
            genomes: vec![],
            top_fitness: 0.0,
            staleness: 0,
        }
    }

    pub fn new_from(species: &Species) -> Self {
        Self { config: species.config, genomes: vec![], top_fitness: species.top_fitness, staleness: species.staleness }
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

    pub fn get_top_genome(&mut self) -> Option<&mut Genome> {
        self.sort_genomes();
        self.get_genome_by_index(0)
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
        let survive_count = (self.genomes.len() / 2) + 1;
        let mut survived = vec![];
        for i in 0..survive_count {
            survived.push(Genome::copy_from(&self.genomes[i], true))
        }
        self.genomes = survived;
    }

    pub fn check_progress(&mut self) {
        let top_fitness = self.top_fitness;
        if let Some(genomo) = self.get_top_genome() {
            if genomo.fitness > top_fitness {
                self.top_fitness = genomo.fitness;
                self.staleness = 0;
            } else {
                self.staleness += 1;
            }
        }
    }

    // assume the genomes are sorted
    pub fn breed_new_child(&mut self) -> Genome {
        let mut rng = rand::thread_rng();
        let length = self.genomes.len();
        self.sort_genomes();
        let mut child = if rng.gen_range(0.0..1.0) < self.config.crossover_chance && length > 1 {
            let parent1 = &self.genomes[0];
            let parent2 = &self.genomes[1];
            Genome::cross_over(parent1, parent2)
        } else {
            Genome::copy_from(&self.genomes[0], true)
        };
        child.mutate();
        child
    }

    pub fn calculate_adjusted_fitness(&mut self) {
        let length = self.genomes.len() as f64;
        self.genomes.iter_mut().for_each(|genome| {
            genome.adjusted_fitness = if genome.fitness > 0.0 {
                genome.fitness / length
            } else {
                0.0
            }
        });
    }

    pub fn total_adjusted_fitness(&self) -> f64 {
        self.genomes.iter().fold(0.0, | acc, genome | acc + genome.adjusted_fitness)
    }

    pub fn debug(&self) {
        for genome in &self.genomes {
            genome.debug();
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use super::super::config::Config;
    #[test]
    fn test_sort_genomes(){
        let config = Config::new(1, 1, 10);
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
    fn test_remove_weak_genomes(){
        let config = Config::new(1, 1, 10);
        let mut species = Species::new(config);
        let mut genome1 = Genome::new(config);
        genome1.fitness = 3.0;
        let mut genome2 = Genome::new(config);
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(config);
        genome3.fitness = 2.0;
        species.genomes = vec![genome1, genome2, genome3];
        species.sort_genomes();

        species.remove_weak_genomes();

        assert_eq!(species.genomes.len(), 1);
        assert_eq!(species.genomes[0].fitness, 3.0);
    }

    #[test]
    fn test_check_progress_making_progress() {
        let config = Config::new(1, 1, 10);
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
        let config = Config::new(1, 1, 10);
        // higher fitness
        let mut species = Species::new(config);
        species.top_fitness = 4.0;
        let mut genome1 = Genome::new(Config::new(1, 1, 10));
        genome1.fitness = 3.0;
        let mut genome2 = Genome::new(Config::new(1, 1, 10));
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(Config::new(1, 1, 10));
        genome3.fitness = 2.0;
        species.genomes = vec![genome1, genome2, genome3];
        species.check_progress();
        assert_eq!(species.staleness, 1);
        assert_eq!(species.top_fitness, 4.0);
    }

    #[test]
    fn test_check_calculate_adjusted_fitness() {
        let config = Config::new(1, 1, 10);
        // higher fitness
        let mut species = Species::new(config);
        let mut genome1 = Genome::new(Config::new(1, 1, 10));
        genome1.fitness = 3.0;
        let mut genome2 = Genome::new(Config::new(1, 1, 10));
        genome2.fitness = 1.0;
        let mut genome3 = Genome::new(Config::new(1, 1, 10));
        genome3.fitness = 2.0;
        species.genomes = vec![genome1, genome2, genome3];
        species.calculate_adjusted_fitness();

        assert_eq!(species.genomes[0].adjusted_fitness, 1.0);
        assert_eq!(species.genomes[1].adjusted_fitness, 1.0 / 3.0);
        assert_eq!(species.genomes[2].adjusted_fitness, 2.0 / 3.0);
        assert_eq!(species.total_adjusted_fitness(), 2.0);
    }
}
