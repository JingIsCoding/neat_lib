use super::{population::*, config::*, genome::*};

#[cfg(test)]
mod tests {
    use crate::{loader::save_load::{FileSaverLoader, Loader}};
    use std::io::{stdout, stdin, Write};
    use super::{Population, Config, Genome};
    const FILE_PATH: &str = "./tmp/xor";
    #[test]
    fn test_xor() {
        let config = Config::new(2, 1, 100);
        let mut population = Population::new(config, "xor");
        while population.generation < 5000 {
            println!("top_fitness at {:?} {:?}", population.generation, population.top_fitness);
            run_double_check(population.get_top_fitness());
            if population.top_fitness > 15.0 {
                run_double_check(population.get_top_fitness());
                println!("succeed at {:?}", population.generation);
                population.save_to(FILE_PATH);
                return
            }
            for genome in population.all_genomes() {
                let mut fitness = 0.0;
                for i in 0..2 {
                    for j in 0..2 {
                        let inputs = vec![i as f64, j as f64];
                        let outputs = genome.evalute_network(&inputs);
                        let expected = (i ^ j) as f64;
                        fitness += 1.0 - (expected - outputs[0]).abs();
                    }
                }
                fitness = fitness * fitness;
                genome.fitness = fitness;
            }
            population.breed_new_generation();
        }
        panic!("failed..");
    }

    fn run_double_check(top_genome: (f64, Option<&mut Genome>)) {
        let (_ ,genome) = top_genome;
        if let Some(genome) = genome {
            println!("top genome is {:?}", genome.ref_id);
            //genome.debug();
            for i in 0..2 {
                for j in 0..2 {
                    let inputs = vec![i as f64, j as f64];
                    let outputs = genome.evalute_network(&inputs);
                    println!("{:?} ^ {:?} is {:?}", i, j, outputs[0])
                }
            }
        }
    }

    #[test]
    fn test_interactive_xor() {
        let mut population: Population = FileSaverLoader::new(FILE_PATH).load().expect("should exist");
        print!("fitness is {:?}", population.top_fitness);
        let (_, top_genomo) = population.get_top_fitness();
        if let Some(genome) = top_genomo {
            let stdin = stdin();
            let mut stdout = stdout();
            loop {
                // Print a prompt to the user
                print!("Enter inputs: q to quit");
                stdout.flush().expect("Failed to flush stdout.");
                // Read a line of input from the user
                let mut input = String::new();
                stdin
                    .read_line(&mut input)
                    .expect("Failed to read line from input.");
                let parts: Vec<&str> = input.split(",").map(| input | input.trim()).collect();
                if parts.len() == 1 {
                    break;
                }
                let inputs: Vec<f64> = parts.iter().map(| input | input.parse().unwrap()).collect();
                let outputs = genome.evalute_network(&inputs);
                println!("outputs are {:?} {:?} ", inputs[0] as usize ^ inputs[1] as usize, outputs);
            }
        }
    }
}
