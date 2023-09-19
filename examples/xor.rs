extern crate neat_lib;

use neat_lib::*;
const FILE_PATH: &str = "./tmp/xor";
const MAX_ROUND: u32 = 1500;

fn main() {
    let mut config = Config::new_from_path("./examples/xor_config.json");
    config.inputs = 2;
    config.outputs = 1;
    config.population = 100;
    let mut pop = Population::new(config, "xor");
    while pop.generation < MAX_ROUND {
        println!("{:?}\n", pop);
        if pop.top_fitness > 15.0 {
            println!("succeed at {:?}", pop.generation);
            pop.save_to(FILE_PATH);
            return
        }
        for genome in pop.all_genomes() {
            let mut fitness = 0.0;
            for i in 0..2 {
                for j in 0..2 {
                    let inputs = vec![i as f64, j as f64];
                    let outputs = genome.evalute(&inputs);
                    if let Ok(outputs) = outputs {
                        let expected = (i ^ j) as f64;
                        fitness += 1.0 - (expected - outputs[0]).abs();
                    }
                }
            }
            fitness = fitness * fitness;
            genome.fitness = fitness;
        }
        pop.breed_new_generation();
    }

    println!("Failed to solve the problem.");
}
