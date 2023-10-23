extern crate neat_lib;

use neat_lib::*;
const FILE_PATH: &str = "./tmp/xor";
const MAX_ROUND: u32 = 1500;

fn evalute(config: &Config, genome: &Genome) -> f64 {
    let mut nn = FeedForwardNetwork::new(config, genome);
    let mut fitness = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            let inputs = vec![i as f64, j as f64];
            let outputs = nn.evalute(inputs);
            let expected = (i ^ j) as f64;
            fitness += 1.0 - (expected - outputs[0]).abs();
        }
    }
    fitness
}

fn main() {
    let mut config = Config::new_from_path("./examples/xor_config.json");
    config.inputs = 2;
    config.outputs = 1;
    config.population = 100;
    let mut pop = Population::new(&config, "xor");
    match pop.run(evalute, EvaluteOption { epochs: 500, fitness_target: Some(3.6) }){
        Ok(genome) => {
            print!("top genome {:?}", genome);
        },
        Err(err) => {
            println!("Failed to solve the problem.");
            print!("err {:?}", err);
        },
    }
}
