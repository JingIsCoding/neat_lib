mod xor;

mod tests {
    use crate::Config;
    use crate::Population;
    const FILE_PATH: &str = "./tmp/xor";

    #[test]
    fn test_xor(){
        let mut config = Config::default();
        config.inputs = 2;
        config.outputs = 1;
        config.population = 100;
        let mut pop = Population::new(config, "xor");
        while pop.generation < 1500 {
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
                        let outputs = genome.evalute_network(&inputs);
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
        if let Ok(outputs) = pop.evaluate(vec![vec![1.0, 1.0]]) {
            println!("{:?}", outputs);
        }
    }
}
