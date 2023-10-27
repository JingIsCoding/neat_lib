# NEAT in rust

Evolving Neural Networks through Augmenting Topologies (NEAT) is a method to dynamically evole a neural network to maximize the performance (fitness) for certain tasks, the original paper can be found [here](https://nn.cs.utexas.edu/downloads/papers/stanley.ec02.pdf). 
This is a rust implementation that is inspired by the original papar as well as [evo-NEAT](https://github.com/vishnugh/evo-NEAT) and the [neat-python](https://github.com/CodeReclaimers/neat-python) project.

This is an experimental project, for anyone who is interesting in ML and Rust language. So feel free to point out errors in the model and contribute.

# Get started

The simplest model to check the effectiveness of the model is xor operation as it is a non linear problem.\

Run the examples:\
```cargo run --example xor```

You might have to try multiple time for the network to find the solution, as it does not alway evolve the necessary node to make it non-learn.

# Concepts
The most unique and facinating characteristic of a NEAT network is that instead of a predefined network topology, the network itself can evolve new struture and change the parameters within the network.

This evolulation happens in the process of mutation, and since the mutation happens slowly and grows in uncertain direction. The network has to run a lot of genomes at the same time so that we will have a higher probability of landing on a solution.

The mutation happens in a process called "crossover", where we randomly take genes and connections from two "parent" genomes to generate a new genome and then mutates it.

As most the mutation are harmful at beginning, this network adopt a concept of species to protect the new genome from being wiped out immediately, which encourages exploration.

Most of the concepts are from the original paper, except for a few tweaks.
1. Gene and Connection gene can be removed.
2. Each gene has a responsiveness and bias attached.

## Gene
A basic unit in a genome, each gene has a unique key so that it will always represent certain features.

## Connection
Connection represent the connection between two genes, it also contains the weight between two genes.

## Innovation number
Each connection will get a incremental innovation number, so that we can distingulish different connections in crossover process and speciation.

## Genome
Genome is a collection of gene and connections.

## Species
Species is a collectoin of genomes that are similar enough to each other.

## Population
Population is a collection of species.
