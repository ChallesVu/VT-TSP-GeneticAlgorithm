# VT-TSP-GeneticAlgorithm
Solve the Traveling Salesman Problem by building a genetic algorithm.
# Traveling Salesman Problem Solver

## Overview

This project provides a solution to the Traveling Salesman Problem (TSP) using a Genetic Algorithm. The solution includes components for generating test data, implementing the genetic algorithm, and solving the TSP.

## Components

1. **DataGenerator:** A utility class for generating test data for the TSP.
2. **GeneticAlgorithm Interface:** Defines basic genetic algorithm operations.
3. **TSPGeneticAlgorithm:** Implements genetic algorithm operations for solving the TSP.
4. **TSPSolver:** A class for solving the TSP using the TSPGeneticAlgorithm.

## Usage

To use the TSP solver:

1. Generate test data using the `DataGenerator`.
2. Read data from the generated file.
3. Set parameters such as the number of iterations, cities, and population size.
4. Create an instance of `TSPGeneticAlgorithm` with random probabilities for crossover and mutation.
5. Create an instance of `TSPSolver` and solve the TSP problem using the `solve` method.

Example:

```java
// Example usage in the main method
DataGenerator data = new DataGenerator();
data.generateTestData("cities.csv");

// Read data from the generated file
String fileName = "cities.csv";
List<String> lines = FileReader.readLines(fileName);

// Set parameters
int iterations = 100;
int cities = lines.size() - 1;
int populationSize = lines.size() - 1;

// Generate random probabilities for crossover and mutation
double crossoverProbability = Math.random();
double mutationProbability = Math.random();

// Create an instance of TSPGeneticAlgorithm
TSPGeneticAlgorithm tspGeneticAlgorithm = new TSPGeneticAlgorithm(crossoverProbability, mutationProbability, lines);

// Create an instance of TSPSolver and solve the TSP problem
TSPSolver tspSolver = new TSPSolver(tspGeneticAlgorithm);
tspSolver.solve(iterations, populationSize, cities);

## Dependencies

No external dependencies are required for this project.

## Contributing

Contributions are welcome! Feel free to open issues or submit pull requests.
