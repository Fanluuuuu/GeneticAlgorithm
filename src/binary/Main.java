package binary;

/**

 主程序类
 */
public class Main {
    public static void main(String[] args) {
        int populationSize = 50;
        double mutationRate = 0.01;
        double crossoverRate = 0.9;
        int elitismCount = 2;
        FitnessFunction fitnessFunction = new SinFitnessFunction();
        SelectionOperator selectionOperator = new RouletteWheelSelection();
        CrossoverOperator crossoverOperator = new SinglePointCrossover();
        MutationOperator mutationOperator = new FlipMutation();

        GeneticAlgorithm geneticAlgorithm = new GeneticAlgorithm(populationSize, mutationRate, crossoverRate,
                elitismCount, fitnessFunction, selectionOperator, crossoverOperator, mutationOperator);

        Population population = new Population(populationSize, 20);

        int generations = 0;
        while (generations < 100 && population.getFittest().getFitness() < 0.99999) {
            population.calculateFitness(fitnessFunction);
            System.out.println("Generation: " + generations + "\tFittest: " + population.getFittest().getFitness());

            population = geneticAlgorithm.evolvePopulation(population);

            generations++;
        }

        System.out.println("\nSolution found in generation " + population.getGeneration());
        System.out.println("Fitness: " + population.getFittest().getFitness());
        System.out.println("Phenotype: " + population.getFittest().getPhenotype());
    }
}
