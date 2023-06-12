package demo1;

public class Main {
    public static void main(String[] args) {
        boolean binaryEncoding = true; // 是否使用二进制编码（true表示使用，false表示使用实数编码）

        GeneticAlgorithm ga = new GeneticAlgorithm(50, 0.01, 0.95, 2, binaryEncoding);
        Population population = new Population(50, binaryEncoding);

        int generationCount = 0;
        int maxGenerations = 1000;

        while (population.getFittest().getFitness() < 10 && generationCount < maxGenerations) {
            population = ga.evolvePopulation(population);
            generationCount++;
            System.out.println("Generation: " + generationCount + ", Fittest: " + population.getFittest().getFitness());
        }
    }
}
