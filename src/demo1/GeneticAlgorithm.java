package demo1;

import java.util.Arrays;
import java.util.Random;

public class GeneticAlgorithm {
    private int populationSize;
    private double mutationRate;
    private double crossoverRate;
    private int elitismCount;
    private boolean binaryEncoding; // 是否使用二进制编码

    public GeneticAlgorithm(int populationSize, double mutationRate, double crossoverRate, int elitismCount, boolean binaryEncoding) {
        this.populationSize = populationSize;
        this.mutationRate = mutationRate;
        this.crossoverRate = crossoverRate;
        this.elitismCount = elitismCount;
        this.binaryEncoding = binaryEncoding;
    }

    public Population evolvePopulation(Population population) {
        Population newPopulation = new Population(population.size(), binaryEncoding);

        for (int i = 0; i < elitismCount; i++) {
            newPopulation.saveIndividual(i, population.getFittest(i));
        }

        for (int i = elitismCount; i < population.size(); i++) {
            Individual parent1 = selectParent(population);
            Individual parent2 = selectParent(population);

            Individual child = crossover(parent1, parent2);

            newPopulation.saveIndividual(i, child);
        }

        for (int i = elitismCount; i < newPopulation.size(); i++) {
            mutate(newPopulation.getIndividual(i));
        }

        return newPopulation;
    }

    private Individual crossover(Individual parent1, Individual parent2) {
        Individual child = new Individual(binaryEncoding);

        if (binaryEncoding) { // 如果使用二进制编码
            for (int i = 0; i < parent1.size(); i++) {
                if (Math.random() <= crossoverRate) {
                    child.setGene(i, parent1.getGene(i));
                } else {
                    child.setGene(i, parent2.getGene(i));
                }
            }
        } else { // 如果使用实数编码
            for (int i = 0; i < parent1.size(); i++) {
                double value1 = (double) parent1.getGene(i);
                double value2 = (double) parent2.getGene(i);
                double minValue = Math.min(value1, value2);
                double maxValue = Math.max(value1, value2);

                double geneValue = minValue + Math.random() * (maxValue - minValue);

                child.setGene(i, geneValue);
            }
        }

        return child;
    }

    private void mutate(Individual individual) {
        if (binaryEncoding) { // 如果使用二进制编码
            for (int i = 0; i < individual.size(); i++) {
                if (Math.random() <= mutationRate) {
                    char gene = (char) individual.getGene(i);

                    if (gene == '0') {
                        gene = '1';
                    } else {
                        gene = '0';
                    }

                    individual.setGene(i, gene);
                }
            }
        } else { // 如果使用实数编码
            for (int i = 0; i < individual.size(); i++) {
                if (Math.random() <= mutationRate) {
                    double gene = (double) individual.getGene(i);
                    gene += Math.random() * 0.2 - 0.1; // 变异范围为[-0.1, 0.1]
                    individual.setGene(i, gene);
                }
            }
        }
    }

    private Individual selectParent(Population population) {
        Random random = new Random();
        int tournamentSize = 5;
        Population tournament = new Population(tournamentSize, binaryEncoding);

        for (int i = 0; i < tournamentSize; i++) {
            int randomIndex = random.nextInt(population.size());
            tournament.saveIndividual(i, population.getIndividual(randomIndex));
        }

        return tournament.getFittest();
    }
}

class Individual {
    private Object[] genes;
    private int fitness;
    private boolean binaryEncoding; // 是否使用二进制编码

    public Individual(boolean binaryEncoding) {
        this.binaryEncoding = binaryEncoding;

        genes = new Object[10];
        if (binaryEncoding) { // 如果使用二进制编码
            char temp = 0;
            for(int i = 0; i < 10; i++){
                genes[i] = temp;
            }
            //this.genes = new char[10]; // 假设每个个体有10个基因
        } else { // 如果使用实数编码
            double temp = 0;
            for(int i = 0; i < 10; i++){
                genes[i] = temp;
            }
            //this.genes = new double[10];
        }
        this.fitness = 0;
        generateIndividual();
    }

    public void generateIndividual() {
        Random random = new Random();

        if (binaryEncoding) { // 如果使用二进制编码
            for (int i = 0; i < size(); i++) {
                //char gene的取值为0或1
                char gene = random.nextInt(2) == 0 ? '0' : '1';
                genes[i] = gene;
            }
        } else { // 如果使用实数编码
            for (int i = 0; i < size(); i++) {
                double gene = random.nextDouble() * 10.0; // 假设取值范围为[0, 10]
                genes[i] = gene;
            }
        }
    }

    public Object getGene(int index) {
        return genes[index];
    }

    public void setGene(int index, Object value) {
        genes[index] = value;
        fitness = 0;
    }

    public int size() {
        return genes.length;
    }

    public int getFitness() {
        if (fitness == 0) {
            if (binaryEncoding) { // 如果使用二进制编码
                for (int i = 0; i < size(); i++) {
                    if ((char) genes[i] == '1') {
                        fitness++;
                    }
                }
            } else { // 如果使用实数编码
                double sum = 0.0;

                for (int i = 0; i < size(); i++) {
                    sum += (double) genes[i];
                }

                fitness = (int) Math.round(Math.abs(sum - 10));
            }
        }

        return fitness;
    }

    @Override
    public String toString() {
        return Arrays.toString(genes);
    }
}

class Population {
    private Individual[] individuals;
    private boolean binaryEncoding; // 是否使用二进制编码

    public Population(int size, boolean binaryEncoding) {
        this.binaryEncoding = binaryEncoding;
        this.individuals = new Individual[size];

        for (int i = 0; i < size(); i++) {
            Individual individual = new Individual(binaryEncoding);
            saveIndividual(i, individual);
        }
    }

    public Individual getIndividual(int index) {
        return individuals[index];
    }

    public void saveIndividual(int index, Individual individual) {
        individuals[index] = individual;
    }

    public Individual getFittest() {
        Individual fittest = individuals[0];

        for (int i = 1; i < size(); i++) {
            if (fittest.getFitness() <= getIndividual(i).getFitness()) {
                fittest = getIndividual(i);
            }
        }

        return fittest;
    }

    public Individual getFittest(int index) {
        Arrays.sort(individuals, (a, b) -> b.getFitness() - a.getFitness());
        return individuals[index];
    }

    public int size() {
        return individuals.length;
    }
}