package binary;

import java.util.*;

/**
 * 适应度函数接口
 */
interface FitnessFunction {
    double calculateFitness(Individual individual);
}

/**
 * sin(x)函数的适应度函数实现类
 */
class SinFitnessFunction implements FitnessFunction {
    public double calculateFitness(Individual individual) {
        double x = individual.getPhenotype();
        return Math.sin(x);
    }
}

/**
 *x*Math.sin(10*Math.PI*2)+2.0
 */
class Sin10PIFitnessFunction implements FitnessFunction {
    public double calculateFitness(Individual individual) {
        double x = individual.getPhenotype();

        //x*sin(10*PI*x)+2.0
        return x*Math.sin(10*Math.PI*2)+2.0;
    }
}

/**
 * 选择操作接口
 */
interface SelectionOperator {
    Individual select(Population population);
}

/**
 * 简单随机选择实现类
 */
class SimpleRandomSelection implements SelectionOperator {
    public Individual select(Population population) {
        Random random = new Random();
        int index = random.nextInt(population.size());
        return population.getIndividual(index);
    }
}

/**
 * 轮盘赌选择实现类
 */
class RouletteWheelSelection implements SelectionOperator {
    public Individual select(Population population) {
        double totalFitness = population.getTotalFitness();
        double value = Math.random() * totalFitness;
        double sum = 0;
        for (int i = 0; i < population.size(); i++) {
            sum += population.getIndividual(i).getFitness();
            if (sum > value) {
                return population.getIndividual(i);
            }
        }
        return null;
    }
}

/**
 * 交叉操作接口
 */
interface CrossoverOperator {
    Individual crossover(Individual parent1, Individual parent2);
}

/**
 * 单点交叉实现类
 */
class SinglePointCrossover implements CrossoverOperator {
    public Individual crossover(Individual parent1, Individual parent2) {
        Individual child = new Individual(parent1.getLength());
        int crossoverPoint = (int) (Math.random() * parent1.getLength());
        for (int i = 0; i < parent1.getLength(); i++) {
            if (i < crossoverPoint) {
                child.setGene(i, parent1.getGene(i));
            } else {
                child.setGene(i, parent2.getGene(i));
            }
        }
        return child;
    }
}

/**
 * 变异操作接口
 */
interface MutationOperator {
    void mutate(Individual individual, double mutationRate);
}

/**
 * 翻转变异实现类
 */
class FlipMutation implements MutationOperator {
    public void mutate(Individual individual, double mutationRate) {
        Random random = new Random();
        for (int i = 0; i < individual.getLength(); i++) {
            if (Math.random() < mutationRate) {
                int gene = individual.getGene(i);
                if (gene == 0) {
                    individual.setGene(i, 1);
                } else {
                    individual.setGene(i, 0);
                }
            }
        }
    }
}

/**
 * 个体类
 */
class Individual {
    private int[] genes; // 基因序列
    private double fitness; // 适应度值

    public Individual(int length) {
        genes = new int[length];
        for (int i = 0; i < length; i++) {
            genes[i] = (int) (Math.random() * 2);
        }
    }

    public Individual(int[] genes) {
        this.genes = genes;
    }

    public void setGenes(int[] genes){
        this.genes = genes;
    }

    public int getGene(int index) {
        return genes[index];
    }

    public void setGene(int index, int value) {
        genes[index] = value;
    }

    public int[] getGenes() {
        return genes;
    }

    public double getPhenotype() {
        int decimalValue = 0;
        for (int i = 0; i < genes.length; i++) {
            decimalValue += genes[i] * Math.pow(2, i);
        }
        double minValue = -10.0;
        double maxValue = 10.0;
        return minValue + (maxValue - minValue) * decimalValue / (Math.pow(2, genes.length) - 1);
    }

    public double getFitness() {
        return fitness;
    }

    public void setFitness(double fitness) {
        this.fitness = fitness;
    }

    public int getLength() {
        return genes.length;
    }

    public void setGeneLength(int length) {
        if (length > genes.length) {
            int[] newGenes = new int[length];
            System.arraycopy(genes, 0, newGenes, 0, genes.length);
            for (int i = genes.length; i < length; i++) {
                newGenes[i] = (int) (Math.random() * 2);
            }
            genes = newGenes;
        }
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < genes.length; i++) {
            sb.append(genes[i]);
        }
        return sb.toString();
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null || !(obj instanceof Individual)) {
            return false;
        }
        Individual other = (Individual) obj;
        return Arrays.equals(genes, other.genes);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(genes);
    }
}

/**
 * 种群类
 */
class Population {
    private Individual[] individuals; // 个体数组
    private double totalFitness; // 总适应度值
    private int previousFittest; // 上一代最优解的适应度值
    private int geneLength; // 基因长度
    private int generation; // 当前世代

    public Population(int size, int geneLength) {
        individuals = new Individual[size];
        for (int i = 0; i < size; i++) {
            individuals[i] = new Individual(geneLength);
        }
        totalFitness = 0;
        previousFittest = 0;
        this.geneLength = geneLength;
        generation = 0;
    }

    public Individual getIndividual(int index) {
        return individuals[index];
    }

    public void saveIndividual(int index, Individual individual) {
        individuals[index] = individual;
    }

    public Individual[] getIndividuals() {
        return individuals;
    }

    public double getTotalFitness() {
        return totalFitness;
    }

    public int getPreviousFittest() {
        return previousFittest;
    }

    public void calculateFitness(FitnessFunction fitnessFunction) {
        totalFitness = 0;
        for (Individual individual : individuals) {
            double fitness = fitnessFunction.calculateFitness(individual);
            individual.setFitness(fitness);
            totalFitness += fitness;
        }
        Arrays.sort(individuals, Collections.reverseOrder(Comparator.comparingDouble(Individual::getFitness)));
        previousFittest = (int) getFittest().getFitness();
    }

    public Individual getFittest() {
        return individuals[0];
    }

    public Individual getFittest(int offset) {
        return individuals[offset];
    }

    public int size() {
        return individuals.length;
    }

    public int getGeneLength() {
        return geneLength;
    }

    public void increaseGeneLength() {
        geneLength++;
        for (Individual individual : individuals) {
            int[] genes = Arrays.copyOf(individual.getGenes(), geneLength);
            for (int i = individual.getLength(); i < geneLength; i++) {
                genes[i] = (int) (Math.random() * 2);
            }
            individual.setGeneLength(geneLength);
            individual.setGenes(genes);
        }
    }

    public int getGeneration() {
        return generation;
    }

    public void incrementGeneration() {
        generation++;
    }
}

/**

 遗传算法类
 */
class GeneticAlgorithm {
    private int populationSize; // 种群大小
    private double mutationRate; // 变异率
    private double crossoverRate; // 交叉率
    private int elitismCount; // 精英数量
    private FitnessFunction fitnessFunction; // 适应度函数
    private SelectionOperator selectionOperator; // 选择操作
    private CrossoverOperator crossoverOperator; // 交叉操作
    private MutationOperator mutationOperator; // 变异操作

    public GeneticAlgorithm(int populationSize, double mutationRate, double crossoverRate, int elitismCount,
                            FitnessFunction fitnessFunction, SelectionOperator selectionOperator,
                            CrossoverOperator crossoverOperator, MutationOperator mutationOperator) {
        this.populationSize = populationSize;
        this.mutationRate = mutationRate;
        this.crossoverRate = crossoverRate;
        this.elitismCount = elitismCount;
        this.fitnessFunction = fitnessFunction;
        this.selectionOperator = selectionOperator;
        this.crossoverOperator = crossoverOperator;
        this.mutationOperator = mutationOperator;
    }

    /**

     进化种群
     */
    public Population evolvePopulation(Population population) {
        Population newPopulation = new Population(populationSize, population.getGeneLength());

// 保留精英
        for (int i = 0; i < elitismCount; i++) {
            newPopulation.saveIndividual(i, population.getFittest(i));
        }

// 繁殖新个体并添加到新种群中
        for (int i = elitismCount; i < population.size(); i++) {
            Individual parent1 = selectionOperator.select(population);
            Individual parent2 = selectionOperator.select(population);
            Individual child = crossoverOperator.crossover(parent1, parent2);
            newPopulation.saveIndividual(i, child);
        }

// 变异除了精英以外的所有个体
        for (int i = elitismCount; i < newPopulation.size(); i++) {
            mutationOperator.mutate(newPopulation.getIndividual(i), mutationRate);
        }

        newPopulation.calculateFitness(fitnessFunction);
        newPopulation.incrementGeneration();

        return newPopulation;
    }
}