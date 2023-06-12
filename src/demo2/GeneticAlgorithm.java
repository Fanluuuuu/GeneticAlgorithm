package demo2;

import java.util.Random;

/**
 * 工厂模式-遗传算法
 */
public class GeneticAlgorithm {
    private int populationSize; // 种群大小
    private double mutationRate; // 变异率
    private double crossoverRate; // 交叉率
    private int elitismCount; // 精英数量
    private int geneLength; // 基因长度
    private FitnessFunction fitnessFunction; // 适应度函数
    private SelectionOperator selectionOperator; // 选择操作
    private CrossoverOperator crossoverOperator; // 交叉操作
    private MutationOperator mutationOperator; // 变异操作

    public GeneticAlgorithm(int populationSize, double mutationRate, double crossoverRate, int elitismCount,
                            int geneLength, FitnessFunction fitnessFunction, SelectionOperator selectionOperator,
                            CrossoverOperator crossoverOperator, MutationOperator mutationOperator) {
        this.populationSize = populationSize;
        this.mutationRate = mutationRate;
        this.crossoverRate = crossoverRate;
        this.elitismCount = elitismCount;
        this.geneLength = geneLength;
        this.fitnessFunction = fitnessFunction;
        this.selectionOperator = selectionOperator;
        this.crossoverOperator = crossoverOperator;
        this.mutationOperator = mutationOperator;
    }

    /**
     * 进化种群
     */
    public Population evolvePopulation(Population population) {
        Population newPopulation = new Population(populationSize, geneLength);

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
            mutationOperator.mutate(newPopulation.getIndividuals()[i], mutationRate);
        }

        // 基因长度随时间变化
        if ((population.getGeneration() % 10) == 0 && population.getFittest().getFitness() == population.getPreviousFittest()) {
            newPopulation.increaseGeneLength();
        }

        return newPopulation;
    }
}

/**
 * 适应度函数工厂
 */
interface FitnessFunction {
    double calculateFitness(Individual individual);
}

class SinFitnessFunction implements FitnessFunction {
    @Override
    public double calculateFitness(Individual individual) {
        double value = individual.getDecimalValue();
        return value * Math.sin(10 * Math.PI * value) + 2.0;
    }
}

/**
 * 选择操作工厂
 */
interface SelectionOperator {
    Individual select(Population population);
}

class RouletteWheelSelection implements SelectionOperator {
    @Override
    public Individual select(Population population) {
        double fitnessSum = 0;
        for (Individual individual : population.getIndividuals()) {
            fitnessSum += individual.getFitness();
        }

        Random random = new Random();
        double roulettePosition = random.nextDouble() * fitnessSum;

        double spinWheel = 0;
        for (Individual individual : population.getIndividuals()) {
            spinWheel += individual.getFitness();
            if (spinWheel >= roulettePosition) {
                return individual;
            }
        }

        return null;
    }
}

/**
 * 交叉操作工厂
 */
interface CrossoverOperator {
    Individual crossover(Individual parent1, Individual parent2);
}

class SinglePointCrossover implements CrossoverOperator {
    @Override
    public Individual crossover(Individual parent1, Individual parent2) {
        Random random = new Random();
        int splitPoint = random.nextInt(parent1.getGenes().length);

        char[] childGenes = new char[parent1.getGenes().length];

        for (int i = 0; i < parent1.getGenes().length; i++) {
            if (i < splitPoint) {
                childGenes[i] = parent1.getGenes()[i];
            } else {
                childGenes[i] = parent2.getGenes()[i];
            }
        }

        return new Individual(childGenes);
    }
}

/**

 变异操作工厂 */
interface MutationOperator {
     void mutate(Individual individual, double mutationRate);
 }
class FlipMutation implements MutationOperator {
    @Override
    public void mutate(Individual individual, double mutationRate) {
        Random random = new Random();
        char[] genes = individual.getGenes();

        for (int i = 0; i < genes.length; i++) {
            if (random.nextDouble() < mutationRate) {
                genes[i] = (genes[i] == '0') ? '1' : '0';
            }
        }

        individual.setGenes(genes);
    }
}

/**

 单个个体
 */
class Individual {
    private char[] genes;
    private double fitness;

    public Individual(char[] genes) {
        this.genes = genes;
        this.fitness = -1;
    }

    /**

     获取基因序列的十进制值 */
    public double getDecimalValue() {
        int decimalValue = 0;
        for (int i = 0; i < genes.length; i++) {
            decimalValue += genes[i] * Math.pow(2, i);
        }
        return decimalValue / (Math.pow(2, genes.length) - 1.0);
    }
    /**

     获取适应度 */
    public double getFitness() {
        if (fitness == -1) {
            fitness = FitnessCalc.getFitness(this);
        } return fitness;
    }
    /**
     设置适应度 */
    public void setFitness(double fitness) {
        this.fitness = fitness;
    }
    public char[] getGenes() {
        return genes;
    }

    public void setGenes(char[] genes) {
        this.genes = genes;
    }
}

/**

 种群
 */
class Population {
    private Individual[] individuals;
    private double previousFittest; // 前一代最优适应度值
    private int geneLength; // 基因长度
    private int generation; // 当前代数

    public Population(int populationSize, int geneLength) {
        individuals = new Individual[populationSize];
        previousFittest = -1;
        this.geneLength = geneLength;
        generation = 0;

        for (int i = 0; i < populationSize; i++) {
            char[] genes = new char[geneLength];
            for (int j = 0; j < geneLength; j++) {
                genes[j] = (Math.random() < 0.5) ? '0' : '1';
            }
            individuals[i] = new Individual(genes);
        }
    }

    /**

     获取最优个体
     */
    public Individual getFittest() {
        Individual fittest = individuals[0];
        for (int i = 1; i < individuals.length; i++) {
            if (individuals[i].getFitness() > fittest.getFitness()) {
                fittest = individuals[i];
            }
        }

// 记录前一代最优值
        if (fittest.getFitness() > previousFittest) {
            previousFittest = fittest.getFitness();
        }

        return fittest;
    }

    /**

     获取前n个最优个体 */
    public Individual getFittest(int n) {
        sortIndividuals();
        return individuals[n];
    }
    /**

     对个体数组进行排序（适应度从高到低） */
    private void sortIndividuals() {
        for (int i = 0; i < individuals.length; i++) {
            for (int j = i + 1; j < individuals.length; j++) {
                if (individuals[i].getFitness() < individuals[j].getFitness()) {
                    Individual temp = individuals[i];
                    individuals[i] = individuals[j];
                    individuals[j] = temp;
                }
            }
        }
    }
    /**
     添加新的个体 */
    public void saveIndividual(int index, Individual individual) {
        individuals[index] = individual;
    }
    /**

     基因长度增加1 */
    public void increaseGeneLength() { geneLength++; for (Individual individual : individuals) { char[] genes = new char[geneLength];
        for (int i = 0; i < individual.getGenes().length; i++) {
            genes[i] = individual.getGenes()[i];
        }
        // 新增的基因位随机设置为'0'或'1'
        genes[individual.getGenes().length] = (Math.random() < 0.5) ? '0' : '1';
        individual.setGenes(genes);
    }
    }

    public int size() {
        return individuals.length;
    }

    public Individual[] getIndividuals() {
        return individuals;
    }

    public double getPreviousFittest() {
        return previousFittest;
    }

    public int getGeneration() {
        return generation;
    }

    public void increaseGeneration() {
        generation++;
    }
}

/**

 适应度计算工具类
 */
class FitnessCalc {
    private static FitnessFunction fitnessFunction = new SinFitnessFunction();

    /**

     获取个体的适应度值 */
    public static double getFitness(Individual individual) {
        return fitnessFunction.calculateFitness(individual);
    }
}
