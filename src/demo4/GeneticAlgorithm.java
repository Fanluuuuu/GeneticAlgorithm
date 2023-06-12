package demo4;

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
                double gene = individual.getGene(i);
                // 对于每个基因，以50%的概率进行增量或减量变异
                if (Math.random() < 0.5){
                    individual.setGene(i, gene + Math.random());
                }else{
                    individual.setGene(i, gene - Math.random());
                }
            }
        }
    }
}

/**
 * 个体类
 */
class Individual {
    private double[] genes; // 基因序列
    private double fitness; // 适应度值

    public Individual(int length) {
        genes = new double[length];
        for (int i = 0; i < length; i++) {
            genes[i] = -10 + Math.random() * 20; // 实数编码，范围在[-10, 10]
        }
    }

    public Individual(double[] genes) {
        this.genes = genes;
    }

    public void setGenes(double[] genes) {
        this.genes = genes;
    }

    public double getGene(int index) {
        return genes[index];
    }

    public void setGene(int index, double value) {
        genes[index] = value;
    }

    public double[] getGenes() {
        return genes;
    }

    public double getPhenotype() {
        // 实数编码中的基因序列直接就是表现型，无需转换
        return genes[0];
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
            double[] newGenes = new double[length];
            System.arraycopy(genes, 0, newGenes, 0, genes.length);
            for (int i = genes.length; i < length; i++) {
                newGenes[i] = -10 + Math.random() * 20;
            }
            genes = newGenes;
        }
    }
}
/**
 * 种群类
 */
class Population {
    private List<Individual> individuals; // 个体列表
    private FitnessFunction fitnessFunction; // 适应度函数

    public Population(int size, int length, FitnessFunction fitnessFunction) {
        individuals = new ArrayList<>();
        for (int i = 0; i < size; i++) {
            Individual individual = new Individual(length);
            individuals.add(individual);
        }
        this.fitnessFunction = fitnessFunction;
    }

    /**
     * 计算种群中所有个体的适应度值之和
     */
    public double getTotalFitness() {
        double totalFitness = 0;
        for (Individual individual : individuals) {
            individual.setFitness(fitnessFunction.calculateFitness(individual));
            totalFitness += individual.getFitness();
        }
        return totalFitness;
    }

    /**
     * 对种群进行选择、交叉和变异操作，生成新一代种群
     */
    public void evolve(double crossoverRate, double mutationRate,
                       SelectionOperator selectionOperator,
                       CrossoverOperator crossoverOperator,
                       MutationOperator mutationOperator) {
        List<Individual> newIndividuals = new ArrayList<>();
        for (int i = 0; i < individuals.size(); i++) {
            Individual parent1 = selectionOperator.select(this);
            Individual parent2 = selectionOperator.select(this);
            // 父母个体交叉产生子代
            if (Math.random() < crossoverRate) {
                Individual child = crossoverOperator.crossover(parent1, parent2);
                // 对子代进行变异
                mutationOperator.mutate(child, mutationRate);
                newIndividuals.add(child);
            } else { // 如果不交叉，则直接复制父母个体到新一代中
                newIndividuals.add(parent1);
                newIndividuals.add(parent2);
            }
        }
        individuals = newIndividuals;
    }

    /**
     * 获取种群中最优秀的个体
     */
    public Individual getFittest() {
        Individual fittest = individuals.get(0);
        double maxFitness = fittest.getFitness();
        for (int i = 1; i < individuals.size(); i++) {
            double fitness = individuals.get(i).getFitness();
            if (fitness > maxFitness) {
                fittest = individuals.get(i);
                maxFitness = fitness;
            }
        }
        return fittest;
    }

    /**
     * 获取种群大小
     */
    public int size() {
        return individuals.size();
    }

    /**
     * 获取指定索引位置的个体
     */
    public Individual getIndividual(int index) {
        return individuals.get(index);
    }
}