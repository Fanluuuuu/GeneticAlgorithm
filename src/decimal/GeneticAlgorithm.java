package decimal;

import org.apache.commons.math3.ml.clustering.CentroidCluster;
import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.Clusterable;
import org.apache.commons.math3.ml.clustering.KMeansPlusPlusClusterer;

import java.util.*;

/**
 * 适应度函数接口
 */
interface FitnessFunction {
    double calculateFitness(Individual individual);
}
/**
 * x*sin(10*Π*x)+2.0函数的适应度函数实现类
 */
class SinFitnessFunction implements FitnessFunction {
    public double calculateFitness(Individual individual) {
        double x = individual.getPhenotype();
        return x*Math.sin(10*Math.PI*2)+2.0;
    }
}
/**
 * 聚类（K均值）适应度函数实现类
 */
class KMeansFitnessFunction implements FitnessFunction {
    public double calculateFitness(Individual individual) {
        double[] genes = individual.getGenes();
        // 将基因序列转换为可用于聚类的数据
        DoubleArrayWrapper[] data = new DoubleArrayWrapper[genes.length];
        for (int i = 0; i < genes.length; i++) {
            data[i] = new DoubleArrayWrapper(new double[]{genes[i]});
        }
        // 使用K均值聚类算法将数据聚类为2个簇，并计算每个簇中心点的距离平方和作为适应度值
        KMeansPlusPlusClusterer<DoubleArrayWrapper> clusterer = new KMeansPlusPlusClusterer<>(2);
        List<CentroidCluster<DoubleArrayWrapper>> clusters = clusterer.cluster(Arrays.asList(data));
        double fitness = 0;
        for (org.apache.commons.math3.ml.clustering.Cluster<DoubleArrayWrapper> cluster : clusters) {
            DoubleArrayWrapper center = getCenter(cluster);
            double sum = 0;
            for (DoubleArrayWrapper point : cluster.getPoints()) {
                sum += Math.pow(point.getPoint()[0] - center.getPoint()[0], 2);
            }
            fitness += sum;
        }
        return fitness;
    }
    /**
     * 计算簇的中心点
     */
    private DoubleArrayWrapper getCenter(Cluster<DoubleArrayWrapper> cluster) {
        double[] center = new double[cluster.getPoints().get(0).getPoint().length];
        for (DoubleArrayWrapper point : cluster.getPoints()) {
            double[] p = point.getPoint();
            for (int i = 0; i < p.length; i++) {
                center[i] += p[i];
            }
        }
        for (int i = 0; i < center.length; i++) {
            center[i] /= cluster.getPoints().size();
        }
        return new DoubleArrayWrapper(center);
    }
}

/**
 * 梯度下降适应度函数实现类
 */
class GradientDescentFitnessFunction implements FitnessFunction {
    public double calculateFitness(Individual individual) {
        double[] genes = individual.getGenes();
        // 解码基因序列为模型参数
        double w = genes[0];
        double b = genes.length > 1 ? genes[1] : 0.0;
        // 构建模型，这里以y=w*x+b为例
        GradientDescentFitnessFunction.Model model = new GradientDescentFitnessFunction.Model(w, b);
        // 训练模型，并计算损失函数作为适应度值
        double[][] xData = new double[100][1];
        double[][] yData = new double[100][1];
        for (int i = 0; i < 100; i++) {
            double x = -1 + 0.02 * i;
            xData[i][0] = x;
            yData[i][0] = Math.sin(x);
        }
        GradientDescentFitnessFunction.Trainer trainer = new GradientDescentFitnessFunction.Trainer(xData, yData);
        double[] params = trainer.train(model, 1000, 0.01);
        double fitness = trainer.loss(model, params);
        return fitness;
    }

    /**
     * 线性模型类
     */
    class Model {
        private double w;
        private double b;

        public Model(double w, double b) {
            this.w = w;
            this.b = b;
        }

        public double predict(double x) {
            return w * x + b;
        }
    }

    /**
     * 模型训练器类
     */
    class Trainer {
        private double[][] xData;
        private double[][] yData;

        public Trainer(double[][] xData, double[][] yData) {
            this.xData = xData;
            this.yData = yData;
        }

        /**
         * 使用梯度下降算法在训练集上训练模型，并返回最优的模型参数
         */
        public double[] train(GradientDescentFitnessFunction.Model model, int maxIter, double learningRate) {
            double[] params = new double[]{model.w, model.b};
            for (int i = 0; i < maxIter; i++) {
                double[] gradient = gradient(model, params);
                params[0] -= learningRate * gradient[0];
                params[1] -= learningRate * gradient[1];
            }
            return params;
        }

        /**
         * 计算模型在训练集上的平均损失函数值
         */
        public double loss(GradientDescentFitnessFunction.Model model, double[] params) {
            double sum = 0;
            for (int i = 0; i < xData.length; i++) {
                double x = xData[i][0];
                double y = yData[i][0];
                double yPred = model.predict(x);
                sum += Math.pow(y - yPred, 2);
            }
            return sum / xData.length;
        }

        /**
         * 计算模型在训练集上的梯度向量
         */
        public double[] gradient(GradientDescentFitnessFunction.Model model, double[] params) {
            double sumW = 0;
            double sumB = 0;
            for (int i = 0; i < xData.length; i++) {
                double x = xData[i][0];
                double y = yData[i][0];
                double yPred = model.predict(x);
                sumW += (yPred - y) * x;
                sumB += (yPred - y);
            }
            double[] gradient = new double[]{sumW / xData.length, sumB / xData.length};
            return gradient;
        }
    }
}

/**
 * Double[] 的包装类，实现 Clusterable 接口
 */
class DoubleArrayWrapper implements Clusterable {
    private final double[] point;

    public DoubleArrayWrapper(double[] point) {
        this.point = point;
    }

    public double[] getPoint() {
        return point;
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
 * 邻居交叉实现类
 */
class NeighborCrossover implements CrossoverOperator {
    private int neighborSize;
    public NeighborCrossover(int neighborSize) {
        this.neighborSize = neighborSize;
    }
    public Individual crossover(Individual parent1, Individual parent2) {
        Individual child = new Individual(parent1.getLength());
        int crossoverPoint = (int) (Math.random() * parent1.getLength());

        // 邻居交叉操作
        for (int i = 0; i < parent1.getLength(); i++) {
            if (i == crossoverPoint) {
                // 取得当前基因的邻居
                List<Integer> neighbors = getNeighbors((int) parent1.getGene(i), parent2, i);
                // 从邻居中随机选择一个
                int selectedNeighbor = neighbors.get((int)(Math.random()*neighbors.size()));
                child.setGene(i, selectedNeighbor);
            } else {
                child.setGene(i, parent1.getGene(i));
            }
        }
        return child;
    }
    /**
     * 获取基因邻居列表
     */
    private List<Integer> getNeighbors(int geneValue, Individual parent, int currentIndex) {
        List<Integer> neighbors = new ArrayList<>();
        int start = Math.max(0, currentIndex - neighborSize);
        int end = Math.min(parent.getLength(), currentIndex + neighborSize + 1);
        for (int i = start; i < end; i++) {
            if (parent.getGene(i) != geneValue) {
                neighbors.add((int) parent.getGene(i));
            }
        }
        return neighbors;
    }
}
/**
 * 加权平均交叉实现类
 */
class WeightedAverageCrossover implements CrossoverOperator {
    private double alpha;
    public WeightedAverageCrossover(double alpha) {
        this.alpha = alpha;
    }
    public Individual crossover(Individual parent1, Individual parent2) {
        Individual child = new Individual(parent1.getLength());
        for (int i = 0; i < parent1.getLength(); i++) {
            double gene1 = parent1.getGene(i);
            double gene2 = parent2.getGene(i);
            double newGene = alpha * gene1 + (1 - alpha) * gene2;
            child.setGene(i, newGene);
        }
        return child;
    }
}
/**
 * 随机平均交叉实现类
 */
class RandomAverageCrossover implements CrossoverOperator {
    private double prob;

    public RandomAverageCrossover(double prob) {
        this.prob = prob;
    }

    public Individual crossover(Individual parent1, Individual parent2) {
        Individual child = new Individual(parent1.getLength());
        for (int i = 0; i < parent1.getLength(); i++) {
            double gene1 = parent1.getGene(i);
            double gene2 = parent2.getGene(i);
            double newGene = Math.random() < prob ? gene1 : gene2;
            child.setGene(i, newGene);
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

    private double phenotype; // 表现型
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

    public Individual(double[] genes, double phenotype) {
        this.genes = genes;
        this.phenotype = phenotype;
    }

    /**
     * 单点交叉
     */
    public Individual crossover(Individual other, double rate) {
        if (Math.random() < rate) {
            int index = (int) (Math.random() * genes.length);
            double[] childGenes = new double[genes.length];
            for (int i = 0; i < genes.length; i++) {
                if (i < index) {
                    childGenes[i] = genes[i];
                } else {
                    childGenes[i] = other.genes[i];
                }
            }
            double childPhenotype = childGenes[0];
            return new Individual(childGenes, childPhenotype);
        } else {
            return Math.random() < 0.5 ? this : other;
        }
    }

    /**
     * 基因翻转变异
     */
    public void mutate(double rate) {
        for (int i = 0; i < genes.length; i++) {
            if (Math.random() < rate) {
                genes[i] = -genes[i];
            }
        }
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

class GeneticAlgorithm {
    private int populationSize; // 种群大小
    private FitnessFunction fitnessFunction; // 适应度函数
    private double crossoverRate; // 交叉概率
    private double mutationRate; // 变异概率

    public GeneticAlgorithm(int populationSize, FitnessFunction fitnessFunction, double crossoverRate, double mutationRate) {
        this.populationSize = populationSize;
        this.fitnessFunction = fitnessFunction;
        this.crossoverRate = crossoverRate;
        this.mutationRate = mutationRate;
    }

    /**
     * 运行遗传算法，并返回适应度最好的个体
     */
    public Individual run(int maxGenerations) {
        List<Individual> population = initializePopulation();
        for (int i = 0; i < maxGenerations; i++) {
            List<Individual> newPopulation = new ArrayList<>();
            while (newPopulation.size() < populationSize) {
                // 轮盘赌选择
                Individual parent1 = rouletteWheelSelection(population);
                Individual parent2 = rouletteWheelSelection(population);

                // 单点交叉
                Individual child = parent1.crossover(parent2, crossoverRate);

                // 基因翻转变异
                child.mutate(mutationRate);

                // 计算适应度
                double fitness = fitnessFunction.calculateFitness(child);
                child.setFitness(fitness);

                newPopulation.add(child);
            }
            population = newPopulation;

            // 输出每代种群中适应度最好的个体
            Individual bestIndividual = getBestIndividual(population);
        }

        // 返回适应度最好的个体
        return getBestIndividual(population);
    }

    /**
     * 初始化种群
     */
    private List<Individual> initializePopulation() {
        List<Individual> population = new ArrayList<>();
        for (int i = 0; i < populationSize; i++) {
            double[] genes = new double[10];
            for (int j = 0; j < genes.length; j++) {
                genes[j] = Math.random() * 10 - 5;
            }
            double phenotype = genes[0];
            population.add(new Individual(genes, phenotype));
        }
        return population;
    }

    /**
     * 轮盘赌选择
     */
    private Individual rouletteWheelSelection(List<Individual> population) {
        double sumFitness = 0;
        for (Individual individual : population) {
            sumFitness += individual.getFitness();
        }
        double rand = Math.random() * sumFitness;
        double currentFitness = 0;
        for (Individual individual : population) {
            currentFitness += individual.getFitness();
            if (currentFitness >= rand) {
                return individual;
            }
        }
        // 如果轮盘赌选择失败，则返回种群中适应度最好的个体
        return getBestIndividual(population);
    }

    /**
     * 获取种群中适应度最好的个体
     */
    private Individual getBestIndividual(List<Individual> population) {
        Individual bestIndividual = population.get(0);
        for (Individual individual : population) {
            if (individual.getFitness() > bestIndividual.getFitness()) {
                bestIndividual = individual;
            }
        }
        return bestIndividual;
    }
}