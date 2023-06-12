package demo2;

/**

 主程序入口
 */
public class Main {
    private static final int POPULATION_SIZE = 50; // 种群大小
    private static final double MUTATION_RATE = 0.1; // 变异率
    private static final double CROSSOVER_RATE = 0.95; // 交叉率
    private static final int ELITISM_COUNT = 2; // 精英数量
    private static final int GENE_LENGTH = 10; // 基因长度
    private static final int MAX_GENERATIONS = 100; // 最大代数

    public static void main(String[] args) {
// 创建遗传算法对象
        GeneticAlgorithm ga = new GeneticAlgorithm(POPULATION_SIZE, MUTATION_RATE, CROSSOVER_RATE, ELITISM_COUNT,
                GENE_LENGTH, FitnessCalc::getFitness, new RouletteWheelSelection(), new SinglePointCrossover(),
                new FlipMutation());

        // 创建初始种群
        Population population = new Population(POPULATION_SIZE, GENE_LENGTH);

        // 进化并输出每一代最优值
        for (int i = 0; i < MAX_GENERATIONS; i++) {
            population = ga.evolvePopulation(population);
            population.increaseGeneration();

            System.out.println("Generation " + population.getGeneration() + ": Best fitness = " +
                    population.getFittest().getFitness());
        }

        System.out.println("Solution found in generation " + population.getGeneration());
        System.out.println("Best individual: " + population.getFittest().getDecimalValue());
    }
}
