package decimal;


/**
 * 主程序类
 */
public class DecimalMain {
    public void decimal(){
        int populationSize = 100;
        int geneLength = 1;
        // 遗传算法参数设置
        double crossoverRate = 0.9;
        double mutationRate = 0.01;
        //选择不同的适应度函数
        Population population = new Population(populationSize, geneLength, new SinFitnessFunction());
        //选择方式
        SelectionOperator selectionOperator = new RouletteWheelSelection();
        //交叉方式
        CrossoverOperator crossoverOperator = new SinglePointCrossover();
        //变异方式
        MutationOperator mutationOperator = new FlipMutation();
        // 记录全局最优个体和适应度值
        Individual globalBestIndividual = null;
        double globalMaxFitness = Double.NEGATIVE_INFINITY;
        // 初始化上一代最优个体
        Individual bestIndividual = population.getFittest();
        double maxFitness = bestIndividual.getFitness();

        // 遗传算法主循环
        int maxGeneration = 0;
        double currentMaxFitness = 0;
        int i = 0;

        while(maxGeneration < 100 && currentMaxFitness <0.99999){
            // 进化种群
            population.evolve(crossoverRate, mutationRate, selectionOperator, crossoverOperator, mutationOperator);

            // 计算当前一代种群适应度值之和
            double totalFitness = population.getTotalFitness();

            // 获取当前一代中最优个体
            Individual currentBestIndividual = population.getFittest();
            currentMaxFitness = currentBestIndividual.getFitness();

            // 判断当前一代中最优个体和上一代最优个体哪个更优，将更优的那个作为下一代最优个体
            if (currentMaxFitness > maxFitness) {
                bestIndividual = currentBestIndividual;
                maxFitness = currentMaxFitness;
            } else {
                currentBestIndividual = bestIndividual;
                currentMaxFitness = maxFitness;
            }

            // 如果当前一代中出现更优的个体，则将其作为全局最优个体
            if (currentMaxFitness > globalMaxFitness) {
                globalBestIndividual = currentBestIndividual;
                globalMaxFitness = currentMaxFitness;
            }

            // 输出当前一代中最优个体和全局最优个体
            System.out.println("Generation " + (i++) + ": Best fitness = " + currentMaxFitness);
            maxGeneration++;
        }

        // 使用K均值聚类算法作为适应度函数进行测试
        GeneticAlgorithm ga2 = new GeneticAlgorithm(100, new KMeansFitnessFunction(), 0.8, 0.1);
        Individual bestIndividual2 = ga2.run(100);

        // 使用梯度下降算法作为适应度函数进行测试
        GeneticAlgorithm ga3 = new GeneticAlgorithm(100, new GradientDescentFitnessFunction(), 0.8, 0.1);
        Individual bestIndividual3 = ga3.run(100);

        System.out.println("Best individual (K-means): " + bestIndividual2.getFitness());
        System.out.println("Best individual (gradient descent): " + bestIndividual3.getFitness());


    }
}
