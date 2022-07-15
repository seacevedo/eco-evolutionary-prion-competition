import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

class Simulation {

    private int numMonomers;
    private int numPolymers;
    private SynchronizedMersenneTwister mtRand = new SynchronizedMersenneTwister();
    private double lambda;
    private double d;
    private double beta;
    private double b;
    private double a;
    private int n;
    private double mu;
    private int simid;

    Simulation(int numMonomers, int numPolymers, double lambda, double d, double beta, double b, double a, int n, double mu, int simid){
        this.numMonomers = numMonomers;
        this.numPolymers = numPolymers;
        this.lambda = lambda;
        this.d = d;
        this.beta = beta;
        this.b = b;
        this.a = a;
        this.n = n;
        this.mu = mu;
        this.simid = simid;
    }

    // This function runs the Gillespie's SSA

    void startSimulation() {

        final int MAX_SIZE = 1000000;
        final double MAX_TIME = 10000;
        int initialPolymerSize = 50;
        double currentTime = 0.0;
        double previousTime = 0.0;

        double[] polymerSizes = new double[MAX_SIZE];
        double[] polymerizationValues = new double[MAX_SIZE];
        double[] polymerizationRates = new double[MAX_SIZE];

        double[] breakageRates = new double[MAX_SIZE];
        double[] breakageValues = new double[MAX_SIZE];


        double[] timeArr = new double[MAX_SIZE];
        int[] monomerArr = new int[MAX_SIZE];
        int[] polymerArr = new int[MAX_SIZE];
        int[] massArr = new int[MAX_SIZE];
        double[] sizeArr = new double[MAX_SIZE];
        double[] meanPolyArr = new double[MAX_SIZE];
        double[] varPolyArr = new double[MAX_SIZE];
        double[] meanBreakageArr = new double[MAX_SIZE];
        double[] varBreakageArr = new double[MAX_SIZE];

        double drawProbability;
        int counter = 0;
        Variance var = new Variance();
        Mean m = new Mean();

        for(int i = 0; i < numPolymers; i++){
            polymerSizes[i] = initialPolymerSize;
            polymerizationValues[i] = beta;
            breakageValues[i] = b;
        }

        while(currentTime < MAX_TIME && numPolymers != 0) {


            double totalPolyRate = setPolyRates(numPolymers, numMonomers,lambda+d*numMonomers+a*numPolymers, polymerizationValues, polymerizationRates);
            double totalBreakageRate = setBreakageRates(numPolymers, totalPolyRate, polymerSizes, breakageRates, breakageValues);

            drawProbability = totalBreakageRate * mtRand.nextDouble();
            currentTime += -Math.log(mtRand.nextDouble())/totalBreakageRate;


            if(drawProbability < lambda){
                numMonomers++;
            }

            else if(drawProbability < lambda+d*numMonomers){
                numMonomers--;
            }

            else if(drawProbability < lambda+d*numMonomers+a*numPolymers){

                int index = mtRand.nextInt(numPolymers);

                polymerizationValues[index] = polymerizationValues[numPolymers-1];
                polymerizationValues[numPolymers-1] = 0;

                breakageValues[index] = breakageValues[numPolymers-1];
                breakageValues[numPolymers-1] = 0;

                polymerSizes[index] = polymerSizes[numPolymers-1];
                polymerSizes[numPolymers-1] = 0;

                numPolymers--;

            }

            else if(drawProbability < totalPolyRate){
                int index = binarySearch(polymerizationRates, 0, numPolymers, drawProbability);
                numMonomers--;
                polymerSizes[index]++;
            }

            else {
                int index = binarySearch(breakageRates, 0, numPolymers, drawProbability);
                breakage(polymerizationValues, breakageValues, polymerSizes, index, polymerizationValues[index], breakageValues[index]);
            }

            mutation(polymerizationValues, breakageValues);

            if(currentTime == 0.0 || currentTime - previousTime > 0.1){

                int sum = 0;

                for(int i = 0; i < numPolymers; i++) {
                    sum += polymerSizes[i];
                }


                timeArr[counter] = currentTime;
                monomerArr[counter] = numMonomers;
                polymerArr[counter] = numPolymers;
                massArr[counter] = sum;
                sizeArr[counter] =  m.evaluate(polymerSizes, 0, numPolymers);
                meanPolyArr[counter] = m.evaluate(polymerizationValues, 0, numPolymers);
                varPolyArr[counter] = var.evaluate(polymerizationValues, 0, numPolymers);
                meanBreakageArr[counter] = m.evaluate(breakageValues, 0, numPolymers);
                varBreakageArr[counter] = var.evaluate(breakageValues, 0, numPolymers);


                previousTime = currentTime;
                counter++;
            }
        }

        outputToFile(timeArr, monomerArr, polymerArr, massArr, sizeArr, meanPolyArr, varPolyArr, meanBreakageArr, varBreakageArr, counter, simid);

    }

    // Outputs time series data from simulation

    private void outputToFile(double[] timeArr, int[] monomerArr, int[] polymerArr, int[] massArr, double[] sizeArr, double[] meanPolyArr, double[] varPolyArr, double[] meanBreakageArr, double[] varBreakageArr, int counter, int simid){

        try{

            PrintWriter pw = new PrintWriter(new File("output_" + simid + ".csv"));
            pw.println("Time,MonomerLevel,PolymerLevel,MassLevel,MeanSize,MeanPolymerization,VariancePolymerization,MeanBreakage,VarianceBreakage");
            for(int i = 0; i < counter; i++){
                pw.println(timeArr[i] + "," + monomerArr[i] + "," + polymerArr[i] + "," + massArr[i] + "," + sizeArr[i] + "," + meanPolyArr[i] + ","  + varPolyArr[i] + "," + meanBreakageArr[i] + ","  + varBreakageArr[i]);
            }
            pw.close();

        }

        catch (IOException e){
            e.printStackTrace();
        }
    }

    // binary search is used to choose event from cumulative sum array


    private int binarySearch(double[] cumDist, int low, int high, double probability){

        while(low < high){

            int mid =(low + high) / 2;
            double distance = cumDist[mid];

            if(distance < probability){
                low = mid + 1;
            }

            else if(distance > probability){
                high = mid;
            }

            else{
                return mid;
            }

        }

        return low;
    }

    // function that handles the fragmentation event

    private void breakage(double[] polymerizationValues, double[] breakageValues, double[] polymerSizes, int index, double polyValue, double breakageValue){

        int polymer_i = mtRand.nextInt((int) polymerSizes[index] - 1) + 1;
        int polymer_j = (int) polymerSizes[index] - polymer_i;

        if(polymer_i >= n && polymer_j >= n){
            polymerizationValues[index] = polyValue;
            breakageValues[index] = breakageValue;
            polymerSizes[index] = polymer_i;


            polymerizationValues[numPolymers] = polyValue;
            breakageValues[numPolymers] = breakageValue;
            polymerSizes[numPolymers] = polymer_j;

            numPolymers += 1;
        }

        else if(polymer_i < n && polymer_j >= n){
            polymerizationValues[index] = polyValue;
            breakageValues[index] = breakageValue;
            polymerSizes[index] = polymer_j;
            numMonomers += polymer_i;
        }

        else if(polymer_i >= n){
            polymerizationValues[index] = polyValue;
            breakageValues[index] = breakageValue;
            polymerSizes[index] = polymer_i;
            numMonomers += polymer_j;
        }

        else {

            polymerizationValues[index] = polymerizationValues[numPolymers-1];
            polymerizationValues[numPolymers-1] = 0;

            polymerSizes[index] = polymerSizes[numPolymers-1];
            polymerSizes[numPolymers-1] = 0;

            breakageValues[index] = breakageValues[numPolymers-1];
            breakageValues[numPolymers-1] = 0;


            numMonomers += polymer_i + polymer_j;
            numPolymers -= 1;
        }

    }

    // constant product mutation of polymerization and fragmentation rates of a polymer
    private void mutation(double[] polymerizationArr, double[] breakageArr){

        if(mtRand.nextDouble() < numPolymers * mu){
            int randInt = mtRand.nextInt(numPolymers);
            mutationPerturbation(breakageArr[randInt], polymerizationArr[randInt], randInt, polymerizationArr, breakageArr);




            if(polymerizationArr[randInt] < 0){
                polymerizationArr[randInt] = 0;
            }

            if(breakageArr[randInt] < 0){
                breakageArr[randInt] = 0;
            }
        }

    }


    // calculates perturbations of polymerization and fragmentation values of a single polymer, where b*beta = k

    private void mutationPerturbation(double b, double beta, int index, double[] polymerizationArr, double[] breakageArr) {

        double delta = (-0.0001 + (0.0001 + 0.0001) * mtRand.nextDouble());
        polymerizationArr[index] = beta + delta;
        breakageArr[index] = b + ((beta*b/polymerizationArr[index]) - b);

    }

    // calculate polymerization reaction rates for each polymer

    private double setPolyRates(int numPolymers, int numMonomers, double cumSum, double[] polymerizationValues, double[] polymerizationRates){

        for (int i = 0; i < numPolymers; i++) {
            cumSum += polymerizationValues[i] * numMonomers;
            polymerizationRates[i] = cumSum;
        }

        return cumSum;
    }
    
    // calculate size-dependent fragmentation reaction rates
    private double setBreakageRates(int numPolymers, double cumSum, double[] polymerSizes, double[] breakageRates, double[] breakageValues) {
        for (int i = 0; i < numPolymers; i++) {
            cumSum += breakageValues[i] * (polymerSizes[i] - 1);
            breakageRates[i] = cumSum;
        }

        return cumSum;
    }

}

