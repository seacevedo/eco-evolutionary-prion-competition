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

        int[] polymerSizes = new int[MAX_SIZE];
        double[] polymerizationValues = new double[MAX_SIZE];
        double[] polymerizationRates = new double[MAX_SIZE];

        double[] breakageRates = new double[MAX_SIZE];

        double[] timeArr = new double[MAX_SIZE];
        int[] monomerArr = new int[MAX_SIZE];
        int[] polymerArr = new int[MAX_SIZE];
        int[] massArr = new int[MAX_SIZE];
        double[] sizeArr = new double[MAX_SIZE];
        double[] meanBreakageArr = new double[MAX_SIZE];
        double[] varBreakageArr = new double[MAX_SIZE];

        double drawProbability;
        int counter = 0;
        Variance var = new Variance();
        Mean m = new Mean();

        for(int i = 0; i < numPolymers; i++){
            polymerSizes[i] = initialPolymerSize;
            polymerizationValues[i] = beta;
        }

        while(currentTime < MAX_TIME && numPolymers != 0) {


            double totalPolyRate = setPolyRates(numPolymers, numMonomers,lambda+d*numMonomers+a*numPolymers, polymerizationValues, polymerizationRates);
            double totalBreakageRate = setBreakageRates(numPolymers, totalPolyRate, polymerSizes, breakageRates);

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
                breakage(polymerizationValues, polymerSizes, index, polymerizationValues[index]);
            }

            mutation(polymerizationValues);

            if(currentTime == 0.0 || currentTime - previousTime > 0.1){

                int sum = 0;

                for(int i = 0; i < numPolymers; i++) {
                    sum += polymerSizes[i];
                }


                timeArr[counter] = currentTime;
                monomerArr[counter] = numMonomers;
                polymerArr[counter] = numPolymers;
                massArr[counter] = sum;
                sizeArr[counter] = (double) sum/numPolymers;
                meanBreakageArr[counter] = m.evaluate(polymerizationValues, 0, numPolymers);
                varBreakageArr[counter] = var.evaluate(polymerizationValues, 0, numPolymers);


                previousTime = currentTime;
                counter++;
            }
        }

        outputToFile(timeArr, monomerArr, polymerArr, massArr, sizeArr, meanBreakageArr, varBreakageArr, counter, simid);

    }

    // Outputs time series data from simulation
    private void outputToFile(double[] timeArr, int[] monomerArr, int[] polymerArr, int[] massArr, double[] sizeArr, double[] meanBreakageArr, double[] varBreakageArr, int counter, int simid){

        try{

            PrintWriter pw = new PrintWriter(new File("output_" + simid + ".csv"));
            pw.println("Time,MonomerLevel,PolymerLevel,MassLevel,MeanSize,MeanPolymerization,VariancePolymerization");
            for(int i = 0; i < counter; i++){
                pw.println(timeArr[i] + "," + monomerArr[i] + "," + polymerArr[i] + "," + massArr[i] + "," + sizeArr[i] + "," + meanBreakageArr[i] + ","  + varBreakageArr[i]);
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

    private void breakage(double[] polymerizationValues, int[] polymerSizes, int index, double breakageValue){

        int polymer_i = mtRand.nextInt(polymerSizes[index] - 1) + 1;
        int polymer_j = polymerSizes[index] - polymer_i;

        if(polymer_i >= n && polymer_j >= n){
            polymerizationValues[index] = breakageValue;
            polymerSizes[index] = polymer_i;


            polymerizationValues[numPolymers] = breakageValue;
            polymerSizes[numPolymers] = polymer_j;

            numPolymers += 1;
        }

        else if(polymer_i < n && polymer_j >= n){
            polymerizationValues[index] = breakageValue;
            polymerSizes[index] = polymer_j;
            numMonomers += polymer_i;
        }

        else if(polymer_i >= n){
            polymerizationValues[index] = breakageValue;
            polymerSizes[index] = polymer_i;
            numMonomers += polymer_j;
        }

        else {

            polymerizationValues[index] = polymerizationValues[numPolymers-1];
            polymerizationValues[numPolymers-1] = 0;

            polymerSizes[index] = polymerSizes[numPolymers-1];
            polymerSizes[numPolymers-1] = 0;

            numMonomers += polymer_i + polymer_j;
            numPolymers -= 1;
        }

    }

    // mutation of polymerization rate of a polymer
    private void mutation(double[] polymerizationArr){

        if(mtRand.nextDouble() < numPolymers * mu){
            int randInt = mtRand.nextInt(numPolymers);
            polymerizationArr[randInt] += (-0.01 + (0.01 + 0.01) * mtRand.nextDouble());

            if(polymerizationArr[randInt] < 0){
                polymerizationArr[randInt] = 0;
            }
        }

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
    private double setBreakageRates(int numPolymers, double cumSum, int[] polymerSizes, double[] breakageRates) {
        for (int i = 0; i < numPolymers; i++) {
            cumSum += b * (polymerSizes[i] - 1);
            breakageRates[i] = cumSum;
        }

        return cumSum;
    }



}
