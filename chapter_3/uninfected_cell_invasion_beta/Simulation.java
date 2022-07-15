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

    Simulation(int numMonomers, int numPolymers, double lambda, double d, double beta, double b, double a, int n){
        this.numMonomers = numMonomers;
        this.numPolymers = numPolymers;
        this.lambda = lambda;
        this.d = d;
        this.beta = beta;
        this.b = b;
        this.a = a;
        this.n = n;
    }

    boolean startSimulation() {

        final int MAX_SIZE = 1000000;
        final double MAX_TIME = 1000;
        int initialPolymerSize = 50;
        double currentTime = 0.0;

        int[] polymerSizes = new int[MAX_SIZE];
        double[] breakageValues = new double[MAX_SIZE];
        double[] breakageRates = new double[MAX_SIZE];


        double drawProbability;

        for(int i = 0; i < numPolymers; i++){
            polymerSizes[i] = initialPolymerSize;
            breakageValues[i] = b;
        }

        while(currentTime < MAX_TIME && numPolymers != 0) {

            double totalRate = setPropensities(numPolymers, lambda+d*numMonomers+a*numPolymers+beta*numMonomers*numPolymers, breakageValues, polymerSizes, breakageRates);

                    drawProbability = totalRate * mtRand.nextDouble();
            currentTime += -Math.log(mtRand.nextDouble())/totalRate;

            if(drawProbability < lambda){
                numMonomers++;
            }

            else if(drawProbability < lambda+d*numMonomers){
                numMonomers--;
            }

            else if(drawProbability < lambda+d*numMonomers+a*numPolymers){

                int index = mtRand.nextInt(numPolymers);

                breakageValues[index] = breakageValues[numPolymers-1];
                breakageValues[numPolymers-1] = 0;

                polymerSizes[index] = polymerSizes[numPolymers-1];
                polymerSizes[numPolymers-1] = 0;

                numPolymers--;

            }

            else if(drawProbability < lambda+d*numMonomers+a*numPolymers+beta*numMonomers*numPolymers){
                int index = mtRand.nextInt(numPolymers);
                numMonomers--;
                polymerSizes[index]++;
            }

            else{
                int index = binarySearch(breakageRates, 0, numPolymers, drawProbability);
                breakage(breakageValues, polymerSizes, index, breakageValues[index]);
            }


        }

        return numPolymers > 0;

    }


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

    private void breakage(double[] breakageValues, int[] polymerSizes, int index, double breakageValue){

        int polymer_i = mtRand.nextInt(polymerSizes[index] - 1) + 1;
        int polymer_j = polymerSizes[index] - polymer_i;

        if(polymer_i >= n && polymer_j >= n){
            breakageValues[index] = breakageValue;
            polymerSizes[index] = polymer_i;

            breakageValues[numPolymers] = breakageValue;
            polymerSizes[numPolymers] = polymer_j;

            numPolymers += 1;
        }

        else if(polymer_i < n && polymer_j >= n){
            breakageValues[index] = breakageValue;
            polymerSizes[index] = polymer_j;
            numMonomers += polymer_i;
        }

        else if(polymer_i >= n){
            breakageValues[index] = breakageValue;
            polymerSizes[index] = polymer_i;
            numMonomers += polymer_j;
        }

        else {

            breakageValues[index] = breakageValues[numPolymers-1];
            breakageValues[numPolymers-1] = 0;

            polymerSizes[index] = polymerSizes[numPolymers-1];
            polymerSizes[numPolymers-1] = 0;

            numMonomers += polymer_i + polymer_j;
            numPolymers -= 1;
        }

    }


    private double setPropensities(int numPolymers, double cumSum, double[] breakageValues, int[] polymerSizes, double[] breakageRates){

        for (int i = 0; i < numPolymers; i++) {
            cumSum += breakageValues[i] * (polymerSizes[i] - 1);
            breakageRates[i] = cumSum;
        }

        return cumSum;
    }

}


