class BurnInSimulation {

    Cell cell;
    private int MAX_TIME;
    private SynchronizedMersenneTwister rand = new SynchronizedMersenneTwister();

    BurnInSimulation(Cell cell, int MAX_TIME) {
        this.cell = cell;
        this.MAX_TIME = MAX_TIME;
    }

    // This function runs the Gillespie's SSA
    void startSimulation() {

        double currentTime = 0.0;
        double drawProbability;



        while (currentTime < MAX_TIME && cell.numPolymers != 0) {

            double totalPolyRate = setPolyRates(cell.numPolymers, cell.numMonomers,cell.lambda+cell.d*cell.numMonomers+cell.a*cell.numPolymers, cell.polymerizationValues, cell.polymerizationRates);
            double totalRate = setBreakageRates(totalPolyRate);

            drawProbability = totalRate * rand.nextDouble();
            currentTime += -Math.log(rand.nextDouble()) / totalRate;

            chooseEvent(drawProbability, totalPolyRate);

        }



    }

    // Event choosing algorithm
    private void chooseEvent(double probability, double totalPolyRate){


        if(probability < cell.lambda){
            cell.numMonomers++;
        }

        else if(probability < cell.lambda+cell.d*cell.numMonomers){
            cell.numMonomers--;
        }

        else if(probability < cell.lambda+cell.d*cell.numMonomers+cell.a*cell.numPolymers){

            int index = rand.nextInt(cell.numPolymers);

            cell.polymerizationValues[index] = cell.polymerizationValues[cell.numPolymers-1];
            cell.polymerizationValues[cell.numPolymers-1] = 0;

            cell.polymerSizes[index] = cell.polymerSizes[cell.numPolymers-1];
            cell.polymerSizes[cell.numPolymers-1] = 0;


            cell.numPolymers--;
        }

        else if(probability < totalPolyRate){
            int index = rand.nextInt(cell.numPolymers);
            cell.numMonomers--;
            cell.polymerSizes[index]++;
        }


        else{
            int index = binarySearch(cell.breakageRates, 0, cell.numPolymers, probability);
            breakage(index, cell.polymerizationValues[index], cell);
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

    private void breakage(int index, double breakageValue, Cell cell){

        int polymer_i = rand.nextInt((cell.polymerSizes[index] - 1) + 1) + 1;
        int polymer_j = cell.polymerSizes[index] - polymer_i;

        if(polymer_i >= cell.n && polymer_j >= cell.n){
            cell.polymerizationValues[index] = breakageValue;
            cell.polymerSizes[index] = polymer_i;

            cell.polymerizationValues[cell.numPolymers] = breakageValue;
            cell.polymerSizes[cell.numPolymers] = polymer_j;

            cell.numPolymers += 1;
        }

        else if(polymer_i < cell.n && polymer_j >= cell.n){
            cell.polymerizationValues[index] = breakageValue;
            cell.polymerSizes[index] = polymer_j;
            cell.numMonomers += polymer_i;
        }

        else if(polymer_i >= cell.n){
            cell.polymerizationValues[index] = breakageValue;
            cell.polymerSizes[index] = polymer_i;
            cell.numMonomers += polymer_j;
        }

        else {

            cell.polymerizationValues[index] = cell.polymerizationValues[cell.numPolymers-1];
            cell.polymerizationValues[cell.numPolymers-1] = 0;

            cell.polymerSizes[index] = cell.polymerSizes[cell.numPolymers-1];
            cell.polymerSizes[cell.numPolymers-1] = 0;

            cell.numMonomers += polymer_i + polymer_j;
            cell.numPolymers -= 1;
        }

    }

    // calculate size-dependent fragmentation reaction rates

    private double setBreakageRates(double cumSum){


        for(int i = 0; i < cell.numPolymers; i++){
            cumSum += cell.b * (cell.polymerSizes[i] - 1);
            cell.breakageRates[i] = cumSum;
        }

        return cumSum;
    }

    // calculate polymerization reaction rates for each polymer

    private double setPolyRates(int numPolymers, int numMonomers, double cumSum, double[] polymerizationValues, double[] polymerizationRates){

        for (int i = 0; i < numPolymers; i++) {
            cumSum += polymerizationValues[i] * numMonomers;
            polymerizationRates[i] = cumSum;
        }

        return cumSum;
    }


}
