class InvasionSimulation {

    Cell cell;
    private int MAX_TIME;
    private double invadingPolyRate;
    private double invadingBreakageRate;
    private SynchronizedMersenneTwister rand = new SynchronizedMersenneTwister();
    int invasionOccurred = 2;

    InvasionSimulation(Cell cell, int MAX_TIME, double invadingPolyRate, double invadingBreakageRate) {
        this.cell = cell;
        this.MAX_TIME = MAX_TIME;
        this.invadingPolyRate = invadingPolyRate;
        this.invadingBreakageRate = invadingBreakageRate;
    }


    int startSimulation() {

        double currentTime = 0.0;
        double drawProbability;



        while(currentTime < MAX_TIME && cell.numPolymers != 0) {


            double totalPolyRate = setPolyRates(cell.numMonomers,cell.lambda+cell.d*cell.numMonomers+cell.a*cell.numPolymers);
            double totalRate = setBreakageRates(totalPolyRate);

            drawProbability = totalRate * rand.nextDouble();
            currentTime += -Math.log(rand.nextDouble())/totalRate;

            chooseEvent(drawProbability, totalPolyRate);

            int invadingStrainCounter = 0;
            int residentStrainCounter = 0;

            for(int i = 0; i < cell.numPolymers; i++){
                if(cell.polymerizationValues[i] == cell.beta && cell.breakageValues[i] == cell.b){
                    residentStrainCounter++;
                }

                else if(cell.polymerizationValues[i] == invadingPolyRate && cell.breakageValues[i] == invadingBreakageRate){
                    invadingStrainCounter++;
                }
            }

            if(invadingStrainCounter < 1){
                invasionOccurred = 0;
                break;
            }

            else if(invadingStrainCounter > residentStrainCounter + residentStrainCounter * 0.3){
                invasionOccurred = 1;
                break;
            }





        }


        return invasionOccurred;


    }



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

            cell.breakageValues[index] = cell.breakageValues[cell.numPolymers-1];
            cell.breakageValues[cell.numPolymers-1] = 0;

            cell.polymerSizes[index] = cell.polymerSizes[cell.numPolymers-1];
            cell.polymerSizes[cell.numPolymers-1] = 0;


            cell.numPolymers--;
        }

        else if(probability < totalPolyRate){
            int index = binarySearch(cell.polymerizationRates, 0, cell.numPolymers, probability);
            cell.numMonomers--;
            cell.polymerSizes[index]++;
        }

        else{
            int index = binarySearch(cell.breakageRates, 0, cell.numPolymers, probability);
            breakage(index, cell.polymerizationValues[index], cell.breakageValues[index], cell);
        }


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

    private void breakage(int index, double polyValue,  double breakageValue, Cell cell) {

        int polymer_i = 0;
        int polymer_j = 0;

        try {
            polymer_i = rand.nextInt((cell.polymerSizes[index] - 1) + 1) + 1;
            polymer_j = cell.polymerSizes[index] - polymer_i;
        } catch (IllegalArgumentException e) {
            e.printStackTrace();
        }

        if(polymer_i >= cell.n && polymer_j >= cell.n){
            cell.polymerizationValues[index] = polyValue;
            cell.breakageValues[index] = breakageValue;
            cell.polymerSizes[index] = polymer_i;

            cell.polymerizationValues[cell.numPolymers] = polyValue;
            cell.breakageValues[cell.numPolymers] = breakageValue;
            cell.polymerSizes[cell.numPolymers] = polymer_j;

            cell.numPolymers += 1;
        }

        else if(polymer_i < cell.n && polymer_j >= cell.n){
            cell.polymerizationValues[index] = polyValue;
            cell.breakageValues[index] = breakageValue;
            cell.polymerSizes[index] = polymer_j;
            cell.numMonomers += polymer_i;
        }

        else if(polymer_i >= cell.n){
            cell.polymerizationValues[index] = polyValue;
            cell.breakageValues[index] = breakageValue;
            cell.polymerSizes[index] = polymer_i;
            cell.numMonomers += polymer_j;
        }

        else {

            cell.polymerizationValues[index] = cell.polymerizationValues[cell.numPolymers-1];
            cell.polymerizationValues[cell.numPolymers-1] = 0;

            cell.breakageValues[index] = cell.breakageValues[cell.numPolymers-1];
            cell.breakageValues[cell.numPolymers-1] = 0;

            cell.polymerSizes[index] = cell.polymerSizes[cell.numPolymers-1];
            cell.polymerSizes[cell.numPolymers-1] = 0;

            cell.numMonomers += polymer_i + polymer_j;
            cell.numPolymers -= 1;
        }
    }

    private double setBreakageRates(double cumSum){


        for(int i = 0; i < cell.numPolymers; i++){
            cumSum += cell.breakageValues[i] * (cell.polymerSizes[i] - 1);
            cell.breakageRates[i] = cumSum;
        }

        return cumSum;
    }

    private double setPolyRates(int numMonomers, double cumSum){

        for (int i = 0; i < cell.numPolymers; i++) {
            cumSum += cell.polymerizationValues[i] * numMonomers;
            cell.polymerizationRates[i] = cumSum;
        }

        return cumSum;
    }



}
