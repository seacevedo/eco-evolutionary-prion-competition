import java.io.Serializable;

class Cell implements Serializable {

    int numPolymers;
    int numMonomers;

    int[] polymerSizes;
    double[] polymerizationValues;
    double[] breakageValues;
    double[] polymerizationRates;
    double[] breakageRates;

    double lambda;
    double d;
    double beta;
    double b;
    double a;
    int n;
    private int startingPolymerSize;


    private int MAX_SIZE;

    Cell(int numPolymers, int numMonomers, double lambda, double d, double beta, double b, double a, int n, int startingPolymerSize, int MAX_SIZE){
        this.numPolymers = numPolymers;
        this.numMonomers = numMonomers;
        this.lambda = lambda;
        this.d = d;
        this.beta = beta;
        this.b = b;
        this.a = a;
        this.n = n;
        this.MAX_SIZE = MAX_SIZE;
        this.polymerSizes = new int[this.MAX_SIZE];
        this.polymerizationValues = new double[this.MAX_SIZE];
        this.breakageValues = new double[this.MAX_SIZE];
        this.startingPolymerSize = startingPolymerSize;
        this.polymerizationRates = new double[this.MAX_SIZE];
        this.breakageRates = new double[this.MAX_SIZE];
    }



    void initPrions(){
        for(int i = 0; i < numPolymers; i++){
            polymerSizes[i] = startingPolymerSize;
            polymerizationValues[i] = beta;
            breakageValues[i] = b;

        }
    }

}
