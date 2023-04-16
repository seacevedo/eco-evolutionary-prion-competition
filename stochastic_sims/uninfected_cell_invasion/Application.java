import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

// This code considers invasion of uninfected cells under dynamics of the 
// nucleated polymerization model (NPM) described by Masel, Jansen, and Nowak, 1999
// Code used to generate figure 2.3 A, C, E

public class Application {

    public static void main(String[] args) {

        int numSim = Integer.parseInt(args[0]); // Number of simulations to run
        int processes = Integer.parseInt(args[1]); // number of simulations to run in parallel
        int id = Integer.parseInt(args[2]); // SLURM BATCH ID
        double b = 0.0;

        // Look for files that contains a fragmentation rate to initilize the simulations; dependent on SLURM ID

        try {
            Path path = Paths.get("/PATH/TO/FRAGMENTATION_RATES/frag_rate_" + id + ".csv");
            b = Double.parseDouble(Files.readAllLines(path).get(0));
        }

        catch (IOException e){
            e.printStackTrace();
        }

        Simulation[] simArray = new Simulation[numSim];
        ExecutorService pool = Executors.newFixedThreadPool(processes);
        ArrayList<CompletableFuture<Boolean>> invasionSimList = new ArrayList<>();

        // parameters of the NPM
        int numMonomers = 500;
        int numPolymers = 1;
        double lambda = 2000;
        double d = 4;
        double beta = 0.03;
        double a = 0.05;
        int n = 6;

        //Initialize simualtions
        for(int i = 0; i < numSim; i++){
            simArray[i] = new Simulation(numMonomers, numPolymers, lambda, d, beta, b, a, n);
        }

        // add to pool of futures, futures will wait for all simulations to finish 
        // before calculating invasion probability (number of simulations were single prion establishes stable population)
        for(Simulation sim : simArray){
            invasionSimList.add(CompletableFuture.supplyAsync(() -> {return sim.startSimulation();}, pool));
        }


        invasionSimList.forEach(CompletableFuture::join);


        pool.shutdown();


        // calculate invasion probability, or the fraction of sims where prion estbalishes a population
        double invasionProbability = 0;
        int numSimsCompleted = 0;
        int numInvasions = 0;

        for(CompletableFuture<Boolean> fut : invasionSimList){
            try{
                if(!fut.isCancelled()){
                    numSimsCompleted++;
                }

                if(fut.get()){
                    invasionProbability++;
                    numInvasions++;
                }
            }

            catch (InterruptedException | ExecutionException e){
                e.printStackTrace();
            }
        }

        invasionProbability = invasionProbability / numSim;


        outputToFile(b, invasionProbability, numInvasions, numSimsCompleted, b + ".csv");

    }

    // output results to csv file
    private static void outputToFile(double breakage, double invasionProbability, int numInvasions, int numSimsCompleted, String fileName){

        try{

            PrintWriter pw = new PrintWriter(new File(fileName));
            pw.println();
            pw.println(breakage + "," + invasionProbability+ "," + numInvasions + "," + numSimsCompleted);
            pw.close();

        }

        catch (IOException e){
            e.printStackTrace();
        }
    }

}
