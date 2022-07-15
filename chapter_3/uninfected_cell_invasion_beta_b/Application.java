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
// Code used to generate figure 3.4A


public class Application {

    public static void main(String[] args) {

        int numSim = Integer.parseInt(args[0]);
        int processes = Integer.parseInt(args[1]);
        int id = Integer.parseInt(args[2]);
        double beta = 0.0;
	    double b = 0.0;

        try {
            Path path1 = Paths.get("PATH/TO/BETA_B_RATES/beta_rate_" + id + ".csv");
	        Path path2 = Paths.get("PATH/TO/BETA_B_RATES/b_rate_" + id + ".csv");
            beta = Double.parseDouble(Files.readAllLines(path1).get(0));
	        b = Double.parseDouble(Files.readAllLines(path2).get(0));
        }

        catch (IOException e){
            e.printStackTrace();
        }

        Simulation[] simArray = new Simulation[numSim];
        ExecutorService pool = Executors.newFixedThreadPool(processes);
        ArrayList<CompletableFuture<Boolean>> invasionSimList = new ArrayList<>();


        int numMonomers = 750;
        int numPolymers = 1;
        double lambda = 3000;
        double d = 4;
        double a = 0.05;
        int n = 6;

        for(int i = 0; i < numSim; i++){
            simArray[i] = new Simulation(numMonomers, numPolymers, lambda, d, beta, b, a, n);
        }

        for(Simulation sim : simArray){
            invasionSimList.add(CompletableFuture.supplyAsync(() -> {return sim.startSimulation();}, pool));
        }


        invasionSimList.forEach(CompletableFuture::join);


        pool.shutdown();

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


        outputToFile(beta, b, invasionProbability, numInvasions, numSimsCompleted, beta + "_" + b + ".csv");

    }

    private static void outputToFile(double beta, double b, double invasionProbability, int numInvasions, int numSimsCompleted, String fileName){

        try{

            PrintWriter pw = new PrintWriter(new File(fileName));
            pw.println();
            pw.println(beta + "," + b + "," + invasionProbability+ "," + numInvasions + "," + numSimsCompleted);
            pw.close();

        }

        catch (IOException e){
            e.printStackTrace();
        }
    }

}
