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
// Code used to generate figure 3.2A

public class Application {

    public static void main(String[] args) {

        int numSim = Integer.parseInt(args[0]);
        int processes = Integer.parseInt(args[1]);
        int id = Integer.parseInt(args[2]);
        double beta = 0.0;

        try {
            Path path = Paths.get("/PATH/TO/POLYMERIZATION_RATES/poly_rate_" + id + ".csv");
            beta = Double.parseDouble(Files.readAllLines(path).get(0));
        }

        catch (IOException e){
            e.printStackTrace();
        }

        Simulation[] simArray = new Simulation[numSim];
        ExecutorService pool = Executors.newFixedThreadPool(processes);
        ArrayList<CompletableFuture<Boolean>> invasionSimList = new ArrayList<>();


        int numMonomers = 500;
        int numPolymers = 1;
        double lambda = 2000;
        double d = 4;
        double b = 0.0009;
        double a = 0.03;
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


        outputToFile(beta, invasionProbability, numInvasions, numSimsCompleted, beta + ".csv");

    }

    private static void outputToFile(double beta, double invasionProbability, int numInvasions, int numSimsCompleted, String fileName){

        try{

            PrintWriter pw = new PrintWriter(new File(fileName));
            pw.println();
            pw.println(beta + "," + invasionProbability+ "," + numInvasions + "," + numSimsCompleted);
            pw.close();

        }

        catch (IOException e){
            e.printStackTrace();
        }
    }

}
