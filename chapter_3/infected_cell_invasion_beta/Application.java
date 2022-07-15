import org.apache.commons.lang3.SerializationUtils;

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

// This code considers invasion of infected cells under dynamics of the 
// nucleated polymerization model (NPM) described by Masel, Jansen, and Nowak, 1999
// Code used to generate figure 3.2 B


public class Application {

    public static void main(String[] args) {

        int numSim = Integer.parseInt(args[0]); // Number of simulations to run
        int numProcesses = Integer.parseInt(args[1]); // number of simulations to run in parallel
        int id = Integer.parseInt(args[2]); // SLURM BATCH ID
        double invadingPolyRate = 0;
        double residentPolyRate = 0;


        // Look for files that contain polymerization rates to initilize the polymerization rates for the invading and resident prions; dependent on SLURM ID

        try {
            Path path1 = Paths.get("PATH/TO/POLYMERIZATION_RATES_INVADER/poly_rate_invader_" + id + ".csv");
            Path path2 = Paths.get("PATH/TO/POLYMERIZATION_RATES_RESIDENT/poly_rate_resident_" + id + ".csv");
            invadingPolyRate = Double.parseDouble(Files.readAllLines(path1).get(0));
            residentPolyRate = Double.parseDouble(Files.readAllLines(path2).get(0));
        }

        catch (IOException e){
            e.printStackTrace();
        }

        // parameters of the NPM

        double lambda = 2000;
        double d = 4;
        double a = 0.05;
        double b = 0.0009;
        int n = 12;
        double c = 0.05;

        int numMonomers = 500;
        int numPolymers = 100;
        int initPolymerSize = 50;

        int MAX_TIME = 5000;
        int MAX_SIZE = 10000;

        ArrayList<CompletableFuture<Void>> burnInList = new ArrayList<>();
        ArrayList<CompletableFuture<Integer>> invasionSimList = new ArrayList<>();
        ExecutorService pool = Executors.newFixedThreadPool(numProcesses);
        InvasionSimulation[] invasionSims = new InvasionSimulation[numSim];

        SynchronizedMersenneTwister rand = new SynchronizedMersenneTwister();
        // simulations where each represents a cell; each one is initilalized with a different polymerization rate (represents two different strains)
        Cell invadingCell = new Cell(numPolymers, numMonomers, lambda, d, invadingPolyRate, b, a, n, initPolymerSize, MAX_SIZE);
        Cell residentCell = new Cell(numPolymers, numMonomers, lambda, d, residentPolyRate, b, a, n, initPolymerSize, MAX_SIZE);

        invadingCell.initPrions();
        residentCell.initPrions();

        BurnInSimulation[] burnInSims = {new BurnInSimulation(invadingCell, MAX_TIME),
                new BurnInSimulation(residentCell, MAX_TIME)};
        
        burnInList.add(CompletableFuture.runAsync(() -> burnInSims[0].startSimulation(), pool));
        burnInList.add(CompletableFuture.runAsync(() -> burnInSims[1].startSimulation(), pool));


        burnInList.forEach(CompletableFuture::join);

        for(int i = 0; i < numSim; i++){


            burnInSims[0] = new BurnInSimulation(invadingCell, 1);
            burnInSims[1] = new BurnInSimulation(residentCell, 1);

            // run prion populations in each cell to equilibrium (Burn-In period)
            burnInList.add(CompletableFuture.runAsync(() -> burnInSims[0].startSimulation(), pool));
            burnInList.add(CompletableFuture.runAsync(() -> burnInSims[1].startSimulation(), pool));

            burnInList.forEach(CompletableFuture::join);        
            // Choose one of the two cells to invade, clone cell state, and inoculate with invading prion
            Cell targetCell = SerializationUtils.clone(residentCell);

            if(invadingCell.numPolymers != 0 && residentCell.numPolymers != 0) {


                int migratedPolymerIndex = rand.nextInt(invadingCell.numPolymers);

                targetCell.polymerSizes[targetCell.numPolymers] = invadingCell.polymerSizes[migratedPolymerIndex];
                targetCell.polymerizationValues[targetCell.numPolymers] = invadingCell.polymerizationValues[migratedPolymerIndex];
                targetCell.numPolymers++;

                invasionSims[i] = new InvasionSimulation(targetCell, MAX_TIME, invadingPolyRate, 10);
            }

            else {
                invasionSims[i] = new InvasionSimulation(targetCell, MAX_TIME, 0, 10);
                invasionSims[i].invasionOccurred = 0;
            }


        }

        for(InvasionSimulation sim : invasionSims){
            invasionSimList.add(CompletableFuture.supplyAsync(() -> {return sim.startSimulation();}, pool));
        }

        invasionSimList.forEach(CompletableFuture::join);

        pool.shutdown();

        // calculate invasion probability, or the fraction of sims where invader replaces the resident population in cell

        double invasionProbability = 0;
        int numSimsCompleted = 0;
        int numSuccessfulInvasions = 0;
        int numRunOutOfTime = 0;

        for(CompletableFuture<Integer> fut : invasionSimList){
            try{
                if(!fut.isCancelled()){
                    numSimsCompleted++;
                }

                if(fut.get() == 1){
                    invasionProbability++;
                    numSuccessfulInvasions++;
                }

                else if(fut.get() == 2){
                    numRunOutOfTime++;
                }
            }

            catch (InterruptedException | ExecutionException e){
                e.printStackTrace();
            }
        }

        invasionProbability = invasionProbability / numSim;


        outputToFile(invadingPolyRate, residentPolyRate, invasionProbability, numSuccessfulInvasions, numRunOutOfTime, numSimsCompleted, invadingPolyRate + "," + residentPolyRate + ".csv");


    }

    // output results to csv file

    private static void outputToFile(double invadingPoly, double hostPoly, double invasionProbability, int numInvasions, int numTimeRunOut, int numSimsCompleted, String fileName){

        try{

            PrintWriter pw = new PrintWriter(new File(fileName));
            pw.println(invadingPoly + "," + hostPoly + "," + invasionProbability + "," + numInvasions + "," + numTimeRunOut + ","  + numSimsCompleted);
            pw.close();

        }

        catch (IOException e){
            e.printStackTrace();
        }
    }

}
