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
// Code used to generate heatmaps in figure 3.4, 3.5, and 3.6


public class Application {

    public static void main(String[] args) {

        int numSim = Integer.parseInt(args[0]);
        int numProcesses = Integer.parseInt(args[1]);
        int id = Integer.parseInt(args[2]);
        double invadingPolyRate = 0;
        double residentPolyRate = 0;
        double invadingBreakageRate = 0;
        double residentBreakageRate = 0;


        try {
            Path path1 = Paths.get("/home/seaceved/invasion_heatmap_beta_b/invader_beta/invader_" + id + ".csv");
            Path path2 = Paths.get("/home/seaceved/invasion_heatmap_beta_b/resident_beta/resident_" + id + ".csv");
            Path path3 = Paths.get("/home/seaceved/invasion_heatmap_beta_b/invader_b/invader_" + id + ".csv");
            Path path4 = Paths.get("/home/seaceved/invasion_heatmap_beta_b/resident_b/resident_" + id + ".csv");
            invadingPolyRate = Double.parseDouble(Files.readAllLines(path1).get(0));
            residentPolyRate = Double.parseDouble(Files.readAllLines(path2).get(0));
            invadingBreakageRate = Double.parseDouble(Files.readAllLines(path3).get(0));
            residentBreakageRate = Double.parseDouble(Files.readAllLines(path4).get(0));
        }

        catch (IOException e){
            e.printStackTrace();
        }

        double lambda = 2000;
        double d = 2;
        double a = 0.05;
        int n = 6;

        int numMonomers = 1000;
        int numPolymers = 10;
        int initPolymerSize = 50;

        int MAX_TIME = 8000;
        int MAX_SIZE = 100000;

        ArrayList<CompletableFuture<Void>> burnInList = new ArrayList<>();
        ArrayList<CompletableFuture<Integer>> invasionSimList = new ArrayList<>();
        ExecutorService pool = Executors.newFixedThreadPool(numProcesses);
        InvasionSimulation[] invasionSims = new InvasionSimulation[numSim];

        SynchronizedMersenneTwister rand = new SynchronizedMersenneTwister();

        Cell invadingCell = new Cell(numPolymers, numMonomers, lambda, d, invadingPolyRate, invadingBreakageRate, a, n, initPolymerSize, MAX_SIZE);
        Cell residentCell = new Cell(numPolymers, numMonomers, lambda, d, residentPolyRate, residentBreakageRate, a, n, initPolymerSize, MAX_SIZE);

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

            burnInList.add(CompletableFuture.runAsync(() -> burnInSims[0].startSimulation(), pool));
            burnInList.add(CompletableFuture.runAsync(() -> burnInSims[1].startSimulation(), pool));

            burnInList.forEach(CompletableFuture::join);

            Cell targetCell = SerializationUtils.clone(residentCell);

            if(invadingCell.numPolymers != 0 && residentCell.numPolymers != 0) {


                int migratedPolymerIndex = rand.nextInt(invadingCell.numPolymers);

                targetCell.polymerSizes[targetCell.numPolymers] = invadingCell.polymerSizes[migratedPolymerIndex];
                targetCell.polymerizationValues[targetCell.numPolymers] = invadingCell.polymerizationValues[migratedPolymerIndex];
                targetCell.breakageValues[targetCell.numPolymers] = invadingCell.breakageValues[migratedPolymerIndex];
                targetCell.numPolymers++;

                invasionSims[i] = new InvasionSimulation(targetCell, MAX_TIME, invadingPolyRate, invadingBreakageRate);
            }

            else {
                invasionSims[i] = new InvasionSimulation(targetCell, MAX_TIME, 0, 0);
                invasionSims[i].invasionOccurred = 0;
            }


        }

        for(InvasionSimulation sim : invasionSims){
            invasionSimList.add(CompletableFuture.supplyAsync(() -> {return sim.startSimulation();}, pool));
        }

        invasionSimList.forEach(CompletableFuture::join);

        pool.shutdown();

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


        outputToFile(invadingPolyRate, residentPolyRate, invadingBreakageRate, residentBreakageRate, invasionProbability, numSuccessfulInvasions, numRunOutOfTime, numSimsCompleted, invadingPolyRate + "," + residentPolyRate + "," + invadingBreakageRate + "," + residentBreakageRate + ".csv");


    }


    private static void outputToFile(double invadingPoly, double hostPoly, double invadingBreakage, double hostBreakage, double invasionProbability, int numInvasions, int numTimeRunOut, int numSimsCompleted, String fileName){

        try{

            PrintWriter pw = new PrintWriter(new File(fileName));
            pw.println(invadingPoly + "," + hostPoly + "," + invadingBreakage + "," + hostBreakage + "," + invasionProbability + "," + numInvasions + "," + numTimeRunOut + ","  + numSimsCompleted);
            pw.close();

        }

        catch (IOException e){
            e.printStackTrace();
        }
    }

}
