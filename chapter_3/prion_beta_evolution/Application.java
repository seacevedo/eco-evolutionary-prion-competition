import java.util.concurrent.*;

// This code considers evolution of polymerization under dynamics of the 
// nucleated polymerization model (NPM) described by Masel, Jansen, and Nowak, 1999
// Code use to generate figure 3.1 A; can set mu = 0 for vanilla NPM

public class Application {

    public static void main(String[] args) {

        int numSim = Integer.parseInt(args[0]); // Number of simulations to run
        int processes = Integer.parseInt(args[1]); // number of simulations to run in parallel

        Simulation[] simArray = new Simulation[numSim];
        ExecutorService pool = Executors.newFixedThreadPool(processes);
        ExecutorCompletionService<Object> completionService = new ExecutorCompletionService<>(pool);

        // parameters of the NPM
        int numMonomers = 500;
        int numPolymers = 10;
        double lambda = 2000;
        double d = 4;
        double beta = 0.015;
        double b = 0.0009;
        double a = 0.05;
        int n = 6;
        double mu = 0.01;

        // Initialize simualtions

        for(int i = 0; i < numSim; i++){
            simArray[i] = new Simulation(numMonomers, numPolymers, lambda, d, beta, b, a, n, mu, i);
        }

        // submit to the pool to execute simulations

        for(final Simulation sim : simArray){
            completionService.submit(Executors.callable(new Runnable() {
                @Override
                public void run() {
                    sim.startSimulation();
                }
            }));
        }

        pool.shutdown();

    }

}

