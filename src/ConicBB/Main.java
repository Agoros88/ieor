/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ConicBB;

import conicSimplex.*;
import ilog.concert.IloException;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

/**
 *
 * @author AGomez
 */
public class Main {

    /**
     * Main method. Solves conic programs with matrix of the form F*Sigma*F'+D. <br>
     * @param args the command line arguments. <br>
     * 0: Size. <br>
     * 1: Beta. <br>
     * 2: Rank of Sigma. <br>
     * 3: Diagonal coefficient. <br>
     * 4: Density of F.
     * 5: Positive. <br>
     * 6: Seed. <br>
     * 7: Method:<br>
     * <\t> 11: Conic. <br>
     * <\t> 21: Outer approximation. <br>
     * <\t> 31: Quadratic. <br>
     * @throws ilog.concert.IloException
     */
    public static void main(String[] args) throws IloException, IOException {
        int size = Integer.parseInt(args[0]), method = Integer.parseInt(args[7]), rank = Integer.parseInt(args[2]);
        boolean positive = Boolean.parseBoolean(args[5]);
        double beta = Double.parseDouble(args[1]), diagonal=Double.parseDouble(args[3]),
                density=Double.parseDouble(args[4]);
        long seed = Long.parseLong(args[6]);
        DagObjective dag=DagObjective.generateGrid(size, rank, diagonal,density, positive, new Random(seed));
        
        int iterations = 0;
        long time=0;
        double sol=Double.POSITIVE_INFINITY, lowerBound=Double.NEGATIVE_INFINITY;

//        System.out.println(dag.arcs.size());
        if (method <=2) {

            ConicSolver solver = new ConicSolver();
            solver.init(dag, beta, method);
//            System.out.println("Beginning optimization");
            time = System.currentTimeMillis();
            solver.cplex.solve();
            
            time = System.currentTimeMillis() - time;
            sol = solver.cplex.getObjValue();
            iterations=solver.cplex.getNnodes();
            lowerBound=sol-solver.cplex.getMIPRelativeGap()*(sol+1e-10);
        }
        
        else if (method>=3)
        {
            BranchAndBound bb=new BranchAndBound(dag, beta, method);
            time = System.currentTimeMillis();
            bb.doBranchBound();
            time = System.currentTimeMillis() - time;
            sol=bb.upperBound;
            iterations=bb.iterations;
            lowerBound=bb.getLowerBound();
            if(method==3)
            {
                System.out.println("Time Primal: "+((QuadraticSolver)bb.solver).timePrimal
                        +"\t Time Dual: "+((QuadraticSolver)bb.solver).timeDual);
            }
        }
//      
        System.out.println("Integer Solution:");
        System.out.println(sol + " \t" +lowerBound+"\t"+ time + "\t" + iterations);

        boolean exists = new File("./results/conicSimplexInteger.csv").exists();

        try (FileWriter out = new FileWriter(new File("./results/conicSimplexInteger.csv"), true)) {
            if (!exists) {
                out.write("Size,Variables,Beta,Rank,Diagonal,Density,Positive,Seed,Method,Sol,LowerBound,Time,Iterations \n");
            }
            out.write(args[0] + "," + dag.arcs.size() +","+beta+","+rank+","+diagonal+","+density+","+positive+","+seed+","+method+","+sol+","+lowerBound+","+ time+ "," + iterations + "\n");
        }

    }

}
