/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lagrangean;

import ilog.concert.IloException;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

/**
 *
 * @author AGomez
 */
public class Lagrangean {

    /**
     * @param args the command line arguments. <br>
     * 0: Size. <br>
     * 1: Beta. <br>
     * 2: Seed. <br>
     * 3: Type: From -5 to -2: diagonal dominant with density 0.02,0.1,0.5,1.
     * -1:Diagonal. 0:Full rank. From 1 to 4: Rank 1 to 4. <br>
     * 4: Abs. <br>
     * 5: Method. 0:conic - 1:QCP - 2:outer - 3:lowRankLP - 4:lowRankLP+outer -
     * 5:lagrangeanLowRank - 6:QP coordinate descent - 7:QP bisection <br>
     * @throws ilog.concert.IloException
     */
    public static void main(String[] args) throws IloException, IOException {
        int size = Integer.parseInt(args[0]), method = Integer.parseInt(args[5]), rank = Integer.parseInt(args[3]);
        boolean abs = Boolean.parseBoolean(args[4]);
        double beta = Double.parseDouble(args[1]), sol = 0;
        long seed = Long.parseLong(args[2]), timeConic, timeLagrangean, timeCoordinate, timeBissection, timeBissectionCoordinate, time = 0;
        Dag dag;
        if (rank >= 0) {
            dag = Dag.generateGridRank(size, rank, abs, new Random(seed));
        } else {
            dag = Dag.generateGrid(size, rank, abs, new Random(seed));
        }
        int iterations = 0;
//        System.out.println(dag.arcs.size());
        if (method <= 5) {

            ConicSolver solver = new ConicSolver();
            solver.init(dag, beta, method);
            time = System.currentTimeMillis();
            solver.cplex.solve();
            time = System.currentTimeMillis() - time;
            sol = solver.cplex.getObjValue();
            iterations = solver.cplex.getNiterations();
//        } else if (method == 6 || method==7) {
//            if (rank <= 0) {
//                return;
//            }
//            LagrangeanRelaxation2 lagrangean = new LagrangeanRelaxation2(dag, beta);
////            System.out.println("");
////        System.out.println("Bissection");
//            time = System.currentTimeMillis();
//            sol = lagrangean.minimizeT();
//            time = System.currentTimeMillis() - time;
        } else if (method == 6 || method==7) {
            LagrangeanRelaxation1 lagrangean = new LagrangeanRelaxation1(dag, beta);
//            System.out.println("");
//            System.out.println("CoordinateDescent:");
            time = System.currentTimeMillis();
            sol = lagrangean.sSearchCoordinate(method);
            time = System.currentTimeMillis() - time;
            iterations = lagrangean.iterations;
        } else if (method == 8 || method==9) {
            LagrangeanRelaxation1 lagrangean = new LagrangeanRelaxation1(dag, beta);
//             System.out.println("");
//        System.out.println("BissectionCoordinate");
            time = System.currentTimeMillis();
            sol = lagrangean.sBissectionCoordinate(method);
            time = System.currentTimeMillis() - time;
            iterations = lagrangean.iterations;
        }
//      
        System.out.println("Solution:");
        System.out.println(sol + " \t" + time + "\t" + iterations);
//        System.out.println("");
//        System.out.println("Golden Section Search");
//        timeLagrangean = System.currentTimeMillis();
//        double solLagrangean = lagrangean.sSearch();
//        timeLagrangean = System.currentTimeMillis() - timeLagrangean;
//        boolean exists = new File("./results/results.csv").exists();
//
//        try (FileWriter out = new FileWriter(new File("./results/results.csv"), true)) {
//            if (!exists) {
//                out.write("Size,Variables,Beta,Seed,Q,Abs,Method,Sol,Time,Iterations \n");
//            }
//            out.write(args[0] + "," + dag.arcs.size() + "," + args[1] + "," + args[2] + "," + args[3]+ "," + args[4] + "," + args[5] + "," + sol + "," + time+ "," + iterations + "\n");
//        }

    }

}
