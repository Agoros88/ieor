/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package conicSimplex;

import ilog.concert.IloException;
import ilog.cplex.IloCplex.BasisStatus;

/**
 * Solves a conic problem using lagrangean relaxation. <br>
 *
 * @author AGomez
 */
public class QuadatricReformulation {

    //--------------------------------------------------------------------------
    // Constants
    //--------------------------------------------------------------------------

    /**
     * Tolerance for the termination criterion.
     */
    private static final double TOLERANCE = 0.01;

    //--------------------------------------------------------------------------
    // Attributes
    //--------------------------------------------------------------------------
    /**
     * Solver for the quadratic problem.
     */
    protected QuadraticSolver solver;

    /**
     * The dag of the problem.
     */
    protected DagObjective dag;

    /**
     * Beta.
     */
    protected double beta;

    /**
     * Current value for t.
     */
    protected double t;

    /**
     * Number of simplex iterations.
     */
    protected int iterations;


    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     *
     * @param dag The DAG.
     * @param beta Beta.
     * @throws ilog.concert.IloException
     */
    public QuadatricReformulation(DagObjective dag, double beta) throws IloException {
        this.dag = dag;
        this.beta = beta;
        solver = new QuadraticSolver();
        iterations = 0;
    }

    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------
   

    /**
     * uses coordinate descent on the quadratic reformulation to solve the conic program. <br>
     *
     * @param method The method to use. <br>
     * @return The optimal S. <br>
     * @throws ilog.concert.IloException
     */
    public double coordinateDescent(int method) throws IloException {
        solver.init(dag, 0, method);

//        solver.useCallback();
        solver.cplex.solve();
        iterations += solver.cplex.getNiterations();
        
        t = Math.sqrt(solver.cplex.getValue(solver.quadraticPart));
//        System.out.println(iterations+"\t"+Double.POSITIVE_INFINITY+"\t"+(beta*t+solver.cplex.getValue(solver.linearPart)));
        double newT,value;
        boolean terminate = false;
        while (!terminate) {
            value = minimizeX(beta / (2 * t));
            newT = Math.sqrt(solver.cplex.getValue(solver.quadraticPart));
            if (Math.abs(newT - t) < TOLERANCE) {
                terminate = true;
            }
//            System.out.println(iterations+"\t"+t+"\t"+(beta*newT+solver.cplex.getValue(solver.linearPart)));
            t = newT;
            
        }
        return beta * t + solver.cplex.getValue(solver.linearPart);
    }

    

   

    /**
     * Finds the best solution for a given lambda. <br>
     *
     * @param lambda Lambda. <br>
     * @return The objective value.
     */
    private double minimizeX(double lambda) throws IloException {

        BasisStatus[] baseVariables = solver.cplex.getBasisStatuses(solver.x), baseConstraints = solver.cplex.getBasisStatuses(solver.constraints);
        solver.changeObjective(lambda);
        solver.cplex.setBasisStatuses(solver.x, baseVariables, solver.constraints, baseConstraints);

        solver.cplex.solve();

        iterations += solver.cplex.getNiterations();

        return solver.cplex.getObjValue() + beta * t - lambda * t * t;
    }

}
