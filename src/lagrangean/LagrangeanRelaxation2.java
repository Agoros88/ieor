/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lagrangean;

import ilog.concert.IloException;
import ilog.cplex.IloCplex;

/**
 * Solves a conic problem using 2 lagrangean relaxations and linear program.
 * <br>
 *
 * @author AGomez
 */
public class LagrangeanRelaxation2 {

    //--------------------------------------------------------------------------
    // Constants
    //--------------------------------------------------------------------------
    /**
     * Tolerance for the termination criterion.
     */
    private static final double TOLERANCE = 0.001;

    /**
     * Tolerance for the termination criterion.
     */
    private static final double S_INTERVAL = 0.01;

    /**
     * Tolerance for coordinate descent.
     */
    private static final double COORDINATE_TOLERANCE = 0.001;
    //--------------------------------------------------------------------------
    // Attributes
    //--------------------------------------------------------------------------
    /**
     * Solver for the quadratic problem.
     */
    protected LPSolver solver;

    /**
     * The dag of the problem.
     */
    protected Dag dag;

    /**
     * Beta.
     */
    protected double omega;

    /**
     * Current value for s.
     */
    protected double t;

    /**
     * Current value for lambda.
     */
    protected double[] lambda;

    /**
     * Time to solve the root problem.
     */
    protected long timeRoot;

    /**
     * Indicates whether it is the first call to the solve.
     */
    protected boolean firstTime = true;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     *
     * @param dag The DAG.
     * @param omega Omega.
     * @throws ilog.concert.IloException
     */
    public LagrangeanRelaxation2(Dag dag, double omega) throws IloException {
        this.dag = dag;
        this.omega = omega;
        solver = new LPSolver();

        solver.init(dag);

    }

    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------
    /**
     * Searches for the best lambda for the i-th coordinate using golden section
     * search. <br>
     *
     * @param i The coordinate to search. <br>
     * @return The value associated with the best lambda. <br>
     * @throws ilog.concert.IloException
     */
    public double minimizeT() throws IloException {
        double tMin = 0.01,
                tMax = dag.arcs.size(),
                tLeft = 0, tRight = 0;
        t = tMin;
        lambda = new double[dag.qVectors.length];

        double valueMin = maximizeLambda();
        t = tMax;

        double valueMax = maximizeLambda();
        double valueLeft = Double.NaN, valueRight = Double.NaN;

        double phi = (Math.sqrt(5) - 1) / 2;
        double innerInterval = 0;
        while (tMax - tMin > TOLERANCE * (innerInterval)) {
            if (Double.isNaN(valueLeft)) {
                tLeft = tMax - phi * (tMax - tMin);
                t = tLeft;
                valueLeft = maximizeLambda();
            }
            if (Double.isNaN(valueRight)) {
                tRight = tMin + phi * (tMax - tMin);
                t = tRight;
                valueRight = maximizeLambda();
            }
            innerInterval = tLeft + tRight;
            if (valueLeft < valueRight) {
                tMax = tRight;
                valueMax = valueRight;
                tRight = tLeft;
                valueRight = valueLeft;
                valueLeft = Double.NaN;
            } else if (valueLeft > valueRight) {

                tMin = tLeft;
                valueMin = valueLeft;
                tLeft = tRight;
                valueLeft = valueRight;
                valueRight = Double.NaN;
            } else {
//                System.out.println("Equal");
                tMin = tLeft;
                valueMin = valueLeft;
                tMax = tRight;
                valueMax = valueRight;
                valueRight = Double.NaN;
                valueLeft = Double.NaN;
            }

        }
        double bestT = tMin, bestValue = valueMin;
        if (valueLeft > bestValue) {
            bestValue = valueLeft;
            bestT = tLeft;
        }
        if (valueRight > bestValue) {
            bestValue = valueRight;
            bestT = tRight;
        }
        if (valueMax > bestValue) {
            bestValue = valueMax;
            bestT = tMax;
        }
        t = bestT;
        return bestValue;
    }

    /**
     * Searches for the best lambda. <br>
     *
     * @return The value associated with the best lambda.
     */
    private double maximizeLambda() throws IloException {
        double previousValue = Double.NEGATIVE_INFINITY, value = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < lambda.length; i++) {
            value = maximizeLambdaI(i);
        }

        while (value - previousValue > TOLERANCE) {
            previousValue = value;
            for (int i = 0; i < lambda.length; i++) {
                value = maximizeLambdaI(i);

            }
//            System.out.println(t+"\t"+previousValue+"\t"+value);
        }
        return value;

    }

    /**
     * Searches for the best lambda for the i-th coordinate using golden section
     * search. <br>
     *
     * @param i The coordinate to search. <br>
     * @return The value associated with the best lambda. <br>
     * @throws ilog.concert.IloException
     */
    private double maximizeLambdaI(int i) throws IloException {
        double lambdaMin = -2 * (Math.sqrt(dag.vertices) - 1) * omega / t,//Negative of the length of a path  * omega /t.
                lambdaMax = 2 * (Math.sqrt(dag.vertices) - 1) * omega / t, //Length of a path *omega / t.
                lambdaLeft = 0, lambdaRight = 0;
        lambda[i] = lambdaMin;

        double valueMin = minimizeX();

        lambda[i] = lambdaMax;
        double valueMax = minimizeX();

        double valueLeft = Double.NaN, valueRight = Double.NaN;

        double phi = (Math.sqrt(5) - 1) / 2;
        double innerInterval = 0;

        while (lambdaMax - lambdaMin > TOLERANCE * innerInterval && (lambdaMax - lambdaMin) > TOLERANCE) {
            if (Double.isNaN(valueLeft)) {
                lambdaLeft = lambdaMax - phi * (lambdaMax - lambdaMin);
                lambda[i] = lambdaLeft;
                valueLeft = minimizeX();
            }
            if (Double.isNaN(valueRight)) {
                lambdaRight = lambdaMin + phi * (lambdaMax - lambdaMin);
                lambda[i] = lambdaRight;
                valueRight = minimizeX();
            }
            innerInterval = lambdaLeft + lambdaRight;
//            System.out.println(t+" "+i+"\t"+valueLeft+"\t"+valueRight);

            if (valueLeft > valueRight) {
                lambdaMax = lambdaRight;
                valueMax = valueRight;
                lambdaRight = lambdaLeft;
                valueRight = valueLeft;
                valueLeft = Double.NaN;
            } else if (valueLeft < valueRight) {

                lambdaMin = lambdaLeft;
                valueMin = valueLeft;
                lambdaLeft = lambdaRight;
                valueLeft = valueRight;
                valueRight = Double.NaN;
            } else {
//                System.out.println("Equal");
                lambdaMin = lambdaLeft;
                valueMin = valueLeft;
                lambdaMax = lambdaRight;
                valueMax = valueRight;
                valueRight = Double.NaN;
                valueLeft = Double.NaN;
            }
//             System.out.println("\t"+lambdaMin+"\t"+lambdaLeft+"\t"+lambdaRight+"\t"+lambdaMax);
//            System.out.println("\t \t"+(lambdaMax-lambdaMin)+"\t"+TOLERANCE+"\t"+((lambdaMax-lambdaMin)<TOLERANCE));

        }
        double bestLambda = lambdaMin, bestValue = valueMin;
        if (valueLeft > bestValue) {
            bestValue = valueLeft;
            bestLambda = lambdaLeft;
        }
        if (valueRight > bestValue) {
            bestValue = valueRight;
            bestLambda = lambdaRight;
        }
        if (valueMax > bestValue) {
            bestValue = valueMax;
            bestLambda = lambdaMax;
        }
        lambda[i] = bestLambda;
        return bestValue;
    }

    /**
     * Finds the best solution for a given lambda. <br>
     *
     * @param lambda Lambda. <br>
     * @return The objective value.
     */
    private double minimizeX() throws IloException {
        if (firstTime) {
            solver.changeObjective(lambda);
//            firstTime = false;
        } else {
            IloCplex.BasisStatus[] baseVariables = solver.cplex.getBasisStatuses(solver.x), baseConstraints = solver.cplex.getBasisStatuses(solver.constraints);
            solver.changeObjective(lambda);
            solver.cplex.setBasisStatuses(solver.x, baseVariables, solver.constraints, baseConstraints);
        }
        solver.cplex.solve();
//        System.out.println(lambda+"\t"+(linearValue+lambda*linearVariance)+"\t"+solver.cplex.getObjValue());
        double resp = solver.cplex.getObjValue() + omega * t / 2;

        for (int i = 0; i < lambda.length; i++) {
            resp -= 0.5 * t * lambda[i] * lambda[i] / omega;
        }
        return resp;
    }

}
