/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lagrangean;

import ilog.concert.IloException;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.BasisStatus;

/**
 * Solves a conic problem using lagrangean relaxation. <br>
 *
 * @author AGomez
 */
public class LagrangeanRelaxation1 {

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
    private static final double COORDINATE_TOLERANCE = 0.1;
    //--------------------------------------------------------------------------
    // Attributes
    //--------------------------------------------------------------------------
    /**
     * Solver for the quadratic problem.
     */
    protected QSolver solver;

    /**
     * The dag of the problem.
     */
    protected Dag dag;

    /**
     * Beta.
     */
    protected double beta;

    /**
     * Value and variance of the shortest path.
     */
    protected double linearValue, linearVariance;

    /**
     * Bounds
     */
    protected double minVariance = 0, minMean = 0, maxMean = Double.POSITIVE_INFINITY, maxVariance = Double.POSITIVE_INFINITY;

    /**
     * Current value for s.
     */
    protected double s;

    /**
     * Time to solve the root problem.
     */
    protected long timeRoot;

    /**
     * Whether to end the simplex loop.
     */
    protected boolean end;

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
    public LagrangeanRelaxation1(Dag dag, double beta) throws IloException {
        this.dag = dag;
        this.beta = beta;
        solver = new QSolver();
        iterations = 0;
//        solver.init(dag, 0);
//        
//        solver.cplex.solve();
//        linearValue=solver.cplex.getObjValue();
//        minMean=linearValue;
//        double[] xl=solver.cplex.getValues(solver.x);
//        linearVariance=0;
//        for (int i = 0; i < xl.length; i++) {
//            for (int j = 0; j < xl.length; j++) {
//                linearVariance+=xl[i]*xl[j]*dag.covariances[i][j];
//            }
//        }
//        maxVariance=linearVariance;
//        minimizeX(Double.POSITIVE_INFINITY);
//        minVariance=solver.cplex.getValue(solver.quadraticPart);
//        maxMean=solver.cplex.getValue(solver.linearPart);

//        System.out.println(linearVariance+" "+solver.cplex.getValue(solver.quadraticPart));
//        
//        System.out.println("");
//        
//        System.out.println("");
    }

    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------
    /**
     * Returns the solution of lambda=beta/2s. <br>
     *
     * @return the corresponding solution.
     */
    private double midLambda() throws IloException {
        return minimizeX(beta / (2 * s));
    }

    /**
     * Searches for the best lambda.
     *
     * @param s The current value for the variance.
     * @throws ilog.concert.IloException
     */
    private double lambdaSearch() throws IloException {
        double lambdaMin = 0, lambdaMax = beta / s, lambdaLeft = 0, lambdaRight = 0;

        double valueMin = linearValue + beta * s, valueMax = minimizeX(lambdaMax), valueLeft = 0, valueRight = 0;

        double initialInterval = lambdaMax - lambdaMin, interval = lambdaMax - lambdaMin;
        double phi = (Math.sqrt(5) - 1) / 2;
        double bestInteger = 0, bestContinuous = 0, temp;

//        bestContinuous=valueMin[0];
        double gap = 1;
        double innerInterval = 0;
//        temp=valueMax[0];
//        bestContinuous=temp>bestContinuous?temp:bestContinuous;
        while (lambdaMax - lambdaMin > TOLERANCE * (innerInterval)) {
            if (valueLeft == 0) {
                lambdaLeft = lambdaMax - phi * (lambdaMax - lambdaMin);
                valueLeft = minimizeX(lambdaLeft);
            }
            if (valueRight == 0) {
                lambdaRight = lambdaMin + phi * (lambdaMax - lambdaMin);
                valueRight = minimizeX(lambdaRight);
            }
            innerInterval = lambdaLeft + lambdaRight;
            //System.out.println(valueLeft+"\t"+valueRight);
            if (valueLeft > valueRight) {
                lambdaMax = lambdaRight;
                valueMax = valueRight;
                lambdaRight = lambdaLeft;
                valueRight = valueLeft;
                valueLeft = 0;
            } else if (valueLeft < valueRight) {

                lambdaMin = lambdaLeft;
                valueMin = valueLeft;
                lambdaLeft = lambdaRight;
                valueLeft = valueRight;
                valueRight = 0;
            } else {
//                System.out.println("Equal");
                lambdaMin = lambdaLeft;
                valueMin = valueLeft;
                lambdaMax = lambdaRight;
                valueMax = valueRight;
                valueRight = 0;
                valueLeft = 0;
            }

            interval = lambdaMax - lambdaMin;
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
        //System.out.println(bestLambda+" in ["+lambdaMin+","+lambdaMax+"]\t"+beta/s);
        return bestValue;
    }

    /**
     * Searches for the optimal S. <br>
     *
     * @return The optimal S.
     */
    public double sSearch() throws IloException {
        minimizeX(Double.POSITIVE_INFINITY);
        minVariance = solver.cplex.getValue(solver.quadraticPart);
        maxMean = solver.cplex.getValue(solver.linearPart);

        double sMin = Math.sqrt(minVariance), sMax = Math.sqrt(linearVariance), sLeft = 0, sRight = 0;
        double phi = (Math.sqrt(5) - 1) / 2;
        s = sMax;
        double valueMin = beta * sMin + maxMean, valueMax = beta * sMax + minMean, valueLeft = 0, valueRight = 0;
        double innerInterval = 0;
        double tempM = 0, tempV = 0;
        double lowerBound = 0, upperBound = Double.POSITIVE_INFINITY;
        //while(sMax-sMin>TOLERANCE*(innerInterval) )
        while ((upperBound - lowerBound) / lowerBound > TOLERANCE) {
            if (valueLeft == 0) {

                sLeft = sMax - phi * (sMax - sMin);
                s = sLeft;
                valueLeft = midLambda();
                tempV = solver.cplex.getValue(solver.quadraticPart);
            }
            if (valueRight == 0) {
                sRight = sMin + phi * (sMax - sMin);
                s = sRight;
                valueRight = midLambda();
                tempM = solver.cplex.getValue(solver.linearPart);
            }
            innerInterval = sLeft + sRight;
            //System.out.println(valueLeft+"\t"+valueRight);
            if (valueLeft > valueRight) {
                minVariance = tempV;
                upperBound = valueRight;
                lowerBound = minMean + beta * Math.sqrt(minVariance);
                //System.out.println(upperBound+"\t"+lowerBound);
                sMin = sLeft;
                valueMin = valueLeft;
                sLeft = sRight;
                valueLeft = valueRight;
                valueRight = 0;

            } else if (valueLeft < valueRight) {
                minMean = tempM;
                upperBound = valueLeft;
                lowerBound = minMean + beta * Math.sqrt(minVariance);
//                System.out.println(upperBound+"\t"+lowerBound);
                sMax = sRight;
                valueMax = valueRight;
                sRight = sLeft;
                valueRight = valueLeft;
                valueLeft = 0;

            } else {
                minMean = tempM;
                minVariance = tempV;
                //System.out.println("Equal");

                upperBound = valueRight;
                lowerBound = minMean + beta * Math.sqrt(minVariance);
//                System.out.println(upperBound+"\t"+lowerBound);
                sMin = sLeft;
                valueMin = valueLeft;
                sMax = sRight;
                valueMax = valueRight;
                valueRight = 0;
                valueLeft = 0;
                minMean = tempM;
                minVariance = tempV;
            }
        }
        double bestS = sMin, bestValue = valueMin;
        if (valueLeft < bestValue && valueLeft > 0) {
            bestValue = valueLeft;
            bestS = sLeft;
        }
        if (valueRight < bestValue && valueRight > 0) {
            bestValue = valueRight;
            bestS = sRight;
        }
        if (valueMax < bestValue) {
            bestValue = valueMax;
            bestS = sMax;
        }
//        System.out.println(bestS+" in ["+sMin+","+sMax+"]\t Original ["+Math.sqrt(this.minVariance)+","+Math.sqrt(linearVariance)+"]");
        return bestValue;
    }

    /**
     * Searches for the optimal S. <br>
     *
     * @param method The method to use (even:normal; odd:lowRank). <br>
     * @return The optimal S.
     */
    public double sSearchCoordinate(int method) throws IloException {
        solver.init(dag, 0, method);

//        solver.useCallback();
        solver.cplex.solve();
        iterations += solver.cplex.getNiterations();
//        s=0;
//        minimizeX(0);
//        for (double xi : solver.cplex.getValues(solver.x)) {
//            System.out.print(xi+" ");
//        }
//        System.out.println("");
//      System.out.println(iterations+"\t"+Double.POSITIVE_INFINITY+"\t"+Double.POSITIVE_INFINITY+"\t"+(Double.POSITIVE_INFINITY)+"\t"+(minMean+beta*Math.sqrt(minVariance)));
        s = Math.sqrt(solver.cplex.getValue(solver.quadraticPart));
        double newS;
//       solver.useCallback();
//        System.out.println(s);
        end = false;
        double value;
        boolean terminate = false;
        while (!terminate) {
            value = midLambda();
            newS = Math.sqrt(solver.cplex.getValue(solver.quadraticPart));
//            if(newS<s)
//            {
//                minMean=solver.cplex.getValue(solver.linearPart);
//            }
//            else if(newS>s)
//            {
//                minVariance=solver.cplex.getValue(solver.quadraticPart);
//            }
            if (Math.abs(newS - s) < COORDINATE_TOLERANCE) {
                terminate = true;
            }
//            if(end)terminate=true;
//            System.out.println(iterations+"\t"+s+"\t"+newS+"\t"+(beta*newS+solver.cplex.getValue(solver.linearPart))+"\t"+(minMean+beta*Math.sqrt(minVariance)));
            s = newS;
        }
        return beta * s + solver.cplex.getValue(solver.linearPart);
    }

    /**
     * Searches for the optimal S using binary search. <br>
     *
     * @return The optimal S.
     */
    public double sBissection() throws IloException {
        timeRoot = System.currentTimeMillis();
        minimizeX(Double.POSITIVE_INFINITY);
        timeRoot = System.currentTimeMillis() - timeRoot;
        minVariance = solver.cplex.getValue(solver.quadraticPart);
        maxMean = solver.cplex.getValue(solver.linearPart);
        double minS = Math.sqrt(Math.max(minVariance, 0)), maxS = Math.sqrt(maxVariance);

        boolean terminate = false;
        double upperBound = Double.POSITIVE_INFINITY, lowerBound = minMean + beta * minS;
        double newS;
        while ((upperBound - lowerBound) / lowerBound > TOLERANCE) {
            s = (minS + maxS) / 2;

            midLambda();
            newS = Math.sqrt(solver.cplex.getValue(solver.quadraticPart));
            if (newS < s) {
                maxS = s;
                minMean = solver.cplex.getValue(solver.linearPart);
            } else if (newS > s) {
                minS = s;
                minVariance = solver.cplex.getValue(solver.quadraticPart);
            }
            upperBound = beta * newS + solver.cplex.getValue(solver.linearPart);
            lowerBound = minMean + beta * Math.sqrt(minVariance);

//            System.out.println(s + "\t" + newS + "\t" + "[" + minS + "," + maxS + "]\t" + upperBound + "\t" + lowerBound + "\t" + (upperBound - lowerBound) / lowerBound);
            //s=newS;
        }
        return upperBound;
    }

    /**
     * Searches for the optimal S using binary search and coordinate descent.
     * <br>
     *
     * @param method The method to use (even:normal; odd:lowRank). <br>
     * @return The optimal S.
     */
    public double sBissectionCoordinate(int method) throws IloException {
        solver.init(dag, 0, method);

        solver.cplex.solve();
        iterations += solver.cplex.getNiterations();
        linearValue = solver.cplex.getObjValue();
        minMean = linearValue;
        double[] xl = solver.cplex.getValues(solver.x);
        linearVariance = 0;
        for (int i = 0; i < xl.length; i++) {
            for (int j = 0; j < xl.length; j++) {
                linearVariance += xl[i] * xl[j] * dag.covariances[i][j];
            }
        }
        maxVariance = linearVariance;

        minimizeX(Double.POSITIVE_INFINITY);
//        timeRoot=System.currentTimeMillis()-timeRoot;
        minVariance = solver.cplex.getValue(solver.quadraticPart);
        maxMean = solver.cplex.getValue(solver.linearPart);
        double minS = Math.sqrt(Math.max(minVariance, 0)), maxS = Math.sqrt(maxVariance);
//        System.out.println(minS +" "+maxS);

        double upperBound = Double.POSITIVE_INFINITY, lowerBound = minMean + beta * minS;
        double newS;
        while ((upperBound - lowerBound) / lowerBound > TOLERANCE) {
            s = (minS + maxS) / 2;

            midLambda();
            newS = Math.sqrt(solver.cplex.getValue(solver.quadraticPart));
            if (newS < s) {
                maxS = newS;
                minMean = solver.cplex.getValue(solver.linearPart);
            } else if (newS > s) {
                minS = newS;
                minVariance = solver.cplex.getValue(solver.quadraticPart);
            }
            upperBound = beta * newS + solver.cplex.getValue(solver.linearPart);
            lowerBound = minMean + beta * Math.sqrt(minVariance);

//            System.out.println(s+"\t"+newS+"\t"+"["+minS+","+maxS+"]\t"+upperBound+"\t"+lowerBound+"\t"+(upperBound-lowerBound)/lowerBound);
            //s=newS;
        }
        return upperBound;
    }

    /**
     * Finds the best solution for a given lambda. <br>
     *
     * @param lambda Lambda. <br>
     * @return The objective value.
     */
    private double minimizeX(double lambda) throws IloException {
//        System.out.println(previousObj);
//        System.out.println(solver.cplex.getObjValue()+beta*s-lambda*s*s);
//        double val=solver.cplex.getValue(solver.linearPart)+lambda*solver.cplex.getValue(solver.quadraticPart)+beta*s-lambda*s*s;
        BasisStatus[] baseVariables = solver.cplex.getBasisStatuses(solver.x), baseConstraints = solver.cplex.getBasisStatuses(solver.constraints);
        solver.changeObjective(lambda);
        solver.cplex.setBasisStatuses(solver.x, baseVariables, solver.constraints, baseConstraints);

        solver.cplex.solve();
//        if(initial==solver.cplex.getObjValue() && solver.cplex.isDualFeasible())
//            end=true;
//        System.out.println(solver.cplex.getObjValue()+beta*s-lambda*s*s);
//        for (double xi : solver.cplex.getValues(solver.x)) {
//            System.out.print(xi+" ");
//        }
//        System.out.println("");
//        System.out.println(lambda+"\t"+(linearValue+lambda*linearVariance)+"\t"+solver.cplex.getObjValue());
        iterations += solver.cplex.getNiterations();
//        System.out.println((solver.cplex.getObjValue()+beta*s-lambda*s*s));
//        if(previousObj-(solver.cplex.getObjValue()+beta*s-lambda*s*s)<0.005)
//            end=true;
//        previousObj=solver.cplex.getObjValue()+beta*s-lambda*s*s;
//        System.out.println(previousObj);
        return solver.cplex.getObjValue() + beta * s - lambda * s * s;
    }

}
