/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package lagrangean;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.concert.IloQuadNumExpr;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Conic solver for the shortest path.
 *
 * @author Andres Gomez.
 */
public class ConicSolver {
    //--------------------------------------------------------------------------
    // Attributes
    //--------------------------------------------------------------------------

    /**
     * Cplex solver.
     */
    protected IloCplex cplex;

    /**
     * X variables
     */
    protected IloNumVar[] x;

    /**
     * Objective.
     */
    protected IloLinearNumExpr obj;

    /**
     * Quadratic part.
     */
    protected IloNumExpr quadraticPart;

    /**
     * Variable modeling the standard deviation.
     */
    protected IloNumVar sd;

    /**
     * Upper bound constraint.
     */
    protected IloRange upperBound;

    /**
     * Redundant integer variable.
     */
    protected IloNumVar redundant;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     *
     * @param dag Directed Acyclic Graph. <br>
     * @param beta Coefficient for the standard deviation.
     * @throws ilog.concert.IloException
     */
//    public ConicSolver(Dag dag, double beta) throws IloException {
//        init(dag, beta);
//        cplex.setOut(null);
//        cplex.setParam(IloCplex.IntParam.Threads, 1);
//    }
    /**
     * Constructor.
     *
     * @throws ilog.concert.IloException
     */
    public ConicSolver() throws IloException {

        cplex = new IloCplex();
        cplex.setOut(null);
        cplex.setParam(IloCplex.IntParam.Threads, 1);
    }

    //--------------------------------------------------------------------------
    // Getters and Setters
    //--------------------------------------------------------------------------
    /**
     * Gets the cplex object corresponding to this model. <br>
     *
     * @return cplex.
     */
    public IloCplex getCplex() {
        return cplex;
    }

    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------   
    /**
     * Creates the CPLEX model. <br>
     *
     * @param dag Directed Acyclic Graph. <br>
     * @param beta Coefficient for the standard deviation. <br>
     * @param configuration Configuration for the solver: <br>
     * 0: barrier. <br>
     * 1: QCP approximation. <br>
     * 2: Outer approximation. <br>
     * 3: Low rank approximation. <br>
     * 4: Low rank approximation+outer approximation.
     * @throws ilog.concert.IloException
     */
    protected void init(Dag dag, double beta, int configuration) throws IloException {

        cplex.clearModel();

        // Inits the variables.
        x = new IloNumVar[dag.arcs.size()];
        for (int i = 0; i < x.length; i++) {
//                x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x" + (int)dag.arcs.get(i)[0] + "," + (int)dag.arcs.get(i)[1]);
            x[i] = cplex.numVar(0, 1, IloNumVarType.Float, "x" + (int) dag.arcs.get(i)[0] + "," + (int) dag.arcs.get(i)[1]);
        }
        sd = cplex.numVar(0, Double.POSITIVE_INFINITY, "z");

        // Objective
//            IloLinearNumExpr obj = cplex.scalProd(x, parameters.getCost());
//            cplex.addObjective(IloObjectiveSense.Maximize, obj);
        // Constraints
            /*Flow conservation*/
        IloLinearNumExpr[] flowConservation = new IloLinearNumExpr[dag.vertices];
        flowConservation[0] = cplex.linearNumExpr(-1);
        flowConservation[flowConservation.length - 1] = cplex.linearNumExpr(1);
        IloNumExpr conic = cplex.prod(sd, sd, -1);
        obj = cplex.linearNumExpr();
        quadraticPart = cplex.quadNumExpr();
        obj.addTerm(sd, beta);

        for (int i = 1; i < flowConservation.length - 1; i++) {
            flowConservation[i] = cplex.linearNumExpr();
        }
        for (int i = 0; i < dag.arcs.size(); i++) {
            double[] arc = dag.arcs.get(i);
            flowConservation[(int) arc[0]].addTerm(x[i], 1);
            flowConservation[(int) arc[1]].addTerm(x[i], -1);
            conic = cplex.sum(conic, cplex.prod(x[i], x[i], dag.covariances[i][i]));
            quadraticPart = cplex.sum(quadraticPart, cplex.prod(x[i], x[i], dag.covariances[i][i]));
            for (int j = i + 1; j < dag.covariances[i].length; j++) {
                if (dag.covariances[i][j] != 0) {
                    conic = cplex.sum(conic, cplex.prod(x[i], x[j], 2 * dag.covariances[i][j]));
                    quadraticPart = cplex.sum(quadraticPart, cplex.prod(x[i], x[j], 2 * dag.covariances[i][j]));
                }
            }
            obj.addTerm(x[i], arc[2]);

        }
        for (int i = 0; i < flowConservation.length; i++) {

            cplex.addEq(flowConservation[i], 0, "Flow conservation at " + i);
        }

        cplex.addMinimize(obj);

        if (configuration == 0) {
            cplex.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Barrier);
            if (beta > 0.1) {
                cplex.addLe(conic, 0.0, "Conic constraint");
            }
        } else if (configuration == 1) {
            cplex.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Barrier);
            if (beta > 0.1) {
                lowRank(dag);
            }
        } else if (configuration == 2) {
            addRedundant();
            cplex.setParam(IloCplex.IntParam.MIQCPStrat, 1);
            if (beta > 0.1) {
                cplex.addLe(conic, 0.0, "Conic constraint");
            }
        } else if (configuration == 3) {
            addRedundant();
            cplex.setParam(IloCplex.IntParam.MIQCPStrat, 1);
            if (beta > 0.1) {
                lowRank(dag);
            }
        } else if (configuration == 4) {
            addRedundant();
            cplex.setParam(IloCplex.IntParam.MIQCPStrat, 2);
            if (beta > 0.1) {
                cplex.addLe(conic, 0.0, "Conic constraint");
            }
        } else if (configuration == 5) {
            addRedundant();
            cplex.setParam(IloCplex.IntParam.MIQCPStrat, 2);
            if (beta > 0.1) {
                lowRank(dag);
            }
//            upperBound=cplex.addGe(constraint, parameters.getThreshold(),"Upper bound");

        }

    }

/**
 * Adds a redundant integer variable to the model.
 *
 * @throws ilog.concert.IloException
 */
private void addRedundant() throws IloException {
        redundant = cplex.boolVar();
        cplex.addLe(redundant, 0);
        cplex.setParam(IloCplex.IntParam.Symmetry, 0);
        cplex.setParam(IloCplex.IntParam.AggFill, 0);
        cplex.setParam(IloCplex.IntParam.AggInd, 0);
        cplex.setParam(IloCplex.IntParam.BndStrenInd, 0);
        cplex.setParam(IloCplex.IntParam.CoeRedInd, 0);
        cplex.setParam(IloCplex.IntParam.DepInd, 0);
        cplex.setParam(IloCplex.IntParam.PreDual, -1);
//        cplex.setParam(IloCplex.BooleanParam.PerInd, false);
        cplex.setParam(IloCplex.BooleanParam.PreInd, false);
        cplex.setParam(IloCplex.IntParam.PreLinear, 0);
        cplex.setParam(IloCplex.IntParam.PrePass, 0);
        cplex.setParam(IloCplex.IntParam.PreslvNd, -1);
        cplex.setParam(IloCplex.IntParam.RelaxPreInd, 0);
        cplex.setParam(IloCplex.IntParam.RepeatPresolve, 0);
        cplex.setParam(IloCplex.IntParam.Reduce, 0);

    }

    /**
     * Reformulates the problem exploiting low rank information. <br>
     *
     * @param dag The graph. <r>
     * @throws IloException
     */
    protected void lowRank(Dag dag) throws IloException {
        IloQuadNumExpr conic = cplex.quadNumExpr();
        conic.addTerm(-1, sd, sd);
        if (dag.qVectors != null) {
            IloNumVar[] y = new IloNumVar[dag.qVectors.length];
            double[] coefficients = new double[y.length];
            for (int i = 0; i < y.length; i++) {
                y[i] = cplex.numVar(0, Double.POSITIVE_INFINITY, IloNumVarType.Float, "|q_" + i + "x|");
                cplex.addGe(y[i], cplex.scalProd(x, dag.qVectors[i]));
                cplex.addGe(y[i], cplex.prod(-1, cplex.scalProd(x, dag.qVectors[i])));
                coefficients[i] = 1;
            }

            conic.addTerms(coefficients, y, y);
            cplex.addLe(conic, 0);
        }
    }

    /**
     *
     * @throws IloException
     */
    private void addOuterApproximation(Dag dag) throws IloException {
        if (dag.qVectors != null) {
            IloNumVar[] y = new IloNumVar[dag.qVectors.length];
            for (int i = 0; i < y.length; i++) {
                y[i] = cplex.numVar(0, Double.POSITIVE_INFINITY, IloNumVarType.Float, "|q_" + i + "x|");
                cplex.addGe(y[i], cplex.scalProd(x, dag.qVectors[i]));
                cplex.addGe(y[i], cplex.prod(-1, cplex.scalProd(x, dag.qVectors[i])));
                cplex.addLe(y[i], sd);
            }

            if (dag.qVectors.length >= 2) {
                for (int i = 0; i < y.length; i++) {
                    for (int j = i + 1; j < y.length; j++) {
                        cplex.addLe(cplex.sum(y[i], y[j]), cplex.prod(Math.sqrt(2), sd));
                    }
                }
            }

            if (dag.qVectors.length >= 3) {
                for (int i = 0; i < y.length; i++) {
                    for (int j = i + 1; j < y.length; j++) {
                        for (int k = j + 1; k < y.length; k++) {
                            cplex.addLe(cplex.sum(y[i], y[j], y[k]), cplex.prod(Math.sqrt(3), sd));
                        }
                    }
                }
            }

            if (dag.qVectors.length >= 4) {
                for (int i = 0; i < y.length; i++) {
                    for (int j = i + 1; j < y.length; j++) {
                        for (int k = j + 1; k < y.length; k++) {
                            for (int l = k + 1; l < y.length; l++) {
                                cplex.addLe(cplex.sum(y[i], y[j], y[k], y[l]), cplex.prod(2, sd));
                            }
                        }
                    }
                }
            }

        }
    }

}
