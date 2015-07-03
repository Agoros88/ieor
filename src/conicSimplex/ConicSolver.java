/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package conicSimplex;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;

/**
 * Conic solver.
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
     */
    protected void init(DagObjective dag, double beta, int configuration) throws IloException {

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
            int[] arc = dag.arcs.get(i);
            flowConservation[arc[0]].addTerm(x[i], 1);
            flowConservation[arc[1]].addTerm(x[i], -1);
            conic = cplex.sum(conic, cplex.prod(x[i], x[i], dag.objective.Matrix(i, i)));
            quadraticPart = cplex.sum(quadraticPart, cplex.prod(x[i], x[i], dag.objective.Matrix(i, i)));
            for (int j = i + 1; j < dag.arcs.size(); j++) {
                if (dag.objective.Matrix(i, j) != 0) {
                    conic = cplex.sum(conic, cplex.prod(x[i], x[j], 2 * dag.objective.Matrix(i, j)));
                    quadraticPart = cplex.sum(quadraticPart, cplex.prod(x[i], x[j], 2 * dag.objective.Matrix(i, j)));
                }
            }
            obj.addTerm(x[i], dag.objective.c[i]);

        }
        for (int i = 0; i < flowConservation.length; i++) {

            cplex.addEq(flowConservation[i], 0, "Flow conservation at " + i);
        }

        cplex.addLe(conic, 0.0, "Conic constraint");
        cplex.addMinimize(obj);
        
        if (configuration == 11) {
            cplex.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Barrier);
        } else if (configuration == 21) {
            cplex.setParam(IloCplex.IntParam.MIQCPStrat, 2);
            addRedundant();
        }
//        else if (configuration == 1) {
//            cplex.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Barrier);
//            if (beta > 0.1) {
//                lowRank(dag);
//            }
//        } else if (configuration == 2) {
//            addRedundant();
//            cplex.setParam(IloCplex.IntParam.MIQCPStrat, 1);
//            if (beta > 0.1) {
//                cplex.addLe(conic, 0.0, "Conic constraint");
//            }
//        } else if (configuration == 3) {
//            addRedundant();
//            cplex.setParam(IloCplex.IntParam.MIQCPStrat, 1);
//            if (beta > 0.1) {
//                lowRank(dag);
//            }
//        } else if (configuration == 4) {
//            addRedundant();
//            cplex.setParam(IloCplex.IntParam.MIQCPStrat, 2);
//            if (beta > 0.1) {
//                cplex.addLe(conic, 0.0, "Conic constraint");
//            }

//            upperBound=cplex.addGe(constraint, parameters.getThreshold(),"Upper bound");
    }

    /**
     * Adds a redundant integer variable to the model.
     *
     * @throws ilog.concert.IloException
     */
    private void addRedundant() throws IloException {
        redundant = cplex.boolVar();
        cplex.addLe(redundant, 0);

        // Turns off presolve
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

        // Turning off cuts
        cplex.setParam(IloCplex.IntParam.EachCutLim, 0);
        cplex.setParam(IloCplex.IntParam.CutPass, -1);
        cplex.setParam(IloCplex.DoubleParam.CutsFactor, 1.0);
        cplex.setParam(IloCplex.IntParam.AggCutLim, 0);
        cplex.setParam(IloCplex.IntParam.Cliques, -1);
        cplex.setParam(IloCplex.IntParam.Covers, -1);
        cplex.setParam(IloCplex.IntParam.DisjCuts, -1);

        cplex.setParam(IloCplex.IntParam.FlowCovers, -1);

        cplex.setParam(IloCplex.IntParam.FlowPaths, -1);
        cplex.setParam(IloCplex.IntParam.FracCuts, -1);
        cplex.setParam(IloCplex.IntParam.GUBCovers, -1);
        cplex.setParam(IloCplex.IntParam.ImplBd, -1);
        cplex.setParam(IloCplex.IntParam.MIRCuts, -1);
        cplex.setParam(IloCplex.IntParam.ZeroHalfCuts, -1);
        cplex.setParam(IloCplex.IntParam.MCFCuts, -1);// Multi commodity cuts.
        cplex.setParam(IloCplex.IntParam.LiftProjCuts, -1); // Lift and project cuts.
        // Turning off heursitics
        cplex.setParam(IloCplex.IntParam.RINSHeur, -1);
        cplex.setParam(IloCplex.IntParam.HeurFreq, -1);

    }


}
