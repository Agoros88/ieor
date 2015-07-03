/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ConicBB;

import conicSimplex.*;
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
    protected IloNumExpr obj;

    /**
     * Value of the linear part.
     */
    protected IloNumExpr linearPart;

    /**
     * Value of the quadratic part.
     */
    protected IloNumExpr quadraticPart;

    /**
     * Variable modeling the standard deviation.
     */
    protected IloNumVar sd;

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
                // Time
        cplex.setParam(IloCplex.DoubleParam.TiLim, 7200);

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
            x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x" + (int) dag.arcs.get(i)[0] + "," + (int) dag.arcs.get(i)[1]);
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
        linearPart = cplex.linearNumExpr();
        quadraticPart = cplex.quadNumExpr();

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
            ((IloLinearNumExpr) linearPart).addTerm(x[i], dag.objective.c[i]);
            ((IloLinearNumExpr) obj).addTerm(x[i], dag.objective.c[i]);

        }
        for (int i = 0; i < flowConservation.length; i++) {

            cplex.addEq(flowConservation[i], 0, "Flow conservation at " + i);
        }
        ((IloLinearNumExpr) obj).addTerm(sd, beta);
        cplex.addLe(conic, 0.0, "Conic constraint");

//        cplex.use(new QuadraticSolver(dag,x));
//        cplex.setParam(IloCplex.IntParam.MIQCPStrat, 1);
        cplex.addMinimize(obj);

        if (configuration == 1) {
            cplex.setParam(IloCplex.IntParam.MIQCPStrat, 1);
            turnOffFeatures();
        } else if (configuration == 2) {
            cplex.setParam(IloCplex.IntParam.MIQCPStrat, 2);
            turnOffFeatures();
        }

    }

    /**
     * Changes the objective. <br>
     *
     * @param lambda The new lambda.
     * @throws ilog.concert.IloException
     */
    public void changeObjective(double lambda) throws IloException {
        //cplex.getObjective().clearExpr();
        if (Double.isInfinite(lambda)) {
            obj = quadraticPart;
        } else {
            obj = cplex.sum(linearPart, cplex.prod(quadraticPart, lambda));
        }
//        for (int i = 0; i < dag.arcs.size(); i++) {
//                double[] arc = dag.arcs.get(i);
//                
//                obj=cplex.sum(obj,cplex.prod(arc[2], x[i]));
//                obj=cplex.sum(obj,cplex.prod(x[i], x[i],dag.covariances[i][i]*lambda));
//                for (int j = i+1; j < dag.covariances[i].length; j++) {
//                    if(dag.covariances[i][j]!=0)
//                    {
//                        obj = cplex.sum(obj, cplex.prod(x[i], x[j],2*dag.covariances[i][j]*lambda));
//                    }
//                }
//
//            }
        cplex.getObjective().setExpr(obj);
    }

    /**
     * Turns off different features of CPLEX
     *
     * @throws IloException
     */
    private void turnOffFeatures() throws IloException {
        // Turns off presolve
        cplex.setParam(IloCplex.IntParam.Symmetry, 0);
        cplex.setParam(IloCplex.IntParam.AggFill, 0);
        cplex.setParam(IloCplex.IntParam.AggInd, 0);
        cplex.setParam(IloCplex.IntParam.BndStrenInd, 0);
        cplex.setParam(IloCplex.IntParam.CoeRedInd, 0);
        cplex.setParam(IloCplex.IntParam.DepInd, 0);
        cplex.setParam(IloCplex.IntParam.PreDual, -1);
        cplex.setParam(IloCplex.BooleanParam.PreInd, false);
        cplex.setParam(IloCplex.IntParam.PreLinear, 0);
        cplex.setParam(IloCplex.IntParam.PrePass, 0);
        cplex.setParam(IloCplex.IntParam.PreslvNd, -1); // No presolve at nodes 
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
        cplex.setParam(IloCplex.IntParam.RINSHeur, -1); // No heuristics for improving feasible solutions
        cplex.setParam(IloCplex.IntParam.HeurFreq, -1); // No periodic heuristics

        // Branch and node selection
        cplex.setParam(IloCplex.IntParam.NodeSel, 1);
        cplex.setParam(IloCplex.IntParam.VarSel, 1);
//        cplex.setParam(IloCplex.IntParam.BrDir, -1);//Processes down branch first
        cplex.setParam(IloCplex.IntParam.DiveType, 1);//Does not dive
        cplex.setParam(IloCplex.IntParam.Probe, -1); // No probing
        cplex.setParam(IloCplex.IntParam.FPHeur, -1); // No feasibility pump
        

    }

}
