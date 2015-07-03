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
import ilog.cplex.IloCplex.BasisStatus;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Quadratic solver for the shortest path.
 *
 * @author Andres Gomez.
 */
public class QCPContinuousSolver extends ContinuousSolver{
    //--------------------------------------------------------------------------
    // Constants
    //--------------------------------------------------------------------------
    /**
     * Tolerance.
     */
    private static final double TOLERANCE=0.001;
    //--------------------------------------------------------------------------
    // Attributes
    //--------------------------------------------------------------------------

    /**
     * Cplex solver.
     */
    protected IloCplex cplex;

    /**
     * X variables.
     */
    protected IloNumVar[] x;
    
    /**
     * Variable modeling the standard deviation.
     */
    protected IloNumVar sd;

    /**
     * Constraints.
     */
    protected IloRange[] constraints;

    /**
     * Objective.
     */
    protected IloNumExpr obj;

    /**
     * Dag.
     */
    protected DagObjective dag;

    /**
     * Value of the linear part.
     */
    protected IloNumExpr linearPart;

    /**
     * Value of the quadratic part.
     */
    protected IloNumExpr quadraticPart;
    
    /**
     * Coefficient of the nonlinear term.
     */
    protected double Omega;
    
   
    
    /**
     * Objective value.
     */
    protected double  objValue;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    /**
     * Constructor. <br>
     *
     * @param dag The dag. <br>
     * @param Omega
     * @throws ilog.concert.IloException
     */
    public QCPContinuousSolver(DagObjective dag, double Omega) throws IloException {

        cplex = new IloCplex();
        cplex.setOut(null);
        cplex.setParam(IloCplex.IntParam.Threads, 1);
        cplex.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Barrier);
        cplex.setParam(IloCplex.IntParam.BarAlg, 1);
        this.Omega = Omega;
        init(dag, 0);
        
    }

    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------   
    /**
     * Creates the CPLEX model. <br>
     *
     * @param dag Directed Acyclic Graph. <br>
     * @param lambda Coefficient for the standard deviation. <br>
     * @param method Method for solving the problem. <br>
     */
    private void init(DagObjective dag, int method) {
        try {
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
        ((IloLinearNumExpr) obj).addTerm(sd, Omega);
        cplex.addLe(conic, 0.0, "Conic constraint");

//        cplex.use(new QuadraticSolver(dag,x));
//        cplex.setParam(IloCplex.IntParam.MIQCPStrat, 1);
        cplex.addMinimize(obj);

        } catch (IloException ex) {
            Logger.getLogger(QCPContinuousSolver.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    /**
     * Updates the bounds on the variables. <br>
     *
     * @param lb Lower bounds on the variables. <br>
     * @param ub Upper bounds on the variables.
     * @throws ilog.concert.IloException
     */
    private void updateBounds(double[] lb, double[] ub) throws IloException {
        for (int i = 0; i < ub.length; i++) {
            x[i].setLB(lb[i]);
            x[i].setUB(ub[i]);

        }
    }
    
     /**
     * Resets the bounds on the variables. <br>
     *
     * @return the lower (0) and upper bounds (1).
     * @throws ilog.concert.IloException
     */
    private double[][] resetBounds() throws IloException {
        double[][] bounds=new double[2][x.length];
        for (int i = 0; i < x.length; i++) {
            x[i].setLB(0);
            x[i].setUB(1);
            bounds[1][i]=1;
        }
        return bounds;
    }

  
    
    @Override
    public NodeInfo solveNode(NodeInfo info) throws IloException
    {
        double t=info.t;
        double[] lowerBounds, upperBounds;
        if(Double.isInfinite(t))
        {
            double[][] bounds=resetBounds();
            lowerBounds=bounds[0];
            upperBounds=bounds[1];
        }
        else
        {
            lowerBounds=info.lowerBounds;
            upperBounds=info.upperBounds;
            updateBounds(lowerBounds, upperBounds);
            
        }
        cplex.solve();
        sol=cplex.getValues(x);
        objValue=cplex.getObjValue();
//        System.out.println(objValue);
        return new NodeInfo(0, null, null, lowerBounds, upperBounds,objValue);
    }

   

}
