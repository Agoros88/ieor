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
public class QuadraticSolver extends ContinuousSolver{
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
    
    protected long timePrimal,timeDual;

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
    public QuadraticSolver(DagObjective dag, double Omega) throws IloException {

        cplex = new IloCplex();
        cplex.setOut(null);
        cplex.setParam(IloCplex.IntParam.Threads, 1);
        cplex.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Primal);
        init(dag, 0, 0);
        timePrimal=0;
        timeDual=0;
        this.Omega = Omega;
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
    private void init(DagObjective dag, double lambda, int method) {
        try {
            this.dag = dag;
//            cplex.clearModel();

            // Inits the variables.
            x = new IloNumVar[dag.arcs.size()];
            for (int i = 0; i < x.length; i++) {
//                x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x(" + (int)dag.arcs.get(i)[0] + "," + (int)dag.arcs.get(i)[1]+")");
                x[i] = cplex.numVar(0, 1, IloNumVarType.Float, ""+i);
            }

            // Objective
//            IloLinearNumExpr obj = cplex.scalProd(x, parameters.getCost());
//            cplex.addObjective(IloObjectiveSense.Maximize, obj);
            // Constraints
            /*Flow conservation*/
            IloLinearNumExpr[] flowConservation = new IloLinearNumExpr[dag.vertices];
            flowConservation[0] = cplex.linearNumExpr(-1);
            flowConservation[flowConservation.length - 1] = cplex.linearNumExpr(1);
            obj = cplex.numExpr();
            linearPart = cplex.linearNumExpr();
            quadraticPart = cplex.quadNumExpr();

            for (int i = 1; i < flowConservation.length - 1; i++) {
                flowConservation[i] = cplex.linearNumExpr();
            }
            for (int i = 0; i < dag.arcs.size(); i++) {
                int[] arc = dag.arcs.get(i);
                flowConservation[arc[0]].addTerm(x[i], 1);
                flowConservation[arc[1]].addTerm(x[i], -1);
                obj = cplex.sum(obj, cplex.prod(dag.objective.c[i], x[i]));
                linearPart = cplex.sum(linearPart, cplex.prod(dag.objective.c[i], x[i]));
                if (lambda > 0) {
                    obj = cplex.sum(obj, cplex.prod(x[i], x[i], dag.objective.Matrix(i, i) * lambda));
                }
                quadraticPart = cplex.sum(quadraticPart, cplex.prod(x[i], x[i], dag.objective.Matrix(i, i)));
                for (int j = i + 1; j < dag.arcs.size(); j++) {
                    if (dag.objective.Matrix(i, j) != 0) {
                        if (lambda > 0) {
                            obj = cplex.sum(obj, cplex.prod(x[i], x[j], 2 * dag.objective.Matrix(i, j) * lambda));
                        }
                        quadraticPart = cplex.sum(quadraticPart, cplex.prod(x[i], x[j], 2 * dag.objective.Matrix(i, j)));
                    }
                }

            }
            constraints = new IloRange[flowConservation.length];
            for (int i = 0; i < flowConservation.length; i++) {

                constraints[i] = cplex.addEq(flowConservation[i], 0, "Flow conservation at " + i);
            }

            cplex.addMinimize(obj);
//            upperBound=cplex.addGe(constraint, parameters.getThreshold(),"Upper bound");

        } catch (IloException ex) {
            Logger.getLogger(QuadraticSolver.class.getName()).log(Level.SEVERE, null, ex);
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

    /**
     * Changes the objective. <br>
     *
     * @param lambda The new lambda.
     * @throws ilog.concert.IloException
     */
    private void changeObjective(double lambda) throws IloException {
        //cplex.getObjective().clearExpr();
       
        obj = cplex.sum(linearPart, cplex.prod(quadraticPart, lambda));
        
        cplex.getObjective().setExpr(obj);
    }
    
    public NodeInfo solveNode(NodeInfo info) throws IloException
    {
        long time;
        double t=info.t;
        double[] lowerBounds=new double[x.length], upperBounds=new double[x.length];
        if(Double.isInfinite(t))
        {
            double[][] bounds=resetBounds();
            
            cplex.solve();
            lowerBounds=bounds[0];
            upperBounds=bounds[1];
        }
        else
        {
            lowerBounds=info.lowerBounds;
            upperBounds=info.upperBounds;
            updateBounds(lowerBounds, upperBounds);
            changeObjective(Omega/(2*t));
            
            cplex.setBasisStatuses(x, info.basisVar, constraints, info.basisRange);
            cplex.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Dual);
            time=System.currentTimeMillis();
            cplex.solve();
            timeDual+=System.currentTimeMillis()-time;
        }
        cplex.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Primal);
        time=System.currentTimeMillis();
        double newT=Math.sqrt(cplex.getValue(quadraticPart));
        BasisStatus[] var,con;
        while (Math.abs(t-newT)>TOLERANCE)
        {
            t=newT;
            var=cplex.getBasisStatuses(x);
            con= cplex.getBasisStatuses(constraints);
            changeObjective(Omega/(2*t));
            cplex.setBasisStatuses(x, var, constraints, con);
            cplex.solve();
            newT=Math.sqrt(cplex.getValue(quadraticPart));
        }
        timePrimal+=System.currentTimeMillis()-time;
        sol=cplex.getValues(x);
        objValue=cplex.getValue(linearPart)+Omega*newT;
//        System.out.println(objValue);
        return new NodeInfo(t, cplex.getBasisStatuses(x), cplex.getBasisStatuses(constraints), lowerBounds, upperBounds,objValue);
    }

   

}
