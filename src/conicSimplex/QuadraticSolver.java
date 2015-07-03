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
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Quadratic solver for the shortest path.
 *
 * @author Andres Gomez.
 */
public class QuadraticSolver {
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
     * Constraints
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
    

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
   

    /**
     * Constructor.
     *
     * @throws ilog.concert.IloException
     */
    public QuadraticSolver() throws IloException {

        cplex = new IloCplex();
        cplex.setOut(null);
        cplex.setParam(IloCplex.IntParam.Threads, 1);
        cplex.setParam(IloCplex.IntParam.RootAlg, IloCplex.Algorithm.Primal);
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
     * @param lambda Coefficient for the standard deviation. <br>
     * @param method Method for solving the problem.
     */
    protected void init(DagObjective dag, double lambda, int method) {
        try {
            this.dag=dag;
            cplex.clearModel();

            // Inits the variables.
            x = new IloNumVar[dag.arcs.size()];
            for (int i = 0; i < x.length; i++) {
//                x[i] = cplex.numVar(0, 1, IloNumVarType.Bool, "x(" + (int)dag.arcs.get(i)[0] + "," + (int)dag.arcs.get(i)[1]+")");
                x[i] = cplex.numVar(0, 1, IloNumVarType.Float, "x" + (int)dag.arcs.get(i)[0] + "," + (int)dag.arcs.get(i)[1]);
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
            linearPart=cplex.linearNumExpr();
            quadraticPart=cplex.quadNumExpr();
            
            for (int i = 1; i < flowConservation.length - 1; i++) {
                flowConservation[i] = cplex.linearNumExpr();
            }
            for (int i = 0; i < dag.arcs.size(); i++) {
                int[] arc = dag.arcs.get(i);
                flowConservation[arc[0]].addTerm(x[i], 1);
                flowConservation[arc[1]].addTerm(x[i], -1);
                obj=cplex.sum(obj,cplex.prod(dag.objective.c[i], x[i]));
                linearPart=cplex.sum(linearPart,cplex.prod(dag.objective.c[i], x[i]));
                if(lambda>0)
                    obj=cplex.sum(obj,cplex.prod(x[i], x[i],dag.objective.Matrix(i, i)*lambda));
                quadraticPart=cplex.sum(quadraticPart,cplex.prod(x[i], x[i],dag.objective.Matrix(i, i)));
                for (int j = i+1; j < dag.arcs.size(); j++) {
                    if(dag.objective.Matrix(i, j)!=0)
                    {
                        if(lambda>0)
                            obj = cplex.sum(obj, cplex.prod(x[i], x[j],2*dag.objective.Matrix(i, j)*lambda));
                        quadraticPart = cplex.sum(quadraticPart, cplex.prod(x[i], x[j],2*dag.objective.Matrix(i, j)));
                    }
                }

            }
            constraints=new IloRange[flowConservation.length];
            for (int i = 0; i < flowConservation.length; i++) {

                constraints[i]=cplex.addEq(flowConservation[i], 0, "Flow conservation at " + i);
            }
           

            cplex.addMinimize(obj);
//            upperBound=cplex.addGe(constraint, parameters.getThreshold(),"Upper bound");
            

        } catch (IloException ex) {
            Logger.getLogger(QuadraticSolver.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
    
    /**
     * Changes the objective. <br>
     * @param lambda The new lambda.
     */
    public void changeObjective(double lambda) throws IloException
    {
        //cplex.getObjective().clearExpr();
        if(Double.isInfinite(lambda))
        {
            obj=quadraticPart;
        }
        else
            obj=cplex.sum(linearPart, cplex.prod(quadraticPart, lambda));
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
    
  
}
