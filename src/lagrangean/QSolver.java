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
import ilog.concert.IloObjective;
import ilog.concert.IloQuadNumExpr;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Quadratic solver for the shortest path.
 *
 * @author Andres Gomez.
 */
public class QSolver {
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
    protected Dag dag;
    
    /**
     * Value of the linear part.
     */
    protected IloNumExpr linearPart;
    
    /**
     * Value of the quadratic part.
     */
    protected IloNumExpr quadraticPart;
    
    /**
     * Current and previous objective.
     */
    protected double currentObj=Double.POSITIVE_INFINITY,improvement=0;
    
    /**
     * Iterations.
     */
    protected int iterations;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
   

    /**
     * Constructor.
     *
     * @throws ilog.concert.IloException
     */
    public QSolver() throws IloException {

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
     * @param lambda Coefficient for the standard deviation.
     */
    protected void init(Dag dag, double lambda, int method) {
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
                double[] arc = dag.arcs.get(i);
                flowConservation[(int)arc[0]].addTerm(x[i], 1);
                flowConservation[(int)arc[1]].addTerm(x[i], -1);
                obj=cplex.sum(obj,cplex.prod(arc[2], x[i]));
                linearPart=cplex.sum(linearPart,cplex.prod(arc[2], x[i]));
                if(lambda>0)
                    obj=cplex.sum(obj,cplex.prod(x[i], x[i],dag.covariances[i][i]*lambda));
                quadraticPart=cplex.sum(quadraticPart,cplex.prod(x[i], x[i],dag.covariances[i][i]));
                for (int j = i+1; j < dag.covariances[i].length; j++) {
                    if(dag.covariances[i][j]!=0)
                    {
                        if(lambda>0)
                            obj = cplex.sum(obj, cplex.prod(x[i], x[j],2*dag.covariances[i][j]*lambda));
                        quadraticPart = cplex.sum(quadraticPart, cplex.prod(x[i], x[j],2*dag.covariances[i][j]));
                    }
                }

            }
            constraints=new IloRange[flowConservation.length];
            for (int i = 0; i < flowConservation.length; i++) {

                constraints[i]=cplex.addEq(flowConservation[i], 0, "Flow conservation at " + i);
            }
           

            cplex.addMinimize(obj);
//            upperBound=cplex.addGe(constraint, parameters.getThreshold(),"Upper bound");
            if((method%2)==1)
            {
                lowRank(dag);
            }

        } catch (IloException ex) {
            Logger.getLogger(QSolver.class.getName()).log(Level.SEVERE, null, ex);
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
    
    
    
    /**
     * Changes the objective. <br>
     * @param lambda The new lambda. <br>
     * @param currentObj The current objective.<br>
     * @param improvement The expected improvement.
     */
    public void changeObjective(double lambda, double currentObj,double improvement) throws IloException
    {
        this.currentObj=currentObj;
        this.improvement=improvement;
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
    
    /**
     * Instructs the solver to use the callback.
     */
    public void useCallback() throws IloException
    {
        cplex.setParam(IloCplex.IntParam.ItLim, 30);
//        cplex.use(new Aborter());
    }

    
    /**
     * Reformulates the problem exploiting low rank information. <br>
     *
     * @param dag The graph. <r>
     * @throws IloException
     */
    protected void lowRank(Dag dag) throws IloException {
        quadraticPart = cplex.quadNumExpr();
        
        if (dag.qVectors != null) {
            IloNumVar[] y = new IloNumVar[x.length+dag.qVectors.length];
           System.arraycopy(x, 0, y, 0, x.length);
           IloRange[] constraints2=new IloRange[constraints.length+dag.qVectors.length];
           System.arraycopy(constraints, 0, constraints2, 0, constraints.length);
            for (int i = 0; i < dag.qVectors.length; i++) {
                y[x.length+i] = cplex.numVar(0, Double.POSITIVE_INFINITY, IloNumVarType.Float, "q_" + (x.length+i) + "x");
                constraints2[constraints.length+i]=(IloRange) cplex.addEq(y[x.length+i], cplex.scalProd(x, dag.qVectors[i]));
                quadraticPart=cplex.sum(quadraticPart,cplex.prod(y[x.length+i], y[x.length+i]));
            }
            x=y;
            constraints=constraints2;
        }
    }
    
    /**
     * Class that aborts the optimization after a single iteration.
     */
    private class Aborter extends IloCplex.SimplexCallback
    {
        //----------------------------------------------------------------------
        // Attributes
        //----------------------------------------------------------------------
       
        /**
         * Whether a feasible solution has been found.
         */
        private int count=0;
        
        
        //----------------------------------------------------------------------
        // Constructor
        //----------------------------------------------------------------------
        @Override
        protected void main() throws IloException {
//            iterations++;
            try
            {
                double obj=getObjValue();
                abort();
//            if(getObjValue()<=currentObj-improvement)
//            {
//                currentObj=obj;
//                iterations++;
//                
////                System.out.println(getDualInfeasibility());
//               
//            }
//            else if(obj<currentObj)
//            {
//                iterations++;
//                abort();
//            }
            }catch(ilog.cplex.CpxException e)
            {
//                System.out.println(count++);
            }
            
//            if(isFeasible()&& getNiterations()>1 )
//                abort();
        }
        
    }

}
