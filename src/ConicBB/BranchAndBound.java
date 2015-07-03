/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ConicBB;

import conicSimplex.DagObjective;
import ilog.concert.IloException;
import java.util.Arrays;
import java.util.TreeSet;

/**
 * Class in charge of branch and bound.
 * @author AGomez
 */
public class BranchAndBound {
    //--------------------------------------------------------------------------
    // Constants
    //--------------------------------------------------------------------------
    /**
     * Precision in numerical computations.
     */
    protected static final double PRECISION=1e-5;
    //--------------------------------------------------------------------------
    // Attributes
    //--------------------------------------------------------------------------
    /**
     * Best bound found in branch and bound tree.
     */
    protected double upperBound;
    
    /**
     * Solver for the quadratic problems.
     */
    protected ContinuousSolver solver;
    
    /**
     * Nodes to process.
     */
    private TreeSet<NodeInfo> nodes; 
    
    /**
     * Number of branch and bound nodes processed.
     */
    protected int iterations;
    
    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    /**
     * Constructor.
     * @param dag The dag. <br>
     * @param Omega The coefficient of the nonlinear term. <br>
     * @param configuration The configuration to be used. <br>
     * @throws ilog.concert.IloException
     */
    public BranchAndBound(DagObjective dag, double Omega, double configuration) throws IloException {
        this.upperBound = Double.POSITIVE_INFINITY;
        if(configuration==3)
            solver=new QuadraticSolver(dag,Omega);
        else if(configuration==4)
            solver=new QCPContinuousSolver(dag,Omega);
        nodes=new TreeSet<>();
    }
    
    
    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------
    /**
     * Branches after solving a node.
     */
    private NodeInfo processNode(NodeInfo info) throws IloException
    {
        if(info.objectiveBound>upperBound)
        {
            nodes.clear();
            return null;
        }
        info=solver.solveNode(info);
        if(info.objectiveBound>upperBound)
        {
//            System.out.println("cutoff"+"\t" +info.objectiveBound+"\t"+upperBound);
            return null;
        }
        double mostInfeasible=0,value;
        int  infeasibleIndex=-1;
        for (int i = 0; i < solver.sol.length; i++) {
            value=solver.sol[i];
            if(Math.min(value, 1-value)>mostInfeasible)
            {
                mostInfeasible=Math.min(value, 1-value);
                infeasibleIndex=i;
            }
        }
        if (mostInfeasible<PRECISION)
        {
            upperBound=info.objectiveBound;
            return null;
//            System.out.println("int"+"\t "+"\t"+upperBound);
        }
        else
        {
            NodeInfo[] infos=info.branchCopy(infeasibleIndex);
//            nodes.addAll(Arrays.asList(infos));
            if(solver.sol[infeasibleIndex]<0.5)
            {
                nodes.add(infos[1]);
                return infos[0];
            }
            else
            {
                nodes.add(infos[0]);
                return infos[1];
            }
        }
        
    }
    
    /**
     * Begins branch and bound.
     */
    public void doBranchBound() throws IloException
    {
        
//        nodes.add(new NodeInfo(Double.POSITIVE_INFINITY, null, null, null, null, Double.NEGATIVE_INFINITY));
        iterations=1;
        System.out.println("Iteration \t Lower bound \t Upper bound \t Gap \t Time");
        double gap=Double.POSITIVE_INFINITY;
        long startingTime=System.currentTimeMillis();
        NodeInfo info=new NodeInfo(Double.POSITIVE_INFINITY, null, null, null, null, Double.NEGATIVE_INFINITY),info2;
        while(!nodes.isEmpty() || info!=null)
        {
            if(info==null)
                info=nodes.pollFirst();
            info=processNode(info);
            
//            if(iterations%Integer.highestOneBit(iterations)==0 || info==null)
//            {
//                info2=nodes.first();
//            gap=Double.isInfinite(upperBound)?Double.POSITIVE_INFINITY:100*(upperBound-info2.objectiveBound)/(upperBound+1e-10);
//                System.out.println(iterations+"\t"+info2.objectiveBound+"\t"+upperBound+"\t"+gap+"%"+"\t"+(System.currentTimeMillis()-startingTime)/1000.0);
//            }
            iterations++;
            if(System.currentTimeMillis()-startingTime>7200000)
            {
//                System.out.println(iterations+"\t"+info.objectiveBound+"\t"+upperBound+"\t"+gap+"%");
                return;
            }
        }
    }
    
    
    /**
     * Gets a lower bound on the objective. <br>
     * @return lowerBound.
     */
    public double getLowerBound()
    {
        if(nodes.isEmpty())
            return upperBound;
        else
            return nodes.first().objectiveBound;
    }
    
    
}
