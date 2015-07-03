/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ConicBB;

import ilog.cplex.IloCplex.BasisStatus;
import ilog.cplex.IloCplex.BranchDirection;

/**
 * Node information to use in branch and bound. <br>
 * @author AGomez
 */
public class NodeInfo implements Comparable<NodeInfo>{
    //--------------------------------------------------------------------------
    // Static ID
    //--------------------------------------------------------------------------
    /**
     * Number of nodes created.
     */
    private static int count=0;;
    //--------------------------------------------------------------------------
    // Attributes
    //--------------------------------------------------------------------------
    /**
     * Id of the node.
     */
    protected int id;
    /**
     * Value of t.
     */
    protected double t;
    
    /**
     * Basis status in the node.
     */
    protected BasisStatus[] basisVar,basisRange;
    
    /**
     * Lower and upper bounds for the variables.
     */
    protected double[] lowerBounds,upperBounds;
    
    /**
     * Bound on the objective.
     */
    protected double objectiveBound;
    
    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    /**
     * Constructor by parameters. <br>
     * @param t The previous value for t. <br>
     * @param basisVar The basis status for the variables. <br>
     * @param basisRange The basis status for the constraints. <br>
     * @param lowerBounds The variable lower bounds. <br>
     * @param upperBounds The variable upper bounds. <br>
     * @param objectiveBound The bound on the objective.
     */
    public NodeInfo(double t, BasisStatus[] basisVar, BasisStatus[] basisRange, double[] lowerBounds, double[] upperBounds, double objectiveBound) {
        this.t = t;
        this.basisVar = basisVar;
        this.basisRange = basisRange;
        this.lowerBounds = lowerBounds;
        this.upperBounds = upperBounds;
        this.objectiveBound = objectiveBound;
        this.id=count++;
    }
    
    
    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------
    
    /**
     * Copies the node info after branching in the argument position.
     * @param pos The position to branch. <br>
     * @return An array containing the two new nodes to process.
     */
    public NodeInfo[] branchCopy(int pos)
    {
        double[] lowerBounds1=new double[lowerBounds.length],
                lowerBounds2=new double[lowerBounds.length],
                upperBounds1=new double[lowerBounds.length],
                upperBounds2=new double[lowerBounds.length];
        System.arraycopy(lowerBounds, 0, lowerBounds1, 0, lowerBounds.length);
        System.arraycopy(lowerBounds, 0, lowerBounds2, 0, lowerBounds.length);
        System.arraycopy(upperBounds, 0, upperBounds1, 0, upperBounds.length);
        System.arraycopy(upperBounds, 0, upperBounds2, 0, upperBounds.length);
        upperBounds1[pos]=0;
        lowerBounds2[pos]=1;
        return new NodeInfo[]{new NodeInfo(t, basisVar, basisRange, lowerBounds1, upperBounds1,objectiveBound),
        new NodeInfo(t, basisVar, basisRange, lowerBounds2, upperBounds2,objectiveBound)};
        
    }

    @Override
    public int compareTo(NodeInfo o) {
        if (Double.compare(objectiveBound, o.objectiveBound)!=0)
            return Double.compare(objectiveBound, o.objectiveBound);
        return Integer.compare(id, o.id);
            
    }
    
}
