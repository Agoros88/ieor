/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ConicBB;

import ilog.concert.IloException;

/**
 * A solver for branch and bound. <br>
 * @author AGomez.
 */
public abstract class ContinuousSolver {
    //--------------------------------------------------------------------------
    // Attributes
    //--------------------------------------------------------------------------
     /**
     * Optimal solution.
     */
    protected double[] sol;
    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------
     /**
     * Solves the node using the warm-start information. <br>
     * @param info The warm start information. <br>
     * @return The updated node info information. <br>
     * @throws ilog.concert.IloException
     */
    public abstract NodeInfo solveNode(NodeInfo info) throws IloException;
}
