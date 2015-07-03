/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package conicSimplex;

import java.util.Random;

/**
 * Models the matrix in the objective. <br>
 *
 * @author AGomez.
 */
public class PSDObjective {

    //--------------------------------------------------------------------------
    // Attributes
    //--------------------------------------------------------------------------

    /**
     * The rank of the problem.
     */
    protected int rank;
    
    /**
     * Linear cost vector.
     */
    public double[] c;

    /**
     * The PSD rxr dense matrix.
     */
    protected double[][] Sigma;

    /**
     * The sparse multiplying Sigma.
     */
    protected double[][] F;

    /**
     * The diagonal matrix.
     */
    protected double[] D;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    /**
     * Constructor by parameters. <br>
     *
     * @param n The dimension of the problem. <br>
     * @param rank The rank of Sigma. <br>
     * @param d The weight of the diagonal matrix. <br>
     * @param rho The density of F. <br>
     * @param positive Whether to use only nonnegative terms. <br>
     * @param random Random number generator.
     */
    public PSDObjective(int n, int rank, double d, double rho, boolean positive, Random random) {
        this.rank = rank;

        double[][] temp = new double[rank][rank];

        // Initializes Sigma
        for (double[] temp1 : temp) {
            for (int j = 0; j < temp1.length; j++) {
                temp1[j] = random.nextDouble() * 2 - 1;
                if (positive) {
                    temp1[j] = Math.abs(temp1[j]);
                }
            }
        }
        Sigma = new double[rank][rank];
        for (int i = 0; i < Sigma.length; i++) {
            for (int j = 0; j < Sigma.length; j++) {
                Sigma[i][j] = 0;
                for (double[] t : temp) {
                    Sigma[i][j] += t[i] * t[j];
                }
            }
        }
        
        // Initializes D and c
        D=new double[n];
        c=new double[n];
        for (int i = 0; i < D.length; i++) {
            D[i]=random.nextDouble()*d;
            c[i]=random.nextDouble();
        }
        
        // Initializes F.
        F=new double[n][rank];
        for (double[] F1 : F) {
            for (int j = 0; j < F1.length; j++) {
                if (random.nextDouble()<rho) {
                    F1[j] = random.nextDouble()*2-1;
                    if (positive) {
                        F1[j] = Math.abs(F1[j]);
                    }
                }
            }   
        }
    }
    
    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------
    /**
     * Gets the coefficient of the matrix of the argument coordinates. <br>
     * @param i The row. <br>
     * @param j The column. <br>
     * @return Matrix[i][j].
     */
    public double Matrix(int i, int j)
    {
        double resp=(i==j)?D[i]:0;
        for (int k = 0; k < rank; k++) {
            for (int l = 0; l < rank; l++) {
                resp+=F[i][l]*Sigma[l][k]*F[j][k];
            }
        }
        return resp;
    }

}
