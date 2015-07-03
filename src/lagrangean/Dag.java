/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lagrangean;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * A sparse directed acyclic graph. <br>
 *
 * @author Andres Gomez.
 */
public class Dag {

    //--------------------------------------------------------------------------
    // Attributes
    //--------------------------------------------------------------------------
    /**
     * Number of vertices.
     */
    public int vertices;

    /**
     * List of arcs. <br>
     * arc: {head,tail,mean}.
     */
    public List<double[]> arcs;

    /**
     * Covariance matrix of the arcs.
     */
    public double[][] covariances;

    /**
     * Vectors defining the pd matrix. <br>
     * dim 1: nank, dim 2:n.
     */
    public double[][] qVectors;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    /**
     * Constructor by parameters. <br>
     *
     * @param vertices The number of vertices. <br>
     * @param arcs The arcs. <br>
     * @param covariances The covariance matrix.
     */
    public Dag(int vertices, List<double[]> arcs, double[][] covariances) {
        this.vertices = vertices;
        this.arcs = arcs;
        this.covariances = covariances;
    }

    /**
     * Constructor by parameters. <br>
     *
     * @param vertices The number of vertices. <br>
     * @param arcs The arcs. <br>
     * @param covariances The covariance matrix. <br>
     * @param qVectors Matrix decomposition.
     */
    public Dag(int vertices, List<double[]> arcs, double[][] covariances, double[][] qVectors) {
        this.vertices = vertices;
        this.arcs = arcs;
        this.covariances = covariances;
        this.qVectors = qVectors;
    }

    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------
    /**
     * Generates a model in a directed acyclic graph. <br>
     *
     * @param size The size of the grid. <br>
     * @param type The type of Q. <br>
     * @param abs Whether to use only nonnegative values.
     * @param r Random number generator. <br>
     * @return The Parameters object.
     */
    public static Dag generateGrid(int size, int type, boolean abs, Random r) {
        int[][] position = new int[size][size];
        int number = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                position[i][j] = number++;
            }
        }
        List<double[]> arcs = new ArrayList<>();
        number = 0;

        // Generates arcs with costs and means.
        double rn1;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i < size - 1) {
                    rn1 = r.nextDouble();
                    arcs.add(new double[]{position[i][j], position[i + 1][j], rn1});
                }
                if (j < size - 1) {
                    rn1 = r.nextDouble();
                    arcs.add(new double[]{position[i][j], position[i][j + 1], rn1});
                }
            }
        }

        double[][] covariances = new double[arcs.size()][arcs.size()];

        if (type == -1)// Diagonal
        {
            for (int i = 0; i < covariances.length; i++) {
                covariances[i][i] = Math.pow(r.nextDouble(), 2);
            }
        } 

        // Rank 1
        //        double[] vector= new double[arcs.size()];
        //        for (int i = 0; i < vector.length; i++) {
        //            vector[i]=-1+2*r.nextDouble();
        //        }
        //        for (int i = 0; i < covariances.length; i++) {
        //            for (int j = 0; j < covariances.length; j++) {
        //                covariances[i][j]=vector[i]*vector[j];
        //               
        //            }     
        //        }
        // Rank n [-1,1]
        //        double[][] matrix= new double[arcs.size()][arcs.size()];
        //        for (int i = 0; i < matrix.length; i++) {
        //            for (int j = 0; j < matrix[i].length; j++) {
        //                matrix[i][j]=-1+2*r.nextDouble();         
        //            }
        //        }
        //        
        //        for (int i = 0; i < covariances.length; i++) {
        //            for (int j = 0; j < covariances.length; j++) {
        //                covariances[i][j]=0;
        //                for (int k = 0; k < matrix[i].length; k++) {
        //                    covariances[i][j]+=matrix[i][k]*matrix[j][k];               
        //                }
        //            }     
        //        }
        // Rank n [0,1]
        //        double[][] matrix= new double[arcs.size()][arcs.size()];
        //        for (double[] matrix1 : matrix) {
        //            for (int j = 0; j < matrix1.length; j++) {
        //                matrix1[j] = r.nextDouble();         
        //            }
        //        }
        //        
        //        for (int i = 0; i < covariances.length; i++) {
        //            for (int j = 0; j < covariances.length; j++) {
        //                covariances[i][j]=0;
        //                for (int k = 0; k < matrix[i].length; k++) {
        //                    covariances[i][j]+=matrix[i][k]*matrix[j][k];               
        //                }
        //            }     
        //        }
        else if (type <= -2 && type >=-5)// Diagonal dominant
        {
            double total;
                double density = 0.0;
                if(type==-2)
                {
                    density=1.1;
                }
                else if (type==-3)
                {
                    density=0.5;
                }
                else if (type==-4)
                {
                    density=0.1;
                }
                else if (type==-5)
                    density=0.02;
                for (int i = 0; i < covariances.length; i++) {
                    for (int j = i + 1; j < covariances.length; j++) {
                        if (r.nextDouble() < density) {
                            rn1 = -1 + 2 * r.nextDouble();
                            if(abs)
                                rn1=Math.abs(rn1);
                            covariances[i][j] = rn1;
                            covariances[j][i] = rn1;
                        }
                    }
                    total = 0;
                    for (int j = 0; j < covariances[i].length; j++) {
                        total += Math.abs(covariances[i][j]);
                    }
                    covariances[i][i] = 2 * total + r.nextDouble();
                    arcs.get(i)[2] = r.nextDouble() * 2 * covariances[i][i];
                }
        }
        
            return new Dag(size * size, arcs, covariances);
        
    }

    /**
     * Generates a model in a directed acyclic graph. <br>
     *
     * @param size The size of the grid. <br>
     * @param rank Rank of the matrix. <br>
     * @param abs Whether to use a matrix with only positive entries or not.
     * <br>
     * @param r Random number generator. <br>
     * @return The Parameters object.
     */
    public static Dag generateGridRank(int size, int rank, boolean abs, Random r) {
        int[][] position = new int[size][size];
        int number = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                position[i][j] = number++;
            }
        }
        List<double[]> arcs = new ArrayList<>();
        number = 0;

        // Generates arcs with costs and means.
        double rn1;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i < size - 1) {
                    rn1 = Math.abs(-1 + 2 * r.nextDouble());
                    arcs.add(new double[]{position[i][j], position[i + 1][j], rn1});
                }
                if (j < size - 1) {
                    rn1 = Math.abs(-1 + 2 * r.nextDouble());
                    arcs.add(new double[]{position[i][j], position[i][j + 1], rn1});
                }
            }
        }

        double[][] covariances = new double[arcs.size()][arcs.size()];

        // Diagonal
//        for (int i = 0; i < covariances.length; i++) {
//            covariances[i][i]=Math.pow(r.nextDouble(), 2);
//        }
        // Rank 1
//        double[] vector= new double[arcs.size()];
//        for (int i = 0; i < vector.length; i++) {
//            vector[i]=-1+2*r.nextDouble();
//        }
//        for (int i = 0; i < covariances.length; i++) {
//            for (int j = 0; j < covariances.length; j++) {
//                covariances[i][j]=vector[i]*vector[j];
//               
//            }     
//        }
        double[][] qVectors = null;
        if (rank > 0) {
            qVectors = new double[rank][arcs.size()];
            for (double[] qVector : qVectors) {
                for (int j = 0; j < qVector.length; j++) {
                    qVector[j] = -1 + 2 * r.nextDouble();
                    if (abs) {
                        qVector[j] = Math.abs(qVector[j]);
                    }
                }
            }

            for (int i = 0; i < covariances.length; i++) {
                for (int j = 0; j < covariances.length; j++) {
                    covariances[i][j] = 0;
                    for (double[] qVector : qVectors) {
                        covariances[i][j] += qVector[i] * qVector[j];
                    }
                }
            }
        } else {
            qVectors = new double[arcs.size()][arcs.size()];
            for (double[] matrix1 : qVectors) {
                for (int j = 0; j < matrix1.length; j++) {
                    matrix1[j] = -1 + 2 * r.nextDouble();
                    if (abs) {
                        matrix1[j] = Math.abs(matrix1[j]);
                    }
                }
            }

            for (int i = 0; i < covariances.length; i++) {
                for (int j = 0; j < covariances.length; j++) {
                    covariances[i][j] = 0;
                    for (double[] qVector : qVectors) {
                        covariances[i][j] += qVector[i] * qVector[j];
                    }
                }
            }
        }

        return new Dag(size * size, arcs, covariances, qVectors);
    }

}
