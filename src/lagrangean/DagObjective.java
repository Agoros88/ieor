package lagrangean;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
import lagrangean.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * A sparse directed acyclic graph. <br>
 *
 * @author Andres Gomez.
 */
public class DagObjective {

    //--------------------------------------------------------------------------
    // Attributes
    //--------------------------------------------------------------------------
    /**
     * Number of vertices.
     */
    public int vertices;

    /**
     * List of arcs. <br>
     * arc: {head,tail}.
     */
    public List<int[]> arcs;

    /**
     * Objective.
     */
    protected PSDObjective objective;

    //--------------------------------------------------------------------------
    // Constructor
    //--------------------------------------------------------------------------
    /**
     * Constructor by parameters. <br>
     *
     * @param vertices The number of vertices. <br>
     * @param arcs The arcs. <br>
     * @param objective The objective function.
     */
    public DagObjective(int vertices, List<int[]> arcs, PSDObjective objective) {
        this.vertices = vertices;
        this.arcs = arcs;
        this.objective = objective;
    }

    //--------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------
    /**
     * Generates a model in a directed acyclic graph, with matrix of the form F*Sigma*F'+D. <br>
     *
     * @param size The size of the grid. <br>
     * @param rank The rank of the matrix. <br>
     * @param diagonal The weight of the diagonal part. <br>
     * @param density The density of the F vector. <br>
     * @param abs Whether to use only nonnegative values.
     * @param r Random number generator. <br>
     * @return The Parameters object.
     */
    public static DagObjective generateGrid(int size, int rank, double diagonal, double density, boolean abs, Random r) {
        int[][] position = new int[size][size];
        int number = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                position[i][j] = number++;
            }
        }
        List<int[]> arcs = new ArrayList<>();
        number = 0;

        // Generates arcs with costs and means.
        double rn1;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i < size - 1) {
                    arcs.add(new int[]{position[i][j], position[i + 1][j]});
                }
                if (j < size - 1) {
                    arcs.add(new int[]{position[i][j], position[i][j + 1]});
                }
            }
        }

        PSDObjective objective2 = new PSDObjective(arcs.size(), rank, diagonal, density, abs, r);

        return new DagObjective(size
                * size, arcs, objective2);

    }

}
