/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ConicBB;

import conicSimplex.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author Andres Gomez.
 */
public class FileGenerator {

    public static void main(String[] args) throws IOException {
        String sContinuous = "java -cp ./dist/LagrangeanConic.jar conicSimplex.ConicSimplex";
        String sInteger = "java -cp ./dist/LagrangeanConic.jar ConicBB.Main";
        
//        int[] sizes={5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
//        int[] sizes={10,30,50};
        int[] sizes = {20};
        double[] betas = {0.1,  1.0};
        int[] ranks = {100};
        double[] diagonals = {0.1, 1.0};
        double[] densities = { 0.1};
        String[] abs = {"false"};
        long[] seeds = {12345};
        int[] methods = {3,4};
        try (FileWriter out = new FileWriter(new File("./runsInteger.bat"))) {
            for (int size : sizes) {
                for (double beta : betas) {
                    for (int rank : ranks) {
                        for (double diagonal : diagonals) {
                            for (double density : densities) {
                                for (String ab : abs) {
                                    for (long seed : seeds) {
                                        for (int method : methods) {
//                                            out.write(sContinuous + " " + size + " " + beta + " "+rank+" "+diagonal+" "+density+ " " + ab+" " + seed   + " " + method + "\n");
                                             out.write(sInteger + " " + size + " " + beta + " "+rank+" "+diagonal+" "+density+ " " + ab+" " + seed   + " " + method + "\n");
                                        }
                                    }

                                }
                            }
                        }
                    }

                }
            }
        }
    }
}
