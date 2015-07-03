/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lagrangean;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author Andres Gomez.
 */
public class FileGenerator {

    public static void main(String[] args) throws IOException {
        String s = "java -jar ./dist/LagrangeanConic.jar";
//        int[] sizes={5,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};
//        int[] sizes={10,30,50};
        int[] sizes = {30};
        double[] betas = {0.5, 1.0, 2.0, 3.0, 5, 10};
        long[] seeds = {12345, 23451, 34512, 45123, 51234};
        int[] ranks = {0};
        String[] abs = {"false", "true"};
        int[] methods = {1, 3,  5, 7, 9};
        //String[] methods={"MP"};
        try (FileWriter out = new FileWriter(new File("./runs.bat"))) {
            for (int size : sizes) {
                for (double beta : betas) {
                    for (long seed : seeds) {
                        for (int rank : ranks) {
                            for (String ab : abs) {
                                for (int method : methods) {
                                    out.write(s + " " + size + " " + beta + " " + seed + " " + rank + " " + ab + " " + method + "\n");
                                }
                            }
                        }
                    }

                }
            }
        }
    }
}

