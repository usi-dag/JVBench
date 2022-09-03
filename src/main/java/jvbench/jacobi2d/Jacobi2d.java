package jvbench.jacobi2d;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorSpecies;

public  class Jacobi2d {

    static private final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_MAX;
    static private final int SPECIES_LENGTH = SPECIES.length();

    public static void serial(int tSteps, int n, double [][] a, double [][] b) {
        int t, i, j;

        for (t = 0; t < tSteps; t++) {
            for (i = 1; i < n - 1; i++) {
                for (j = 1; j < n - 1; j++)
                    b[i][j] = (0.2) * (a[i][j] + a[i][j - 1] + a[i][1 + j] + a[1 + i][j] + a[i - 1][j]);
            }
            for (i = 1; i < n - 1; i++) {
                for (j = 1; j < n - 1; j++)
                    a[i][j] = (0.2) * (b[i][j] + b[i][j-1] + b[i][1+j] + b[1+i][j] + b[i-1][j]);
            }
        }

    }

    public static void vector(int tSteps, int n, double [][] a, double [][] b) {
        int t;

        for (t = 0; t < tSteps; t++)
        {
            kernelJacobi2dVector(tSteps,n, a,b);
            kernelJacobi2dVector(tSteps,n, b, a);
        }
    }


    private static void kernelJacobi2dVector(int tSteps, int n, double [][] a, double [][] b) {
        int i, j;
        DoubleVector xU;
        DoubleVector xUtmp;
        DoubleVector xUleft;
        DoubleVector xUright;
        DoubleVector xUtop;
        DoubleVector xUbottom;


        int size_y = n - 1;
        int size_x = n - 2;

        int limit = SPECIES.loopBound(size_x);

        for (i = 1; i < size_y; i++) {
            for (j = 1; j < limit; j += SPECIES_LENGTH) {
                xU = DoubleVector.fromArray(SPECIES, a[i], j);
                xUtop = DoubleVector.fromArray(SPECIES, a[i-1], j);
                xUbottom = DoubleVector.fromArray(SPECIES, a[i+1], j);
                xUleft = DoubleVector.fromArray(SPECIES, a[i], j-1);
                xUright = DoubleVector.fromArray(SPECIES, a[i], j+1);
                xUtmp = xUleft.add(xUright);
                xUtmp = xUtmp.add(xUtop);
                xUtmp = xUtmp.add(xUbottom);
                xUtmp = xUtmp.add(xU);
                xUtmp = xUtmp.mul(0.2);
                xUtmp.intoArray(b[i], j);
            }
        }

//        for (j = 1; j < limit; j += SPECIES_LENGTH) {
//
//            xU = DoubleVector.fromArray(SPECIES, first[1], j);
//            xUtop = DoubleVector.fromArray(SPECIES, first[0], j);
//            xUbottom = DoubleVector.fromArray(SPECIES, first[2], j);
//
//            for (i = 1; i <= size_y; i++) {
//                if (i != 1) {
//                    xUtop = xU;
//                    xU = xUbottom;
//                    xUbottom = DoubleVector.fromArray(SPECIES, first[i + 1], j);
//                }
//
//                xUleft = DoubleVector.fromArray(SPECIES, first[i], j-1);
//                xUright = DoubleVector.fromArray(SPECIES, first[i], j+1);
//                xUtmp = xUleft.add(xUright);
//                xUtmp = xUtmp.add(xUtop);
//                xUtmp = xUtmp.add(xUbottom);
//                xUtmp = xUtmp.add(xU);
//                xUtmp = xUtmp.mul(0.2);
//                xUtmp.intoArray(second[i], j);
//
//            }
//        }


        for (i = 1; i < n - 1; i++) {
            for (j = limit + 1; j < n - 1; j++)
                b[i][j] = (0.2) * (a[i][j] + a[i][j - 1] + a[i][1 + j] + a[1 + i][j] + a[i - 1][j]);
        }
    }
}
