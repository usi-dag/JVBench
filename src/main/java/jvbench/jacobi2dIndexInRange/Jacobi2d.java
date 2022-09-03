package jvbench.jacobi2dIndexInRange;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorSpecies;

public  class Jacobi2d {

    static private final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_MAX;
    static private final int SPECIES_LENGTH = SPECIES.length();

    public static void main(String[] args) {
        int size = 10000;
        int steps = 14;

        final double [][] a = new double[size][size];
        final double [][] b = new double[size][size];

        int i, j;
        for (i = 0; i < size; i++) {
            for (j = 0; j < size; j++) {
                a[i][j] = ((double) i * (j + 2) + 2) / size;
                b[i][j] = ((double) i * (j + 3) + 3) / size;
            }
        }

        vector(steps, size, a, b);
    }

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
        int size_x = n - 1;


        for (i = 1; i < size_y; i++) {
            for (j = 1; j < size_x; j += SPECIES_LENGTH) {
                VectorMask<Double> mask = SPECIES.indexInRange(j, size_x);
                xU = DoubleVector.fromArray(SPECIES, a[i], j, mask);
                xUtop = DoubleVector.fromArray(SPECIES, a[i-1], j, mask);
                xUbottom = DoubleVector.fromArray(SPECIES, a[i+1], j, mask);
                xUleft = DoubleVector.fromArray(SPECIES, a[i], j-1, mask);
                xUright = DoubleVector.fromArray(SPECIES, a[i], j+1, mask);
                xUtmp = xUleft.add(xUright);
                xUtmp = xUtmp.add(xUtop);
                xUtmp = xUtmp.add(xUbottom);
                xUtmp = xUtmp.add(xU);
                xUtmp = xUtmp.mul(0.2);
                xUtmp.intoArray(b[i], j, mask);
            }
        }
    }
}
