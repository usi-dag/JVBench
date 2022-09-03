package jvbench.jacobi2d;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class Jacobi2dTest {

    private static final int size = 20;

    private static final double [][] a = new double[size][size];
    private static final double [][] b = new double[size][size];

    private static final double [][] v_a = new double[size][size];
    private static final double [][] v_b = new double[size][size];


    private void initArray(int n, double [][] a, double [][] b) {
        int i, j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                a[i][j] = ((double) i * (j + 2) + 2) / n;
                b[i][j] = ((double) i * (j + 3) + 3) / n;
            }
        }
    }

    @BeforeEach
    void setup() {
        initArray(size, a, b);
        initArray(size, v_a, v_b);
    }

    @Test
    void standard() {

        Jacobi2d.serial(1, size, a, b);

        Jacobi2d.vector(1, size, v_a, v_b);


        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                assertTrue((a[i][j] - 0.01) < v_a[i][j] && v_a[i][j] < (a[i][j] + 0.01));
                assertTrue((b[i][j] - 0.01) < v_b[i][j] && v_b[i][j] < (b[i][j] + 0.01));
            }
        }
    }


}