package jvbench.somierIndexInRange;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class SomierTest {
    int ntSteps = 1;
    int n = 32;

    @Test
    void indexInRange() {
        double [] vS = new double[3];
        double [] xS = new double[3];
        double [] vV = new double[3];
        double [] xV = new double[3];

        Somier.init(ntSteps, n);
        Somier.scalar();
        vS[0] = Somier.v[0][n / 2][n / 2][n / 2];
        vS[1] = Somier.v[1][n / 2][n / 2][n / 2];
        vS[2] = Somier.v[2][n / 2][n / 2][n / 2];
        xS[0] = Somier.x[0][n / 2][n / 2][n / 2];
        xS[1] = Somier.x[1][n / 2][n / 2][n / 2];
        xS[2] = Somier.x[2][n / 2][n / 2][n / 2];

        Somier.init(ntSteps, n);
        Somier.vector();
        vV[0] = Somier.v[0][n / 2][n / 2][n / 2];
        vV[1] = Somier.v[1][n / 2][n / 2][n / 2];
        vV[2] = Somier.v[2][n / 2][n / 2][n / 2];
        xV[0] = Somier.x[0][n / 2][n / 2][n / 2];
        xV[1] = Somier.x[1][n / 2][n / 2][n / 2];
        xV[2] = Somier.x[2][n / 2][n / 2][n / 2];

        assertArrayEquals(vS, vV);
        assertArrayEquals(xS, xV);
    }
}