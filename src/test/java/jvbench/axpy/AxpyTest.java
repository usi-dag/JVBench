package jvbench.axpy;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class AxpyTest {

    private static final double [] x = new double[10];
    private static final double [] y = new double[10];

    private static void initVector(double [] pv, int n, double value) {
        for (int i = 0; i < n; i++) {
            pv[i] = value;
        }
    }

    @BeforeEach
    public void setup() {
        initVector(x, 10, 1.0);
        initVector(y, 10, 2.0);
    }


    @Test
    void axpyRef() {
        Axpy.scalar(1.0, x, y, 10);
        for (double value : y) {
            assertEquals(3.0, value);
        }


    }

    @Test
    void axpyVector() {
        Axpy.vector(1.0, x, y, 10);
        for (double value : y) {
            assertEquals(3.0, value);

        }
    }

    @Test
    void axpyScalarFMA() {
        Axpy.scalarFMA(1.0, x, y, 10);
        for (double value : y) {
            assertEquals(3.0, value);

        }
    }


    @Test
    void axpyVectorFMA() {
        Axpy.vectorFMA(1.0, x, y, 10);
        for (double value : y) {
            assertEquals(3.0, value);

        }
    }

    @Test
    void axpyBroadcastExternal() {
        Axpy.vectorBroadcastExternal(1.0, x, y, 10);
        for (double value : y) {
            assertEquals(3.0, value);

        }
    }


    @Test
    void axpyIndexInRange() {
        Axpy.vectorIndexInRange(1.0, x, y, 10);
        for (double value : y) {
            assertEquals(3.0, value);

        }
    }

}