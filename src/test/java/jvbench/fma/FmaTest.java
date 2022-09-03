package jvbench.fma;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class FmaTest {


    @Test
    public void simpleTest() {
        float res1 = Math.fma(1, 2, 3);
        float res2 = (float) 1 * 2 + 3;

        assertEquals(res1, res2);
    }

    @Test
    public void infPosTest() {
        float res1 = Math.fma(Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY);
        float res2 = Float.POSITIVE_INFINITY * Float.POSITIVE_INFINITY + Float.POSITIVE_INFINITY;

        assertEquals(res1, res2);
    }

    @Test
    public void infNegTest() {
        float res1 = Math.fma(Float.NEGATIVE_INFINITY, Float.NEGATIVE_INFINITY, Float.NEGATIVE_INFINITY);
        float res2 = Float.NEGATIVE_INFINITY * Float.NEGATIVE_INFINITY + Float.NEGATIVE_INFINITY;

        assertEquals(res1, res2);
    }

    @Test
    public void onlyMul() {
        float res1 = Math.fma(3, 4, Float.NaN);
        float res2 =  3 * 4;
    }

    @Test
    public void onlyAdd() {
        float res1 = Math.fma(3, Float.NaN, 4);
        float res2 =  3 + 4;
    }
}
