package jvbench.axpy;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorSpecies;


public class Axpy {

    static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_MAX;
    static final int SPECIES_LENGTH = SPECIES.length();

    public static void scalar(double a, double [] dx, double [] dy, int n) {
        for (int i = 0; i < n; i++) {
            dy[i] += a * dx[i];
        }

    }

    public static void scalarFMA(double a, double [] dx, double [] dy, int n) {
        for (int i = 0; i < n; i++) {
//            dy[i] += a * dx[i];
            dy[i] = Math.fma(a, dx[i], dy[i]);
        }

    }

    public static void vector(double a, double [] dx, double [] dy, int n) {
        int i;

//        DoubleVector v_a = DoubleVector.broadcast(SPECIES, a);

        int limit = SPECIES.loopBound(n);
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector v_dx = DoubleVector.fromArray(SPECIES, dx, i);
            DoubleVector v_dy = DoubleVector.fromArray(SPECIES, dy, i);
            DoubleVector v_res = v_dx.mul(a).add(v_dy); // NOTE it calls automatically broadcast (passing a scalar or vector should not result in any changes of performance)
            v_res.intoArray(dy, i); // jfr
        }

        for (; i < n; i++) {
            dy[i] += a * dx[i];
        }

    }

    public static void vectorIndexInRange(double a, double [] dx, double [] dy, int n) {
        int i;

//        DoubleVector v_a = DoubleVector.broadcast(SPECIES, a);

//        int limit = SPECIES.loopBound(n);
        for (i = 0; i < n; i += SPECIES_LENGTH) {
            VectorMask<Double> mask = SPECIES.indexInRange(i, n);
            DoubleVector v_dx = DoubleVector.fromArray(SPECIES, dx, i, mask);
            DoubleVector v_dy = DoubleVector.fromArray(SPECIES, dy, i, mask);
            DoubleVector v_res = v_dx.mul(a).add(v_dy);
            v_res.intoArray(dy, i, mask); // jfr
        }

    }


    public static void vectorBroadcastExternal(double a, double [] dx, double [] dy, int n) {
        int i;

        // Note test if create bradcast out the loop and not every time in mul() affects the performance
        DoubleVector v_a = DoubleVector.broadcast(SPECIES, a);

        int limit = SPECIES.loopBound(n);
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector v_dx = DoubleVector.fromArray(SPECIES, dx, i);
            DoubleVector v_dy = DoubleVector.fromArray(SPECIES, dy, i);
            DoubleVector v_res = v_dx.mul(v_a).add(v_dy); // NOTE it calls automatically broadcast (passing a scalar or vector should not result in any changes of performance)
            v_res.intoArray(dy, i); // jfr
        }

        for (; i < n; i++) {
            dy[i] += a * dx[i];
        }

    }

    public static void vectorFMA(double a, double [] dx, double [] dy, int n) {
        int i;

        DoubleVector v_a = DoubleVector.broadcast(SPECIES, a);

        int limit = SPECIES.loopBound(n);
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector v_dx = DoubleVector.fromArray(SPECIES, dx, i);
            DoubleVector v_dy = DoubleVector.fromArray(SPECIES, dy, i);
            // test if using fma is better than manual mul and add
            DoubleVector v_res = v_dx.fma(v_a, v_dy);
            v_res.intoArray(dy, i); // jfr
        }

        for (; i < n; i++) {
            dy[i] += a * dx[i];
        }

    }

}
