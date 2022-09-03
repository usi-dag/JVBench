package jvbench.swaptionsIndexInRange;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.LongVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorSpecies;


public class RanUnif {

    private static final VectorSpecies<Double> DOUBLE_SPECIES = DoubleVector.SPECIES_MAX;
    private static final int DOUBLE_SPECIES_LENGTH = DOUBLE_SPECIES.length();

    private static final VectorSpecies<Long> LONG_SPECIES = LongVector.SPECIES_MAX;
    private static final int LONG_SPECIES_LENGTH = DOUBLE_SPECIES.length();

    public static double ranUnif(double [] seeds) {
        // uniform random number generator

        double ix, k1;
        double dRes;

        ix = seeds[0];
        k1 = ix / 127773L;
        ix = 16807L * (ix - k1 * 127773L) - k1 * 2836L;
        if (ix < 0) ix = ix + 2147483647L;
        seeds[0] = ix;
        dRes = (ix * 4.656612875e-10);
        return (dRes);

    }

    public static void ranUnifVector(double[] s, int iFactors, int iN, int BLOCKSIZE, DoubleArray[] randZ) {

        // uniform random number generator
//        unsigned long int gvl = __builtin_epi_vsetvl(BLOCKSIZE, __epi_e64, __epi_m1);
//        int limit = DOUBLE_SPECIES.loopBound(BLOCKSIZE);
        DoubleVector    k1;
//        LongVector      zero;
        VectorMask<Double> mask1;
        DoubleVector    dRes;

//        LongVector    cons1     = LongVector.broadcast(LONG_SPECIES, 127773L);// _MM_SET_i64(127773,gvl);
        DoubleVector    cons2     = DoubleVector.broadcast(DOUBLE_SPECIES, 16807L);//_MM_SET_i64(16807,gvl);
//        LongVector    cons3     = LongVector.broadcast(LONG_SPECIES, 2836L);//_MM_SET_i64(2836,gvl);
//        LongVector    cons4     = LongVector.broadcast(LONG_SPECIES, 2147483647L);//_MM_SET_i64(2147483647,gvl);
//        DoubleVector    cons5     = DoubleVector.broadcast(DOUBLE_SPECIES, 4.656612875E-10); //_MM_SET_f64(4.656612875E-10,gvl);
        DoubleVector    xSeed       = DoubleVector.fromArray(DOUBLE_SPECIES, s, 0); // _MM_LOAD_i64(s,gvl);

        for (int l=0;l<=iFactors-1;++l){
            for (int j=1;j<=iN-1;++j){
                for (int b = 0; b < BLOCKSIZE; b += DOUBLE_SPECIES_LENGTH) {
                    VectorMask<Double> mask = DOUBLE_SPECIES.indexInRange(b, BLOCKSIZE);
                    k1 = xSeed.div(127773L); // _MM_DIV_i64(xSeed,cons1,gvl);
                    xSeed = (cons2.mul(xSeed.sub(k1.mul(127773L)))).sub(k1.mul(2836L));  // _MM_SUB_i64(_MM_MUL_i64(cons2,_MM_SUB_i64(xSeed,_MM_MUL_i64(k1,cons1,gvl),gvl),gvl), _MM_MUL_i64(k1,cons3,gvl) , gvl);
                    // zero    = LongVector.zero(LONG_SPECIES); // _MM_SET_i64(0,gvl);
                    mask1 = xSeed.lt(0);//  _MM_VMSLT_i64(xSeed,zero,gvl);
                    xSeed = xSeed.add(2147483647L, mask1); // _MM_ADD_i64_MASK(xSeed,xSeed,cons4,mask1,gvl);FF

                    dRes = xSeed.mul(4.656612875E-10);//  _MM_MUL_f64( cons5,_MM_VFCVT_F_X_f64(xSeed,gvl),gvl);

//                _MM_STORE_f64(&randZ[l][BLOCKSIZE*j], dRes,gvl);
                    dRes.intoArray(randZ[l].array, (BLOCKSIZE * j) + b + randZ[l].getIterator(), mask);
                }

            }
        }
    }
}
