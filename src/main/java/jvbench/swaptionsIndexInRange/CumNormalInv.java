package jvbench.swaptionsIndexInRange;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

public class CumNormalInv {

    static double[] a  = {
        2.50662823884,
                -18.61500062529,
                41.39119773534,
                -25.44106049637
    };

    static double [] b = {
        -8.47351093090,
                23.08336743743,
                -21.06224101826,
                3.13082909833
    };

    static double [] c = {
        0.3374754822726147,
                0.9761690190917186,
                0.1607979714918209,
                0.0276438810333863,
                0.0038405729373609,
                0.0003951896511919,
                0.0000321767881768,
                0.0000002888167364,
                0.0000003960315187
    };
    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_MAX;


    public static double cumNormalInv(double u) {
        // Returns the inverse of cumulative normal distribution function.
        // Reference: Moro, B., 1995, "The Full Monte," RISK (February), 57-58.

        double x, r;

        x = u - 0.5;
        if( Math.abs (x) < 0.42 )
        {
            r = x * x;
            r = x * ((( a[3]*r + a[2]) * r + a[1]) * r + a[0])/
                    ((((b[3] * r+ b[2]) * r + b[1]) * r + b[0]) * r + 1.0);
            return (r);
        }

        r = u;
        if( x > 0.0 ) r = 1.0 - u;
        r = Math.log(-Math.log(r));
        r = c[0] + r * (c[1] + r *
                (c[2] + r * (c[3] + r *
                        (c[4] + r * (c[5] + r * (c[6] + r * (c[7] + r*c[8])))))));
        if( x < 0.0 ) r = -r;

        return (r);

    }

    public static void cumNormalInvVector(DoubleArray u, int iu, DoubleArray output, int iOutput, VectorMask<Double> mask) {
        // Returns the inverse of cumulative normal distribution function.
        // Reference: Moro, B., 1995, "The Full Monte," RISK (February), 57-58.

        DoubleVector   x;
        DoubleVector   r1;
        DoubleVector   r;

        DoubleVector   zero    = DoubleVector.zero(SPECIES);
        DoubleVector   one     = DoubleVector.broadcast(SPECIES, 1.0);
        DoubleVector   Cons1   = DoubleVector.broadcast(SPECIES, 0.5);
        DoubleVector   Cons2   = DoubleVector.broadcast(SPECIES, 0.42);
        DoubleVector   vU      = DoubleVector.fromArray(SPECIES, u.array, iu, mask);

        DoubleVector   a0      = DoubleVector.broadcast(SPECIES, a[0]);
        DoubleVector   a1      = DoubleVector.broadcast(SPECIES, a[1]);
        DoubleVector   a2      = DoubleVector.broadcast(SPECIES, a[2]);
        DoubleVector   a3      = DoubleVector.broadcast(SPECIES, a[3]);

        DoubleVector   b0      = DoubleVector.broadcast(SPECIES, b[0]);
        DoubleVector   b1      = DoubleVector.broadcast(SPECIES, b[1]);
        DoubleVector   b2      = DoubleVector.broadcast(SPECIES, b[2]);
        DoubleVector   b3      = DoubleVector.broadcast(SPECIES, b[3]);

        DoubleVector   c0      = DoubleVector.broadcast(SPECIES, c[0]);
        DoubleVector   c1      = DoubleVector.broadcast(SPECIES, c[1]);
        DoubleVector   c2      = DoubleVector.broadcast(SPECIES, c[2]);
        DoubleVector   c3      = DoubleVector.broadcast(SPECIES, c[3]);
        DoubleVector   c4      = DoubleVector.broadcast(SPECIES, c[4]);
        DoubleVector   c5      = DoubleVector.broadcast(SPECIES, c[5]);
        DoubleVector   c6      = DoubleVector.broadcast(SPECIES, c[6]);
        DoubleVector   c7      = DoubleVector.broadcast(SPECIES, c[7]);
        DoubleVector   c8      = DoubleVector.broadcast(SPECIES, c[8]);

        VectorMask<Double> mask1;
        VectorMask<Double>  mask2;
        VectorMask<Double>  mask3;

        x = vU.sub(0.5); //  _MM_SUB_f64(vU,Cons1 ,gvl);


        r = x.mul(x); // _MM_MUL_f64(x,x ,gvl);

//        r = _MM_DIV_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(a3,r ,gvl),a2,gvl),r,gvl),a1,gvl),r,gvl),a0,gvl),x,gvl),_MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(b3,r ,gvl),b2,gvl),r,gvl),b1,gvl),r,gvl),b0,gvl),r,gvl),one,gvl),gvl);
        r = x.mul(((a3.mul(r).add(a2)).mul(r).add(a1)).mul(r).add(a0)).div((((b3.mul(r).add(b2)).mul(r).add(b1)).mul(r).add(b0)).mul(r).add(1.0), mask);

        // SECOND PART
        r1 = vU;
        mask2  = x.compare(VectorOperators.GT, 0.0); // _MM_VFGT_f64(x,zero,gvl);
        DoubleVector merge = one.sub(vU);
        r1 = DoubleVector.zero(SPECIES).add(r1, mask2.not()).add(merge, mask2); //  _MM_SUB_f64_MASK(r1,one,vU,mask2,gvl); //sub(vs2,vs1)
        Cons1 = r1.lanewise(VectorOperators.LOG, mask); // _MM_LOG_f64(r1,gvl);
        r1 = Cons1.lanewise(VectorOperators.NEG);// _MM_VFSGNJN_f64(Cons1,Cons1,gvl);
        r1 = r1.lanewise(VectorOperators.LOG, mask); // _MM_LOG_f64(r1,gvl);

//        r1 = _MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(_MM_ADD_f64(_MM_MUL_f64(c8,r1 ,gvl),c7,gvl),r1,gvl),c6,gvl),r1,gvl),c5,gvl),r1,gvl),c4,gvl),r1,gvl),c3,gvl),r1,gvl),c2,gvl),r1,gvl),c1,gvl),r1,gvl),c0,gvl);
        r1 = (((((((c8.mul(r1).add(c7)).mul(r1).add(c6)).mul(r1).add(c5)).mul(r1).add(c4)).mul(r1).add(c3)).mul(r1).add(c2)).mul(r1).add(c1)).mul(r1).add(c0);
        mask3  = x.lt(0); //  _MM_VFLT_f64(x,zero,gvl);
        r1 = DoubleVector.zero(SPECIES).add(r1, mask3).add(r1.lanewise(VectorOperators.NEG), mask3.not());

        mask1  = x.abs().lt(Cons2); // _MM_VFLT_f64(_MM_VFSGNJX_f64(x,x,gvl),Cons2,gvl);
        r = DoubleVector.zero(SPECIES).add(r1, mask1).add(r, mask1.not()); // _MM_MERGE_f64(r1,r, mask1,gvl);

        r.intoArray(output.array, iOutput, mask); // _MM_STORE_f64(output,r,gvl);
    }
}
