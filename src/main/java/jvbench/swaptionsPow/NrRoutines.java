package jvbench.swaptionsPow;

public interface NrRoutines {
//    void nrerror(char[] error_text);
//
//    int choldc(double[][] a, int n);
//
//    void gaussj(double[][] a, int n, double[][] b, int m);
//
//    int[] ivector(long nl, long nh);
//
//    void free_ivector(int[] v, long nl, long nh);

    static DoubleArray dvector(int nl, int nh) {
        DoubleArray v;

        v = new DoubleArray(nh-nl+2, -nl + 1); // (FTYPE *) malloc((size_t)((nh - nl + 2) * sizeof(FTYPE)));
//        if (!v) nrerror("allocation failure in dvector()");
        return v; // - nl + 1;
    }


    static DoubleArray [] dmatrix(int nrl, int nrh, int ncl, int nch) {
        int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;

//        double[][] m = new double[nrow + 1][];

        DoubleArray [] m = new DoubleArray[nrow + 1];


        // allocate rows and set pointers to them
        m[nrl] = new DoubleArray(nrow * ncol + 1); // new double[nrow * ncol + 1]; // (FTYPE *) malloc((size_t)((nrow*ncol+1)*sizeof(FTYPE)));
//        if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
//        m[nrl] += 1;
//        m[nrl] -= ncl;
//        m[nrl].setIterator(1 - ncl);

        for (i = nrl + 1; i <= nrh; i++) {
            m[i] = new DoubleArray(nrow*ncol+1, ncol); // m[i-1]; //m[i - 1]/* +ncol */; // TODO memory leak 395 GiB
        }

        // return pointer to array of pointers to rows
        return m;
    }

//    void free_dmatrix(double[][] m, long nrl, long nrh, long ncl, long nch);
}
