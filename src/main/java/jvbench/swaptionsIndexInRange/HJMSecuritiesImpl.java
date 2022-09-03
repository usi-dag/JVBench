package jvbench.swaptionsIndexInRange;


import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorSpecies;

public class HJMSecuritiesImpl {

    static int BLOCK_SIZE = 64;
    static int NUM_TRIALS = 128; // -sm
    static int nThreads = 1; // -nt
    static int nSwaptions = 1; // -ns
    static int iN = 11;
    static int iFactors = 3;
    static double [] seed = new double[]{1979}; // new Seed(1979); // -sd
    static long swaptionSeed;
    static double [] swaptionSeedVector;
    static Parm[] swaptions;

    public static boolean useVectorAPI = true;
    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_MAX;
    private static final int SPECIES_LENGTH = SPECIES.length();

    public static void init(int numTrials, int nSwaption) {
        NUM_TRIALS = numTrials;
        nSwaptions = nSwaption;

        int i, j;

        DoubleArray[] factors;


        swaptionSeed = (long) (2147483647L * RanUnif.ranUnif(seed));


        // initialize input dataset
//        NrRoutines.dmatrix(0, iFactors - 1, 0, iN - 2);
        factors = NrRoutines.dmatrix(0, iFactors - 1, 0, iN - 2);
//the three rows store vol data for the three factors
        factors[0].set(0, .01);
        factors[0].set(1, .01);
        factors[0].set(2, .01);
        factors[0].set(3, .01);
        factors[0].set(4, .01);
        factors[0].set(5, .01);
        factors[0].set(6, .01);
        factors[0].set(7, .01);
        factors[0].set(8, .01);
        factors[0].set(9, .01);

        factors[1].set(0, .009048);
        factors[1].set(1, .008187);
        factors[1].set(2, .007408);
        factors[1].set(3, .006703);
        factors[1].set(4, .006065);
        factors[1].set(5, .005488);
        factors[1].set(6, .004966);
        factors[1].set(7, .004493);
        factors[1].set(8, .004066);
        factors[1].set(9, .003679);

        factors[2].set(0, .001000);
        factors[2].set(1, .000750);
        factors[2].set(2, .000500);
        factors[2].set(3, .000250);
        factors[2].set(4, .000000);
        factors[2].set(5, -.000250);
        factors[2].set(6, -.000500);
        factors[2].set(7, -.000750);
        factors[2].set(8, -.001000);
        factors[2].set(9, -.001250);


        swaptions = new Parm[nSwaptions];

        int k;
        for (i = 0; i < nSwaptions; i++) {
            swaptions[i] = new Parm();
            swaptions[i].Id = i;
            swaptions[i].iN = iN;
            swaptions[i].iFactors = iFactors;
            swaptions[i].dYears = 5.0 + ((int) (60 * RanUnif.ranUnif(seed))) * 0.25; //5 to 20 years in 3 month intervals

            swaptions[i].dStrike =
                    0.1 + ((int) (49 * RanUnif.ranUnif(seed))) * 0.1; //strikes ranging from 0.1 to 5.0 in steps of 0.1
            swaptions[i].dCompounding = 0;
            swaptions[i].dMaturity = 1.0;
            swaptions[i].dTenor = 2.0;
            swaptions[i].dPaymentInterval = 1.0;

            swaptions[i].pdYield = NrRoutines.dvector(0, iN - 1);
            swaptions[i].pdYield.set(0, .1);
            for (j = 1; j <= swaptions[i].iN - 1; ++j)
                swaptions[i].pdYield.set(j, swaptions[i].pdYield.get(j - 1) + .005);

            swaptions[i].ppdFactors = NrRoutines.dmatrix(0, swaptions[i].iFactors - 1, 0, swaptions[i].iN - 2);
            for (k = 0; k <= swaptions[i].iFactors - 1; ++k)
                for (j = 0; j <= swaptions[i].iN - 2; ++j) {
                    swaptions[i].ppdFactors[k].set(j, factors[k].get(j));

                }
        }

        HJM.pdPayoffDiscountFactors = NrRoutines.dvector(0, iN * BLOCK_SIZE - 1);
        HJM.pdDiscountingRatePath = NrRoutines.dvector(0, iN * BLOCK_SIZE - 1);
        HJM.ppdHJMPath = NrRoutines.dmatrix(0, iN - 1, 0, iN * BLOCK_SIZE - 1);    // **** per Trial data **** //
        HJM.pdForward = NrRoutines.dvector(0, iN - 1);
        HJM.ppdDrifts = NrRoutines.dmatrix(0, iFactors - 1, 0, iN - 2);
        HJM.pdTotalDrift = NrRoutines.dvector(0, iN - 2);


    }

    public static void benchmark(boolean vectorize) {
        useVectorAPI = vectorize;
        worker(0);
    }


    public static void main(String[] args) {
        int iSuccess = 0;
        int i, j;

        useVectorAPI = true;
       init(64, 16384);


        // start time
        int threadId = 0;
        worker(threadId);


        for (i = 0; i < nSwaptions; i++) {
            System.out.printf("Swaption %d: [SwaptionPrice: %.10f StdError: %.10f] \n",
                    i, swaptions[i].dSimSwaptionMeanPrice, swaptions[i].dSimSwaptionStdError);
        }
    }

    private static void worker(int threadId) {
        int tid = threadId;
        double[] pdSwaptionPrice = new double[2];

        int beg, end, chunksize;

        if (tid < (nSwaptions % nThreads)) {
            chunksize = nSwaptions / nThreads + 1;
            beg = tid * chunksize;
            end = (tid + 1) * chunksize;
        } else {
            chunksize = nSwaptions / nThreads;
            int offsetThread = nSwaptions % nThreads;
            int offset = offsetThread * (chunksize + 1);
            beg = offset + (tid - offsetThread) * chunksize;
            end = offset + (tid - offsetThread + 1) * chunksize;

            if (tid == nThreads - 1)
                end = nSwaptions;
            int BLOCK_SIZE_AUX;

            for (int i = beg; i < end; i++) {

                if (useVectorAPI) {
                    int limit = SPECIES.loopBound(NUM_TRIALS);
                    swaptionSeedVector = new double[limit];
                    for (int j = 0; j < limit; j++) {
                        swaptionSeedVector[j] = swaptionSeed + i; // + j + ((long) i * limit);
                    }

                    BLOCK_SIZE_AUX = BLOCK_SIZE;

                } else {
                    //                swaptionSeedVector = new long[1]; // (long *) malloc(1 * sizeof(long));
                    swaptionSeedVector =  new double[1];
                    swaptionSeedVector[0] = swaptionSeed + i;// new Seed(swaptionSeed + i);
                    BLOCK_SIZE_AUX = BLOCK_SIZE;
                }



                int iSuccess = HJM.HJM_Swaption_Blocking(pdSwaptionPrice, swaptions[i].dStrike,
                        swaptions[i].dCompounding, swaptions[i].dMaturity,
                        swaptions[i].dTenor, swaptions[i].dPaymentInterval,
                        swaptions[i].iN, swaptions[i].iFactors, swaptions[i].dYears,
                        swaptions[i].pdYield, swaptions[i].ppdFactors,
                        swaptionSeedVector, NUM_TRIALS, BLOCK_SIZE_AUX, 0);
                assert (iSuccess == 1);
                swaptions[i].dSimSwaptionMeanPrice = pdSwaptionPrice[0];
                swaptions[i].dSimSwaptionStdError = pdSwaptionPrice[1];

            }




        }


    }
}

class Seed {
    private long seed;

    public Seed(long seed) {
        this.seed = seed;
    }

    public void add(long value) {
        this.seed += value;
    }

    public long getSeed() {
        return seed;
    }

    public void setSeed(long seed) {
        this.seed = seed;
    }
}


class DoubleArray {
    public final double[] array;
    private int iterator = 0;

    public DoubleArray(int size) {
        array = new double[size];
    }

    public DoubleArray(int size, int start) {
        array = new double[size];
        iterator = start;
    }

    public int getIterator() {
        return iterator;
    }

    public void setIterator(int iterator) {
        this.iterator = iterator;
    }

    public void set(int index, double value) {
        array[iterator + index] = value;
    }

    public double get(int index) {
        return array[iterator + index];
    }

}

class Parm {
    int Id;
    double dSimSwaptionMeanPrice;
    double dSimSwaptionStdError;
    double dStrike;
    double dCompounding;
    double dMaturity;
    double dTenor;
    double dPaymentInterval;
    int iN;
    double dYears;
    int iFactors;
    DoubleArray pdYield;
    DoubleArray[] ppdFactors;

    public Parm(int id,
                double dSimSwaptionMeanPrice,
                double dSimSwaptionStdError,
                double dStrike,
                double dCompounding,
                double dMaturity,
                double dTenor,
                double dPaymentInterval,
                int iN,
                double dYears,
                int iFactors,
                DoubleArray pdYield,
                DoubleArray[] ppdFactors) {
        Id = id;
        this.dSimSwaptionMeanPrice = dSimSwaptionMeanPrice;
        this.dSimSwaptionStdError = dSimSwaptionStdError;
        this.dStrike = dStrike;
        this.dCompounding = dCompounding;
        this.dMaturity = dMaturity;
        this.dTenor = dTenor;
        this.dPaymentInterval = dPaymentInterval;
        this.iN = iN;
        this.dYears = dYears;
        this.iFactors = iFactors;
        this.pdYield = pdYield;
        this.ppdFactors = ppdFactors;
    }

    public Parm() {

    }

    public int getId() {
        return Id;
    }

    public void setId(int id) {
        Id = id;
    }

    public double getdSimSwaptionMeanPrice() {
        return dSimSwaptionMeanPrice;
    }

    public void setdSimSwaptionMeanPrice(double dSimSwaptionMeanPrice) {
        this.dSimSwaptionMeanPrice = dSimSwaptionMeanPrice;
    }

    public double getdSimSwaptionStdError() {
        return dSimSwaptionStdError;
    }

    public void setdSimSwaptionStdError(double dSimSwaptionStdError) {
        this.dSimSwaptionStdError = dSimSwaptionStdError;
    }

    public double getdStrike() {
        return dStrike;
    }

    public void setdStrike(double dStrike) {
        this.dStrike = dStrike;
    }

    public double getdCompounding() {
        return dCompounding;
    }

    public void setdCompounding(double dCompounding) {
        this.dCompounding = dCompounding;
    }

    public double getdMaturity() {
        return dMaturity;
    }

    public void setdMaturity(double dMaturity) {
        this.dMaturity = dMaturity;
    }

    public double getdTenor() {
        return dTenor;
    }

    public void setdTenor(double dTenor) {
        this.dTenor = dTenor;
    }

    public double getdPaymentInterval() {
        return dPaymentInterval;
    }

    public void setdPaymentInterval(double dPaymentInterval) {
        this.dPaymentInterval = dPaymentInterval;
    }

    public int getiN() {
        return iN;
    }

    public void setiN(int iN) {
        this.iN = iN;
    }

    public double getdYears() {
        return dYears;
    }

    public void setdYears(double dYears) {
        this.dYears = dYears;
    }

    public int getiFactors() {
        return iFactors;
    }

    public void setiFactors(int iFactors) {
        this.iFactors = iFactors;
    }

    public DoubleArray getPdYield() {
        return pdYield;
    }

    public void setPdYield(DoubleArray pdYield) {
        this.pdYield = pdYield;
    }

    public DoubleArray[] getPpdFactors() {
        return ppdFactors;
    }

    public void setPpdFactors(DoubleArray[] ppdFactors) {
        this.ppdFactors = ppdFactors;
    }
}
