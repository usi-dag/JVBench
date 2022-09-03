package jvbench.blackscholesMergeWihtoutBlend;

import jdk.incubator.vector.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class Blackscholes {

    static private final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_MAX;
    static private final int SPECIES_LENGTH = SPECIES.length();

    static final int nThreads = 1;
    static final String outputFileName = "";
    static int numOptions;
    static OptionData[] data;
    static float [] prices;

    private static final int PAD =  256;
    private static final int LINESIZE =  64;
    private static final int NUM_RUNS = 100;

    private static float [] buffer;
    private static float [] sptprice;
    private static float [] strike;
    private static float [] rate;
    private static float [] volatility;
    private static float [] otime;


    private static float[] buffer2;
    static int [] otype;

    private static final float inv_sqrt_2xPI = (float) 0.39894228040143270286;

    public static void init(String inputFileName) {
        try (BufferedReader reader = new BufferedReader(new FileReader(inputFileName))) {

            numOptions = Integer.parseInt(reader.readLine());

            buffer = new float[5 * numOptions + PAD];
            sptprice = new float[5 * numOptions + PAD];
            strike = new float[5 * numOptions + PAD];
            rate = new float[5 * numOptions + PAD];
            volatility = new float[5 * numOptions + PAD];
            otime = new float[5 * numOptions + PAD];

            data = new OptionData[numOptions];
            prices = new float[numOptions];

            buffer2 = new float[5 * numOptions];
            otype = new int[5 * numOptions];

            for (int loopNum = 0; loopNum < numOptions; loopNum++) {
                String line = reader.readLine();
                String [] options = line.split(" ");
                if (options.length != 9) {
                    throw new IllegalArgumentException("Invalid number of options: " + options.length + " should be 9");
                }
                data[loopNum] = new OptionData(
                        Float.parseFloat(options[0]),
                        Float.parseFloat(options[1]),
                        Float.parseFloat(options[2]),
                        Float.parseFloat(options[3]),
                        Float.parseFloat(options[4]),
                        Float.parseFloat(options[5]),
                        options[6].charAt(0),
                        Float.parseFloat(options[7]),
                        Float.parseFloat(options[8])
                );

            }

        } catch (IOException e) {
            System.err.println("ERROR: Unable to read file " + inputFileName + "\n" + e.getMessage());
        }

        for (int i = 0; i < numOptions; i++) {
            otype[i] = (data[i].getOptionType() == 'P') ? 1 : 0;
            sptprice[i] = data[i].getS();
            strike[i] = data[i].getStrike();
            rate[i] = data[i].getR();
            volatility[i] = data[i].getV();
            otime[i] = data[i].getT();
        }
    }

    public static float [] getPrices() {
        return prices;
    }

    public static OptionData[] getData() {return data;}

    public static void scalar() {
        for (int j=0; j < NUM_RUNS; j++) {
            for (int i=0; i<numOptions; i++) {
                float price =  blkSchlsEqEuroNoDiv(sptprice[i], strike[i], rate[i], volatility[i], otime[i], otype[i], 0);
                prices[i] = price;
            }
        }
    }

    public static void vector() {
        for (int j=0; j < NUM_RUNS; j++) {
            int limit = SPECIES.loopBound(numOptions);
            int i;
            for (i=0; i<limit; i += SPECIES_LENGTH) {
                blkSchlsEqEuroNoDivVector(i);
            }

            for (; i<numOptions; i++) {
                float price =  blkSchlsEqEuroNoDiv(sptprice[i], strike[i], rate[i], volatility[i], otime[i], otype[i], 0);
                prices[i] = price;
            }
        }
    }

    private static void blkSchlsEqEuroNoDivVector(int i) {

        FloatVector xStockPrice;
        FloatVector xStrikePrice;
        FloatVector xRiskFreeRate;
        FloatVector xVolatility;
        FloatVector xTime;
        FloatVector xSqrtTime;

        FloatVector xLogTerm;
        FloatVector xD1, xD2;
        FloatVector xPowerTerm;
        FloatVector xDen;

        FloatVector xRatexTime;
        FloatVector xFutureValueX;

        VectorMask<Integer> xMask;
        IntVector xOtype;
        IntVector  xZero;

        FloatVector xOptionPrice;
        FloatVector xOptionPrice1;
        FloatVector xOptionPrice2;
        FloatVector xfXd1;
        FloatVector xfXd2;

        xStockPrice = FloatVector.fromArray(SPECIES, sptprice, i);
        xStrikePrice = FloatVector.fromArray(SPECIES, strike, i);
        xRiskFreeRate = FloatVector.fromArray(SPECIES, rate, i);
        xVolatility = FloatVector.fromArray(SPECIES, volatility, i);

        xTime = FloatVector.fromArray(SPECIES, otime, i);
        xSqrtTime = xTime.sqrt();

        xLogTerm = (xStockPrice.div(xStrikePrice)).lanewise(VectorOperators.LOG);

        xPowerTerm = xVolatility.mul(xVolatility);
        xPowerTerm = xPowerTerm.mul(0.5f);

        xD1 = xRiskFreeRate.add(xPowerTerm);
        xD1 = xD1.mul(xTime);
        xD1 = xD1.add(xLogTerm);


        xDen = xVolatility.mul(xSqrtTime);
        xD1 = xD1.div(xDen);
        xD2 = xD1.sub(xDen);

        xfXd1 = cndfSIMD(xD1);
        xfXd2 = cndfSIMD(xD2);


        xRatexTime = xRiskFreeRate.mul(xTime);
        xRatexTime = xRatexTime.lanewise(VectorOperators.NEG);
        xFutureValueX = xRatexTime.lanewise(VectorOperators.EXP);
        xStrikePrice = FloatVector.fromArray(SPECIES, strike, i);
        xFutureValueX = xFutureValueX.mul(xStrikePrice);


        xOtype = IntVector.fromArray(IntVector.SPECIES_MAX, otype, i);
        xZero = IntVector.zero(IntVector.SPECIES_MAX);
        xMask = xZero.eq(xOtype);


        xOptionPrice1 = xStockPrice.mul(xfXd1);
        xOptionPrice1 = xOptionPrice1.sub(xFutureValueX.mul(xfXd2));

        xfXd1 = (FloatVector.broadcast(SPECIES, 1.0f).sub(xfXd1));
        xfXd2 = (FloatVector.broadcast(SPECIES, 1.0f).sub(xfXd2));

        xOptionPrice2 = xFutureValueX.mul(xfXd2);
        xOptionPrice2 = xOptionPrice2.sub(xStockPrice.mul(xfXd1));
        VectorMask<Float> vMask = xMask.cast(FloatVector.SPECIES_MAX);
//        xOptionPrice = xOptionPrice2.blend(xOptionPrice1, xMask.cast(FloatVector.SPECIES_MAX));
        xOptionPrice = FloatVector.zero(SPECIES).add(xOptionPrice2, vMask.not()).add(xOptionPrice1, vMask);

        xOptionPrice.intoArray(prices, i);
    }

    private static FloatVector cndfSIMD(FloatVector xInput) {
        FloatVector xNPrimeofX;
        FloatVector xK2;
        FloatVector xK2_2;
        FloatVector xK2_3;
        FloatVector xK2_4;
        FloatVector xK2_5;
        FloatVector xLocal;
        FloatVector xLocal_1;
        FloatVector xLocal_2;
        FloatVector xLocal_3;

        VectorMask<Float> xMask;

        FloatVector expValues;

        FloatVector xOne = FloatVector.broadcast(SPECIES, 1);


        xMask = xInput.lt(0.0f);


        xInput = xInput.lanewise(VectorOperators.NEG, xMask);

        expValues = xInput.mul(xInput);
        expValues = expValues.mul(-0.5f);
        expValues = expValues.lanewise(VectorOperators.EXP);

        xNPrimeofX = expValues;
        xNPrimeofX = xNPrimeofX.mul(inv_sqrt_2xPI);


        xK2 = xInput.mul(0.2316419f).add(1.0f);
        xK2 = xOne.div(xK2);

        xK2_2 = xK2.mul(xK2);
        xK2_3 = xK2_2.mul(xK2);
        xK2_4 = xK2_3.mul(xK2);
        xK2_5 = xK2_4.mul(xK2);

        xLocal_1 = xK2.mul(0.319381530f);
        xLocal_2 = xK2_2.mul(-0.356563782f);
        xLocal_3 = xK2_3.mul(1.781477937f);
        xLocal_2 = xLocal_2.add(xLocal_3);
        xLocal_3 = xK2_4.mul(-1.821255978f);
        xLocal_2 = xLocal_2.add(xLocal_3);
        xLocal_3 = xK2_5.mul(1.330274429f);
        xLocal_2 = xLocal_2.add(xLocal_3);


        xLocal_1 = xLocal_2.add(xLocal_1);

        xLocal = xLocal_1.mul(xNPrimeofX);
        xLocal = xOne.sub(xLocal);

//        xLocal = xLocal.blend(xOne.sub(xLocal), xMask);
        xLocal = FloatVector.zero(SPECIES).add(xOne.sub(xLocal), xMask).add(xLocal, xMask.not());

        return xLocal;




    }

    private static float blkSchlsEqEuroNoDiv(float sptprice, float strike, float rate, float volatility, float time, int otype, float timet) {
        float OptionPrice;

        // local private working variables for the calculation
        float xStockPrice;
        float xStrikePrice;
        float xRiskFreeRate;
        float xVolatility;
        float xTime;
        float xSqrtTime;

        float logValues;
        float xLogTerm;
        float xD1;
        float xD2;
        float xPowerTerm;
        float xDen;
        float d1;
        float d2;
        float FutureValueX;
        float NofXd1;
        float NofXd2;
        float NegNofXd1;
        float NegNofXd2;

        xStockPrice = sptprice;
        xStrikePrice = strike;
        xRiskFreeRate = rate;
        xVolatility = volatility;

        xTime = time;
        xSqrtTime = (float) Math.sqrt(xTime);

        logValues = (float) Math.log( sptprice / strike );

        xLogTerm = logValues;

        xPowerTerm = xVolatility * xVolatility;
        xPowerTerm = (xPowerTerm * 0.5f);


        xD1 = xRiskFreeRate + xPowerTerm;
        xD1 = xD1 * xTime;
        xD1 = xD1 + xLogTerm;

        xDen = xVolatility * xSqrtTime;
        xD1 = xD1 / xDen;
        xD2 = xD1 -  xDen;

        d1 = xD1;
        d2 = xD2;

        NofXd1 = CNDF( d1 );
        NofXd2 = CNDF( d2 );

        FutureValueX = (float) (strike * ( Math.exp( -(rate)*(time) ) ));
        if (otype == 0) {
            OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
        } else {
            NegNofXd1 = (float) (1.0 - NofXd1);
            NegNofXd2 = (float) (1.0 - NofXd2);
            OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
        }

        return OptionPrice;
    }

    private static float CNDF(float inputX) {
        int sign;

        float outputX;
        float xInput;
        float xNPrimeofX;
        float expValues;
        float xK2;
        float xK2_2, xK2_3;
        float xK2_4, xK2_5;
        float xLocal, xLocal_1;
        float xLocal_2, xLocal_3;

        // Check for negative value of InputX
        if (inputX < 0.0f) {
            inputX = -inputX;
            sign = 1;
        } else
            sign = 0;

        xInput = inputX;

        // Compute NPrimeX term common to both four & six decimal accuracy calcs
        expValues = (float) Math.exp(-0.5f * inputX * inputX);
        xNPrimeofX = expValues;
        xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

        xK2 =  (0.2316419f * xInput);
        xK2 =  (1.0f + xK2);
        xK2 =  (1.0f / xK2);
        xK2_2 = xK2 * xK2;
        xK2_3 = xK2_2 * xK2;
        xK2_4 = xK2_3 * xK2;
        xK2_5 = xK2_4 * xK2;

        xLocal_1 = (xK2 * 0.319381530f);
        xLocal_2 = (xK2_2 * (-0.356563782f));
        xLocal_3 = (xK2_3 * 1.781477937f);
        xLocal_2 = xLocal_2 + xLocal_3;
        xLocal_3 = (xK2_4 * (-1.821255978f));
        xLocal_2 = xLocal_2 + xLocal_3;
        xLocal_3 = (xK2_5 * 1.330274429f);
        xLocal_2 = xLocal_2 + xLocal_3;

        xLocal_1 = xLocal_2 + xLocal_1;
        xLocal   = xLocal_1 * xNPrimeofX;
        xLocal   = (1.0f - xLocal);

        outputX  = xLocal;

        if (sign == 1) {
            outputX = (1.0f - outputX);
        }

        return outputX;
    }
}
