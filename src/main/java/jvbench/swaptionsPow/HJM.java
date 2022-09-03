package jvbench.swaptionsPow;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;


public class HJM {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_MAX;
    private static final int SPECIES_LENGTH = SPECIES.length();//vector to store discount factors for the rate path along which the swaption
    protected static DoubleArray pdPayoffDiscountFactors;//vector to store rate path along which the swaption payoff will be discounted
    protected static DoubleArray pdDiscountingRatePath;// **** per Trial data **** //
    protected static DoubleArray[] ppdHJMPath;
    protected static DoubleArray pdForward;
    protected static DoubleArray[] ppdDrifts;
    protected static DoubleArray pdTotalDrift;

    static int HJM_Yield_to_Forward(DoubleArray pdForward, int iN, DoubleArray pdYield) {
        //This function computes forward rates from supplied yield rates.

        int iSuccess = 0;
        int i;

        //forward curve computation
        pdForward.set(0, pdYield.get(0));
        for (i = 1; i <= iN - 1; ++i) {
            pdForward.set(i, (i + 1) * pdYield.get(i) - i * pdYield.get(i - 1));    //as per formula
            //printf("pdForward: %f = (%d+1)*%f - %d*%f \n", pdForward[i], i, pdYield[i], i, pdYield[i-1]);
        }
        iSuccess = 1;
        return iSuccess;

    }


    static int HJM_Drifts(DoubleArray pdTotalDrift, DoubleArray[] ppdDrifts, int iN, int iFactors, double dYears, DoubleArray[] ppdFactors) {
        //This function computes drift corrections required for each factor for each maturity based on given factor volatilities

        int iSuccess = 0;
        int i, j, l; //looping variables
        double ddelt = (double) (dYears / iN);
        double dSumVol;

        //computation of factor drifts for shortest maturity
        for (i = 0; i <= iFactors - 1; ++i)
            ppdDrifts[i].set(0, 0.5 * ddelt * (ppdFactors[i].get(0)) * (ppdFactors[i].get(0)));

        //computation of factor drifts for other maturities
        for (i = 0; i <= iFactors - 1; ++i)
            for (j = 1; j <= iN - 2; ++j) {
                ppdDrifts[i].set(j, 0);
                for (l = 0; l <= j - 1; ++l)
                    ppdDrifts[i].set(j, ppdDrifts[i].get(j) - ppdDrifts[i].get(l));
                dSumVol = 0;
                for (l = 0; l <= j; ++l)
                    dSumVol += ppdFactors[i].get(l);
                ppdDrifts[i].set(j, ppdDrifts[i].get(j) + 0.5 * ddelt * (dSumVol) * (dSumVol));
            }

        //computation of total drifts for all maturities
        for (i = 0; i <= iN - 2; ++i) {
            pdTotalDrift.set(i, 0);
            for (j = 0; j <= iFactors - 1; ++j)
                pdTotalDrift.set(i, pdTotalDrift.get(i) + ppdDrifts[j].get(i));
        }

        iSuccess = 1;
        return iSuccess;

    }


    static int HJM_Swaption_Blocking(double[] pdSwaptionPrice,
                                     double dStrike,
                                     double dCompounding,
                                     double dMaturity,
                                     double dTenor,
                                     double dPaymentInterval,
                                     int iN,
                                     int iFactors,
                                     double dYears,
                                     DoubleArray pdYield,
                                     DoubleArray[] ppdFactors,
                                     long[] iRndSeed,
                                     int lTrials,
                                     int BLOCKSIZE,
                                     int tid) {

        int iSuccess = 0;
        int i;
        int b; //block looping variable
        int l; //looping variables

        double ddelt = (double) (dYears / iN);        //ddelt = HJM matrix time-step width. e.g. if dYears = 5yrs and
        //iN = no. of time points = 10, then ddelt = step length = 0.5yrs
        int iFreqRatio = (int) (dPaymentInterval / ddelt + 0.5);    // = ratio of time gap between swap payments and HJM step-width.
        //e.g. dPaymentInterval = 1 year. ddelt = 0.5year. This implies that a swap
        //payment will be made after every 2 HJM time steps.

        double dStrikeCont;        //Strike quoted in continuous compounding convention.
        //As HJM rates are continuous, the K in max(R-K,0) will be dStrikeCont and not dStrike.
        if (dCompounding == 0) {
            dStrikeCont = dStrike;    //by convention, dCompounding = 0 means that the strike entered by user has been quoted
            //using continuous compounding convention
        } else {
            //converting quoted strike to continuously compounded strike
            dStrikeCont = (1 / dCompounding) * Math.log(1 + dStrike * dCompounding);
        }
        //e.g., let k be strike quoted in semi-annual convention. Therefore, 1$ at the end of
        //half a year would earn = (1+k/2). For converting to continuous compounding,
        //(1+0.5*k) = exp(K*0.5)
        // => K = (1/0.5)*ln(1+0.5*k)

        //HJM Framework vectors and matrices
        int iSwapVectorLength;  // Length of the HJM rate path at the time index corresponding to swaption maturity.

        // *******************************
        // ppdHJMPath = dmatrix(0,iN-1,0,iN-1);
//        ppdHJMPath = NrRoutines.dmatrix(0, iN - 1, 0, iN * BLOCKSIZE - 1);    // **** per Trial data **** //
//        pdForward = NrRoutines.dvector(0, iN - 1);
//        ppdDrifts = NrRoutines.dmatrix(0, iFactors - 1, 0, iN - 2);
//        pdTotalDrift = NrRoutines.dvector(0, iN - 2);

        //==================================
        // **** per Trial data **** //
        //payoff will be discounted
        DoubleArray pdSwapRatePath;        //vector to store the rate path along which the swap payments made will be discounted
        DoubleArray pdSwapDiscountFactors;    //vector to store discount factors for the rate path along which the swap
        //payments made will be discounted
        DoubleArray pdSwapPayoffs;        //vector to store swap payoffs


        int iSwapStartTimeIndex;
        int iSwapTimePoints;
        double dSwapVectorYears;

        double dSwaptionPayoff;
        double dDiscSwaptionPayoff;
        double dFixedLegValue;

        // Accumulators
        double dSumSimSwaptionPrice;
        double dSumSquareSimSwaptionPrice;

        // Final returned results
        double dSimSwaptionMeanPrice;
        double dSimSwaptionStdError;

        // *******************************
//        pdPayoffDiscountFactors = NrRoutines.dvector(0, iN * BLOCKSIZE - 1);
//        pdDiscountingRatePath = NrRoutines.dvector(0, iN * BLOCKSIZE - 1);
        // *******************************

        iSwapVectorLength = (int) (iN - dMaturity / ddelt + 0.5);  //This is the length of the HJM rate path at the time index
        //printf("iSwapVectorLength = %d\n",iSwapVectorLength);
        //corresponding to swaption maturity.
        // *******************************
        pdSwapRatePath = NrRoutines.dvector(0, iSwapVectorLength * BLOCKSIZE - 1);
        pdSwapDiscountFactors = NrRoutines.dvector(0, iSwapVectorLength * BLOCKSIZE - 1);
        // *******************************
        pdSwapPayoffs = NrRoutines.dvector(0, iSwapVectorLength - 1);


        iSwapStartTimeIndex = (int) (dMaturity / ddelt + 0.5);  //Swap starts at swaption maturity
        iSwapTimePoints = (int) (dTenor / ddelt + 0.5);      //Total HJM time points corresponding to the swap's tenor
        dSwapVectorYears = (double) (iSwapVectorLength * ddelt);


        //now we store the swap payoffs in the swap payoff vector
        for (i = 0; i <= iSwapVectorLength - 1; ++i)
            pdSwapPayoffs.set(i, 0.0); //initializing to zero
        for (i = iFreqRatio; i <= iSwapTimePoints; i += iFreqRatio) {
            if (i != iSwapTimePoints)
                pdSwapPayoffs.set(i, Math.exp(dStrikeCont * dPaymentInterval) - 1); //the bond pays coupon equal to this amount
            if (i == iSwapTimePoints)
                pdSwapPayoffs.set(i, Math.exp(dStrikeCont * dPaymentInterval)); //at terminal time point, bond pays coupon plus par amount
        }

        //generating forward curve at t=0 from supplied yield curve
        iSuccess = HJM_Yield_to_Forward(pdForward, iN, pdYield);
        if (iSuccess != 1) return iSuccess;

        //computation of drifts from factor volatilities
        iSuccess = HJM_Drifts(pdTotalDrift, ppdDrifts, iN, iFactors, dYears, ppdFactors);
        if (iSuccess != 1) return iSuccess;

        dSumSimSwaptionPrice = 0.0;
        dSumSquareSimSwaptionPrice = 0.0;

        //printf("lTrials = %d\n",lTrials);
        //printf("BLOCKSIZE = %d\n",BLOCKSIZE);
        //Simulations begin:
        for (l = 0; l <= lTrials - 1; l += BLOCKSIZE) {

            int BLOCKSIZE_AUX = Math.min(BLOCKSIZE, (lTrials - l));
            //For each trial a new HJM Path is generated
            iSuccess = HJM_SimPath_Forward_Blocking(ppdHJMPath, iN, iFactors, dYears, pdForward, pdTotalDrift, ppdFactors, iRndSeed, BLOCKSIZE_AUX); /* GC: 51% of the time goes here */
            if (iSuccess != 1) return iSuccess;

            //now we compute the discount factor vector

            for (i = 0; i <= iN - 1; ++i) {
                for (b = 0; b <= BLOCKSIZE_AUX - 1; b++) {
                    pdDiscountingRatePath.set(BLOCKSIZE_AUX * i + b, ppdHJMPath[i].get(0 + b)); // why they do 0 + b??
                }
            }


            if (HJMSecuritiesImpl.useVectorAPI) {
                iSuccess = Discount_Factors_Blocking_Vector(pdPayoffDiscountFactors, iN, dYears, pdDiscountingRatePath, BLOCKSIZE_AUX);
//        iSuccess = Discount_Factors_Blocking_vector(pdPayoffDiscountFactors, iN, dYears, pdDiscountingRatePath, BLOCKSIZE_AUX); /* 15% of the time goes here */

            } else {
                iSuccess = Discount_Factors_Blocking(pdPayoffDiscountFactors, iN, dYears, pdDiscountingRatePath, BLOCKSIZE_AUX); /* 15% of the time goes here */
            }


            if (iSuccess != 1) return iSuccess;

            //now we compute discount factors along the swap path
            for (i = 0; i <= iSwapVectorLength - 1; ++i) {
                for (b = 0; b < BLOCKSIZE_AUX; b++) {
                    pdSwapRatePath.set(i * BLOCKSIZE_AUX + b, ppdHJMPath[iSwapStartTimeIndex].get(i * BLOCKSIZE_AUX + b));
                }
            }

            if (HJMSecuritiesImpl.useVectorAPI) {
                //        iSuccess = Discount_Factors_Blocking_vector(pdSwapDiscountFactors, iSwapVectorLength, dSwapVectorYears, pdSwapRatePath, BLOCKSIZE_AUX);
                iSuccess = Discount_Factors_Blocking_Vector(pdSwapDiscountFactors, iSwapVectorLength, dSwapVectorYears, pdSwapRatePath, BLOCKSIZE_AUX);
            } else {
                iSuccess = Discount_Factors_Blocking(pdSwapDiscountFactors, iSwapVectorLength, dSwapVectorYears, pdSwapRatePath, BLOCKSIZE_AUX);
            }


            if (iSuccess != 1) return iSuccess;

            // ========================
            // Simulation

            if (HJMSecuritiesImpl.useVectorAPI) {
//                unsigned long int gvl = __builtin_epi_vsetvl(BLOCKSIZE_AUX, __epi_e64, __epi_m1);
                int limit = SPECIES.loopBound(BLOCKSIZE_AUX);

                DoubleVector xpdSwapDiscountFactors;
                DoubleVector xpdSwapPayoffs;
                DoubleVector xdFixedLegValue = DoubleVector.zero(SPECIES); //_MM_SET_f64(0.0,gvl);
//                DoubleVector zero = DoubleVector.broadcast(SPECIES, 0.0); //_MM_SET_f64(0.0,gvl);
                DoubleVector oNE = DoubleVector.broadcast(SPECIES, 0.0); //_MM_SET_f64(1.0,gvl);
                DoubleVector xdSumSimSwaptionPrice;
                DoubleVector xdSumSquareSimSwaptionPrice;

                for (b = 0; b < limit; b += SPECIES_LENGTH) {
                    xdFixedLegValue = DoubleVector.zero(SPECIES);
                    for (i = 0; i <= iSwapVectorLength - 1; ++i) {
                        xpdSwapDiscountFactors = DoubleVector.fromArray(SPECIES, pdSwapDiscountFactors.array, (i * BLOCKSIZE_AUX) + pdSwapDiscountFactors.getIterator() + b); // _MM_LOAD_f64(&pdSwapDiscountFactors[i*BLOCKSIZE_AUX],gvl);
                        xpdSwapPayoffs = DoubleVector.broadcast(SPECIES, pdSwapPayoffs.get(i)); // _MM_SET_f64(pdSwapPayoffs[i],gvl);
                        xdFixedLegValue = xpdSwapPayoffs.mul(xpdSwapDiscountFactors).add(xdFixedLegValue); // _MM_MACC_f64(xdFixedLegValue,xpdSwapPayoffs,xpdSwapDiscountFactors,gvl);
                    }

                    xdFixedLegValue = (xdFixedLegValue.sub(1)).max(0); // _MM_MAX_f64(_MM_SUB_f64(xdFixedLegValue,oNE,gvl), zero,gvl);
                    xdFixedLegValue = xdFixedLegValue.mul(DoubleVector.broadcast(SPECIES, pdPayoffDiscountFactors.get(iSwapStartTimeIndex * BLOCKSIZE_AUX))); // _MM_MUL_f64(xdFixedLegValue,_MM_SET_f64(pdPayoffDiscountFactors[iSwapStartTimeIndex*BLOCKSIZE_AUX],gvl),gvl);

                    // ========= end simulation ======================================
                    //xdSumSimSwaptionPrice       = _MM_LOAD_f64(&dSumSimSwaptionPrice,1);
                    //xdSumSquareSimSwaptionPrice = _MM_LOAD_f64(&dSumSquareSimSwaptionPrice,1);

//                xdSumSimSwaptionPrice       = DoubleVector.broadcast(SPECIES, dSumSimSwaptionPrice); // _MM_SET_f64(dSumSimSwaptionPrice,1);
//                xdSumSquareSimSwaptionPrice = DoubleVector.broadcast(SPECIES, dSumSquareSimSwaptionPrice); // _MM_SET_f64(dSumSquareSimSwaptionPrice,1);

                    // accumulate into the aggregating variables =====================
                    dSumSimSwaptionPrice += xdFixedLegValue.reduceLanes(VectorOperators.ADD); //_MM_REDSUM_f64(xdFixedLegValue,xdSumSimSwaptionPrice,gvl);
//                    dSumSquareSimSwaptionPrice += xdFixedLegValue.mul(xdFixedLegValue).reduceLanes(VectorOperators.ADD);// _MM_REDSUM_f64(_MM_MUL_f64(xdFixedLegValue,xdFixedLegValue,gvl),xdSumSquareSimSwaptionPrice,gvl);
                    dSumSquareSimSwaptionPrice += xdFixedLegValue.pow(2).reduceLanes(VectorOperators.ADD);// _MM_REDSUM_f64(_MM_MUL_f64(xdFixedLegValue,xdFixedLegValue,gvl),xdSumSquareSimSwaptionPrice,gvl);


//                _MM_STORE_f64(&dSumSimSwaptionPrice,xdSumSimSwaptionPrice,1);
//                _MM_STORE_f64(&dSumSquareSimSwaptionPrice,xdSumSquareSimSwaptionPrice,1);
//                FENCE();
                }
            } else {
                for (b = 0; b < BLOCKSIZE_AUX; b++) {
                    dFixedLegValue = 0.0;
                    for (i = 0; i <= iSwapVectorLength - 1; ++i) {
                        dFixedLegValue += pdSwapPayoffs.get(i) * pdSwapDiscountFactors.get(i * BLOCKSIZE_AUX + b);
                    }
                    dSwaptionPayoff = Math.max(dFixedLegValue - 1.0, 0);
                    dDiscSwaptionPayoff = dSwaptionPayoff * pdPayoffDiscountFactors.get(iSwapStartTimeIndex * BLOCKSIZE_AUX + b);
                    // ========= end simulation ======================================
                    // accumulate into the aggregating variables =====================
                    dSumSimSwaptionPrice += dDiscSwaptionPayoff;
                    dSumSquareSimSwaptionPrice += dDiscSwaptionPayoff * dDiscSwaptionPayoff;
                } // END BLOCK simulation
            }
        }

        // Simulation Results Stored
        dSimSwaptionMeanPrice = dSumSimSwaptionPrice / lTrials;
        dSimSwaptionStdError = Math.sqrt((dSumSquareSimSwaptionPrice - dSumSimSwaptionPrice * dSumSimSwaptionPrice / lTrials) / (lTrials - 1.0)) / Math.sqrt((double) lTrials);

        //results returned
        pdSwaptionPrice[0] = dSimSwaptionMeanPrice;
        pdSwaptionPrice[1] = dSimSwaptionStdError;

        iSuccess = 1;
        return iSuccess;
    }

    private static int Discount_Factors_Blocking_Vector(DoubleArray pdDiscountFactors, int iN, double dYears, DoubleArray pdRatePath, int BLOCKSIZE) {
        int i, j, b;                //looping variables
        int iSuccess;            //return variable

        double ddelt;            //HJM time-step length
        ddelt = (double) (dYears / iN);

        // unsigned long int gvl = __builtin_epi_vsetvl(BLOCKSIZE, __epi_e64, __epi_m1);
        int limit = SPECIES.loopBound(BLOCKSIZE);
        DoubleVector xDdelt;
        DoubleVector xpdRatePath;


        DoubleArray pdexpRes;
        pdexpRes = NrRoutines.dvector(0, (iN - 1) * BLOCKSIZE - 1);

        //precompute the exponientials
        for (j = 0; j <= (iN - 1) * BLOCKSIZE - 1; j += BLOCKSIZE) {
            for (b = 0; b < limit; b += SPECIES_LENGTH) {
                xpdRatePath = DoubleVector.fromArray(SPECIES, pdRatePath.array, j + pdRatePath.getIterator() + b); // _MM_LOAD_f64(&pdRatePath[j],gvl);
                xpdRatePath = xpdRatePath.mul(ddelt); // _MM_MUL_f64(xpdRatePath,xDdelt,gvl);
                xpdRatePath = xpdRatePath.lanewise(VectorOperators.NEG).lanewise(VectorOperators.EXP); // _MM_EXP_f64(_MM_VFSGNJN_f64(xpdRatePath,xpdRatePath,gvl) ,gvl);
                xpdRatePath.intoArray(pdexpRes.array, j + pdexpRes.getIterator() + b); //_MM_STORE_f64(&pdexpRes[j],xpdRatePath,gvl);
            }
        }
        // TODO scalar part

        //initializing the discount factor vector
        for (i = 0; i < (iN) * BLOCKSIZE; i += BLOCKSIZE) {
            for (b = 0; b < limit; b += SPECIES_LENGTH) {
                DoubleVector.broadcast(SPECIES, 1.0).intoArray(pdDiscountFactors.array, i + pdDiscountFactors.getIterator() + b);// _MM_STORE_f64(&pdDiscountFactors[i],_MM_SET_f64(1.0,gvl),gvl);
            }
        }
        // TODO scalar part

        for (i = 1; i <= iN - 1; ++i) {
            for (j = 0; j <= i - 1; ++j) {
                for (b = 0; b < limit; b += SPECIES_LENGTH) {

                    xpdRatePath = DoubleVector.fromArray(SPECIES, pdDiscountFactors.array, (i * BLOCKSIZE) + pdDiscountFactors.getIterator() + b); // _MM_LOAD_f64(&pdDiscountFactors[i*BLOCKSIZE],gvl);
                    xDdelt = DoubleVector.fromArray(SPECIES, pdexpRes.array, (j * BLOCKSIZE) + pdexpRes.getIterator() + b); // _MM_LOAD_f64(&pdexpRes[j*BLOCKSIZE],gvl);
                    xpdRatePath = xpdRatePath.mul(xDdelt); // _MM_MUL_f64(xpdRatePath,xDdelt,gvl);
                    xpdRatePath.intoArray(pdDiscountFactors.array, (i * BLOCKSIZE) + pdDiscountFactors.getIterator() + b); // _MM_STORE_f64(&pdDiscountFactors[i*BLOCKSIZE],xpdRatePath,gvl);
                }
            }
        }
//        FENCE();
//        free_dvector(pdexpRes, 0,(iN-1)*BLOCKSIZE-1);
        iSuccess = 1;
        return iSuccess;
    }

    static int Discount_Factors_Blocking(DoubleArray pdDiscountFactors, int iN, double dYears, DoubleArray pdRatePath, int BLOCKSIZE) {

        int i, j, b;                //looping variables
        int iSuccess;            //return variable

        double ddelt;            //HJM time-step length
        ddelt = (double) (dYears / iN);
        DoubleArray pdexpRes;

        pdexpRes = NrRoutines.dvector(0, (iN - 1) * BLOCKSIZE - 1);

        //precompute the exponientials
        for (j = 0; j <= (iN - 1) * BLOCKSIZE - 1; ++j) {
            pdexpRes.set(j, pdexpRes.get(j) - pdRatePath.get(j) * ddelt);
        }
        for (j = 0; j <= (iN - 1) * BLOCKSIZE - 1; ++j) {
            pdexpRes.set(j, Math.exp(pdexpRes.get(j)));
        }


        //initializing the discount factor vector
        for (i = 0; i < (iN) * BLOCKSIZE; ++i)
            pdDiscountFactors.set(i, 1.0);

        for (i = 1; i <= iN - 1; ++i) {
            //printf("\nVisiting timestep %d : ",i);
            for (b = 0; b < BLOCKSIZE; b++) {
                //printf("\n");
                for (j = 0; j <= i - 1; ++j) {
                    pdDiscountFactors.set(i * BLOCKSIZE + b, pdDiscountFactors.get(i * BLOCKSIZE + b) * pdexpRes.get(j * BLOCKSIZE + b));
                    //printf("(%f) ",pdexpRes[j*BLOCKSIZE + b]);
                }
            } // end Block loop
        }

        iSuccess = 1;
        return iSuccess;
    }

    static int HJM_SimPath_Forward_Blocking(DoubleArray[] ppdHJMPath,
                                            int iN,
                                            int iFactors,
                                            double dYears,
                                            DoubleArray pdForward,
                                            DoubleArray pdTotalDrift,
                                            DoubleArray[] ppdFactors,
                                            long[] lRndSeed,
                                            int BLOCKSIZE) {
        int iSuccess = 0;
        int i, j, l; //looping variables
        DoubleArray[] pdZ; //vector to store random normals
        DoubleArray[] randZ; //vector to store random normals
        double dTotalShock; //total shock by which the forward curve is hit at (t, T-t)
        double ddelt, sqrt_ddelt; //length of time steps

        ddelt = (double) (dYears / iN);
        sqrt_ddelt = Math.sqrt(ddelt);

        pdZ = NrRoutines.dmatrix(0, iFactors - 1, 0, iN * BLOCKSIZE - 1); //assigning memory
        randZ = NrRoutines.dmatrix(0, iFactors - 1, 0, iN * BLOCKSIZE - 1); //assigning memory

        // =====================================================
        // t=0 forward curve stored iN first row of ppdHJMPath
        // At time step 0: insert expected drift
        // rest reset to 0

        if (HJMSecuritiesImpl.useVectorAPI) {
            DoubleVector xZero = DoubleVector.zero(SPECIES);
            //for(int b=0; b<BLOCKSIZE; b++){
            int limit = SPECIES.loopBound(BLOCKSIZE);
            for (j = 0; j <= iN - 1; j++) {
//                _MM_STORE_f64(&ppdHJMPath[0][BLOCKSIZE*j],_MM_SET_f64(pdForward[j],gvl),gvl);
                for (int b = 0; b < limit; b += SPECIES_LENGTH) {
                    // TODO SCALAR
                    DoubleVector.broadcast(SPECIES, pdForward.get(j)).intoArray(ppdHJMPath[0].array, (BLOCKSIZE * j) + b + ppdHJMPath[0].getIterator());

                }
                for (i = 1; i <= iN - 1; ++i) {
//                    _MM_STORE_f64(&ppdHJMPath[i][BLOCKSIZE*j],xZero,gvl);
                    for (int b = 0; b < limit; b += SPECIES_LENGTH) {
                        // TODO SCALAR
                        xZero.intoArray(ppdHJMPath[i].array, (BLOCKSIZE * j) + b + ppdHJMPath[i].getIterator());
                    }

                } //initializing HJMPath to zero
            }
            //}

        } else {
            // TODO  INVERT LOOP for cache hit
            for (int b = 0; b < BLOCKSIZE; b++) {
                for (j = 0; j <= iN - 1; j++) {
                    ppdHJMPath[0].set(BLOCKSIZE * j + b, pdForward.get(j));

                    for (i = 1; i <= iN - 1; ++i) {
                        ppdHJMPath[i].set(BLOCKSIZE * j + b, 0);
                    } //initializing HJMPath to zero
                }
            }
        }

        if (HJMSecuritiesImpl.useVectorAPI) {
            RanUnif.ranUnifVector(lRndSeed, iFactors, iN, BLOCKSIZE, randZ);
            // TODO REM SCALAR PART
        } else {

            // =====================================================
            // sequentially generating random numbers

            for (int b = 0; b < BLOCKSIZE; b++) {
                for (int s = 0; s < 1; s++) {
                    for (j = 1; j <= iN - 1; ++j) {
                        for (l = 0; l <= iFactors - 1; ++l) {
                            //compute random number in exact same sequence
                            randZ[l].set(BLOCKSIZE * j + b + s, RanUnif.ranUnif(lRndSeed));  /* 10% of the total executition time */
                        }
                    }
                }
            }
        }
        serialB(pdZ, randZ, BLOCKSIZE, iN, iFactors);


        if (HJMSecuritiesImpl.useVectorAPI) {
            double pdDriftxddelt;
            // =====================================================
            // Generation of HJM Path1 Vector
//            gvl = __builtin_epi_vsetvl(BLOCKSIZE, __epi_e64, __epi_m1);
            int limit = SPECIES.loopBound(BLOCKSIZE);
            DoubleVector xdTotalShock;


            //for(int b=0; b<BLOCKSIZE; b++){ // b is the blocks
            for (j = 1; j <= iN - 1; ++j) {// j is the timestep

                for (l = 0; l <= iN - (j + 1); ++l) { // l is the future steps
                    // _MM_SET_f64(0.0,gvl);
                    pdDriftxddelt = pdTotalDrift.get(l) * ddelt;

//                    _MM_STORE_f64(&(ppdHJMPath[j][BLOCKSIZE*l]),
//                            _MM_ADD_f64(_MM_LOAD_f64(&ppdHJMPath[j-1][BLOCKSIZE*(l+1)],gvl),
//                                _MM_ADD_f64(_MM_SET_f64(pdDriftxddelt,gvl),
//                                    _MM_MUL_f64(_MM_SET_f64(sqrt_ddelt,gvl),xdTotalShock,gvl),gvl),gvl),gvl);
                    for (int b = 0; b < limit; b += SPECIES_LENGTH) {
                        xdTotalShock = DoubleVector.zero(SPECIES);
                        for (i = 0; i <= iFactors - 1; ++i) {// i steps through the stochastic factors
//                        xdTotalShock = _MM_ADD_f64(xdTotalShock, _MM_MUL_f64(_MM_SET_f64(ppdFactors[i][l],gvl), _MM_LOAD_f64(&pdZ[i][BLOCKSIZE*j],gvl),gvl),gvl);
                            // TODO find a way to reduce lanewise and not across vectors
                            xdTotalShock = xdTotalShock.add(ppdFactors[i].get(l)).mul(DoubleVector.fromArray(SPECIES, pdZ[i].array, (BLOCKSIZE * j) + b + pdZ[i].getIterator()));
                        }
                        DoubleVector.fromArray(SPECIES, ppdHJMPath[j - 1].array, (BLOCKSIZE * (l + 1) + b + ppdHJMPath[j - 1].getIterator()))
                                .add((DoubleVector.broadcast(SPECIES, pdDriftxddelt)
                                        .add((DoubleVector.broadcast(SPECIES, sqrt_ddelt)
                                                .mul(xdTotalShock))))).intoArray(ppdHJMPath[j].array, (BLOCKSIZE * l) + b + ppdHJMPath[j].getIterator());
                        //as per formula
                    }
                }
            }
        } else {

            // =====================================================
            // Generation of HJM Path1
            for (int b = 0; b < BLOCKSIZE; b++) { // b is the blocks
                for (j = 1; j <= iN - 1; ++j) {// j is the timestep

                    for (l = 0; l <= iN - (j + 1); ++l) { // l is the future steps
                        dTotalShock = 0;

                        for (i = 0; i <= iFactors - 1; ++i) {// i steps through the stochastic factors
                            dTotalShock += ppdFactors[i].get(l) * pdZ[i].get(BLOCKSIZE * j + b);
                        }

                        ppdHJMPath[j].set(BLOCKSIZE * l + b, ppdHJMPath[j - 1].get(BLOCKSIZE * (l + 1) + b) + pdTotalDrift.get(l) * ddelt + sqrt_ddelt * dTotalShock);
                        //as per formula
                    }
                }
            } // end Blocks
            // -----------------------------------------------------
        }
        iSuccess = 1;
        return iSuccess;
    }

    static void serialB(DoubleArray[] pdZ, DoubleArray[] randZ, int BLOCKSIZE, int iN, int iFactors) {
        int limit = SPECIES.loopBound(BLOCKSIZE);
        if (HJMSecuritiesImpl.useVectorAPI) {
            for (int l = 0; l <= iFactors - 1; ++l) {
                for (int j = 1; j <= iN - 1; ++j) {
                    for (int b = 0; b < limit; b += SPECIES_LENGTH) {
//                    unsigned long int gvl = __builtin_epi_vsetvl(BLOCKSIZE, __epi_e64, __epi_m1);
                        CumNormalInv.cumNormalInvVector(randZ[l], BLOCKSIZE * j + randZ[l].getIterator() + b, pdZ[l], BLOCKSIZE * j + pdZ[l].getIterator() + b);
                        //FENCE();
                    }
                }
            }
        } else {

            for (int l = 0; l <= iFactors - 1; ++l) {
                for (int j = 1; j <= iN - 1; ++j) {
                    for (int b = 0; b < BLOCKSIZE; b++) {
                        pdZ[l].set(BLOCKSIZE * j + b, CumNormalInv.cumNormalInv(randZ[l].get(BLOCKSIZE * j + b)));  /* 18% of the total executition time */
//    			printf("CumNormalInv output: %f, input: %f\n",pdZ[l][BLOCKSIZE*j + b],randZ[l][BLOCKSIZE*j + b]);
                    }
                }
            }
        }
    }


}