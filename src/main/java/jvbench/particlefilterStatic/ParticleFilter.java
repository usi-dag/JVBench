package jvbench.particlefilterStatic;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

public class ParticleFilter {

    private static final  VectorSpecies<Double> DOUBLE_SPECIES = DoubleVector.SPECIES_MAX;
    public static final int DOUBLE_SPECIES_LENGTH = DOUBLE_SPECIES.length();
//    private static  VectorSpecies<Integer> INT_SPECIES = IntVector.SPECIES_MAX;
//    private static final int INT_SPECIES_LENGTH = INT_SPECIES.length();
    /**
     * @value M value for Linear Congruential Generator (LCG)
     */
    private static final long M = Integer.MAX_VALUE;

    /**
     * @value A value for LCG
     */
    private static final int A = 1103515245;

    /**
     * @value C value for LCG
     */
    private static final int C = 12345;

    private static final double PI = 3.1415926535897932;
    static double[] resXV;
    static double[] resYV;
    static double[] resWV;

    static double[] resX;
    static double[] resY;
    static double[] resW;

    static final boolean [] mask = new boolean[DOUBLE_SPECIES_LENGTH];
    static final VectorMask<Double> xMask = VectorMask.fromArray(DOUBLE_SPECIES, mask, 0);


    /**
     * Takes in a double and returns an integer that approximates to that double
     * @return if the mantissa < .5 => return value < input value; else return value > input value
     */
    private static double roundDouble(double value) {
        int newValue = (int) value;
        if (value - newValue < .5)  return newValue;
        else return newValue + 1;
    }

    /**
     * Set values of the 3D array to a newValue if that value is equal to the testValue
     * @param testValue The value to be replaced
     * @param newValue The value to replace testValue with
     * @param array3D The image vector
     * @param dimX The x dimension of the frame
     * @param dimY The y dimension of the frame
     * @param dimZ The number of frames
     */
    private static void setIf(int testValue, int newValue, int [] array3D, int dimX, int dimY, int dimZ) {
        int x, y, z;
        for (x = 0; x < dimX; x++) {
            for (y = 0; y < dimY; y++) {
                for (z = 0; z < dimZ; z++) {
                    if (array3D[x * dimY * dimZ+y * dimZ+z] == testValue) {
                        array3D[x * dimY * dimZ+y * dimZ+z] = newValue;
                    }
                }
            }
        }
    }

    /**
     * Generates a uniformly distributed random number using the provided seed and GCC's settings for the Linear Congruential Generator (LCG)
     * @see <a href="http://en.wikipedia.org/wiki/Linear_congruential_generator">http://en.wikipedia.org/wiki/Linear_congruential_generator</a>
     * This function is thread-safe
     * @param seed The seed array
     * @param index The specific index of the seed to be advanced
     * @return a uniformly distributed number [0, 1)
     */
    private static double randu(double[] seed, int index) {
        int num = (A* (int) seed[index] + C); // TODO cast to int to compare result with C version
        seed[index] = (num % M);
        return Math.abs(seed[index]/((double) M));
    }


    private static DoubleVector randuVector(double[] seed, int index, double [] result, double [] num) {

        for (int x = index; x < index+DOUBLE_SPECIES_LENGTH; x++) {
            num[x-index] = (A* (int) seed[x] + C);
            seed[x] = (num[x-index] % M);
            result[x-index] = Math.abs(seed[x]/( M));
//            System.out.printf("%d: fabs(seed[%d] / (M)) = fabs(%f / %d) = %f\n", x, x, seed[x], M, Math.abs(seed[x]/M));
        }
        DoubleVector xResult;

        xResult = DoubleVector.fromArray(DOUBLE_SPECIES, result, 0);
//        VectorMask<Double> gvl = DOUBLE_SPECIES.indexInRange(index, seed.length);
//        System.out.println(index + " - " +gvl);
//        DoubleVector xSeed = DoubleVector.fromArray(DOUBLE_SPECIES, seed, index, gvl);
//
//        xSeed = xSeed.mul(A);
//        xSeed = xSeed.add(C);
//        xSeed = xSeed.sub((xSeed.div(M)).mul(M)); // TODO how to roundDown intermadiate result
//        xSeed.intoArray(seed, index, gvl);
//
//        DoubleVector xResult;
        // cast IntVector (convert) to Vector, then cast Vector to DoubleVector
        // otherswise only able to apply division with vector and not double.
//        xResult = (xSeed.convert(VectorOperators.I2D, 0).reinterpretAsDoubles()).div((double) M).abs();
//        xResult = (xSeed).div((double) M).abs();
//        System.out.println("xResult: vec " + xResult);
        return xResult;
    }

    /**
     * Generates a normally distributed random number using the Box-Muller transformation
     * This function is thread-safe
     * @param seed The seed array
     * @param index The specific index of the seed to be advanced
     * @return a double representing random number generated using the Box-Muller algorithm
     * @see <a href="http://en.wikipedia.org/wiki/Normal_distribution">section computing value for normal random distribution</a>
     */
    private static double randn(double[] seed, int index) {
        /*Box-Muller algorithm*/
        double u = randu(seed, index);
        double v = randu(seed, index);
        double cosine = Math.cos(2*PI*v);
        double rt = -2*Math.log(u);
        return Math.sqrt(rt)*cosine;
    }


    private static DoubleVector randnVector(double[] seed, int index, double [] randuVectorResult, double [] randuVectorNum) {
        DoubleVector xU = randuVector(seed, index, randuVectorResult, randuVectorNum);
        DoubleVector xV = randuVector(seed, index, randuVectorResult, randuVectorNum);
        DoubleVector xCosine;
        DoubleVector xRt;

        xV = xV.mul(Math.PI*2.0);
        xCosine = xV.lanewise(VectorOperators.COS);
//        System.out.println("XU before log Numerical: " + xU);
        xU = xU.lanewise(VectorOperators.LOG);
//        System.out.println("XU after log Numerical: " + xU);
        xRt = xU.mul(-2.0);

        return (xRt.sqrt()).mul(xCosine);
    }

    /**
     * Sets values of 3D matrix using randomly generated numbers from a normal distribution
     * @param array3D The video to be modified
     * @param dimX The x dimension of the frame
     * @param dimY The y dimension of the frame
     * @param dimZ The number of frames
     * @param seed The seed array
     */

    static void addNoise(int[] array3D, int dimX, int dimY, int dimZ, double[] seed){
        int x, y, z;
        for(x = 0; x < dimX; x++){
            for(y = 0; y < dimY; y++){
                for(z = 0; z < dimZ; z++){
                    array3D[x * dimY  * dimZ + y * dimZ + z] = array3D[x * dimY * dimZ + y * dimZ + z] + (int)(5*randn(seed, 0));
                }
            }
        }
    }

    /**
     * Fills a radius x radius matrix representing the disk
     * @param disk The pointer to the disk to be made
     * @param radius  The radius of the disk to be made
     */
    static void strelDisk(int[] disk, int radius)
    {
        int diameter = radius*2 - 1;
        int x, y;
        for(x = 0; x < diameter; x++){
            for(y = 0; y < diameter; y++){
                double distance = Math.sqrt(Math.pow((double)(x-radius+1),2) + Math.pow((double)(y-radius+1),2));
                if(distance < radius)
                    disk[x*diameter + y] = 1;
            }
        }
    }

    /**
     * Dilates the provided video
     * @param matrix The video to be dilated
     * @param posX The x location of the pixel to be dilated
     * @param posY The y location of the pixel to be dilated
     * @param posZ The z location of the pixel to be dilated
     * @param dimX The x dimension of the frame
     * @param dimY The y dimension of the frame
     * @param dimZ The number of frames
     * @param error The error radius
     */
    static void dilate_matrix(int[] matrix, int posX, int posY, int posZ, int dimX, int dimY, int dimZ, int error)
    {
        int startX = posX - error;
        while(startX < 0)
            startX++;
        int startY = posY - error;
        while(startY < 0)
            startY++;
        int endX = posX + error;
        while(endX > dimX)
            endX--;
        int endY = posY + error;
        while(endY > dimY)
            endY--;
        int x,y;
        for(x = startX; x < endX; x++){
            for(y = startY; y < endY; y++){
                double distance = Math.sqrt( Math.pow((double)(x-posX),2) + Math.pow((double)(y-posY),2) );
                if(distance < error)
                    matrix[x*dimY*dimZ + y*dimZ + posZ] = 1;
            }
        }
    }

    /**
     * Dilates the target matrix using the radius as a guide
     * @param matrix The reference matrix
     * @param dimX The x dimension of the video
     * @param dimY The y dimension of the video
     * @param dimZ The z dimension of the video
     * @param error The error radius to be dilated
     * @param newMatrix The target matrix
     */
    static void imdilate_disk(int[] matrix, int dimX, int dimY, int dimZ, int error, int[] newMatrix)
    {
        int x, y, z;
        for(z = 0; z < dimZ; z++){
            for(x = 0; x < dimX; x++){
                for(y = 0; y < dimY; y++){
                    if(matrix[x*dimY*dimZ + y*dimZ + z] == 1){
                        dilate_matrix(newMatrix, x, y, z, dimX, dimY, dimZ, error);
                    }
                }
            }
        }
    }

    /**
     * Fills a 2D array describing the offsets of the disk object
     * @param se The disk object
     * @param numOnes The number of ones in the disk
     * @param neighbors The array that will contain the offsets
     * @param radius The radius used for dilation
     */
    static void getneighbors(int[] se, int numOnes, double[] neighbors, int radius){
        int x, y;
        int neighY = 0;
        int center = radius - 1;
        int diameter = radius*2 -1;
        for(x = 0; x < diameter; x++){
            for(y = 0; y < diameter; y++){
                if(se[x*diameter + y] != 0){
                    neighbors[neighY*2] = (int)(y - center);
                    neighbors[neighY*2 + 1] = (int)(x - center);
                    neighY++;
                }
            }
        }
    }

    /**
     * The synthetic video sequence we will work with here is composed of a
     * single moving object, circular in shape (fixed radius)
     * The motion here is a linear motion
     * the foreground intensity and the backgrounf intensity is known
     * the image is corrupted with zero mean Gaussian noise
     * @param I The video itself
     * @param IszX The x dimension of the video
     * @param IszY The y dimension of the video
     * @param Nfr The number of frames of the video
     * @param seed The seed array used for number generation
     */
    public static void videoSequence(int[] I, int IszX, int IszY, int Nfr, double[] seed){
        int k;
        int max_size = IszX*IszY*Nfr;
        /*get object centers*/
        int x0 = (int)roundDouble(IszY/2.0);
        int y0 = (int)roundDouble(IszX/2.0);
        I[x0 *IszY *Nfr + y0 * Nfr  + 0] = 1;

        /*move point*/
        int xk, yk, pos;
        for(k = 1; k < Nfr; k++){
            xk = Math.abs(x0 + (k-1));
            yk = Math.abs(y0 - 2*(k-1));
            pos = yk * IszY * Nfr + xk *Nfr + k;
            if(pos >= max_size)
                pos = 0;
            I[pos] = 1;
        }

        /*dilate matrix*/
        int [] newMatrix = new int [IszX * IszY * Nfr];
        imdilate_disk(I, IszX, IszY, Nfr, 5, newMatrix);
        int x, y;
        for(x = 0; x < IszX; x++){
            for(y = 0; y < IszY; y++){
                for(k = 0; k < Nfr; k++){
                    I[x*IszY*Nfr + y*Nfr + k] = newMatrix[x*IszY*Nfr + y*Nfr + k];
                }
            }
        }

        /*define background, add noise*/
        setIf(0, 100, I, IszX, IszY, Nfr);
        setIf(1, 228, I, IszX, IszY, Nfr);
        /*add noise*/
        addNoise(I, IszX, IszY, Nfr, seed);
    }

    /**
     * Determines the likelihood sum based on the formula: SUM( (IK[IND] - 100)^2 - (IK[IND] - 228)^2)/ 100
     * @param I The 3D matrix
     * @param ind The current ind array
     * @param numOnes The length of ind array
     * @return A double representing the sum
     */
    double calcLikelihoodSum(int [] I, int [] ind, int numOnes){
        double likelihoodSum = 0.0;
        int y;
        for(y = 0; y < numOnes; y++)
            likelihoodSum += (Math.pow((I[ind[y]] - 100),2) - Math.pow((I[ind[y]]-228),2))/50.0;
        return likelihoodSum;
    }

    /**
     * Finds the first element in the CDF that is greater than or equal to the provided value and returns that index
     * This function uses sequential search
     * @param CDF The CDF
     * @param lengthCDF The length of CDF
     * @param value The value to be found
     * @return The index of value in the CDF; if value is never found, returns the last index
     */
    static int findIndex(double[] CDF, int lengthCDF, double value){
        int index = -1;
        int x;

        // for(int a = 0; a < lengthCDF; a++)
        // {
        // System.out.println("%f ",CDF[a]);
        // }
        // System.out.println("\n");

        // System.out.println("CDF[x] >= value ,%f >= %f \n",CDF[0],value);

        for(x = 0; x < lengthCDF; x++){
            if(CDF[x] >= value){
                index = x;
                break;
            }
        }
        if(index == -1){
            return lengthCDF-1;
        }
        return index;
    }

    /**
     * Finds the first element in the CDF that is greater than or equal to the provided value and returns that index
     * This function uses binary search before switching to sequential search
     * @param CDF The CDF
     * @param beginIndex The index to start searching from
     * @param endIndex The index to stop searching
     * @param value The value to find
     * @return The index of value in the CDF; if value is never found, returns the last index
     * Use at your own risk; not fully tested
     */
    int findIndexBin(double [] CDF, int beginIndex, int endIndex, double value){
        if(endIndex < beginIndex)
            return -1;
        int middleIndex = beginIndex + ((endIndex - beginIndex)/2);
        /*check the value*/
        if(CDF[middleIndex] >= value)
        {
            /*check that it's good*/
            if(middleIndex == 0)
                return middleIndex;
            else if(CDF[middleIndex-1] < value)
                return middleIndex;
            else if(CDF[middleIndex-1] == value)
            {
                while(middleIndex > 0 && CDF[middleIndex-1] == value)
                    middleIndex--;
                return middleIndex;
            }
        }
        if(CDF[middleIndex] > value)
            return findIndexBin(CDF, beginIndex, middleIndex+1, value);
        return findIndexBin(CDF, middleIndex-1, endIndex, value);
    }


    /**
     * The implementation of the particle filter using OpenMP for many frames
     * @see <a href="http://openmp.org/wp/">openmp</a>
     * This function is designed to work with a video of several frames. In addition, it references a provided MATLAB function which takes the video, the objxy matrix and the x and y arrays as arguments and returns the likelihoods
     * @param I The video to be run
     * @param IszX The x dimension of the video
     * @param IszY The y dimension of the video
     * @param Nfr The number of frames
     * @param seed The seed array used for random number generation
     * @param Nparticles The number of particles to be used
     */
    static void  particleFilter(int[] I, int IszX, int IszY, int Nfr, double[] seed, int Nparticles){

        int max_size = IszX*IszY*Nfr;

        //original particle centroid
        double xe = roundDouble(IszY/2.0);
        double ye = roundDouble(IszX/2.0);

        //expected object locations, compared to center
        int radius = 5;
        int diameter = radius*2 - 1;
        int [] disk = new int [diameter*diameter];
        strelDisk(disk, radius);
        int countOnes = 0;
        int x, y;
        for(x = 0; x < diameter; x++){
            for(y = 0; y < diameter; y++){
                if(disk[x*diameter + y] == 1)
                    countOnes++;
            }
        }



        double [] objxy = new double [countOnes*2];
        getneighbors(disk, countOnes, objxy, radius);

        //initial weights are all equal (1/Nparticles)
        double [] weights = new double  [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);

        for(x = 0; x < Nparticles; x++){
            weights[x] = 1/((double)(Nparticles));
        }

        //initial likelihood to 0.0
        double [] likelihood = new double [Nparticles];
        double [] arrayX = new double [Nparticles];
        double [] arrayY = new double [Nparticles];
        double [] xj = new double [Nparticles];
        double [] yj = new double [Nparticles];
        double [] CDF = new double [Nparticles];
        double [] u = new double [Nparticles];
        int [] ind = new int [countOnes * Nparticles];

        for(x = 0; x < Nparticles; x++){
            arrayX[x] = xe;
            arrayY[x] = ye;
        }
        int k;

        int indX, indY;
        for(k = 1; k < Nfr; k++){
            //apply motion model
            //draws sample from motion model (random walk). The only prior information
            //is that the object moves 2x as fast as in the y direction
            for(x = 0; x < Nparticles; x++){
                arrayX[x] += 1 + 5*randn(seed, x);
                arrayY[x] += -2 + 2*randn(seed, x);
            }
            //particle filter likelihood
            for(x = 0; x < Nparticles; x++){
                //compute the likelihood: remember our assumption is that you know
                // foreground and the background image intensity distribution.
                // Notice that we consider here a likelihood ratio, instead of
                // p(z|x). It is possible in this case. why? a hometask for you.
                //calc ind
                for(y = 0; y < countOnes; y++){
                    indX = (int) (roundDouble(arrayX[x]) + objxy[y*2 + 1]);
                    indY = (int) (roundDouble(arrayY[x]) + objxy[y*2]);
                    ind[x*countOnes + y] = Math.abs(indX*IszY*Nfr + indY*Nfr + k);
                    if(ind[x*countOnes + y] >= max_size)
                        ind[x*countOnes + y] = 0;
                }
                likelihood[x] = 0;
                for(y = 0; y < countOnes; y++)
                    likelihood[x] += (Math.pow((I[ind[x*countOnes + y]] - 100),2) - Math.pow((I[ind[x*countOnes + y]]-228),2))/50.0;
                likelihood[x] = likelihood[x]/((double) countOnes);
            }
            // update & normalize weights
            // using equation (63) of Arulampalam Tutorial
            for(x = 0; x < Nparticles; x++){
                weights[x] = weights[x] * Math.exp(likelihood[x]);
            }
            double sumWeights = 0;
            //#pragma omp parallel for private(x) reduction(+:sumWeights)
            for(x = 0; x < Nparticles; x++){
                sumWeights += weights[x];
            }
            for(x = 0; x < Nparticles; x++){
                weights[x] = weights[x]/sumWeights;
            }
            xe = 0;
            ye = 0;
            for(x = 0; x < Nparticles; x++){
                xe += arrayX[x] * weights[x];
                ye += arrayY[x] * weights[x];
            }


            CDF[0] = weights[0];
            for(x = 1; x < Nparticles; x++){
                CDF[x] = weights[x] + CDF[x-1];
            }
            double u1 = (1/((double)(Nparticles)))*randu(seed, 0);
            for(x = 0; x < Nparticles; x++){
                u[x] = u1 + x/((double)(Nparticles));
            }
            int j, i;

            for(j = 0; j < Nparticles; j++){
                i = findIndex(CDF, Nparticles, u[j]);
                if(i == -1)
                    i = Nparticles-1;
                xj[j] = arrayX[i];
                yj[j] = arrayY[i];

            }

            for(x = 0; x < Nparticles; x++){
                //reassign arrayX and arrayY
                arrayX[x] = xj[x];
                arrayY[x] = yj[x];
                weights[x] = 1/((double)(Nparticles));
            }

            resX = arrayX;
            resY = arrayY;
            resW = weights;
        }
    }

    public static void particleFilterVector(int[] I,
                                            int IszX,
                                            int IszY,
                                            int Nfr,
                                            double[] seed,
                                            double[] randuVectorResult,
                                            double[] randuVectorNum,
                                            int Nparticles){


        int max_size = IszX*IszY*Nfr;
        //original particle centroid
        double xe = roundDouble(IszY/2.0);
        double ye = roundDouble(IszX/2.0);

        //expected object locations, compared to center
        int radius = 5;
        int diameter = radius*2 - 1;
        int [] disk = new int [diameter*diameter];
        strelDisk(disk, radius);
        int countOnes = 0;
        int x, y;
        for(x = 0; x < diameter; x++){
            for(y = 0; y < diameter; y++){
                if(disk[x*diameter + y] == 1)
                    countOnes++;
            }
        }

        double [] objxy = new double[countOnes*2];
        getneighbors(disk, countOnes, objxy, radius);

        //initial weights are all equal (1/Nparticles)
        double [] weights = new double [Nparticles];

        int limit = DOUBLE_SPECIES.loopBound(Nparticles);
        DoubleVector xweights = DoubleVector.broadcast(DOUBLE_SPECIES, 1.0/(double)(Nparticles));

        for(x = 0; x < limit; x=x+DOUBLE_SPECIES_LENGTH){
            xweights.intoArray(weights, x);
        }

        for (; x < Nparticles; x++) {
            weights[x] = 1.0/(double)(Nparticles);
        }

        //initial likelihood to 0.0
        double [] likelihood = new double [Nparticles];
        double [] arrayX = new double [Nparticles];
        double [] arrayY = new double [Nparticles];
        double [] xj = new double [Nparticles];
        double [] yj = new double [Nparticles];
        double [] CDF = new double [Nparticles];
        double [] u = new double [Nparticles];
        int [] ind = new int [countOnes * Nparticles];
        double [] locations = new double [Nparticles];

        DoubleVector xArrayX = DoubleVector.broadcast(DOUBLE_SPECIES, xe);
        DoubleVector xArrayY = DoubleVector.broadcast(DOUBLE_SPECIES, ye);

        int i;
        for(i = 0; i < limit; i=i+DOUBLE_SPECIES_LENGTH){
            xArrayX.intoArray(arrayX, i);
            xArrayY.intoArray(arrayY, i);
        }

        for (; i < Nparticles; i++) {
            arrayX[i] = xe;
            arrayY[i] = ye;
        }


        DoubleVector    xAux;

        int k;
        int indX, indY;


        for(k = 1; k < Nfr; k++){
            //apply motion model
            //draws sample from motion model (random walk). The only prior information
            //is that the object moves 2x as fast as in the y direction
            for(x = 0; x < limit; x=x+DOUBLE_SPECIES_LENGTH){
                xArrayX = DoubleVector.fromArray(DOUBLE_SPECIES, arrayX, x); // _MM_LOAD_f64(&arrayX[x],gvl);
                xAux = randnVector(seed, x,randuVectorResult,randuVectorNum);
                xAux = xAux.mul(5.0); // _MM_MUL_f64(xAux, _MM_SET_f64(5.0,gvl),gvl);
                xAux = xAux.add(1.0); // _MM_ADD_f64(xAux, _MM_SET_f64(1.0,gvl),gvl);
                xArrayX = xAux.add(xArrayX); // _MM_ADD_f64(xAux, xArrayX ,gvl);
                xArrayX.intoArray(arrayX, x); // _MM_STORE_f64(&arrayX[x],xArrayX,gvl);

                xArrayY = DoubleVector.fromArray(DOUBLE_SPECIES, arrayY, x); // _MM_LOAD_f64(&arrayY[x],gvl);
                xAux = randnVector(seed, x,randuVectorResult,randuVectorNum);
                xAux = xAux.mul(2.0); // _MM_MUL_f64(xAux, _MM_SET_f64(2.0,gvl),gvl);
                xAux = xAux.add(-2.0); // _MM_ADD_f64(xAux, _MM_SET_f64(-2.0,gvl),gvl);
                xArrayY = xAux.add(xArrayY); // _MM_ADD_f64(xAux, xArrayY ,gvl);
                xArrayY.intoArray(arrayY, x); // _MM_STORE_f64(&arrayY[x],xArrayY,gvl);
            }

            for (; x < Nparticles; x++) {
                arrayX[x] += 1 + 5*randn(seed, x);
                arrayY[x] += -2 + 2*randn(seed, x);
            }

            //particle filter likelihood
            for(x = 0; x < Nparticles; x++){
                //compute the likelihood: remember our assumption is that you know
                // foreground and the background image intensity distribution.
                // Notice that we consider here a likelihood ratio, instead of
                // p(z|x). It is possible in this case. why? a hometask for you.
                //calc ind
                for(y = 0; y < countOnes; y++){
                    indX = (int) (roundDouble(arrayX[x]) + objxy[y*2 + 1]);
                    indY = (int) (roundDouble(arrayY[x]) + objxy[y*2]);
                    ind[x*countOnes + y] = Math.abs(indX*IszY*Nfr + indY*Nfr + k);
                    if(ind[x*countOnes + y] >= max_size)
                        ind[x*countOnes + y] = 0;
                }
                likelihood[x] = 0;
                for(y = 0; y < countOnes; y++)
                    likelihood[x] += (Math.pow((I[ind[x*countOnes + y]] - 100),2) - Math.pow((I[ind[x*countOnes + y]]-228),2))/50.0;
                likelihood[x] = likelihood[x]/((double) countOnes);
            }
            // update & normalize weights
            // using equation (63) of Arulampalam Tutorial
            for(x = 0; x < Nparticles; x++){
                weights[x] = weights[x] * Math.exp(likelihood[x]);
            }
            double sumWeights = 0;
            for(x = 0; x < Nparticles; x++){
                sumWeights += weights[x];
            }
            for(x = 0; x < Nparticles; x++){
                weights[x] = weights[x]/sumWeights;
            }
            xe = 0;
            ye = 0;
            // estimate the object location by expected values
            for(x = 0; x < Nparticles; x++){
                xe += arrayX[x] * weights[x];
                ye += arrayY[x] * weights[x];
            }

            CDF[0] = weights[0];
            for(x = 1; x < Nparticles; x++){
                CDF[x] = weights[x] + CDF[x-1];
            }
            double u1 = (1/((double)(Nparticles)))*randu(seed, 0);
            for(x = 0; x < Nparticles; x++){
                u[x] = u1 + x/((double)(Nparticles));
            }

            int j;

            VectorMask<Double> xComp;


            DoubleVector xCDF;
            DoubleVector xU;
            DoubleVector xArray;

            int vector_complete;
            int valid;

            for(i = 0; i < limit; i=i+DOUBLE_SPECIES_LENGTH){
                vector_complete = 0;
                xArray  = DoubleVector.broadcast(DOUBLE_SPECIES, Nparticles-1);
                xU      = DoubleVector.fromArray(DOUBLE_SPECIES, u, i/*, gvl*/);
                for(j = 0; j < Nparticles; j++){
                    xCDF = DoubleVector.broadcast(DOUBLE_SPECIES, CDF[j]);
                    xComp = xCDF.compare(VectorOperators.GE, xU/*, gvl*/);
                    xComp = xComp.not().eq(xMask);
			        valid = xComp.firstTrue();
                    if(valid != DOUBLE_SPECIES_LENGTH && valid + i < Nparticles)
                    {
//                        xArray = xArray.blend(j, xComp);
                        xArray = DoubleVector.zero(DOUBLE_SPECIES).add(j, xComp/*.and(gvl)*/).add(xArray, xComp.not()/*.and(gvl)*/);
                        vector_complete = xComp.or(xMask).trueCount();


                    }
                    if(vector_complete == DOUBLE_SPECIES_LENGTH){ break; }
                }
                xArray.intoArray(locations, i/*, gvl*/);

            }
            for(j = limit; j < Nparticles; j++){
                i = findIndex(CDF, Nparticles, u[j]);
                if(i == -1)
                    i = Nparticles-1;
                xj[j] = arrayX[i];
                yj[j] = arrayY[i];

            }


            for(j = 0; j < limit; j++){
                i = (int) locations[j];
                xj[j] = arrayX[i];
                yj[j] = arrayY[i];
            }

            for(x = 0; x < Nparticles; x++){
                //reassign arrayX and arrayY
                arrayX[x] = xj[x];
                arrayY[x] = yj[x];
                weights[x] = 1/((double)(Nparticles));

            }

            resXV = arrayX;
            resYV = arrayY;
            resWV = weights;
        }
    }


    public static void main(String[] args) {

        final boolean vectorize = true;
        int IszX, IszY, Nfr, Nparticles;

        IszX = 16;
        IszY = 16;
        Nfr = 16;
        Nparticles = 20;

        double[] seed = new double[Nparticles];
        int i;

        for (i = 0; i < Nparticles; i++) {
            // in c -> time(0)*i; TODO Check if correct
            seed[i] = i; // (int) /*System.currentTimeMillis() * i; // */ Instant.now().getEpochSecond() * i;
        }


        // allocate Matrix
        int [] I = new int [IszX * IszY * Nfr];


        videoSequence(I, IszX, IszY, Nfr, seed);



        double [] randuVectorResult = new double [DOUBLE_SPECIES_LENGTH];
        double [] randuVectorNum = new double[DOUBLE_SPECIES_LENGTH];


        if (vectorize) {
            particleFilterVector(I, IszX, IszY, Nfr, seed, randuVectorResult, randuVectorNum, Nparticles);
        } else {
            particleFilter(I, IszX, IszY, Nfr, seed, Nparticles);
        }
    }


}
