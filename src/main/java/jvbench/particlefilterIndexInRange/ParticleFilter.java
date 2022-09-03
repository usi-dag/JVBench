package jvbench.particlefilterIndexInRange;

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
//        long long start = get_time();
        //original particle centroid
        double xe = roundDouble(IszY/2.0);
        double ye = roundDouble(IszX/2.0);

        //expected object locations, compared to center
        int radius = 5;
        int diameter = radius*2 - 1;
        int [] disk = new int [diameter*diameter]; // (int *)malloc(diameter*diameter*sizeof(int));
        strelDisk(disk, radius);
        int countOnes = 0;
        int x, y;
        for(x = 0; x < diameter; x++){
            for(y = 0; y < diameter; y++){
                if(disk[x*diameter + y] == 1)
                    countOnes++;
            }
        }

        //System.out.println("countOnes = %d \n",countOnes); // 69

        double [] objxy = new double [countOnes*2]; // (double *)malloc(countOnes*2*sizeof(double));
        getneighbors(disk, countOnes, objxy, radius);

//        long long get_neighbors = get_time();
//        System.out.println("TIME TO GET NEIGHBORS TOOK: %f\n", elapsed_time(start, get_neighbors));
        //initial weights are all equal (1/Nparticles)
        double [] weights = new double  [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        //#pragma omp parallel for shared(weights, Nparticles) private(x)
        for(x = 0; x < Nparticles; x++){
            weights[x] = 1/((double)(Nparticles));
        }
//        long long get_weights = get_time();
//        System.out.println("TIME TO GET WEIGHTSTOOK: %f\n", elapsed_time(get_neighbors, get_weights));
        //initial likelihood to 0.0
        double [] likelihood = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        double [] arrayX = new double [Nparticles];  // (double *)malloc(sizeof(double)*Nparticles);
        double [] arrayY = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        double [] xj = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        double [] yj = new double [Nparticles]; //(double *)malloc(sizeof(double)*Nparticles);
        double [] CDF = new double [Nparticles]; //(double *)malloc(sizeof(double)*Nparticles);
        double [] u = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        int [] ind = new int [countOnes * Nparticles]; // (int*)malloc(sizeof(int)*countOnes*Nparticles);
        //#pragma omp parallel for shared(arrayX, arrayY, xe, ye) private(x)
        for(x = 0; x < Nparticles; x++){
            arrayX[x] = xe;
            arrayY[x] = ye;
        }
        int k;

//        System.out.println("TIME TO SET ARRAYS TOOK: %f\n", elapsed_time(get_weights, get_time()));
        int indX, indY;
        for(k = 1; k < Nfr; k++){
//            long long set_arrays = get_time();
            //apply motion model
            //draws sample from motion model (random walk). The only prior information
            //is that the object moves 2x as fast as in the y direction
            //#pragma omp parallel for shared(arrayX, arrayY, Nparticles, seed) private(x)
            for(x = 0; x < Nparticles; x++){
                arrayX[x] += 1 + 5*randn(seed, x);
                arrayY[x] += -2 + 2*randn(seed, x);
            }
//            long long error = get_time();
//            System.out.println("TIME TO SET ERROR TOOK: %f\n", elapsed_time(set_arrays, error));
            //particle filter likelihood
            //#pragma omp parallel for shared(likelihood, I, arrayX, arrayY, objxy, ind) private(x, y, indX, indY)
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
//            long long likelihood_time = get_time();
//            System.out.println("TIME TO GET LIKELIHOODS TOOK: %f\n", elapsed_time(error, likelihood_time));
            // update & normalize weights
            // using equation (63) of Arulampalam Tutorial
            //#pragma omp parallel for shared(Nparticles, weights, likelihood) private(x)
            for(x = 0; x < Nparticles; x++){
                weights[x] = weights[x] * Math.exp(likelihood[x]);
            }
//            long long exponential = get_time();
//            System.out.println("TIME TO GET EXP TOOK: %f\n", elapsed_time(likelihood_time, exponential));
            double sumWeights = 0;
            //#pragma omp parallel for private(x) reduction(+:sumWeights)
            for(x = 0; x < Nparticles; x++){
                sumWeights += weights[x];
            }
//            long long sum_time = get_time();
//            System.out.println("TIME TO SUM WEIGHTS TOOK: %f\n", elapsed_time(exponential, sum_time));
            //#pragma omp parallel for shared(sumWeights, weights) private(x)
            for(x = 0; x < Nparticles; x++){
                weights[x] = weights[x]/sumWeights;
            }
//            long long normalize = get_time();
//            System.out.println("TIME TO NORMALIZE WEIGHTS TOOK: %f\n", elapsed_time(sum_time, normalize));
            xe = 0;
            ye = 0;
            // estimate the object location by expected values
            //#pragma omp parallel for private(x) reduction(+:xe, ye)
            for(x = 0; x < Nparticles; x++){
                xe += arrayX[x] * weights[x];
                ye += arrayY[x] * weights[x];
            }
//            long long move_time = get_time();
//            System.out.println("TIME TO MOVE OBJECT TOOK: %f\n", elapsed_time(normalize, move_time));
//            System.out.println("XE: %lf\n", xe);
//            System.out.println("YE: %lf\n", ye);
            double distance = Math.sqrt( Math.pow((double)(xe-(int)roundDouble(IszY/2.0)),2) + Math.pow((double)(ye-(int)roundDouble(IszX/2.0)),2) );
//            System.out.println("%lf\n", distance);
            //display(hold off for now)

            //pause(hold off for now)

            //resampling


            CDF[0] = weights[0];
            for(x = 1; x < Nparticles; x++){
                CDF[x] = weights[x] + CDF[x-1];
            }
//            long long cum_sum = get_time();
//            System.out.println("TIME TO CALC CUM SUM TOOK: %f\n", elapsed_time(move_time, cum_sum));
            double u1 = (1/((double)(Nparticles)))*randu(seed, 0);
            //#pragma omp parallel for shared(u, u1, Nparticles) private(x)
            for(x = 0; x < Nparticles; x++){
                u[x] = u1 + x/((double)(Nparticles));
            }
//            long long u_time = get_time();
//            System.out.println("TIME TO CALC U TOOK: %f\n", elapsed_time(cum_sum, u_time));
            int j, i;

            //#pragma omp parallel for shared(CDF, Nparticles, xj, yj, u, arrayX, arrayY) private(i, j)
            for(j = 0; j < Nparticles; j++){
                i = findIndex(CDF, Nparticles, u[j]);
                if(i == -1)
                    i = Nparticles-1;
                //System.out.println("%d ", i);
                xj[j] = arrayX[i];
                yj[j] = arrayY[i];

            }
            //System.out.println("\n");

//            long long xyj_time = get_time();
//            System.out.println("TIME TO CALC NEW ARRAY X AND Y TOOK: %f\n", elapsed_time(u_time, xyj_time));

            //#pragma omp parallel for shared(weights, Nparticles) private(x)
            for(x = 0; x < Nparticles; x++){
                //reassign arrayX and arrayY
                arrayX[x] = xj[x];
                arrayY[x] = yj[x];
                weights[x] = 1/((double)(Nparticles));

//                System.out.println(x + " x: " + arrayX[x] + ", y: " + arrayY[x] + ", w: " + weights[x]);

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
//        long long start = get_time();
        //original particle centroid
        double xe = roundDouble(IszY/2.0);
        double ye = roundDouble(IszX/2.0);

        //expected object locations, compared to center
        int radius = 5;
        int diameter = radius*2 - 1;
        int [] disk = new int [diameter*diameter]; // (int *)malloc(diameter*diameter*sizeof(int));
        strelDisk(disk, radius);
        int countOnes = 0;
        int x, y;
        for(x = 0; x < diameter; x++){
            for(y = 0; y < diameter; y++){
                if(disk[x*diameter + y] == 1)
                    countOnes++;
            }
        }

        //printf("countOnes = %d \n",countOnes); // 69

        double [] objxy = new double[countOnes*2]; // (double *)malloc(countOnes*2*sizeof(double));
        getneighbors(disk, countOnes, objxy, radius);

//        long long get_neighbors = get_time();
//        printf("TIME TO GET NEIGHBORS TOOK: %f\n", elapsed_time(start, get_neighbors));
        //initial weights are all equal (1/Nparticles)
        double [] weights = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        //#pragma omp parallel for shared(weights, Nparticles) private(x)
    /*
    for(x = 0; x < Nparticles; x++){
        weights[x] = 1/((double)(Nparticles));
    }*/
//        unsigned long int gvl = __builtin_epi_vsetvl(Nparticles, __epi_e64, __epi_m1);

        DoubleVector xweights = DoubleVector.broadcast(DOUBLE_SPECIES, 1.0/(double)(Nparticles));
        for(x = 0; x < Nparticles; x=x+DOUBLE_SPECIES_LENGTH){
            VectorMask<Double> mask = DOUBLE_SPECIES.indexInRange(x, Nparticles);
            xweights.intoArray(weights, x, mask);
        }

        //FENCE();

//        long long get_weights = get_time();
//        printf("TIME TO GET WEIGHTSTOOK: %f\n", elapsed_time(get_neighbors, get_weights));
        //initial likelihood to 0.0
        double [] likelihood = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        double [] arrayX = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        double [] arrayY = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        double [] xj = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        double [] yj = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        double [] CDF = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        double [] u = new double [Nparticles]; // (double *)malloc(sizeof(double)*Nparticles);
        int [] ind = new int [countOnes * Nparticles]; // (int*)malloc(sizeof(int)*countOnes*Nparticles);
        // Se usa adentro del for, aqui para no repetir
        double [] locations = new double [Nparticles]; // (long int *)malloc(sizeof(long int)*Nparticles);

    /*
    //#pragma omp parallel for shared(arrayX, arrayY, xe, ye) private(x)
    for(x = 0; x < Nparticles; x++){
        arrayX[x] = xe;
        arrayY[x] = ye;
    }
    */
        DoubleVector xArrayX = DoubleVector.broadcast(DOUBLE_SPECIES, xe);
        DoubleVector xArrayY = DoubleVector.broadcast(DOUBLE_SPECIES, ye);
        int i;
        for(i = 0; i < Nparticles; i=i+DOUBLE_SPECIES_LENGTH){
            VectorMask<Double> mask = DOUBLE_SPECIES.indexInRange(i, Nparticles);
            xArrayX.intoArray(arrayX, i, mask);
            xArrayY.intoArray(arrayY, i, mask);
        }


        DoubleVector    xAux;

        int k;
//        printf("TIME TO SET ARRAYS TOOK: %f\n", elapsed_time(get_weights, get_time()));
        int indX, indY;

        boolean [] mask2 = new boolean[DOUBLE_SPECIES_LENGTH]; // TODO check memory LEAK

        for(k = 1; k < Nfr; k++){
//            long long set_arrays = get_time();
            //apply motion model
            //draws sample from motion model (random walk). The only prior information
            //is that the object moves 2x as fast as in the y direction
            for(x = 0; x < Nparticles; x=x+DOUBLE_SPECIES_LENGTH){
                VectorMask<Double> mask = DOUBLE_SPECIES.indexInRange(x, Nparticles);
                xArrayX = DoubleVector.fromArray(DOUBLE_SPECIES, arrayX, x, mask); // _MM_LOAD_f64(&arrayX[x],gvl);
                xAux = randnVector(seed, x,randuVectorResult,randuVectorNum);
                xAux = xAux.mul(5.0); // _MM_MUL_f64(xAux, _MM_SET_f64(5.0,gvl),gvl);
                xAux = xAux.add(1.0); // _MM_ADD_f64(xAux, _MM_SET_f64(1.0,gvl),gvl);
                xArrayX = xAux.add(xArrayX); // _MM_ADD_f64(xAux, xArrayX ,gvl);
                xArrayX.intoArray(arrayX, x, mask); // _MM_STORE_f64(&arrayX[x],xArrayX,gvl);

                xArrayY = DoubleVector.fromArray(DOUBLE_SPECIES, arrayY, x, mask); // _MM_LOAD_f64(&arrayY[x],gvl);
                xAux = randnVector(seed, x,randuVectorResult,randuVectorNum);
                xAux = xAux.mul(2.0); // _MM_MUL_f64(xAux, _MM_SET_f64(2.0,gvl),gvl);
                xAux = xAux.add(-2.0); // _MM_ADD_f64(xAux, _MM_SET_f64(-2.0,gvl),gvl);
                xArrayY = xAux.add(xArrayY); // _MM_ADD_f64(xAux, xArrayY ,gvl);
                xArrayY.intoArray(arrayY, x, mask); // _MM_STORE_f64(&arrayY[x],xArrayY,gvl);
            }

//            FENCE();
        /*
        //#pragma omp parallel for shared(arrayX, arrayY, Nparticles, seed) private(x)
        for(x = 0; x < Nparticles; x++){
            arrayX[x] += 1 + 5*randn(seed, x);
            arrayY[x] += -2 + 2*randn(seed, x);
        }
        */
//            long long error = get_time();
//            printf("TIME TO SET ERROR TOOK: %f\n", elapsed_time(set_arrays, error));
            //particle filter likelihood
            //#pragma omp parallel for shared(likelihood, I, arrayX, arrayY, objxy, ind) private(x, y, indX, indY)
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
//            long long likelihood_time = get_time();
//            printf("TIME TO GET LIKELIHOODS TOOK: %f\n", elapsed_time(error, likelihood_time));
            // update & normalize weights
            // using equation (63) of Arulampalam Tutorial
            //#pragma omp parallel for shared(Nparticles, weights, likelihood) private(x)
            for(x = 0; x < Nparticles; x++){
                weights[x] = weights[x] * Math.exp(likelihood[x]);
            }
//            long long exponential = get_time();
//            printf("TIME TO GET EXP TOOK: %f\n", elapsed_time(likelihood_time, exponential));
            double sumWeights = 0;
            //#pragma omp parallel for private(x) reduction(+:sumWeights)
            for(x = 0; x < Nparticles; x++){
                sumWeights += weights[x];
            }
//            long long sum_time = get_time();
//            printf("TIME TO SUM WEIGHTS TOOK: %f\n", elapsed_time(exponential, sum_time));
            //#pragma omp parallel for shared(sumWeights, weights) private(x)
            for(x = 0; x < Nparticles; x++){
                weights[x] = weights[x]/sumWeights;
            }
//            long long normalize = get_time();
//            printf("TIME TO NORMALIZE WEIGHTS TOOK: %f\n", elapsed_time(sum_time, normalize));
            xe = 0;
            ye = 0;
            // estimate the object location by expected values
            //#pragma omp parallel for private(x) reduction(+:xe, ye)
            for(x = 0; x < Nparticles; x++){
                xe += arrayX[x] * weights[x];
                ye += arrayY[x] * weights[x];
            }
//            long long move_time = get_time();
//            printf("TIME TO MOVE OBJECT TOOK: %f\n", elapsed_time(normalize, move_time));
//            printf("XE: %lf\n", xe);
//            printf("YE: %lf\n", ye);
            double distance = Math.sqrt( Math.pow((double)(xe-(int)roundDouble(IszY/2.0)),2) + Math.pow((double)(ye-(int)roundDouble(IszX/2.0)),2) );
//            printf("%lf\n", distance);
            //display(hold off for now)

            //pause(hold off for now)

            //resampling

            CDF[0] = weights[0];
            for(x = 1; x < Nparticles; x++){
                CDF[x] = weights[x] + CDF[x-1];
            }
//            long long cum_sum = get_time();
//            printf("TIME TO CALC CUM SUM TOOK: %f\n", elapsed_time(move_time, cum_sum));
            double u1 = (1/((double)(Nparticles)))*randu(seed, 0);
            //#pragma omp parallel for shared(u, u1, Nparticles) private(x)
            for(x = 0; x < Nparticles; x++){
                u[x] = u1 + x/((double)(Nparticles));
            }
//            long long u_time = get_time();
//            printf("TIME TO CALC U TOOK: %f\n", elapsed_time(cum_sum, u_time));

            int j;

            VectorMask<Double> xComp; //_MMR_MASK_i64           xComp;
            VectorMask<Double> xMask; // _MMR_MASK_i64           xMask;

            DoubleVector xCDF; //_MMR_f64          xCDF;
            DoubleVector xU; //_MMR_f64          xU;
            DoubleVector xArray; // _MMR_i64          xArray;

            int vector_complete;
            int valid;
//            gvl     = __builtin_epi_vsetvl(Nparticles, __epi_e64, __epi_m1);
            xMask   = VectorMask.fromArray(DOUBLE_SPECIES, mask2, 0);
            for(i = 0; i < Nparticles; i=i+DOUBLE_SPECIES_LENGTH){
                VectorMask<Double> mask = DOUBLE_SPECIES.indexInRange(i, Nparticles);
//                gvl     = __builtin_epi_vsetvl(Nparticles-i, __epi_e64, __epi_m1);
//                VectorMask<Double> gvl = DOUBLE_SPECIES.indexInRange(i, Nparticles);
                vector_complete = 0;
                // xMask   = VectorMask.fromArray(DOUBLE_SPECIES, mask, 0); //  _MM_CAST_i1_i64(__builtin_epi_vbroadcast_1xi64(0,gvl));
                xArray  = DoubleVector.broadcast(DOUBLE_SPECIES, Nparticles-1); // _MM_SET_i64(Nparticles-1,gvl);
                xU      = DoubleVector.fromArray(DOUBLE_SPECIES, u, i , mask); // _MM_LOAD_f64(&u[i],gvl);
                for(j = 0; j < Nparticles; j++){
                    xCDF = DoubleVector.broadcast(DOUBLE_SPECIES, CDF[j]); // _MM_SET_f64(CDF[j],gvl);
                    xComp = xCDF.compare(VectorOperators.GE, xU, mask); // _MM_VFGE_f64(xCDF,xU,gvl);
//                    xComp = (xComp.or(xMask)).andNot(xComp.and(xMask)); // _MM_VMXOR_i64(xComp,xMask,gvl);
                    xComp = xComp.not().eq(xMask);
			        valid = xComp.firstTrue(); //_MM_VMFIRST_i64(xComp,gvl);
                    if(valid != DOUBLE_SPECIES_LENGTH && valid + i < Nparticles)
                    {

//                        xArray = xArray.blend(j, xComp);
                        xArray = DoubleVector.zero(DOUBLE_SPECIES).add(j, xComp.and(mask)).add(xArray, xComp.not().and(mask)); // _MM_MERGE_i64(xArray,_MM_SET_i64(j,gvl),xComp,gvl);
//                        xMask = xComp.or(xMask); // _MM_VMOR_i64(xComp,xMask,gvl);


                        vector_complete = xComp.or(xMask).trueCount(); // xMask.trueCount(); // _MM_VMPOPC_i64(xMask,gvl);


                    }
                    if(vector_complete == DOUBLE_SPECIES_LENGTH){ break; }
                }
                xArray.intoArray(locations, i/*, gvl*/); // _MM_STORE_i64(&locations[i],xArray,gvl);
                // ---------------------------------------
                //xArray = _MM_MUL_i64(xArray,_MM_SET_i64(8,gvl),gvl); // Position in elements to position in bytes
                //xarrayX = _MM_LOAD_INDEX_f64(&arrayX[i],xArray,gvl);
                //xarrayY = _MM_LOAD_INDEX_f64(&arrayY[i],xArray,gvl);
                //_MM_STORE_f64(&xj[i],xarrayX,gvl);
                //_MM_STORE_f64(&yj[i],xarrayY,gvl);
                // This commented lines corresponds to the scalar code below
            }

//            FENCE();

            //#pragma omp parallel for shared(CDF, Nparticles, xj, yj, u, arrayX, arrayY) private(i, j)
            for(j = 0; j < Nparticles; j++){
                i = (int) locations[j];
                xj[j] = arrayX[i];
                yj[j] = arrayY[i];
            }
            // for(j = 0; j < Nparticles; j++){ printf("%lf ", xj[i]); } printf("\n");
            // for(j = 0; j < Nparticles; j++){ printf("%lf ", yj[i]); } printf("\n");

//            long long xyj_time = get_time();
//            printf("TIME TO CALC NEW ARRAY X AND Y TOOK: %f\n", elapsed_time(u_time, xyj_time));

            //#pragma omp parallel for shared(weights, Nparticles) private(x)
            for(x = 0; x < Nparticles; x++){
                //reassign arrayX and arrayY
                arrayX[x] = xj[x];
                arrayY[x] = yj[x];
                weights[x] = 1/((double)(Nparticles));

//                System.out.println(x + " x: " + arrayX[x] + ", y: " + arrayY[x] + ", w: " + weights[x]);
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
