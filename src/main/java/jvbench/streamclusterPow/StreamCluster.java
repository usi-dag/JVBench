package jvbench.streamclusterPow;

import jdk.incubator.vector.FloatVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

public class StreamCluster {

    private static int nProc = 1;
    private static long SEED = 1;
    private static final int SP = 1;
    private static final int ITER = 3;
    private static final int CACHE_LINE = 64;
    private static boolean[] switchMembership;
    private static boolean[] isCenterArray;
    private static int[] centerTable;

    static float[] block;
    static float[] centerBlock;



    private static final Random random = new Random(SEED);
    private static boolean useVectorAPI;
    private static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_MAX;
    private static final int SPECIES_LENGTH = SPECIES.length();

    private static double srand48(long seedval) {
        SEED = seedval & 0xFFFFFFFF;
        return (double) SEED / (1L << 48);
    }

    static long lrand48() {
        return  random.nextInt() & Integer.MAX_VALUE;
    }

    public static void main(String[] args) {

        useVectorAPI = true;

        String outputFileName;
        String inputFileName;

        long kMin;
        long kMax;
        long n;
        int chunkSize;
        long clusterSize;
        int dim;


        kMin = 3;
        kMax = 10;
        dim = 3;
        n = 128;
        chunkSize = 128;
        clusterSize = 10;
        inputFileName = "input.txt";
        outputFileName = "output.txt";
        nProc = 1;
        PStream stream;
        srand48(SEED);

        if ( n > 0) {
            stream = new SimStream(n);
        } else {
            stream = new FileStream(inputFileName);
        }

        streamCluster(stream, kMin, kMax, dim, chunkSize, clusterSize, outputFileName, useVectorAPI);
    }

    public static PStream init(String inputFileName) {

        nProc = 1;
        PStream stream;
        srand48(SEED);

        stream = new FileStream(inputFileName);

        return stream;

    }

    public static void streamCluster(PStream stream, long kMin, long kMax, int dim, int chunkSize, long centerSize, String outputFileName, boolean vectorize) {

        useVectorAPI = vectorize;

        block       = new float[(int) (chunkSize * dim)];
        centerBlock = new float[(int) (chunkSize * dim)];
        long[] centerIDs = new long[(int) (centerSize * dim)];


        Points points = new Points(dim, chunkSize, chunkSize);

        Point[] pointArray = points.getPoints();
        for (int i = 0; i < chunkSize; i++) {
//            pointArray[i].setCoord(block[i*dim]);
            pointArray[i] = new Point();
            pointArray[i].setCoordIndex(i*dim);
        }


        Points centers = new Points(dim, 0, (int) centerSize);
        centers.setDim(dim);

        Point[] pointCenterArray = centers.getPoints();
        for (int i = 0; i < centerSize; i++) {
            pointCenterArray[i] = new Point();
            pointCenterArray[i].setCoordIndex(i*dim);
            pointCenterArray[i].setWeight(1);
        }

        centers.setPoints(pointCenterArray);

        long idOffset = 0;
        long[] kFinal = new long[1];

        while (true) {
            int numRead = stream.read(block, dim, chunkSize);

//            System.out.println(Arrays.toString(block));

            points.setNum(numRead);

            for (int i = 0; i < points.getNum(); i++) {
                points.getPoints()[i].setWeight(1);
            }


            switchMembership = new boolean[(int) points.getNum()];
            isCenterArray = new boolean[(int) points.getNum()];
            centerTable = new int[(int) points.getNum()];

            localSearch(points, kMin, kMax, kFinal, false);
            contCenters(points, false);
            if (kFinal[0] + centers.getNum() > centerSize) {
                System.err.println("oops! no more space for centers");
                System.exit(1);
            }

            copyCenters(points, centers, centerIDs, idOffset);
            idOffset += numRead;

            if (stream.feof()) {
                break;
            }

        }

        //finally cluster all temp centers

        switchMembership = new boolean[(int) points.getNum()];
        isCenterArray = new boolean[(int) points.getNum()];
        centerTable = new int[(int) points.getNum()];

        localSearch(centers, kMin, kMax, kFinal, true);
        contCenters(centers, true);
//        outCentersIDs(centers, centerIDs, outputFileName);


    }

    private static void outCentersIDs(Points centers, long[] centerIDs, String outputFileName) {
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName))) {
            int [] is_a_median = new int[(int) centers.getNum()]; // (int*)calloc( sizeof(int), centers->num );
            for( int i =0 ; i< centers.getNum(); i++ ) {
                is_a_median[centers.getPoints()[i].getAssign()] = 1;
            }

            for( int i = 0; i < centers.getNum(); i++ ) {
                if( is_a_median[i] != 0 ) {
                    writer.write(centerIDs[i] + ""); // fprintf(fp, "%u\n", centerIDs[i]);
//                    System.out.println(centerIDs[i]);
                    writer.newLine();
                    writer.write(centers.getPoints()[i].getWeight() + "");// fprintf(fp, "%lf\n", centers->p[i].weight);
//                    System.out.println(centers.getPoints()[i].getWeight());
                    writer.newLine();
                    for( int k = 0; k < centers.getDim(); k++ ) {
                        writer.write(centers.getPoints()[i].getCoord(k, false) + " ");// fprintf(fp, "%lf ", centers->p[i].coord[k]);
//                        System.out.print(centers.getPoints()[i].getCoord(k, false) + " ");
                    }
//                    System.out.println("\n");
                    writer.newLine();
                    writer.newLine();
//                    fprintf(fp,"\n\n");
                }
            }
            

        } catch (IOException e) {
            System.err.println("Unable to open file = " + outputFileName + " " + e);
        }
    }

    private static void copyCenters(Points points, Points centers, long[] centerIDs, long offset) {
        int i;
        int k;

        boolean [] is_a_median = new boolean[(int) points.getNum()]; // (bool *) calloc(points->num, sizeof(bool));

        /* mark the centers */
        for ( i = 0; i < points.getNum(); i++ ) {
            is_a_median[points.getPoints()[i].getAssign()] = true;
        }

        k= centers.getNum();

        /* count how many  */
        for ( i = 0; i < points.getNum(); i++ ) {
            if ( is_a_median[i] ) {
//                System.arraycopy(centers.getPoints()[k].getCoord(), 0, points.getPoints()[i].getCoord(), 0, points.getDim());
//                memcpy( centers->p[k].coord, points->p[i].coord, points->dim * sizeof(float));
                for (int j = 0; j < points.getDim(); j++) {
                    float coord = points.getPoints()[i].getCoord(j, false);
                    centers.getPoints()[k].setCoord(coord, j, true);
                }
                centers.getPoints()[k].setWeight(points.getPoints()[i].getWeight());
                centerIDs[k] = i + offset;
                k++;
            }
        }

        centers.setNum(k);

    }

    private static void contCenters(Points points, boolean isCenter) {
        int i, ii;
        float relweight;

        for (i=0;i<points.getNum();i++) {
            /* compute relative weight of this point to the cluster */
            if (points.getPoints()[i].getAssign() != i) {
                relweight = points.getPoints()[points.getPoints()[i].getAssign()].getWeight() + points.getPoints()[i].getWeight();
                relweight = points.getPoints()[i].getWeight()/relweight;
                for (ii=0;ii< points.getDim();ii++) {
                    float coord = points.getPoints()[points.getPoints()[i].getAssign()].getCoord(ii, isCenter);
                    coord *= 1.0 - relweight;
                    coord += points.getPoints()[i].getCoord(ii, isCenter)*relweight;
                    points.getPoints()[points.getPoints()[i].getAssign()].setCoord(coord, ii, isCenter);
                }
                int weight = (int) points.getPoints()[points.getPoints()[i].getAssign()].getWeight();
                weight += points.getPoints()[i].getWeight();
                points.getPoints()[points.getPoints()[i].getAssign()].setWeight(weight);
            }
        }
    }

    private static void localSearch(Points points, long kMin, long kMax, long[] kFinal, boolean isCenter) {

        Median[] arg = new Median[nProc];

        for (int i = 0; i < nProc; i++) {
            arg[i] = new Median(points, kMin, kMax, kFinal, i);
        }

        localSearchSub(arg[0], isCenter);
    }

    private static void localSearchSub(Median arg, boolean isCenter) {
        pkmedian(arg.getPoints(), arg.getKmin(), arg.getKmax(), arg.getKfinal(), arg.getPid(), isCenter);
    }

    private static float pkmedian(Points points, long kmin, long kmax, long[] kfinal, int pid, boolean isCenter) {

        int i;
        double cost;
        double lastcost;
        double hiz, loz, z;

        long[] k = new long[1];
        int[] feasible = new int[(int) points.getNum()];
        int numfeasible = 0;
        double[] hizs = new double[nProc];



        if (pid == 0) hizs = new double[nProc]; // (double*)calloc(nproc,sizeof(double));
        hiz = loz = 0.0;
        long numberOfPoints = points.getNum();
        long ptDimension = points.getDim();

        //my block
        int bsize = points.getNum() / nProc;
        int k1 = bsize * pid;
        long k2 = k1 + bsize;
        if (pid == nProc - 1) k2 = points.getNum();


        double myhiz = 0;
        for (long kk = k1; kk < k2; kk++) {
            myhiz += dist(points.getPoints()[(int) kk], points.getPoints()[0], (int) ptDimension, isCenter) * points.getPoints()[(int) kk].getWeight();
        }
        hizs[pid] = myhiz;


        for (i = 0; i < nProc; i++) {
            hiz += hizs[i];
        }

        loz = 0.0;
        z = (hiz + loz) / 2.0;
        /* NEW: Check whether more centers than points! */
        if (points.getNum() <= kmax) {
            /* just return all points as facilities */
            for (int kk = k1; kk < k2; kk++) {
                points.getPoints()[kk].setAssign(kk);
                points.getPoints()[kk].setCost(0);
            }
            cost = 0;
            if (pid == 0) {
                kfinal[0] = k[0];
            }
            return (float) cost;
        }
        if (pid == 0) shuffle(points);

        cost = pspeedy(points, z, k, pid);

        i = 0;
        /* give speedy SP chances to get at least kmin/2 facilities */
        while ((k[0] < kmin) && (i < SP)) {
            cost = pspeedy(points, z, k, pid);
            i++;
        }



        /* if still not enough facilities, assume z is too high */
        while (k[0] < kmin) {
            if (i >= SP) {
                hiz = z;
                z = (hiz + loz) / 2.0;
                i = 0;
            }
            if (pid == 0) {
                shuffle(points);
            }
            cost = pspeedy(points, z, k, pid);
            i++;
        }



        /* now we begin the binary search for real */
        /* must designate some points as feasible centers */
        /* this creates more consistancy between FL runs */
        /* helps to guarantee correct # of centers at the end */

        if (pid == 0) {
            numfeasible = selectfeasibleFast(points, feasible, kmin, pid);
            for (i = 0; i < points.getNum(); i++) {
                isCenterArray[points.getPoints()[i].getAssign()] = true;
            }
        }

        while (true) {
            /* first get a rough estimate on the FL solution */
            lastcost = cost;
            cost = pFL(points, feasible, numfeasible,
                    z, k, cost, (long) (ITER * kmax * Math.log((double) kmax)), 0.1, pid);

            /* if number of centers seems good, try a more accurate FL */
            if (((k[0] <= (1.1) * kmax) && (k[0] >= (0.9) * kmin)) ||
                    ((k[0] <= kmax + 2) && (k[0] >= kmin - 2))) {

      /* may need to run a little longer here before halting without
   improvement */
                cost = pFL(points, feasible, numfeasible,
                        z, k, cost, (long) (ITER * kmax * Math.log((double) kmax)), 0.001, pid);
            }

            if (k[0] > kmax) {
                /* facilities too cheap */
                /* increase facility cost and up the cost accordingly */
                loz = z;
                z = (hiz + loz) / 2.0;
                cost += (z - loz) * k[0];
            }
            if (k[0] < kmin) {
                /* facilities too expensive */
                /* decrease facility cost and reduce the cost accordingly */
                hiz = z;
                z = (hiz + loz) / 2.0;
                cost += (z - hiz) * k[0];
            }

            /* if k is good, return the result */
            /* if we're stuck, just give up and return what we have */

            if (((k[0] <= kmax) && (k[0] >= kmin)) || ((loz >= (0.999) * hiz))) {
                break;
            }
        }

        //clean up...
        if (pid == 0) {
            kfinal[0] = k[0];
        }

        return (float) cost;
    }

    private static double pFL(Points points, int[] feasible, int numfeasible, double z, long[] k, double cost, long iter, double e, int pid) {
        long i;
        long x;
        double change;
        long numberOfPoints;

        change = cost;
        /* continue until we run iter iterations without improvement */
        /* stop instead if improvement is less than e */
        while (change/cost > e) {
            change = 0.0;
            numberOfPoints = points.getNum();
            /* randomize order in which centers are considered */

            if( pid == 0 ) {
                intshuffle(feasible, numfeasible);
            }
            for (i=0;i<iter;i++) {
                x = i%numfeasible;
                change += pgain(feasible[(int) x], points, z, k, pid, false);
            }
            cost -= change;

        }
        return(cost);
    }

    private static double pgain(int x, Points points, double z, long[] numcenters, int pid, boolean isCenter) {

        //my block
        long bsize = points.getNum()/nProc;
        long k1 = bsize * pid;
        long k2 = k1 + bsize;
        if( pid == nProc-1 ) k2 = points.getNum();

        int i;
        int number_of_centers_to_close = 0;

        double [] work_mem;
        double gl_cost_of_opening_x;
        int gl_number_of_centers_to_close;

        //each thread takes a block of working_mem.
        int stride = (int) (numcenters[0]+2);
        //make stride a multiple of CACHE_LINE
        int cl = CACHE_LINE/ (Double.SIZE/Byte.SIZE); //sizeof(double);
        if( stride % cl != 0 ) {
            stride = cl * ( stride / cl + 1);
        }
        int K = stride -2 ; // K==*numcenters

        //my own cost of opening x
        double cost_of_opening_x = 0;

//        if( pid==0 )    {
            work_mem = new double[stride*(nProc+1)];
            gl_cost_of_opening_x = 0;
            gl_number_of_centers_to_close = 0;
//        }

  /*For each center, we have a *lower* field that indicates
    how much we will save by closing the center.
    Each thread has its own copy of the *lower* fields as an array.
    We first build a table to index the positions of the *lower* fields.
  */

        int count = 0;
        for(i = (int) k1; i < k2; i++ ) {
            if( isCenterArray[i] ) {
                centerTable[i] = count++;
            }
        }
        work_mem[pid*stride] = count;


        if( pid == 0 ) {
            int accum = 0;
            for( int p = 0; p < nProc; p++ ) {
                int tmp = (int)work_mem[p*stride];
                work_mem[p*stride] = accum;
                accum += tmp;
            }
        }


        for(i = (int) k1; i < k2; i++ ) {
            if( isCenterArray[i] ) {
                centerTable[i] += (int)work_mem[pid*stride];
            }
        }

        //now we finish building the table. clear the working memory.
        for (i = (int) k1; i < k2-k1; i++) {
            switchMembership[i] = false; //memset(switchMembership + k1, 0, (k2-k1)*sizeof(bool));

        }

        for (i = pid *stride; i < stride; i++) {
            work_mem[i] = 0; // memset(work_mem+pid*stride, 0, stride*sizeof(double));

        }


        if( pid== 0 ) {
            for (i = nProc * stride; i < stride; i++) {
                work_mem[i] = 0; // memset(work_mem + nproc * stride, 0, stride * sizeof( double));
            }
        }


        //my *lower* fields
        double [] lower = new double[work_mem.length]; // = &work_mem[pid*stride];
        System.arraycopy(work_mem, pid*stride, lower, 0, work_mem.length - pid*stride);
        //global *lower* fields
        double [] gl_lower = new double[work_mem.length]; // = &work_mem[nProc*stride];
        System.arraycopy(work_mem, 0, gl_lower, nProc*stride, work_mem.length - nProc*stride);

//printf("----------------------------------------------------\n");
        for (i = (int) k1; i < k2; i++ ) {
//    printf("dim = %d \n" , points->dim);
            double x_cost = dist(points.getPoints()[i], points.getPoints()[x], points.getDim(), isCenter ) * points.getPoints()[i].getWeight();
            double current_cost = points.getPoints()[i].getCost();

            if ( x_cost < current_cost ) {

                // point i would save cost just by switching to x
                // (note that i cannot be a median,
                // or else dist(p[i], p[x]) would be 0)

                switchMembership[i] = true;
                cost_of_opening_x += x_cost - current_cost;

            } else {

                // cost of assigning i to x is at least current assignment cost of i

                // consider the savings that i's **current** median would realize
                // if we reassigned that median and all its members to x;
                // note we've already accounted for the fact that the median
                // would save z by closing; now we have to subtract from the savings
                // the extra cost of reassigning that median and its members
                int assign = points.getPoints()[i].getAssign();
                lower[centerTable[assign]] += current_cost - x_cost;
            }
        }


        // at this time, we can calculate the cost of opening a center
        // at x; if it is negative, we'll go through with opening it

        for (i = (int) k1; i < k2; i++ ) {
            if( isCenterArray[i] ) {
                double low = z;
                //aggregate from all threads
                for( int p = 0; p < nProc; p++ ) {
                    low += work_mem[centerTable[i]+p*stride];
                }
                gl_lower[centerTable[i]] = low;
                if ( low > 0 ) {
                    // i is a median, and
                    // if we were to open x (which we still may not) we'd close i

                    // note, we'll ignore the following quantity unless we do open x
                    ++number_of_centers_to_close;
                    cost_of_opening_x -= low;
                }
            }
        }
        //use the rest of working memory to store the following
        work_mem[pid*stride + K] = number_of_centers_to_close;
        work_mem[pid*stride + K+1] = cost_of_opening_x;


        //  printf("thread %d cost complete\n",pid);

        if( pid==0 ) {
            gl_cost_of_opening_x = z;
            //aggregate
            for( int p = 0; p < nProc; p++ ) {
                gl_number_of_centers_to_close += (int)work_mem[p*stride + K];
                gl_cost_of_opening_x += work_mem[p*stride+K+1];
            }
        }

        // Now, check whether opening x would save cost; if so, do it, and
        // otherwise do nothing

        if ( gl_cost_of_opening_x < 0 ) {
            //  we'd save money by opening x; we'll do it
            for (i = (int) k1; i < k2; i++ ) {
                boolean close_center = gl_lower[centerTable[points.getPoints()[i].getAssign()]] > 0 ;
                if ( switchMembership[i] || close_center ) {
                    // Either i's median (which may be i itself) is closing,
                    // or i is closer to x than to its current median
                    points.getPoints()[i].setCost(
                             (points.getPoints()[i].getWeight() * dist(points.getPoints()[i], points.getPoints()[x], points.getDim(), isCenter))
                    );
                    points.getPoints()[i].setAssign(x);
                }
            }
            for(i = (int) k1; i < k2; i++ ) {
                if( isCenterArray[i] && gl_lower[centerTable[i]] > 0 ) {
                    isCenterArray[i] = false;
                }
            }
            if( x >= k1 && x < k2 ) {
                isCenterArray[x] = true;
            }

            if( pid==0 ) {
                numcenters[0] = numcenters[0] + 1 - gl_number_of_centers_to_close;
            }
        }
        else {
            if( pid==0 )
                gl_cost_of_opening_x = 0;  // the value we'll return
        }

        return -gl_cost_of_opening_x;
    }

    private static void intshuffle(int[] intarray, int length) {
        int i, j;
        int temp;
        for (i=0;i<length;i++) {
            j= (int) ((lrand48()%(length - i))+i);
            temp = intarray[i];
            intarray[i]=intarray[j];
            intarray[j]=temp;
        }
    }

    private static int selectfeasibleFast(Points points, int[] feasible, long kmin, int pid) {
        int numfeasible = points.getNum();
        if (numfeasible > (ITER*kmin*Math.log((double)kmin)))
            numfeasible = (int)(ITER*kmin*Math.log((double)kmin));

//        feasible = new int[numfeasible];

        float [] accumweight;
        float totalweight;

  /*
     Calcuate my block.
     For now this routine does not seem to be the bottleneck, so it is not parallelized.
     When necessary, this can be parallelized by setting k1 and k2 to
     proper values and calling this routine from all threads ( it is called only
     by thread 0 for now ).
     Note that when parallelized, the randomization might not be the same and it might
     not be difficult to measure the parallel speed-up for the whole program.
   */
        //  long bsize = numfeasible;
        long k1 = 0;
        long k2 = numfeasible;

        float w;
        int l,r,k;

        /* not many points, all will be feasible */
        if (numfeasible == points.getNum()) {
            for (int i = (int) k1; i<k2; i++)
                feasible[i] = i;
            return numfeasible;
        }
        accumweight= new float[(int) points.getNum()];


        accumweight[0] = points.getPoints()[0].getWeight();
        totalweight=0;
        for( int i = 1; i < points.getNum(); i++ ) {
            accumweight[i] = accumweight[i-1] + points.getPoints()[i].getWeight();
        }
        totalweight=accumweight[points.getNum() -1];

        for(int i = (int) k1; i<k2; i++ ) {
            w = (lrand48()/(float)Integer.MAX_VALUE)*totalweight;
            //binary search
            l=0;
            r= points.getNum()-1;
            if( accumweight[0] > w )  {
                feasible[i]=0;
                continue;
            }
            while( l+1 < r ) {
                k = (l+r)/2;
                if( accumweight[k] > w ) {
                    r = k;
                }
                else {
                    l=k;
                }
            }
            feasible[i]=r;
        }


        return numfeasible;
    }

    private static double pspeedy(Points points, double z, long[] kcenter, int pid) {
        //my block
        long bsize = points.getNum() / nProc;
        long k1 = bsize * pid;
        long k2 = k1 + bsize;
        if (pid == nProc - 1) k2 = points.getNum();

        double totalcost = 0;
        boolean open = false;
        double[] costs = new double[nProc]; //cost for each thread.
        int i = 0;



        /* create center at first point, send it to itself */
        for (int k = (int) k1; k < k2; k++) {
            double distance = dist(points.getPoints()[k], points.getPoints()[0], points.getDim(), false);
            points.getPoints()[k].setCost((float) (distance * points.getPoints()[k].getWeight()));
            points.getPoints()[k].setAssign(0);
        }

        if (pid == 0) {
            kcenter[0] = 1;
            costs = new double[nProc];
        }


        if (pid != 0) { // we are not the master threads. we wait until a center is opened.
            while (true) {
                if (i >= points.getNum()) {
                    break;
                }

                for (int k = (int) k1; k < k2; k++) {
                    double distance = dist(points.getPoints()[i], points.getPoints()[k], points.getDim(), false);
                    if (distance * points.getPoints()[k].getWeight() < points.getPoints()[k].getCost()) {
                        points.getPoints()[k].setCost((float) distance * points.getPoints()[k].getWeight());
                        points.getPoints()[k].setAssign(i);
                    }
                }
            }
        } else { // I am the master thread. I decide whether to open a center and notify others if so.
            for (i = 1; i < points.getNum(); i++) {
                boolean to_open = ((float) lrand48() / (float) Integer.MAX_VALUE) < (points.getPoints()[i].getCost() / z);
                if (to_open) {
                    kcenter[0] += 1;

                    open = true;

                    for (int k = (int) k1; k < k2; k++) {
                        double distance = dist(points.getPoints()[i], points.getPoints()[k], points.getDim(), false);
                        if (distance * points.getPoints()[k].getWeight() < points.getPoints()[k].getCost()) {
                            points.getPoints()[k].setCost((float) (distance * points.getPoints()[k].getWeight()));
                            points.getPoints()[k].setAssign(i);
                        }
                    }

                    open = false;
                }
            }
            open = true;
        }

        open = false;
        double mytotal = 0;
        for (int k = (int) k1; k < k2; k++) {
            mytotal += points.getPoints()[k].getCost();
        }
        costs[pid] = mytotal;

        // aggregate costs from each thread
        if (pid == 0) {
            totalcost = z * (kcenter[0]);
            for (i = 0; i < nProc; i++) {
                totalcost += costs[i];
            }
        }

        return (totalcost);
    }

    private static void shuffle(Points points) {
        int i, j;
        Point temp;
        for (i = 0; i < points.getNum() - 1; i++) {
            j = (int) ((lrand48() % (points.getNum() - i)) + i);
            temp = points.getPoints()[i];
            points.getPoints()[i] = points.getPoints()[j];
            points.getPoints()[j] = temp;
        }
    }

    private static double dist(Point p1, Point p2, int dim, boolean isCenter) {

        int i;
        float result = 0.0F;

        if (useVectorAPI) {

            float [] coord;
            if (isCenter) {
                coord = centerBlock;
            } else {
                coord = block;
            }

            int limit = SPECIES.loopBound(dim);

            FloatVector result1,result2, _aux, _diff, _coord1, _coord2;

            result1 = FloatVector.zero(SPECIES); // _MM_SET_f32(0.0,gvl);
//            result2 = FloatVector.zero(SPECIES); // _MM_SET_f32(0.0,gvl);

            int startIndex1 = p1.getCoordIndex();
            int startIndex2 = p2.getCoordIndex();

            for (i=0;i<limit;i=i+SPECIES_LENGTH) {

//                gvl = __builtin_epi_vsetvl(dim-i, __epi_e32, __epi_m1);
                _coord1 = FloatVector.fromArray(SPECIES, coord, startIndex1 + i); // _MM_LOAD_f32(&(p1.coord[i]),gvl);
                _coord2 = FloatVector.fromArray(SPECIES, coord, startIndex2 + i); // _MM_LOAD_f32(&(p2.coord[i]),gvl);

                _diff = _coord2.sub(_coord1); // _MM_SUB_f32(_coord2,_coord1,gvl);
                result1  = _diff.pow(2); // // _MM_MACC_f32(result1,_diff,_diff,gvl);
            }
            // result2 = _MM_REDSUM_f32(result1,result2,gvl);
            result = result1.reduceLanes(VectorOperators.ADD); // _MM_VGETFIRST_f32(result2);

            for (; i < dim; i++)
                result += (p1.getCoord(i, isCenter) - p2.getCoord(i, isCenter)) * (p1.getCoord(i, isCenter) - p2.getCoord(i, isCenter));



            return (result);

        } else {

            for (i = 0; i < dim; i++)
                result += (p1.getCoord(i, isCenter) - p2.getCoord(i, isCenter)) * (p1.getCoord(i, isCenter) - p2.getCoord(i, isCenter));
            //printf("result = %f \n",result);
            return (result);
        }
    }

}
