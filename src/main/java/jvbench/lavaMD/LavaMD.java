package jvbench.lavaMD;


import jdk.incubator.vector.FloatVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

import java.io.*;
import java.util.Random;


public class LavaMD {

    private static int NUMBER_PAR_PER_BOX = 96;
    private static final Random random = new Random(42);
    private static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_MAX;
    private static final int SPECIES_LENGTH = SPECIES.length();
    private static final ParStr parCPU = new ParStr();
    private static final DimStr dimCPU = new DimStr();
    private static boolean USE_VECTOR_API = false;
    private static BoxStr[] boxCPU;
    private static FourVector[] rvCPU;
    private static float[] qvCPU;
    protected static FourVector[] fvCPU;
    protected static FourVector[] fvCPUV;


    public static void init(String inputFilename) {
        // counters
        int i, j, k, l, m, n;

        int nh;

        try (BufferedReader reader = new BufferedReader(new FileReader(inputFilename))) {

            String [] params = reader.readLine().split(" "); // read param
            dimCPU.cores_arg = Integer.parseInt(params[0]); // not used
            dimCPU.boxes1d_arg = Integer.parseInt(params[1]);
            NUMBER_PAR_PER_BOX = Integer.parseInt(params[2]);

            parCPU.alpha = 0.5f;

            // total number of boxes
            dimCPU.number_boxes = dimCPU.boxes1d_arg * dimCPU.boxes1d_arg * dimCPU.boxes1d_arg;

            // how many particles space has in each direction
            dimCPU.space_elem = dimCPU.number_boxes * NUMBER_PAR_PER_BOX;
            dimCPU.space_mem =  dimCPU.space_elem; // * sizeof(FOUR_VECTOR);
            dimCPU.space_mem2 = dimCPU.space_elem; // * sizeof(fp);

            // box array
            dimCPU.box_mem = dimCPU.number_boxes; // * sizeof(box_str);

            //======================================================================================================================================================150
            //	SYSTEM MEMORY
            //======================================================================================================================================================150

            //====================================================================================================100
            //	BOX
            //====================================================================================================100

            // allocate boxes
            boxCPU = new BoxStr[dimCPU.box_mem];

            // initialize number of home boxes
            nh = 0;


            // home boxes in z direction
            for(i=0; i< dimCPU.boxes1d_arg; i++){
                // home boxes in y direction
                for(j=0; j< dimCPU.boxes1d_arg; j++){
                    // home boxes in x direction
                    for(k=0; k< dimCPU.boxes1d_arg; k++){

                        // current home box
                        boxCPU[nh] = new BoxStr();
                        boxCPU[nh].x = k;
                        boxCPU[nh].y = j;
                        boxCPU[nh].z = i;
                        boxCPU[nh].number = nh;
                        boxCPU[nh].offset = nh * NUMBER_PAR_PER_BOX;

                        // initialize number of neighbor boxes
                        boxCPU[nh].nn = 0;

                        // neighbor boxes in z direction
                        for(l=-1; l<2; l++){
                            // neighbor boxes in y direction
                            for(m=-1; m<2; m++){
                                // neighbor boxes in x direction
                                for(n=-1; n<2; n++){

                                    // check if (this neighbor exists) and (it is not the same as home box)
                                    if(		(((i + l) >= 0 && (j + m) >= 0 && (k + n) >= 0) && ((i + l) < dimCPU.boxes1d_arg && (j + m) < dimCPU.boxes1d_arg && (k + n) < dimCPU.boxes1d_arg))	&&
                                            !(l == 0 && m == 0 && n == 0)){ // Simplified booleans

                                        // current neighbor box
                                        boxCPU[nh].nei[boxCPU[nh].nn] = new NeiStr();
                                        boxCPU[nh].nei[boxCPU[nh].nn].x = (k+n);
                                        boxCPU[nh].nei[boxCPU[nh].nn].y = (j+m);
                                        boxCPU[nh].nei[boxCPU[nh].nn].z = (i+l);
                                        boxCPU[nh].nei[boxCPU[nh].nn].number =	(boxCPU[nh].nei[boxCPU[nh].nn].z * dimCPU.boxes1d_arg * dimCPU.boxes1d_arg) +
                                                (boxCPU[nh].nei[boxCPU[nh].nn].y * dimCPU.boxes1d_arg) +
                                                boxCPU[nh].nei[boxCPU[nh].nn].x;
                                        boxCPU[nh].nei[boxCPU[nh].nn].offset = boxCPU[nh].nei[boxCPU[nh].nn].number * NUMBER_PAR_PER_BOX;

                                        // increment neighbor box
                                        boxCPU[nh].nn = boxCPU[nh].nn + 1;

                                    }

                                } // neighbor boxes in x direction
                            } // neighbor boxes in y direction
                        } // neighbor boxes in z direction

                        // increment home box
                        nh = nh + 1;

                    } // home boxes in x direction
                } // home boxes in y direction
            } // home boxes in z direction


            //====================================================================================================100
            //	PARAMETERS, DISTANCE, CHARGE AND FORCE
            //====================================================================================================100

            // random generator seed set to random value - time in this case
//        srand(time(NULL));

            // input (distances)
            rvCPU = new FourVector[dimCPU.space_mem]; // (FOUR_VECTOR*)malloc(dim_cpu.space_mem);
            qvCPU = new float[dimCPU.space_mem2];
            for(i=0; i< dimCPU.space_elem; i=i+1){
                String[] values = reader.readLine().split(" ");
                rvCPU[i] = new FourVector();
                rvCPU[i].v = Float.parseFloat(values[0]);			// get a number in the range 0.1 - 1.0
                rvCPU[i].x = Float.parseFloat(values[1]);			// get a number in the range 0.1 - 1.0
                rvCPU[i].y = Float.parseFloat(values[2]);			// get a number in the range 0.1 - 1.0
                rvCPU[i].z = Float.parseFloat(values[3]);			// get a number in the range 0.1 - 1.0
                qvCPU[i] = Float.parseFloat(values[4]);
            }


            // output (forces)
            fvCPU = new FourVector[dimCPU.space_mem]; // (FOUR_VECTOR*)malloc(dim_cpu.space_mem);
            fvCPUV = new FourVector[dimCPU.space_mem]; // (FOUR_VECTOR*)malloc(dim_cpu.space_mem);
            for(i=0; i< dimCPU.space_elem; i=i+1){
                fvCPU[i] = new FourVector();
                fvCPU[i].v = 0;								// set to 0, because kernels keeps adding to initial value
                fvCPU[i].x = 0;								// set to 0, because kernels keeps adding to initial value
                fvCPU[i].y = 0;								// set to 0, because kernels keeps adding to initial value
                fvCPU[i].z = 0;								// set to 0, because kernels keeps adding to initial value
                fvCPUV[i] = new FourVector();
                fvCPUV[i].v = 0;
                fvCPUV[i].x = 0;
                fvCPUV[i].y = 0;
                fvCPUV[i].z = 0;
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }




    }

    public static void serial() {
        kernelCPU(parCPU, dimCPU, boxCPU, rvCPU, qvCPU, fvCPU);
    }

    public static void vector() {
        kernelCPUVector(parCPU, dimCPU, boxCPU, rvCPU, qvCPU, fvCPUV);
    }


    private static void kernelCPU(ParStr par, DimStr dim, BoxStr[] box, FourVector[] rv, float[] qv, FourVector[] fv) {
        // parameters
        float alpha;
        float a2;

        // counters
        int i, j, k, l;

        // home box
        long first_i;
        int rAindex;
        int fAindex;
        FourVector [] rA; // FOUR_VECTOR* rA;
        FourVector [] fA; // FOUR_VECTOR* fA;

        // neighbor box
        int pointer;
        long first_j;
        FourVector [] rB; // FOUR_VECTOR* rB;
        float [] qB;

        // common
        float r2;
        float u2;
        float fs;
        float vij;
        float fxij,fyij,fzij;
        ThreeVector d = new ThreeVector(); //THREE_VECTOR d;

//        time1 = get_time();

        //======================================================================================================================================================150
        //	MCPU SETUP
        //======================================================================================================================================================150

        //omp_set_num_threads(dim.cores_arg);

//        time2 = get_time();

        //======================================================================================================================================================150
        //	INPUTS
        //======================================================================================================================================================150

        alpha = par.alpha;
        a2 = 2.0f*alpha*alpha;

//        time3 = get_time();

        //======================================================================================================================================================150
        //	PROCESS INTERACTIONS
        //======================================================================================================================================================150

        //#pragma omp	parallel for \
        //			private(i, j, k) \
        //			private(first_i, rA, fA) \
        //			private(pointer, first_j, rB, qB) \
        //			private(r2, u2, fs, vij, fxij, fyij, fzij, d)
        for(l=0; l<dim.number_boxes; l=l+1){

            //------------------------------------------------------------------------------------------100
            //	home box - box parameters
            //------------------------------------------------------------------------------------------100

            first_i = box[l].offset;												// offset to common arrays

            //------------------------------------------------------------------------------------------100
            //	home box - distance, force, charge and type parameters from common arrays
            //------------------------------------------------------------------------------------------100

//            rA = &rv[first_i];
//            fA = &fv[first_i];
            rAindex = (int) first_i;
            fAindex = (int) first_i;

            //------------------------------------------------------------------------------------------100
            //	Do for the # of (home+neighbor) boxes
            //------------------------------------------------------------------------------------------100
            for (k=0; k<(1+box[l].nn); k++)
            {

                //----------------------------------------50
                //	neighbor box - get pointer to the right box
                //----------------------------------------50

                if(k==0){
                    pointer = l;													// set first box to be processed to home box
                }
                else{
                    pointer = box[l].nei[k-1].number;							// remaining boxes are neighbor boxes
                }

                //----------------------------------------50
                //	neighbor box - box parameters
                //----------------------------------------50

                first_j = box[pointer].offset;

                //----------------------------------------50
                //	neighbor box - distance, force, charge and type parameters
                //----------------------------------------50

//                rB = &rv[first_j];
//                qB = &qv[first_j];
                int rBindex = (int) first_j;
                int qBindex = (int) first_j;

                //----------------------------------------50
                //	Do for the # of particles in home box
                //----------------------------------------50
                for (i=0; i<NUMBER_PAR_PER_BOX; i=i+1){

                    // do for the # of particles in current (home or neighbor) box
                    for (j=0; j<NUMBER_PAR_PER_BOX; j=j+1){

                        // // coefficients
//                        r2 = rA[i].v + rB[j].v - DOT(rA[i],rB[j]);
                        r2 = rv[i + rAindex].v + rv[j + rBindex].v - DOT(rv[i + rAindex],rv[j + rBindex]);
                        u2 = a2*r2;
                        vij= (float) Math.exp(-u2);
                        fs = 2.F*vij;
//                        d.x = rA[i].x  - rB[j].x;
                        d.x = rv[i + rAindex].x  - rv[j + rBindex].x;
//                        d.y = rA[i].y  - rB[j].y;
                        d.y = rv[i + rAindex].y  - rv[j + rBindex].y;
//                        d.z = rA[i].z  - rB[j].z;
                        d.z = rv[i + rAindex].z  - rv[j + rBindex].z;
                        fxij=fs*d.x;
                        fyij=fs*d.y;
                        fzij=fs*d.z;

                        // forces
//                        fA[i].v +=  qB[j]*vij;
                        fv[i + fAindex].v +=  qv[j + qBindex]*vij;
                        fv[i + fAindex].x +=  qv[j+ qBindex]*fxij;
                        fv[i + fAindex].y +=  qv[j+ qBindex]*fyij;
                        fv[i + fAindex].z +=  qv[j+ qBindex]*fzij;

                    } // for j

                } // for i

            } // for k

        } // for l
    }

    private static void kernelCPUVector(ParStr par, DimStr dim, BoxStr[] box, FourVector[] rv, float[] qv, FourVector[] fv) {

        // parameters
        float alpha;
        float a2;

        // counters
        int i, j, k, l;

        // home box
        int first_i;
        FourVector[] rA;
        FourVector[] fA;
        int rAindex;
        int fAindex;

        // neighbor box
        int pointer;
        int first_j;
        FourVector[] rB;


        int rBindex;
        float [] qB;
        int qBindex;

        // common
        float r2;
        float u2;
        float fs;
        float vij;
        float fxij,fyij,fzij;
        ThreeVector d = new ThreeVector();

//        time1 = get_time();

        //======================================================================================================================================================150
        //	MCPU SETUP
        //======================================================================================================================================================150

        //omp_set_num_threads(dim.cores_arg);

//        time2 = get_time();

        //======================================================================================================================================================150
        //	INPUTS
        //======================================================================================================================================================150

        alpha = par.alpha;
        a2 = 2.0f*alpha*alpha;

//        time3 = get_time();

        //======================================================================================================================================================150
        //	PROCESS INTERACTIONS
        //======================================================================================================================================================150

        //#pragma omp	parallel for \
        //			private(i, j, k) \
        //			private(first_i, rA, fA) \
        //			private(pointer, first_j, rB, qB) \
        //			private(r2, u2, fs, vij, fxij, fyij, fzij, d)
        for(l=0; l<dim.number_boxes; l=l+1){

            //------------------------------------------------------------------------------------------100
            //	home box - box parameters
            //------------------------------------------------------------------------------------------100

            first_i = box[l].offset;												// offset to common arrays

            //------------------------------------------------------------------------------------------100
            //	home box - distance, force, charge and type parameters from common arrays
            //------------------------------------------------------------------------------------------100

            //rA = &rv[first_i];
            //fA = &fv[first_i];
            rAindex = first_i;
            fAindex = first_i;
            //------------------------------------------------------------------------------------------100
            //	Do for the # of (home+neighbor) boxes
            //------------------------------------------------------------------------------------------100
            for (k=0; k<(1+box[l].nn); k++)
            {

                //----------------------------------------50
                //	neighbor box - get pointer to the right box
                //----------------------------------------50

                if(k==0){
                    pointer = l;													// set first box to be processed to home box
                }
                else{
                    pointer = box[l].nei[k-1].number;							// remaining boxes are neighbor boxes
                }

                //----------------------------------------50
                //	neighbor box - box parameters
                //----------------------------------------50

                first_j = box[pointer].offset;

                //----------------------------------------50
                //	neighbor box - distance, force, charge and type parameters
                //----------------------------------------50

//                rB = &rv[first_j];
                rBindex = first_j;
//                qB = &qv[first_j];
                qBindex = first_j;

                float [] rBx = new float[NUMBER_PAR_PER_BOX];
                float [] rBy = new float[NUMBER_PAR_PER_BOX];
                float [] rBz = new float[NUMBER_PAR_PER_BOX];
                float [] rBv = new float[NUMBER_PAR_PER_BOX];

                for (i = 0; i < NUMBER_PAR_PER_BOX; i++) {
                    rBv[i] = rv[i + rBindex].v;
                    rBx[i] = rv[i + rBindex].x;
                    rBy[i] = rv[i + rBindex].y;
                    rBz[i] = rv[i + rBindex].z;
                }

                //----------------------------------------50
                //	Do for the # of particles in home box
                //----------------------------------------50
                for (i=0; i<NUMBER_PAR_PER_BOX; i=i+1){

//                    unsigned long int gvl = __builtin_epi_vsetvl(NUMBER_PAR_PER_BOX, __epi_e32, __epi_m1);
                    int limit = SPECIES.loopBound(NUMBER_PAR_PER_BOX);

                    FloatVector xr2;
                    FloatVector xDOT;
                    FloatVector xu2;
                    FloatVector xa2 = FloatVector.broadcast(SPECIES, a2);// _MM_SET_f32(a2,gvl);
                    FloatVector xvij;
                    FloatVector xrA_v = FloatVector.broadcast(SPECIES, rv[i + rAindex].v); //_MM_SET_f32(rA[i].v,gvl);
                    FloatVector xrA_x = FloatVector.broadcast(SPECIES, rv[i + rAindex].x); //_MM_SET_f32(rA[i].x,gvl);
                    FloatVector xrA_y = FloatVector.broadcast(SPECIES, rv[i + rAindex].y); //_MM_SET_f32(rA[i].y,gvl);
                    FloatVector xrA_z = FloatVector.broadcast(SPECIES, rv[i + rAindex].z ); //_MM_SET_f32(rA[i].z,gvl);
                    FloatVector xrB_v;
                    FloatVector xrB_x;
                    FloatVector xrB_y;
                    FloatVector xrB_z;
                    FloatVector xd_x;
                    FloatVector xd_y;
                    FloatVector xd_z;
                    FloatVector xfxij;
                    FloatVector xfyij;
                    FloatVector xfzij;
                    FloatVector xfs;
                    FloatVector xqB;
                    FloatVector xfA_v = FloatVector.zero(SPECIES); // _MM_SET_f32(0.0,gvl);
                    FloatVector xfA_x = FloatVector.zero(SPECIES); // _MM_SET_f32(0.0,gvl);
                    FloatVector xfA_y = FloatVector.zero(SPECIES); // _MM_SET_f32(0.0,gvl);
                    FloatVector xfA_z = FloatVector.zero(SPECIES); // _MM_SET_f32(0.0,gvl);
                    float xfA_1_v = 0.0f; // _MM_SET_f32(0.0,1);
                    float xfA_1_x = 0.0f; // _MM_SET_f32(0.0,1);
                    float xfA_1_y = 0.0f; // _MM_SET_f32(0.0,1);
                    float xfA_1_z = 0.0f; // _MM_SET_f32(0.0,1);

                    // do for the # of particles in current (home or neighbor) box
                    for (j=0; j < limit; j+=SPECIES_LENGTH){
//                        gvl = __builtin_epi_vsetvl(NUMBER_PAR_PER_BOX-j, __epi_e32, __epi_m1);
                        // coefficients
//                        System.out.println("rBv length " + rBv.length + " index j + rBindex = " + j + " + " + rBindex + " = " + (j + rBindex));
                        xrB_v = FloatVector.fromArray(SPECIES, rBv, j); // _MM_LOAD_STRIDE_f32(&rB[j].v,4,gvl);
                        xrB_x = FloatVector.fromArray(SPECIES, rBx, j); // _MM_LOAD_STRIDE_f32(&rB[j].x,4,gvl);
                        xrB_y = FloatVector.fromArray(SPECIES, rBy, j); // _MM_LOAD_STRIDE_f32(&rB[j].y,4,gvl);
                        xrB_z = FloatVector.fromArray(SPECIES, rBz, j); // _MM_LOAD_STRIDE_f32(&rB[j].z,4,gvl);
                        //r2 = rA[i].v + rB[j].v - DOT(rA[i],rB[j]);
                        xr2    = xrA_v.add(xrB_v); // _MM_ADD_f32(xrA_v, xrB_v,gvl);
                        xDOT   = xrA_x.mul(xrB_x); //_MM_MUL_f32(xrA_x, xrB_x,gvl);
                        xDOT   = xrA_y.mul(xrB_y).add(xDOT); // _MM_MACC_f32(xDOT,xrA_y,xrB_y,gvl);
                        xDOT   = xrA_z.mul(xrB_z).add(xDOT); // _MM_MACC_f32(xDOT,xrA_z,xrB_z,gvl);
                        xr2    = xr2.sub(xDOT); // _MM_SUB_f32(xr2, xDOT,gvl);
                        //u2 = a2*r2;
                        xu2    = xa2.mul(xr2); // _MM_MUL_f32(xa2, xr2,gvl);
                        //vij= exp(-u2);
                        xvij   = xu2.lanewise(VectorOperators.NEG).lanewise(VectorOperators.EXP); //  _MM_EXP_f32(_MM_VFSGNJN_f32(xu2,xu2,gvl),gvl);
                        //fs = 2.*vij;
                        xfs    = xvij.mul(2.0f); // _MM_MUL_f32(_MM_SET_f32(2.0f,gvl), xvij,gvl);
                        //d.x = rA[i].x  - rB[j].x;
                        xd_x   = xrA_x.sub(xrB_x); // _MM_SUB_f32(xrA_x, xrB_x,gvl);
                        //d.y = rA[i].y  - rB[j].y;
                        xd_y   = xrA_y.sub(xrB_y); // _MM_SUB_f32(xrA_y, xrB_y,gvl);
                        //d.z = rA[i].z  - rB[j].z;
                        xd_z   = xrA_z.sub(xrB_z); // _MM_SUB_f32(xrA_z, xrB_z,gvl);
                        //fxij=fs*d.x;
                        xfxij  = xfs.mul(xd_x); // _MM_MUL_f32(xfs, xd_x,gvl);
                        //fyij=fs*d.y;
                        xfyij  = xfs.mul(xd_y); // _MM_MUL_f32(xfs, xd_y,gvl);
                        //fzij=fs*d.z;
                        xfzij  = xfs.mul(xd_z); //_MM_MUL_f32(xfs, xd_z,gvl);

                        // forces
                        //fA[i].v +=  qB[j]*vij;
                        //fA[i].x +=  qB[j]*fxij;
                        //fA[i].y +=  qB[j]*fyij;
                        //fA[i].z +=  qB[j]*fzij;
//                        gvl = __builtin_epi_vsetvl(NUMBER_PAR_PER_BOX, __epi_e32, __epi_m1);
                        xqB = FloatVector.fromArray(SPECIES, qv, j + qBindex); // _MM_LOAD_f32(&qB[j],gvl);
                        xfA_v   = xqB.mul(xvij).add(xfA_v); // _MM_MACC_f32(xfA_v,xqB,xvij,gvl);
                        xfA_x   = xqB.mul(xfxij).add(xfA_x); // _MM_MACC_f32(xfA_x,xqB,xfxij,gvl);
                        xfA_y   = xqB.mul(xfyij).add(xfA_y); // _MM_MACC_f32(xfA_y,xqB,xfyij,gvl);
                        xfA_z   = xqB.mul(xfzij).add(xfA_z); // _MM_MACC_f32(xfA_z,xqB,xfzij,gvl);
                    } // for j

                    for (; j < NUMBER_PAR_PER_BOX; j++) {

                        // // coefficients
//                        r2 = rA[i].v + rB[j].v - DOT(rA[i],rB[j]);
                        r2 = rv[i + rAindex].v + rv[j+rBindex].v - DOT(rv[i + rAindex],rv[j+rBindex]);
                        u2 = a2*r2;
                        vij= (float) Math.exp(-u2);
                        fs = 2.f*vij;
//                        d.x = rA[i].x  - rB[j].x;
                        d.x = rv[i + rAindex].x  - rv[j+rBindex].x;
//                        d.y = rA[i].y  - rB[j].y;
                        d.y = rv[i + rAindex].y  - rv[j+rBindex].y;
//                        d.z = rA[i].z  - rB[j].z;
                        d.z = rv[i + rAindex].z  - rv[j+rBindex].z;
                        fxij=fs*d.x;
                        fyij=fs*d.y;
                        fzij=fs*d.z;

                        // forces
//                        fA[i].v +=  qB[j]*vij;
                        xfA_1_v +=  qv[j+qBindex]*vij;
                        xfA_1_x +=  qv[j+qBindex]*fxij;
                        xfA_1_y +=  qv[j+qBindex]*fyij;
                        xfA_1_z +=  qv[j+qBindex]*fzij;

                    }

//                    gvl = __builtin_epi_vsetvl(NUMBER_PAR_PER_BOX, __epi_e32, __epi_m1);

                    xfA_1_v  += xfA_v.reduceLanes(VectorOperators.ADD); // _MM_REDSUM_f32(xfA_v,xfA_1_v,gvl);
                    xfA_1_x  += xfA_x.reduceLanes(VectorOperators.ADD); // _MM_REDSUM_f32(xfA_x,xfA_1_x,gvl);
                    xfA_1_y  += xfA_y.reduceLanes(VectorOperators.ADD); // _MM_REDSUM_f32(xfA_y,xfA_1_y,gvl);
                    xfA_1_z  += xfA_z.reduceLanes(VectorOperators.ADD); // _MM_REDSUM_f32(xfA_z,xfA_1_z,gvl);
                    fv[i + fAindex].v += xfA_1_v; // _MM_STORE_f32(&fA[i].v, xfA_1_v,1);
                    fv[i + fAindex].x += xfA_1_x; // _MM_STORE_f32(&fA[i].x, xfA_1_x,1);
                    fv[i + fAindex].y += xfA_1_y; // _MM_STORE_f32(&fA[i].y, xfA_1_y,1);
                    fv[i + fAindex].z += xfA_1_z; // _MM_STORE_f32(&fA[i].z, xfA_1_z,1);
                } // for i

            } // for k

        } // for l

    }

    private static float DOT(FourVector A, FourVector B) {
        return ((A.x) * (B.x)) + ((A.y)*(B.y)) + ((A.z) * (B.z));
    }

}
