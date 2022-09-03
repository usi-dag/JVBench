package jvbench.particlefilterLoadMaskInsideLoop;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class ParticleFilterTest {

    int IszX =8;
    int IszY = 8;
    int Nfr= 3;
    int Nparticles= 4;
    double[] seed;
    int [] I;

    double [] randuVectorResult = new double [ParticleFilter.DOUBLE_SPECIES_LENGTH];
    double [] randuVectorNum = new double[ParticleFilter.DOUBLE_SPECIES_LENGTH];

    public void setUp() {
        int i;
        I = new int [IszX * IszY * Nfr];
        seed = new double[Nparticles];
        for (i = 0; i < Nparticles; i++) {
            seed[i] = i;
        }

        ParticleFilter.videoSequence(I, IszX, IszY, Nfr, seed);
    }

    @Test
    void loadMaskInsideLoop() {
        setUp();
        ParticleFilter.particleFilter(
                I,
                IszX,
                IszY,
                Nfr,
                seed,
                Nparticles
        );
        setUp();
       ParticleFilter.particleFilterVector(
                I,
                IszX,
                IszY,
                Nfr,
                seed,
                randuVectorResult,
                randuVectorNum,
                Nparticles
        );


        for(int i = 0; i < Nparticles; i++){
            assertEquals(ParticleFilter.resX[i], ParticleFilter.resXV[i]);
            assertEquals(ParticleFilter.resY[i], ParticleFilter.resYV[i]);
            assertEquals(ParticleFilter.resW[i], ParticleFilter.resWV[i]);

        }
    }

}