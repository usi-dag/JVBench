package jvbench.particlefilter;

import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

import static jvbench.particlefilter.ParticleFilter.DOUBLE_SPECIES_LENGTH;

public class ParticlefilterBenchmark {


    @State(Scope.Thread)
    public static class MyState {

        int IszX = Integer.parseInt(System.getProperty("x", "128"));
        int IszY = Integer.parseInt(System.getProperty("y", "128"));
        int Nfr= Integer.parseInt(System.getProperty("z", "24"));
        int Nparticles= Integer.parseInt(System.getProperty("np", "32768"));
        double[] seed;
        int [] I;

        double [] randuVectorResult = new double [DOUBLE_SPECIES_LENGTH];
        double [] randuVectorNum = new double[DOUBLE_SPECIES_LENGTH];


        @Setup(Level.Trial)
        public void setup() {
            int i;
            I = new int [IszX * IszY * Nfr];
            seed = new double[Nparticles];
            for (i = 0; i < Nparticles; i++) {
                seed[i] = i; // (int) /*System.currentTimeMillis() * i; // */ Instant.now().getEpochSecond() * i;
            }

            ParticleFilter.videoSequence(I, IszX, IszY, Nfr, seed);
        }

        @TearDown(Level.Iteration)
        public void reset() {
            int i;
            for (i = 0; i < Nparticles; i++) {
                seed[i] = i; // (int) /*System.currentTimeMillis() * i; // */ Instant.now().getEpochSecond() * i;
            }
        }
     }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void serial(MyState state) {
        ParticleFilter.particleFilter(
                state.I,
                state.IszX,
                state.IszY,
                state.Nfr,
                state.seed,
                state.Nparticles
        );
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void autoVec(MyState state) {
        ParticleFilter.particleFilter(
                state.I,
                state.IszX,
                state.IszY,
                state.Nfr,
                state.seed,
                state.Nparticles
        );
    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void explicitVec(MyState state) {
        ParticleFilter.particleFilterVector(
                state.I,
                state.IszX,
                state.IszY,
                state.Nfr,
                state.seed,
                state.randuVectorResult,
                state.randuVectorNum,
                state.Nparticles
        );
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void fullVec(MyState state) {
        ParticleFilter.particleFilterVector(
                state.I,
                state.IszX,
                state.IszY,
                state.Nfr,
                state.seed,
                state.randuVectorResult,
                state.randuVectorNum,
                state.Nparticles
        );
    }

}
