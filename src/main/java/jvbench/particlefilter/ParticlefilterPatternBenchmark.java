package jvbench.particlefilter;

import jvbench.particlefilterUpdateMaskInsideLoop.ParticleFilter;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

import static jvbench.particlefilter.ParticleFilter.DOUBLE_SPECIES_LENGTH;


public class ParticlefilterPatternBenchmark {


    @State(Scope.Thread)
    public static class XorState {

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

            jvbench.particlefilterXor.ParticleFilter.videoSequence(I, IszX, IszY, Nfr, seed);
        }

        @TearDown(Level.Iteration)
        public void reset() {
            int i;
            for (i = 0; i < Nparticles; i++) {
                seed[i] = i; // (int) /*System.currentTimeMillis() * i; // */ Instant.now().getEpochSecond() * i;
            }
        }
    }

    @State(Scope.Thread)
    public static class MergeWithBlendState {

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

            jvbench.particlefilterMergeWithBlend.ParticleFilter.videoSequence(I, IszX, IszY, Nfr, seed);
        }

        @TearDown(Level.Iteration)
        public void reset() {
            int i;
            for (i = 0; i < Nparticles; i++) {
                seed[i] = i; // (int) /*System.currentTimeMillis() * i; // */ Instant.now().getEpochSecond() * i;
            }
        }
    }

    @State(Scope.Thread)
    public static class UpdateMaskInsideLoopState {

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

    @State(Scope.Thread)
    public static class LoadMaskInsideLoopState {

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

            jvbench.particlefilterLoadMaskInsideLoop.ParticleFilter.videoSequence(I, IszX, IszY, Nfr, seed);
        }

        @TearDown(Level.Iteration)
        public void reset() {
            int i;
            for (i = 0; i < Nparticles; i++) {
                seed[i] = i; // (int) /*System.currentTimeMillis() * i; // */ Instant.now().getEpochSecond() * i;
            }
        }
    }

    @State(Scope.Thread)
    public static class IndexInRangeState {

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

            jvbench.particlefilterIndexInRange.ParticleFilter.videoSequence(I, IszX, IszY, Nfr, seed);
        }

        @TearDown(Level.Iteration)
        public void reset() {
            int i;
            for (i = 0; i < Nparticles; i++) {
                seed[i] = i; // (int) /*System.currentTimeMillis() * i; // */ Instant.now().getEpochSecond() * i;
            }
        }
    }

    @State(Scope.Thread)
    public static class NoSecondMaskState {

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

            jvbench.particlefilterNoSecondMask.ParticleFilter.videoSequence(I, IszX, IszY, Nfr, seed);
        }

        @TearDown(Level.Iteration)
        public void reset() {
            int i;
            for (i = 0; i < Nparticles; i++) {
                seed[i] = i; // (int) /*System.currentTimeMillis() * i; // */ Instant.now().getEpochSecond() * i;
            }
        }
    }

    @State(Scope.Thread)
    public static class StaticMaskState {

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

            jvbench.particlefilterStatic.ParticleFilter.videoSequence(I, IszX, IszY, Nfr, seed);
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
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void xorExtended(XorState state) { // xor extended
        jvbench.particlefilterXor.ParticleFilter.particleFilterVector(
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
    public void mergeWithBlend(MergeWithBlendState state) {
        jvbench.particlefilterMergeWithBlend.ParticleFilter.particleFilterVector(
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
    public void updateMaskInsideLoop(UpdateMaskInsideLoopState state) {
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
    public void loadMaskInsideLoop(LoadMaskInsideLoopState state) {
        jvbench.particlefilterLoadMaskInsideLoop.ParticleFilter.particleFilterVector(
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
    public void indexInRange(IndexInRangeState state) {
        jvbench.particlefilterIndexInRange.ParticleFilter.particleFilterVector(
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
    public void noSecondMask(NoSecondMaskState state) {
        jvbench.particlefilterNoSecondMask.ParticleFilter.particleFilterVector(
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
    public void staticMask(StaticMaskState state) {
        jvbench.particlefilterStatic.ParticleFilter.particleFilterVector(
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
