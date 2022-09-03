package jvbench.canneal;

import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.infra.Blackhole;

import java.util.concurrent.TimeUnit;

public class CannealPatternBenchmark {

    @State(Scope.Thread)
    public static class MinSizeCheckState {

       jvbench.cannealMinSizeCheck.Canneal canneal; // size

        @Setup(Level.Trial)
        public void setup() {
                    canneal = new jvbench.cannealMinSizeCheck.Canneal(
                            Integer.parseInt(System.getProperty("N_THREADS","1")),
                            Integer.parseInt(System.getProperty("N_SWAPS","10000")),
                            Integer.parseInt(System.getProperty("TEMP","2000")),
                            System.getProperty("NETLIST","src/main/resources/canneal/input/2500000.nets"),
                            Integer.parseInt(System.getProperty("N_STEPS","300"))
                    );

                    canneal.init();
            }

    }


    @State(Scope.Thread)
    public static class ReductionExternalState {

        Canneal canneal; // reduction

        @Setup(Level.Trial)
        public void setup() {
            canneal = new Canneal(
                    Integer.parseInt(System.getProperty("N_THREADS","1")),
                    Integer.parseInt(System.getProperty("N_SWAPS","10000")),
                    Integer.parseInt(System.getProperty("TEMP","2000")),
                    System.getProperty("NETLIST","src/main/resources/canneal/input/2500000.nets"),
                    Integer.parseInt(System.getProperty("N_STEPS","300"))
            );

            canneal.init();
        }

    }

    @State(Scope.Thread)
    public static class ReductionInternalState {

        jvbench.cannealReduction.Canneal canneal; // reduction

        @Setup(Level.Trial)
        public void setup() {
            canneal = new jvbench.cannealReduction.Canneal(
                    Integer.parseInt(System.getProperty("N_THREADS","1")),
                    Integer.parseInt(System.getProperty("N_SWAPS","10000")),
                    Integer.parseInt(System.getProperty("TEMP","2000")),
                    System.getProperty("NETLIST","src/main/resources/canneal/input/2500000.nets"),
                    Integer.parseInt(System.getProperty("N_STEPS","300"))
            );

            canneal.init();
        }

    }

    @State(Scope.Thread)
    public static class NoMinSizeCheckState {

        jvbench.cannealNoMinSizeCheck.Canneal canneal;

        @Setup(Level.Trial)
        public void setup() {
            canneal = new jvbench.cannealNoMinSizeCheck.Canneal(
                    Integer.parseInt(System.getProperty("N_THREADS","1")),
                    Integer.parseInt(System.getProperty("N_SWAPS","10000")),
                    Integer.parseInt(System.getProperty("TEMP","2000")),
                    System.getProperty("NETLIST","src/main/resources/canneal/input/2500000.nets"),
                    Integer.parseInt(System.getProperty("N_STEPS","300"))
            );

            canneal.init();
        }
    }

    @State(Scope.Thread)
    public static class IndexInRangeState {

        jvbench.cannealIndexInRange.Canneal canneal; // reduction

        @Setup(Level.Trial)
        public void setup() {
            canneal = new jvbench.cannealIndexInRange.Canneal(
                    Integer.parseInt(System.getProperty("N_THREADS","1")),
                    Integer.parseInt(System.getProperty("N_SWAPS","10000")),
                    Integer.parseInt(System.getProperty("TEMP","2000")),
                    System.getProperty("NETLIST","src/main/resources/canneal/input/2500000.nets"),
                    Integer.parseInt(System.getProperty("N_STEPS","300"))
            );

            canneal.init();
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public static void minSizeCheck(MinSizeCheckState state, Blackhole blackhole) {
        state.canneal.vector();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public static void reductionExternal(ReductionExternalState state, Blackhole blackhole) {
        state.canneal.vector();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public static void reductionInternal(ReductionInternalState state, Blackhole blackhole) {
        state.canneal.vector();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public static void noMinSizeCheck(NoMinSizeCheckState state, Blackhole blackhole) {
        state.canneal.vector();
    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public static void indexInRange(IndexInRangeState state, Blackhole blackhole) {
        state.canneal.vector();
    }


}
