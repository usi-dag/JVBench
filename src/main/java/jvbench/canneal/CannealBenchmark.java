package jvbench.canneal;

import jvbench.cannealReduction.Canneal;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.infra.Blackhole;

import java.util.concurrent.TimeUnit;

public class CannealBenchmark {

    @State(Scope.Thread)
    public static class MyState {

        Canneal canneal;

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



    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public static void serial(MyState state, Blackhole blackhole) {
        state.canneal.scalar();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public static void autoVec(MyState state, Blackhole blackhole) {
        state.canneal.scalar();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public static void explicitVec(MyState state, Blackhole blackhole) {
        state.canneal.vector();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public static void fullVec(MyState state, Blackhole blackhole) {
        state.canneal.vector();
    }


}
