package jvbench.lavaMD;

import jvbench.lavaMDFma.LavaMD;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

public class LavaMDPatternBenchmark {


    @State(Scope.Thread)
    public static class FmaState {

        static String inputFilename = System.getProperty("input", "/lavaMD/input/lavaMD_127776.input");


        @Setup(Level.Trial)
        public  void setup() {
            LavaMD.init(inputFilename);
        }
    }

    @State(Scope.Thread)
    public static class ReductionState {

        static String inputFilename = System.getProperty("input", "/lavaMD/input/lavaMD_127776.input");


        @Setup(Level.Trial)
        public  void setup() {
            jvbench.lavaMDReduction.LavaMD.init(inputFilename);
        }
    }

    @State(Scope.Thread)
    public static class IndexInRangeState {

        static String inputFilename = System.getProperty("input", "/lavaMD/input/lavaMD_127776.input");


        @Setup(Level.Trial)
        public  void setup() {
            LavaMD.init(inputFilename);
        }
    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void fma(FmaState state) {
        LavaMD.vector();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void fmaScalar(FmaState state) {
        LavaMD.serial();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void reduction(ReductionState state) {
        jvbench.lavaMDReduction.LavaMD.vector();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void indexInRange(IndexInRangeState state) {
        LavaMD.vector();
    }

}
