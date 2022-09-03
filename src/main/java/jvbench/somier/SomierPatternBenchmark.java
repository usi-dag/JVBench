package jvbench.somier;

import jvbench.somierIndexInRange.Somier;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

public class SomierPatternBenchmark {


    @State(Scope.Thread)
    public static class FmaState {

        int ntSteps = Integer.parseInt(System.getProperty("steps", "10"));
        int n = Integer.parseInt(System.getProperty("n", "128"));

        @Setup(Level.Trial)
        public void setup() {
            jvbench.somierFma.Somier.init(ntSteps, n);
        }
    }

    @State(Scope.Thread)
    public static class PowState {

        int ntSteps = Integer.parseInt(System.getProperty("steps", "10"));
        int n = Integer.parseInt(System.getProperty("n", "128"));

        @Setup(Level.Trial)
        public void setup() {
            jvbench.somierPow.Somier.init(ntSteps, n);
        }
    }

    @State(Scope.Thread)
    public static class IndexInRangeState {

        int ntSteps = Integer.parseInt(System.getProperty("steps", "10"));
        int n = Integer.parseInt(System.getProperty("n", "128"));

        @Setup(Level.Trial)
        public void setup() {
            Somier.init(ntSteps, n);
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void fma(FmaState state) {
        jvbench.somierFma.Somier.vector();
    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void fmaScalar(FmaState state) {
        jvbench.somierFma.Somier.scalar();
    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void pow(PowState state) {
        jvbench.somierPow.Somier.vector();
    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void indexInRange(IndexInRangeState state) {
        Somier.vector();
    }
}
