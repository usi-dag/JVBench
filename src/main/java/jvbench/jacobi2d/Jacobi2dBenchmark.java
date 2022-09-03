package jvbench.jacobi2d;

import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.infra.Blackhole;

import java.util.concurrent.TimeUnit;

public class Jacobi2dBenchmark {

    @State(Scope.Thread)
    public static class MyState {

        int size = Integer.parseInt(System.getProperty("size","10000"));
        int tSteps = Integer.parseInt(System.getProperty("tsteps","14"));

        public final double [][] a = new double[size][size];
        public final double [][] b = new double[size][size];


        private void initArray(int n, double [][] a, double [][] b) {
            int i, j;
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    a[i][j] = ((double) i * (j + 2) + 2) / n;
                    b[i][j] = ((double) i * (j + 3) + 3) / n;
                }
            }
        }

        @Setup(Level.Trial)
        public void setup() {
            initArray(size, a, b);
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void serial(MyState state, Blackhole blackhole) {
        Jacobi2d.serial(state.tSteps, state.size, state.a, state.b);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void autoVec(MyState state, Blackhole blackhole) {
        Jacobi2d.serial(state.tSteps, state.size, state.a, state.b);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void explicitVec(MyState state, Blackhole blackhole) {
        Jacobi2d.vector(state.tSteps, state.size, state.a, state.b);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void fullVec(MyState state, Blackhole blackhole) {
        Jacobi2d.vector(state.tSteps, state.size, state.a, state.b);
    }

}
