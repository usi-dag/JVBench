package jvbench.axpy;

import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.infra.Blackhole;

import java.util.concurrent.TimeUnit;

public class AxpyPatternBenchmark {

    @State(Scope.Thread)
    public static class MyState {

        public int size;
        public double [] dy;
        public double [] dx;

        public double [] dyCopy;



        @Setup(Level.Trial)
        public void setup() {

            size = Integer.parseInt(System.getProperty("size", "70000")) * 1024;
            dx = new double[size];
            dy = new double[size];
            dyCopy = new double[size];

            for (int i = 0; i < size; i++) {
                dx[i] = 1.0;
                dy[i] = 2.0;
            }
            System.arraycopy(dy, 0, dyCopy, 0, size);
        }

        @TearDown(Level.Iteration)
        public void reset() {
            // after each iteration should reset dy
            System.arraycopy(dyCopy, 0, dy, 0, size);

        }



    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void fmaScalar(MyState state, Blackhole blackhole) {
        Axpy.scalarFMA(1.0, state.dx, state.dy, state.size);
        blackhole.consume(state.dx);
        blackhole.consume(state.dy);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void indexInRange(MyState state, Blackhole blackhole) {
        Axpy.vectorIndexInRange(1.0, state.dx, state.dy, state.size);
        blackhole.consume(state.dx);
        blackhole.consume(state.dy);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void broadcastExternal(MyState state, Blackhole blackhole) {
        Axpy.vectorBroadcastExternal(1.0, state.dx, state.dy, state.size);
        blackhole.consume(state.dx);
        blackhole.consume(state.dy);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void fma(MyState state, Blackhole blackhole) {
        Axpy.vectorFMA(1.0, state.dx, state.dy, state.size);
        blackhole.consume(state.dx);
        blackhole.consume(state.dy);
    }
}
