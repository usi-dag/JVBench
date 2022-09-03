package jvbench.blackscholes;

import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.infra.Blackhole;

import java.util.concurrent.TimeUnit;

public class BlackscholesBenchmark {

    @State(Scope.Thread)
    public static class MyState {

        @Setup(Level.Trial)
        public static void setup() {
            Blackscholes.init(System.getProperty("input","./src/main/resources/blackscholes/input/in_512K.input"));
        }

    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void serial(Blackhole blackhole, MyState myState) {
        Blackscholes.scalar();
        blackhole.consume(Blackscholes.getPrices());
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void autoVec(Blackhole blackhole, MyState myState) {
        Blackscholes.scalar();
        blackhole.consume(Blackscholes.getPrices());
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void explicitVec(Blackhole blackhole, MyState myState) {
        Blackscholes.vector();
        blackhole.consume(Blackscholes.getPrices());
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void fullVec(Blackhole blackhole, MyState myState) {
        Blackscholes.vector();
        blackhole.consume(Blackscholes.getPrices());
    }
}
