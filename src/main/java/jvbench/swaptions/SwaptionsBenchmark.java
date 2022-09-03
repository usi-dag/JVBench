package jvbench.swaptions;

import jvbench.swaptionsConversion.HJMSecuritiesImpl;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

public class SwaptionsBenchmark {



    @State(Scope.Thread)
    public static class MyState {

        int ms = Integer.parseInt(System.getProperty("ms", "16384"));

        int ns = Integer.parseInt(System.getProperty("ns", "64"));

        @Setup(Level.Trial)
        public void setup() {

            HJMSecuritiesImpl.init(ms, ns);
        }

    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void serial(MyState state) {
        HJMSecuritiesImpl.benchmark(false);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void autoVec(MyState state) {
        HJMSecuritiesImpl.benchmark(false);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void explicitVec(MyState state) {
        HJMSecuritiesImpl.benchmark(true);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void fullVec(MyState state) {
        HJMSecuritiesImpl.benchmark(true);
    }



}
