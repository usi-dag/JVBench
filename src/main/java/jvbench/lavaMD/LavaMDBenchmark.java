package jvbench.lavaMD;

import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

public class LavaMDBenchmark {

    @State(Scope.Thread)
    public static class MyState {

        static String inputFilename = System.getProperty("input", "/lavaMD/input/lavaMD_127776.input");


        @Setup(Level.Trial)
        public  void setup() {
            LavaMD.init(inputFilename);
        }
    }



    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void serial(MyState state) {
        LavaMD.serial();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void autoVec(MyState state) {
        LavaMD.serial();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void explicitVec(MyState state) {
        LavaMD.vector();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void fullVec(MyState state) {
        LavaMD.vector();
    }
}
