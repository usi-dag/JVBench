package jvbench.pathfinder;

import jvbench.pathfinderConditionInsideLoop.PathFinder;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

public class PathfinderBenchmark {



    @State(Scope.Thread)
    public static class MyState {
        @Setup(Level.Trial)
        public void setup() {
            jvbench.pathfinderConditionInsideLoop.PathFinder.init(System.getProperty("input", "src/main/resources/pathfinder/input/pathfinder_5000_5000.input"));
        }
    }




    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void serial(MyState state) {
        jvbench.pathfinderConditionInsideLoop.PathFinder.scalar();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void autoVec(MyState state) {
        jvbench.pathfinderConditionInsideLoop.PathFinder.scalar();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5, jvmArgsAppend = {"-XX:-UseSuperWord"})
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void explicitVec(MyState state) {
        jvbench.pathfinderConditionInsideLoop.PathFinder.vector();
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void fullVec(MyState state) {
        PathFinder.vector();
    }
}
