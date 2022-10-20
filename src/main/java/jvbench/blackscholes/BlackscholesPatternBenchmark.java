package jvbench.blackscholes;

import jvbench.blackscholesMergeWihtoutBlend.Blackscholes;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.infra.Blackhole;

import java.util.concurrent.TimeUnit;

public class BlackscholesPatternBenchmark {


    @State(Scope.Thread)
    public static class mergeWithoutBlendState {

        @Setup(Level.Trial)
        public static void setup() {
            jvbench.blackscholesMergeWihtoutBlend.Blackscholes.init(System.getProperty("input","/blackscholes/input/in_512K.input"));
        }

    }

    @State(Scope.Thread)
    public static class PowState {

        @Setup(Level.Trial)
        public static void setup() {
            jvbench.blackscholesPow.Blackscholes.init(System.getProperty("input","/blackscholes/input/in_512K.input"));
        }

    }

    @State(Scope.Thread)
    public static class IndexInRangeState {

        @Setup(Level.Trial)
        public static void setup() {
            jvbench.blackscholesIndexInRange.Blackscholes.init(System.getProperty("input","/blackscholes/input/in_512K.input"));
        }

    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    // blend -> zero.add(mask).add(mask.not) with mask
    public void mergeWithoutBlend(Blackhole blackhole, mergeWithoutBlendState myState) {
        jvbench.blackscholesMergeWihtoutBlend.Blackscholes.vector();
        blackhole.consume(Blackscholes.getPrices());
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void pow(Blackhole blackhole, PowState myState) {
        jvbench.blackscholesPow.Blackscholes.vector();
        blackhole.consume(jvbench.blackscholesPow.Blackscholes.getPrices());
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void indexInRange(Blackhole blackhole, IndexInRangeState myState) {
        jvbench.blackscholesIndexInRange.Blackscholes.vector();
        blackhole.consume(jvbench.blackscholesIndexInRange.Blackscholes.getPrices());
    }
}
