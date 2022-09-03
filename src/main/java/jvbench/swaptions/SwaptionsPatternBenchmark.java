package jvbench.swaptions;

import jvbench.swaptionsIndexInRange.HJMSecuritiesImpl;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

public class SwaptionsPatternBenchmark {


    @State(Scope.Thread)
    public static class PowState {

        int ms = Integer.parseInt(System.getProperty("ms", "16384"));

        int ns = Integer.parseInt(System.getProperty("ns", "64"));

        @Setup(Level.Trial)
        public void setup() {
            jvbench.swaptionsPow.HJMSecuritiesImpl.init(ms, ns);
        }

    }

    @State(Scope.Thread)
    public static class ConversionState {

        int ms = Integer.parseInt(System.getProperty("ms", "16384"));

        int ns = Integer.parseInt(System.getProperty("ns", "64"));

        @Setup(Level.Trial)
        public void setup() {
            jvbench.swaptionsConversion.HJMSecuritiesImpl.init(ms, ns);
        }

    }

    @State(Scope.Thread)
    public static class IndexInRangeState {

        int ms = Integer.parseInt(System.getProperty("ms", "16384"));

        int ns = Integer.parseInt(System.getProperty("ns", "64"));

        @Setup(Level.Trial)
        public void setup() {
            HJMSecuritiesImpl.init(ms, ns);
        }

    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void pow(PowState state) {
        jvbench.swaptionsPow.HJMSecuritiesImpl.benchmark(true);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void conversion(ConversionState state) {
        jvbench.swaptionsConversion.HJMSecuritiesImpl.benchmark(true);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void indexInRange(IndexInRangeState state) {
        HJMSecuritiesImpl.benchmark(true);
    }



}
