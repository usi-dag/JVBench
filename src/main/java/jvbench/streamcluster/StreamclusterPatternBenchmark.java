package jvbench.streamcluster;

import jvbench.streamclusterReduction.PStream;
import jvbench.streamclusterReduction.StreamCluster;
import org.openjdk.jmh.annotations.*;

import java.util.concurrent.TimeUnit;

public class StreamclusterPatternBenchmark {

    @State(Scope.Thread)
    public static class ReductionState {

        long kMin = Long.parseLong(System.getProperty("k1", "3"));
        long kMax = Long.parseLong(System.getProperty("k2", "10"));
        int dim = Integer.parseInt(System.getProperty("dim", "128"));

        int chunkSize = Integer.parseInt(System.getProperty("chunksize", "128"));
        int clusterSize = Integer.parseInt(System.getProperty("clustersize", "10"));;
        String inputFileName = System.getProperty("infile", "src/main/resources/streamcluster/input/streamcluster_128_128.input");
        String outputFileName = System.getProperty("outfile", "output.txt");
        PStream stream;


        @Setup(Level.Trial)
        public void setup() {
            stream = StreamCluster.init(inputFileName);
        }
    }

    @State(Scope.Thread)
    public static class MyStatePow {

        long kMin = Long.parseLong(System.getProperty("k1", "3"));
        long kMax = Long.parseLong(System.getProperty("k2", "10"));
        int dim = Integer.parseInt(System.getProperty("dim", "128"));

        int chunkSize = Integer.parseInt(System.getProperty("chunksize", "128"));
        int clusterSize = Integer.parseInt(System.getProperty("clustersize", "10"));;
        String inputFileName = System.getProperty("infile", "src/main/resources/streamcluster/input/streamcluster_128_128.input");
        String outputFileName = System.getProperty("outfile", "output.txt");
        jvbench.streamclusterPow.PStream stream;


        @Setup(Level.Trial)
        public void setup() {
            stream = jvbench.streamclusterPow.StreamCluster.init(inputFileName);
        }
    }

    @State(Scope.Thread)
    public static class IndexInRangeState {

        long kMin = Long.parseLong(System.getProperty("k1", "3"));
        long kMax = Long.parseLong(System.getProperty("k2", "10"));
        int dim = Integer.parseInt(System.getProperty("dim", "128"));

        int chunkSize = Integer.parseInt(System.getProperty("chunksize", "128"));
        int clusterSize = Integer.parseInt(System.getProperty("clustersize", "10"));;
        String inputFileName = System.getProperty("infile", "src/main/resources/streamcluster/input/streamcluster_128_128.input");
        String outputFileName = System.getProperty("outfile", "output.txt");
        jvbench.streamclusterIndexInRange.PStream stream;


        @Setup(Level.Trial)
        public void setup() {
            stream = jvbench.streamclusterIndexInRange.StreamCluster.init(inputFileName);
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void reduction(ReductionState state) {
        StreamCluster.streamCluster(state.stream, state.kMin, state.kMax, state.dim, state.chunkSize, state.clusterSize, state.outputFileName, true);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void pow(MyStatePow state) {
        jvbench.streamclusterPow.StreamCluster.streamCluster(state.stream, state.kMin, state.kMax, state.dim, state.chunkSize, state.clusterSize, state.outputFileName, true);
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    @BenchmarkMode(Mode.SingleShotTime)
    @Fork(value = 5)
    @Warmup(iterations = 10)
    @Measurement(iterations = 10)
    public void indexInRange(IndexInRangeState state) {
        jvbench.streamclusterIndexInRange.StreamCluster.streamCluster(state.stream, state.kMin, state.kMax, state.dim, state.chunkSize, state.clusterSize, state.outputFileName, true);
    }
}
