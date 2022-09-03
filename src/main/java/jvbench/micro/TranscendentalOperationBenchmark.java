package jvbench.micro;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.infra.Blackhole;

import java.util.Random;
import java.util.concurrent.TimeUnit;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 VECTOR_OP_TAN   = 101,
 VECTOR_OP_TANH  = 102,
 VECTOR_OP_SIN   = 103,
 VECTOR_OP_SINH  = 104,
 VECTOR_OP_COS   = 105,
 VECTOR_OP_COSH  = 106,
 VECTOR_OP_ASIN  = 107,
 VECTOR_OP_ACOS  = 108,
 VECTOR_OP_ATAN  = 109,
 VECTOR_OP_ATAN2 = 110,
 VECTOR_OP_CBRT  = 111,
 VECTOR_OP_LOG   = 112,
 VECTOR_OP_LOG10 = 113,
 VECTOR_OP_LOG1P = 114,
 VECTOR_OP_POW   = 115,
 VECTOR_OP_EXP   = 116,
 VECTOR_OP_EXPM1 = 117,
 VECTOR_OP_HYPOT = 118,
 */
@BenchmarkMode(Mode.SingleShotTime)
@Fork(value = 5)
@Warmup(iterations = 10)
@Measurement(iterations = 10)
public class TranscendentalOperationBenchmark {

    private static final VectorSpecies<Double> DOUBLE_SPECIES = DoubleVector.SPECIES_PREFERRED;
    private static final int SPECIES_LENGTH = DOUBLE_SPECIES.length();

    private static final int ITERATIONS_LIMIT_VECTOR = 10_000_000;
    private static final int ITERATIONS_LIMIT = ITERATIONS_LIMIT_VECTOR * SPECIES_LENGTH;

    private static final int SEED = 42;
    private static final Random random = new Random(SEED);

    @State(Scope.Thread)
    public static class UnaryVectorState {


        @Param({"TAN"  ,
                "TANH" ,
                "SIN"  ,
                "SINH" ,
                "COS"  ,
                "COSH" ,
                "ASIN" ,
                "ACOS" ,
                "ATAN" ,
                "CBRT" ,
                "LOG"  ,
                "LOG10",
                "LOG1P",
                "EXP"  ,
                "EXPM1"
        }
        )
        public String operation;
        public VectorOperators.Unary type;
        public Function<Double, Double> reference;

        public DoubleVector vec = DoubleVector.broadcast(DOUBLE_SPECIES, 1.0);
        public double x = 1;
        public double y = 2;



        @Setup(Level.Trial)
        public void setup() {
            switch (operation) {
                case "TAN" -> {
                    type = VectorOperators.TAN;
                    reference = Math::tan;
                }
                case "TANH" -> {
                    type = VectorOperators.TANH;
                    reference = Math::tanh;
                }
                case "SIN" -> {
                    type = VectorOperators.SIN;
                    reference = Math::sin;
                }
                case "SINH" -> {
                    type = VectorOperators.SINH;
                    reference = Math::sinh;
                }
                case "COS" -> {
                    type = VectorOperators.COS;
                    reference = Math::cos;
                }
                case "COSH" -> {
                    type = VectorOperators.COSH;
                    reference = Math::cosh;
                }
                case "ASIN" -> {
                    type = VectorOperators.ASIN;
                    reference = Math::asin;
                }
                case "ACOS" -> {
                    type = VectorOperators.ACOS;
                    reference = Math::acos;
                }
                case "ATAN" -> {
                    type = VectorOperators.ATAN;
                    reference = Math::atan;
                }
                case "CBRT" -> {
                    type = VectorOperators.CBRT;
                    reference = Math::cbrt;
                }
                case "LOG" -> {
                    type = VectorOperators.LOG;
                    reference = Math::log;
                }
                case "LOG10" -> {
                    type = VectorOperators.LOG10;
                    reference = Math::log10;
                }
                case "LOG1P" -> {
                    type = VectorOperators.LOG1P;
                    reference = Math::log1p;
                }
                case "EXP" -> {
                    type = VectorOperators.EXP;
                    reference = Math::exp;
                }
                case "EXPM1" -> {
                    type = VectorOperators.EXPM1;
                    reference = Math::expm1;
                }
                default -> {
                    System.out.println("ERROR NO FUNCTION FOR CASE: " + operation);
                }
            }
        }

        @TearDown(Level.Trial)
        public void reset() {

        }
    }

    @State(Scope.Thread)
    public static class BinaryVectorState {
        @Param({
                "ATAN2",
                "POW"  ,
                "HYPOT"
        }
        )
        public String operation;
        public VectorOperators.Binary type;
        public BiFunction<Double, Double, Double> reference;

        public DoubleVector vec = DoubleVector.broadcast(DOUBLE_SPECIES, 1.0);
        public double x = 1;
        public double y = 2;


        @Setup(Level.Trial)
        public void setup() {
            switch (operation) {
                case "ATAN2" -> {
                    type = VectorOperators.ATAN2;
                    reference = Math::atan2;
                }
                case "POW" -> {
                    type = VectorOperators.POW;
                    reference = Math::pow;
                }
                case "HYPOT" -> {
                    type = VectorOperators.HYPOT;
                    reference = Math::hypot;
                }
                default -> {
                    System.out.println("ERROR NO FUNCTION FOR CASE: " + operation);
                }
            }
        }

        @TearDown(Level.Trial)
        public void reset() {

        }
    }

    @State(Scope.Thread)
    public static class ExecutionPlan {

        public double [] x = new double[ITERATIONS_LIMIT];
        public double [] y = new double[ITERATIONS_LIMIT];

        @Setup(Level.Trial)
        public void setup() {
            for (int i = 0; i < ITERATIONS_LIMIT; i++) {
                x[i] = random.nextDouble();
                y[i] = random.nextDouble();
            }
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarTAN(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.tan(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorTAN(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.TAN));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.tan(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarTANH(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.tanh(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorTANH(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.TANH));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.tanh(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarSIN(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.sin(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorSIN(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.SIN));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.sin(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarSINH(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.sinh(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorSINH(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.SINH));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.sinh(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarCOS(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.cos(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorCOS(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.COS));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.cos(state.x[i]));
        }
    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarCOSH(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.cosh(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorCOSH(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.COSH));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.cosh(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarASIN(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.asin(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorASIN(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.ASIN));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.asin(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarACOS(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.acos(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorACOS(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.ACOS));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.acos(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarATAN(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.atan(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorATAN(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.ATAN));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.atan(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarCBRT(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.cbrt(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorCBRT(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.CBRT));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.cbrt(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarLOG(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.log(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorLOG(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.LOG));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.log(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarLOG10(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.log10(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorLOG10(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.LOG10));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.log10(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarLOG1P(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.log1p(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorLOG1P(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.LOG1P));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.log1p(state.x[i]));
        }
    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarEXP(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.exp(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorEXP(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.EXP));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.exp(state.x[i]));
        }
    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarEXPM1(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.expm1(state.x[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorEXPM1(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.EXPM1));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.expm1(state.x[i]));
        }
    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarATAN2(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.atan2(state.x[i], state.y[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorATAN2(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            DoubleVector yV = DoubleVector.fromArray(DOUBLE_SPECIES, state.y, i);
            blackhole.consume(xV.lanewise(VectorOperators.ATAN2, yV));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.atan2(state.x[i], state.y[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarPOW(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.pow(state.x[i], state.y[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorPOW(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            DoubleVector yV = DoubleVector.fromArray(DOUBLE_SPECIES, state.y, i);
            blackhole.consume(xV.lanewise(VectorOperators.POW, yV));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.pow(state.x[i], state.y[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarPOW0(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.pow(state.x[i], 0));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorPOW0(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.POW, 0));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.pow(state.x[i], 0));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarPOW1(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.pow(state.x[i], 1));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorPOW1(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.POW, 1));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.pow(state.x[i], 1));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarPOW2(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.pow(state.x[i], 2));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorPOW2(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.POW, 2));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.pow(state.x[i], 2));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarMUL0(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(state.x[i] * 0);
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorMUL0(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.MUL, 0));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(state.x[i] * 0);
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarMUL1(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(state.x[i] * 1);
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorMUL1(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.MUL, 1));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(state.x[i] * 1);
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarMUL2(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(state.x[i] * state.x[i]);
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorMUL2(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            blackhole.consume(xV.lanewise(VectorOperators.MUL, xV));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(state.x[i] * state.x[i]);
        }
    }


    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void scalarHYPOT(ExecutionPlan state, Blackhole blackhole) {
        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.hypot(state.x[i], state.y[i]));
        }
    }

    @Benchmark
    @OutputTimeUnit(TimeUnit.SECONDS)
    public void vectorHYPOT(ExecutionPlan state, Blackhole blackhole) {
        final int limit = DOUBLE_SPECIES.loopBound(state.x.length);
        int i;
        for (i = 0; i < limit; i += SPECIES_LENGTH) {
            DoubleVector xV = DoubleVector.fromArray(DOUBLE_SPECIES, state.x, i);
            DoubleVector yV = DoubleVector.fromArray(DOUBLE_SPECIES, state.y, i);
            blackhole.consume(xV.lanewise(VectorOperators.HYPOT, yV));
        }

        for (; i < ITERATIONS_LIMIT; i++) {
            blackhole.consume(Math.hypot(state.x[i], state.y[i]));
        }
    }




//    @Benchmark
//    @OutputTimeUnit(TimeUnit.SECONDS)
//    @Fork(jvmArgsAppend = {"-XX:-UseSuperWord"})
//    public void unaryScalar(UnaryVectorState state, Blackhole sink) {
//        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
//            state.reference.apply(state.x);
//        }
//    }
//
//    @Benchmark
//    @OutputTimeUnit(TimeUnit.SECONDS)
//    public void unary(UnaryVectorState state, Blackhole sink) {
//        for (int i = 0; i < ITERATIONS_LIMIT_VECTOR; i++) {
//            state.vec.lanewise(state.type);
//        }
//    }
//
//    @Benchmark
//    @OutputTimeUnit(TimeUnit.SECONDS)
//    @Fork(jvmArgsAppend = {"-XX:-UseSuperWord"})
//    public void binaryScalar(BinaryVectorState state, Blackhole sink) {
//        for (int i = 0; i < ITERATIONS_LIMIT; i++) {
//            state.reference.apply(state.x, state.y);
//        }
//    }
//
//    @Benchmark
//    @OutputTimeUnit(TimeUnit.SECONDS)
//    public void binary(BinaryVectorState state, Blackhole sink) {
//        for (int i = 0; i < ITERATIONS_LIMIT_VECTOR; i++) {
//            state.vec.lanewise(state.type, state.y);
//        }
//    }


}
