package jvbench.somierIndexInRange;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorSpecies;

public class Somier {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_MAX;
    private static final int SPECIES_LENGTH = SPECIES.length();

    private static final double[] xCenter = new double[3];
    private static final double[][] ref = new double[1000][3];
    private static final double dt = 0.001;
    private static final double M = 1.0;
    private static final double spring_K = 10.0;
    private static int err;

    private static int size = 10;
    private static int steps = 4;

    static double[][][][] x;
    static double[][][][] v;
    private static double[][][][] a;
    private static double[][][][] f;

    public static void init(int ntSteps, int n) {
        steps = ntSteps;
        size = n;


        x = new double[3][n][n][n];
        v = new double[3][n][n][n];
        a = new double[3][n][n][n];
        f = new double[3][n][n][n];
        double[][][][] fRef = new double[3][n][n][n];

        clear4D(n, f);
        clear4D(n, a);
        clear4D(n, v);
        clear4D(n, fRef);

        initX(n, x);

        v[0][n / 2][n / 2][n / 2] = 0.1;
        v[1][n / 2][n / 2][n / 2] = 0.1;
        v[2][n / 2][n / 2][n / 2] = 0.1;
    }


    public static void scalar() {
        for (int nt = 0; nt < steps - 1; nt++) {

            // reset aggregate states
            xCenter[0] = 0;
            xCenter[1] = 0;
            xCenter[2] = 0;
            clear4D(size, f);

            Forces.computeForces(size, x, f);
            accelaration(size, a, f, M);
            velocities(size, v, a, dt);
            positions(size, x, v, dt);

            computeStats(size, x, xCenter);
        }
    }

    public static void vector() {
        for (int nt = 0; nt < steps - 1; nt++) {

            // reset aggregate states
            xCenter[0] = 0;
            xCenter[1] = 0;
            xCenter[2] = 0;
            clear4D(size, f);

            Forces.computeForcesPrevec(size, x, f);    // printf ("Computed forces\n"); print_4D(size, "F", F); printf ("\n");
            accelerationVector(size, a, f, M);         // printf ("Computed Accelerations\n"); print_4D(size, "A", A); printf ("\n");
            velocitiesVector(size, v, a, dt);          // printf ("Computed Velocities\n"); print_4D(size, "V", V); printf ("\n");
            positionsVector(size, x, v, dt);           // printf ("Computed Positions\n"); print_4D(size, "X", X); printf ("\n");


            computeStats(size, x, xCenter);
        }
    }

    public static void main(String[] args) {
       int steps = 4;
       int N = 128;

        init(steps, N);
        scalar();
//        vector();


        System.out.println("\tV=" + v[0][N / 2][N / 2][N / 2] + ", " + v[1][N / 2][N / 2][N / 2] + ", " + v[2][N / 2][N / 2][N / 2] +
                "\t\t X= " + x[0][N / 2][N / 2][N / 2] + ", " + x[1][N / 2][N / 2][N / 2] + ", " + x[2][N / 2][N / 2][N / 2]);
    }

    private static void positionsVector(int n, double[][][][] x, double[][][][] v, double dt) {
        int i, j, k;

        DoubleVector vV0, vV1, vV2, vX0, vX1, vX2;

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k += SPECIES_LENGTH) {
                    VectorMask<Double> mask = SPECIES.indexInRange(k, n);
                    vX0 = DoubleVector.fromArray(SPECIES, x[0][i][j], k, mask);
                    vV0 = DoubleVector.fromArray(SPECIES, v[0][i][j], k, mask);
                    vX0 = vV0.mul(dt).add(vX0);
                    vX0.intoArray(x[0][i][j], k, mask);

                    vX1 = DoubleVector.fromArray(SPECIES, x[1][i][j], k, mask);
                    vV1 = DoubleVector.fromArray(SPECIES, v[1][i][j], k, mask);
                    vX1 = vV1.mul(dt).add(vX1);
                    vX1.intoArray(x[1][i][j], k, mask);

                    vX2 = DoubleVector.fromArray(SPECIES, x[2][i][j], k, mask);
                    vV2 = DoubleVector.fromArray(SPECIES, v[2][i][j], k, mask);
                    vX2 = vV2.mul(dt).add(vX2);
                    vX2.intoArray(x[2][i][j], k, mask);
                }
            }
    }

    private static void velocitiesVector(int n, double[][][][] v, double[][][][] a, double dt) {
        int i, j, k;
        //#dear compiler: please fuse next two loops if you can

        DoubleVector vV0, vV1, vV2, vA0, vA1, vA2;


        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k += SPECIES_LENGTH) {
                    VectorMask<Double> mask = SPECIES.indexInRange(k, n);
                    vV0 = DoubleVector.fromArray(SPECIES, v[0][i][j], k, mask);
                    vA0 = DoubleVector.fromArray(SPECIES, a[0][i][j], k, mask);
                    vV0 = vA0.mul(dt).add(vV0);
                    vV0.intoArray(v[0][i][j], k, mask);

                    vV1 = DoubleVector.fromArray(SPECIES, v[1][i][j], k, mask);
                    vA1 = DoubleVector.fromArray(SPECIES, a[1][i][j], k, mask);
                    vV1 = vA1.mul(dt).add(vV1);
                    vV1.intoArray(v[1][i][j], k, mask);

                    vV2 = DoubleVector.fromArray(SPECIES, v[2][i][j], k, mask);
                    vA2 = DoubleVector.fromArray(SPECIES, a[2][i][j], k, mask);
                    vV2 = vA2.mul(dt).add(vV2);
                    vV2.intoArray(v[2][i][j], k, mask);
                }
            }

    }

    private static void accelerationVector(int n, double[][][][] a, double[][][][] f, double m) {

        int i, j, k;


        DoubleVector vF0, vF1, vF2, vA0, vA1, vA2;
        double invM = 1 / M;

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k += SPECIES_LENGTH) {
                    VectorMask<Double> mask = SPECIES.indexInRange(k, n);
                    vF0 = DoubleVector.fromArray(SPECIES, f[0][i][j], k, mask);
                    vA0 = vF0.mul(invM);
                    vA0.intoArray(a[0][i][j], k, mask);

                    vF1 = DoubleVector.fromArray(SPECIES, f[1][i][j], k, mask);
                    vA1 = vF1.mul(invM);
                    vA1.intoArray(a[1][i][j], k, mask);

                    vF2 = DoubleVector.fromArray(SPECIES, f[2][i][j], k, mask);
                    vA2 = vF2.mul(invM);
                    vA2.intoArray(a[2][i][j], k, mask);
                }
            }
        }
    }

    private static void computeStats(int n, double[][][][] x, double[] xCenter) {
        // TODO possible vectorization
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    xCenter[0] += x[0][i][j][k];
                    xCenter[1] += x[1][i][j][k];
                    xCenter[2] += x[2][i][j][k];
                }
            }
        }
        xCenter[0] /= (n * n * n);
        xCenter[1] /= (n * n * n);
        xCenter[2] /= (n * n * n);
    }

    private static void positions(int n, double[][][][] x, double[][][][] v, double dt) {
        int i, j, k;
        //#dear compiler: please fuse next two loops if you can
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                for (k = 0; k < n; k++) {
                    x[0][i][j][k] += v[0][i][j][k] * dt;
                    x[1][i][j][k] += v[1][i][j][k] * dt;
                    x[2][i][j][k] += v[2][i][j][k] * dt;
                }
    }

    private static void velocities(int n, double[][][][] v, double[][][][] a, double dt) {
        int i, j, k;
        //#dear compiler: please fuse next two loops if you can
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    v[0][i][j][k] += a[0][i][j][k] * dt;
                    v[1][i][j][k] += a[1][i][j][k] * dt;
                    v[2][i][j][k] += a[2][i][j][k] * dt;
                }
            }
    }

    private static void accelaration(int n, double[][][][] a, double[][][][] f, double m) {
        int i, j, k;
        //#dear compiler: please fuse next two loops if you can
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                for (k = 0; k < n; k++) {
                    a[0][i][j][k] = f[0][i][j][k] / M;
                    a[1][i][j][k] = f[1][i][j][k] / M;
                    a[2][i][j][k] = f[2][i][j][k] / M;
                }
    }


    private static void printState(int n, double[][][][] x, double[] xCenter, int nt) {
        System.out.print("t=" + nt + "\t");
        System.out.println(
                "XC=" + xCenter[0] + "," + xCenter[1] + "," + xCenter[2] +
                        " X[n/2-1] = " + x[0][n / 2 - 1][n / 2][n / 2] + "," + x[1][n / 2 - 1][n / 2][n / 2] + "," + x[2][n / 2 - 1][n / 2][n / 2] +
                        " X[n/2] = " + x[0][n / 2][n / 2][n / 2] + "," + x[1][n / 2][n / 2][n / 2] + "," + x[2][n / 2][n / 2][n / 2] +
                        " V[n/2+1] = " + x[0][n / 2 + 1][n / 2][n / 2] + "," + x[1][n / 2 + 1][n / 2][n / 2] + "," + x[2][n / 2 + 1][n / 2][n / 2]
        );
//            F[0][n/2][n/2][n/2], F[1][n/2][n/2][n/2], F[2][n/2][n/2][n/2],
//            V[0][n/2][n/2][n/2], V[1][n/2][n/2][n/2], V[2][n/2][n/2][n/2]);
//        printf ("\n ");
    }

    private static void initX(int n, double[][][][] x) {
        int i, j, k;

        xCenter[0] = 0;
        xCenter[1] = 0;
        xCenter[2] = 0;

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                for (k = 0; k < n; k++) {
                    x[0][i][j][k] = i;
                    x[1][i][j][k] = j;
                    x[2][i][j][k] = k;

                    xCenter[0] += x[0][i][j][k];
                    xCenter[1] += x[1][i][j][k];
                    xCenter[2] += x[2][i][j][k];
                }

        xCenter[0] /= (n * n * n);
        xCenter[1] /= (n * n * n);
        xCenter[2] /= (n * n * n);


//   X[n/2][n/2][n/2][0] += 0.5; X[n/2][n/2][n/2][1] += 0.5; X[n/2][n/2][n/2][2] += 0.5;
//   X[n/2][n/2][n/2][0] += 0.5; X[n/2][n/2][n/2][1] += 0.5; 
//   X[n/2][n/2][n/2][0] += 0.5;  
    }

    private static void clear4D(int n, double[][][][] x) {
        int i, j, k;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                for (k = 0; k < n; k++) {
                    x[0][i][j][k] = 0.0;
                    x[1][i][j][k] = 0.0;
                    x[2][i][j][k] = 0.0;
                }
    }


    public static double getSpringK() {
        return spring_K;
    }


}
