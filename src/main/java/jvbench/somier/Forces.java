package jvbench.somier;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorSpecies;

public class Forces {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_MAX;
    private static final int SPECIES_LENGTH = SPECIES.length();

    public static void computeForces(int n, double[][][][] x, double[][][][] f) {
        for (int i = 1; i < n - 1; i++) {
            for (int j = 1; j < n - 1; j++) {
                for (int k = 1; k < n - 1; k++) {
                    forceContribution(n, x, f, i, j, k, i - 1, j, k);
                    forceContribution(n, x, f, i, j, k, i + 1, j, k);
                    forceContribution(n, x, f, i, j, k, i, j - 1, k);
                    forceContribution(n, x, f, i, j, k, i, j + 1, k);
                    forceContribution(n, x, f, i, j, k, i, j, k - 1);
                    forceContribution(n, x, f, i, j, k, i, j, k + 1);
                }
            }
        }
    }

    private static void forceContribution(int n, double[][][][] x, double[][][][] f, int i, int j, int k, int neig_i, int neig_j, int neig_k) {
        double dx, dy, dz, dl, spring_F, FX, FY, FZ;

        assert (i >= 1);
        assert (j >= 1);
        assert (k >= 1);
        assert (i < n - 1);
        assert (j < n - 1);
        assert (k < n - 1);
        assert (neig_i >= 0);
        assert (neig_j >= 0);
        assert (neig_k >= 0);
        assert (neig_i < n);
        assert (neig_j < n);
        assert (neig_k < n);

        dx = x[0][neig_i][neig_j][neig_k] - x[0][i][j][k];
        dy = x[1][neig_i][neig_j][neig_k] - x[1][i][j][k];
        dz = x[2][neig_i][neig_j][neig_k] - x[2][i][j][k];
        dl = Math.sqrt(dx * dx + dy * dy + dz * dz);
        spring_F = 0.25 * Somier.getSpringK() * (dl - 1);
        FX = spring_F * dx / dl;
        FY = spring_F * dy / dl;
        FZ = spring_F * dz / dl;
        f[0][i][j][k] += FX;
        f[1][i][j][k] += FY;
        f[2][i][j][k] += FZ;
    }

    public static void computeForcesPrevec(int n, double[][][][] x, double[][][][] f) {
        int limit = SPECIES.loopBound(n - 2);
        for (int i = 1; i < n - 1; i++) {
            for (int j = 1; j < n - 1; j++) {
                forceContrVec(n, x, f, i, j, i, j + 1);
                forceContrVec(n, x, f, i, j, i - 1, j);
                forceContrVec(n, x, f, i, j, i + 1, j);
                forceContrVec(n, x, f, i, j, i, j - 1);
//            forceContrPrevec (n, x, f, i, j, i,   j-1);   //fails if force_contr_vec
                kForceContrPrevec(n, x, f, i, j);

                // scalar part
                for (int k = limit; k < n - 1; k++) {
                    forceContribution(n, x, f, i, j, k, i - 1, j, k);
                    forceContribution(n, x, f, i, j, k, i + 1, j, k);
                    forceContribution(n, x, f, i, j, k, i, j - 1, k);
                    forceContribution(n, x, f, i, j, k, i, j + 1, k);
                    forceContribution(n, x, f, i, j, k, i, j, k - 1);
                    forceContribution(n, x, f, i, j, k, i, j, k + 1);
                }
            }
        }
    }

    private static void kForceContrPrevec(int n, double[][][][] x, double[][][][] f, int i, int j) {
        double dx, dy, dz, dl, spring_F, FX, FY,FZ;

        for (int k=1; k<n-1; k++) {
            dx=x[0][i][j][k-1]-x[0][i][j][k];
            dy=x[1][i][j][k-1]-x[1][i][j][k];
            dz=x[2][i][j][k-1]-x[2][i][j][k];
            dl = Math.sqrt(dx*dx + dy*dy + dz*dz);
            spring_F = 0.25 * Somier.getSpringK()*(dl-1);
            FX = spring_F * dx/dl;
            FY = spring_F * dy/dl;
            FZ = spring_F * dz/dl;
            f[0][i][j][k] += FX;
            f[1][i][j][k] += FY;
            f[2][i][j][k] += FZ;
            dx=x[0][i][j][k+1]-x[0][i][j][k];
            dy=x[1][i][j][k+1]-x[1][i][j][k];
            dz=x[2][i][j][k+1]-x[2][i][j][k];
            dl = Math.sqrt(dx*dx + dy*dy + dz*dz);
            spring_F = 0.25 * Somier.getSpringK()*(dl-1);
            FX = spring_F * dx/dl;
            FY = spring_F * dy/dl;
            FZ = spring_F * dz/dl;
            f[0][i][j][k] += FX;
            f[1][i][j][k] += FY;
            f[2][i][j][k] += FZ;
        }
    }

    private static void forceContrVec(int n, double[][][][] x, double[][][][] f, int i, int j, int neig_i, int neig_j) {
        int limit = SPECIES.loopBound(n - 2);

        DoubleVector v1 = DoubleVector.broadcast(SPECIES, 1.0);
        DoubleVector vSprK = DoubleVector.broadcast(SPECIES, 0.25 * Somier.getSpringK());
        int k;
        for (k = 1; k < limit; k += SPECIES_LENGTH) {

            DoubleVector vX1 = DoubleVector.fromArray(SPECIES, x[0][neig_i][neig_j], k);
            DoubleVector vX2 = DoubleVector.fromArray(SPECIES, x[0][i][j], k);
            DoubleVector vDx = vX1.sub(vX2);
            DoubleVector vDx2 = vDx.mul(vDx);

            DoubleVector vX3 = DoubleVector.fromArray(SPECIES, x[1][neig_i][neig_j], k);
            DoubleVector vX4 = DoubleVector.fromArray(SPECIES, x[1][i][j], k);
            DoubleVector vDy = vX3.sub(vX4);
            DoubleVector vDy2 = vDy.mul(vDy);

            DoubleVector vX5 = DoubleVector.fromArray(SPECIES, x[2][neig_i][neig_j], k);
            DoubleVector vX6 = DoubleVector.fromArray(SPECIES, x[2][i][j], k);
            DoubleVector vDz = vX5.sub(vX6);
            DoubleVector vDz2 = vDz.mul(vDz);

            DoubleVector vDl = vDx2.add(vDy2).add(vDz2).sqrt();

            DoubleVector vDl1 = vDl.sub(v1);
            DoubleVector vSprF = vSprK.mul(vDl1);
            DoubleVector vDFX = vDx.div(vDl);
            DoubleVector vDFY = vDy.div(vDl);
            DoubleVector vDFZ = vDz.div(vDl);

            DoubleVector vFX = DoubleVector.fromArray(SPECIES, f[0][i][j], k);
            DoubleVector vFY = DoubleVector.fromArray(SPECIES, f[1][i][j], k);
            DoubleVector vFZ = DoubleVector.fromArray(SPECIES, f[2][i][j], k);

            vFX = vSprF.mul(vDFX).add(vFX);
            vFY = vSprF.mul(vDFY).add(vFY);
            vFZ = vSprF.mul(vDFZ).add(vFZ);

            vFX.intoArray(f[0][i][j], k);
            vFY.intoArray(f[1][i][j], k);
            vFZ.intoArray(f[2][i][j], k);


        }


    }
}
