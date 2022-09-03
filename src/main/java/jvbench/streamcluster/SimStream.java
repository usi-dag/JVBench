package jvbench.streamcluster;

import static jvbench.streamcluster.StreamCluster.lrand48;

public class SimStream implements PStream {

    private long n;

    public SimStream(long n) {
        this.n = n;
    }

    @Override
    public int read(float[] dest, int dim, int num) {
        int count = 0;
        for (int i = 0; i < num && n > 0; i++) {
            for (int k = 0; k < dim; k++) {
//                dest[i*dim + k] = lrand48()/(float) Integer.MAX_VALUE;
                dest[i*dim + k] = lrand48()/(float) Integer.MAX_VALUE;
            }
            n--;
            count++;
        }

        return count;
    }

    @Override
    public boolean feof() {
        return true;
    }
}
