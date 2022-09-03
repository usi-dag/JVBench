package jvbench.streamcluster;

public interface PStream {

    int read(float [] dest, int dim, int num);


    boolean feof();
}
