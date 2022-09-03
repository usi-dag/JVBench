package jvbench.canneal;

public class Rng {

    private final MTRand rng;
    private static int seed = 0;


    public Rng() {
        rng = new MTRand(seed++);
    }

    public long rand() {
        return rng.randInt();
    }

    public long rand(int max) {
        return rng.randInt(max - 1);
    }

    public double drand() {
        return rng.rand();
    }
}
