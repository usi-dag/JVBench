package jvbench.streamclusterIndexInRange;

public class Median {

    Points points;
    long kmin;
    long kmax;
    long[] kfinal;
    int pid;

    public Median(Points points, long kmin, long kmax, long[] kfinal, int pid) {
        this.points = points;
        this.kmin = kmin;
        this.kmax = kmax;
        this.kfinal = kfinal;
        this.pid = pid;
    }

    public Points getPoints() {
        return points;
    }

    public long getKmin() {
        return kmin;
    }

    public long getKmax() {
        return kmax;
    }

    public long[] getKfinal() {
        return kfinal;
    }

    public int getPid() {
        return pid;
    }
}
