package jvbench.streamclusterPow;


public class Points {

    private int num;
    private int dim;
    private Point[] points;

    public Points(int dim, int chunkSize, int size) {
        this.num = chunkSize;
        this.dim = dim;
        points = new Point[size];
    }


    public int getNum() {
        return num;
    }

    public void setNum(int num) {
        this.num = num;
    }

    public int getDim() {
        return dim;
    }

    public void setDim(int dim) {
        this.dim = dim;
    }

    public Point[] getPoints() {
        return points;
    }

    public void setPoints(Point[] points) {
        this.points = points;
    }
}
