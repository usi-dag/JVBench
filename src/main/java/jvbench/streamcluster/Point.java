package jvbench.streamcluster;

public class Point {


    private float weight;
//    private float [] coord;
    private int assign;
    private double cost;
    private int coordIndex;

    public Point() {

    }

    public Point(float weight, int assign, int cost) {
        this.weight = weight;
//        this.coord = coord;
        this.assign = assign;
        this.cost = cost;
    }

    public void setCoordIndex(int index) {
        this.coordIndex = index;
    }

    public int getCoordIndex() {
        return this.coordIndex;
    }

    public float getWeight() {
        return weight;
    }

    public float getCoord(int index, boolean isCenter) {
        if (isCenter) {
            return StreamCluster.centerBlock[coordIndex + index];
        }
        return StreamCluster.block[coordIndex + index];
    }

    public int getAssign() {
        return assign;
    }

    public double getCost() {
        return cost;
    }

    public void setWeight(float weight) {
        this.weight = weight;
    }

//    public void setCoord(float [] coord) {
//        this.coord = coord;
//    }

    public void setCoord(float value, int index, boolean isCenter) {
        if (isCenter) {
            StreamCluster.centerBlock[coordIndex + index] = value;
        } else {
            StreamCluster.block[coordIndex + index] = value;
        }
    }

    public void setAssign(int assign) {
        this.assign = assign;
    }

    public void setCost(double cost) {
        this.cost = cost;
    }
}
