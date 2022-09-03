package jvbench.blackscholesPow;

public class OptionData {

    private float s;
    private float strike;
    private float r;
    private float divq;
    private float v;
    private float t;

    private char OptionType; // "P"=PUT, "C"=CALL
    private float divs;
    private float DGrefval;

    public OptionData(float s, float strike, float r, float divq, float v, float t, char optionType, float divs, float DGrefval) {
        this.s = s;
        this.strike = strike;
        this.r = r;
        this.divq = divq;
        this.v = v;
        this.t = t;
        OptionType = optionType;
        this.divs = divs;
        this.DGrefval = DGrefval;
    }

    public float getS() {
        return s;
    }

    public void setS(float s) {
        this.s = s;
    }

    public float getStrike() {
        return strike;
    }

    public void setStrike(float strike) {
        this.strike = strike;
    }

    public float getR() {
        return r;
    }

    public void setR(float r) {
        this.r = r;
    }

    public float getDivq() {
        return divq;
    }

    public void setDivq(float divq) {
        this.divq = divq;
    }

    public float getV() {
        return v;
    }

    public void setV(float v) {
        this.v = v;
    }

    public float getT() {
        return t;
    }

    public void setT(float t) {
        this.t = t;
    }

    public char getOptionType() {
        return OptionType;
    }

    public void setOptionType(char optionType) {
        OptionType = optionType;
    }

    public float getDivs() {
        return divs;
    }

    public void setDivs(float divs) {
        this.divs = divs;
    }

    public float getDGrefval() {
        return DGrefval;
    }

    public void setDGrefval(float DGrefval) {
        this.DGrefval = DGrefval;
    }
}
