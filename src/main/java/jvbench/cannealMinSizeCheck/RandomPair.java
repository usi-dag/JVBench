package jvbench.cannealMinSizeCheck;

public record RandomPair<T>(T first, T second) {

    public T getFirst() {
        return first;
    }

    public T getSecond() {
        return second;
    }
}
