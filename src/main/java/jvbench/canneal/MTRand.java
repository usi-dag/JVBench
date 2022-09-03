package jvbench.canneal;

import java.util.Optional;

public class MTRand {

    // typodef uint32 == unsigned long
    public static int N = 624; // length of first state vector
    public static int SAVE = N + 1; // length of array for save()
    protected static int M = 397;
    protected static int left;  // number of values left before reload needed
    protected long [] state = new long[SAVE]; // internal state
    protected int pNext; // next value to get from state


    // initialize with first simple uint32
    public MTRand(long oneSeed) {
        seed(oneSeed);
    }


    // auto-initialize with /dev/urandom or time() and clock()
//    public MTRand() {
//        seed();
//    }

    // or an array
    public MTRand(long [] bigSeed, Optional<Long> seedLength) {
        seed(bigSeed, seedLength);
    }

    public double rand() {
        return ((double) randInt()) * (1.0/4294967295.0);
    }

    public double rand(double n) {
        return rand() * n;
    }

    public double randExc() {
        return randInt() * (1.0/4294967296.0);
    }

    public double randExc(double n) {
        return randExc() * n;
    }

    double randDblExc() {
        return (randInt() + 0.5 ) * (1.0/4294967296.0);
    }

    double randDblExc(double n) {
        return  randDblExc() * n;
    }

    public long randInt() {

        if (left == 0) reload();

        --left;

        long s1;
        s1 = state[pNext++];
        s1 ^= (s1 >> 11);
        s1 ^= (s1 <<  7) & 0x9d2c5680L;
        s1 ^= (s1 << 15) & 0xefc60000L;
        return ( s1 ^ (s1 >> 18) );
    }

    public long randInt(long n) {
        long used = n;
        used |= used >> 1;
        used |= used >> 2;
        used |= used >> 4;
        used |= used >> 8;
        used |= used >> 16;

        // Draw numbers until one is found in [0,n]
        long i;
        do {
            i = randInt() & used;  // toss unused bits to shorten search
        }

        while( i > n );
        return i;
    }

    public double rand53() {
        long a = randInt() >> 5;
        long b = randInt() >> 6;

        return ( a * 67108864.0 + b ) * (1.0/9007199254740992.0);
    }

    public double randNorm(Optional<Double> mean, Optional<Double> variance) {
        double r = Math.sqrt(-2.0 * Math.log(1.0-randDblExc())) * variance.orElse(0.0);
        double phi = 2.0 * 3.14159265358979323846264338328 * randExc();
        return mean.orElse(0.0) + r * Math.cos(phi);
    }

    public void seed(long oneSeed) {
        initialize(oneSeed);
        reload();
    }

    public void seed(long  [] bigSeed, Optional<Long> seedLength) {
        initialize(1965021L);
        int i = 1;
        int j = 0;
        int k = (int) (N < seedLength.orElse(0L) ? N : seedLength);
        for( ; k > 0; --k ) {
            state[i] = state[i] ^ ( (state[i-1] ^ (state[i-1] >> 30)) * 1664525L );
            state[i] += ( bigSeed[j] & 0xffffffffL ) + j;
            state[i] &= 0xffffffffL;
            ++i;  ++j;
            if( i >= N ) { state[0] = state[N-1];  i = 1; }
            if( j >= seedLength.orElse(0L) ) j = 0;
        }
        for( k = N - 1; k > 0; --k ) {
            state[i] = state[i] ^ ( (state[i-1] ^ (state[i-1] >> 30)) * 1566083941L );
            state[i] -= i;
            state[i] &= 0xffffffffL;
            ++i;
            if( i >= N ) { state[0] = state[N-1];  i = 1; }
        }
        state[0] = 0x80000000L;  // MSB is 1, assuring non-zero initial array
        reload();
    }

    public void save(long [] saveArray) {
        for (int i = 0; i < N; i++) {
            saveArray[i] = state[i];
        }
        saveArray[N] = left;
    }

    public void load(long [] loadArray) {
        for (int i = 0; i < N; i++) {
            state[i] = loadArray[i];
        }

        left = (int) loadArray[N];
        pNext = N-left;
    }

    protected void initialize(long oneSeed) {
        state[0] =  oneSeed & 0xffffffffL;
        int i = 1;
        for (; i < N; i++) {
            state[i] =  ( 1812433253L * ( state[i-1] ^ (state[i-1] >> 30) ) + i ) & 0xffffffffL;
        }
    }

    protected void reload() {
        int i;
        int j = 0;

        for (i = 0; i < N - M; i++) {
            state[j] = twist(state[j+M], state[j], state[j + 1]);
            j++;
        }

        for (i = 0; i < M - 1; i++) {
            state[j] = twist(state[j + (M-N)], state[j], state[j+1]);
            j++;
        }

        state[j] = twist(state[j + (M-N)], state[j], state[j+1]);
        left = N;
        pNext = 0;
    }

    protected long hiBit(long u) {
        return u & 0x80000000L;
    }

    protected long loBit(long u) {
        return u & 0x00000001L;
    }

    protected long loBits(long u) {
        return u & 0x7fffffffL;
    }

    protected long mixBits(long u, long v) {
        return hiBit(u) | loBits(v);
    }

    protected long twist(long m, long s0, long s1) {
        return m ^ (mixBits(s0,s1)>>1) ^ (-loBit(s1) & 0x9908b0dfL);
    }

//    static long hash( time_t t, clock_t c );

}
