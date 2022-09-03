package jvbench.swaptions;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class HJMSecuritiesImplTest {

    private static int ms = 16384;
    private static int ns = 1;

    @BeforeEach
    public void setUp() {
        HJMSecuritiesImpl.init(ms, ns);
    }

    @Test
    void scalar() {
        HJMSecuritiesImpl.benchmark(false);

        for (int i = 0; i < ns; i++) {
            assertTrue(HJMSecuritiesImpl.swaptions[i].dSimSwaptionStdError < 0.001);
//            System.out.printf("Swaption %d: [SwaptionPrice: %.10f StdError: %.10f] \n",
//                    i, HJMSecuritiesImpl.swaptions[i].dSimSwaptionMeanPrice, HJMSecuritiesImpl.swaptions[i].dSimSwaptionStdError);
        }
    }

    @Test
    void vector() {
        HJMSecuritiesImpl.benchmark(true);
        for (int i = 0; i < ns; i++) {
            assertTrue(HJMSecuritiesImpl.swaptions[i].dSimSwaptionStdError < 0.001);
//            System.out.printf("Swaption %d: [SwaptionPrice: %.10f StdError: %.10f] \n",
//                    i, HJMSecuritiesImpl.swaptions[i].dSimSwaptionMeanPrice, HJMSecuritiesImpl.swaptions[i].dSimSwaptionStdError);
        }
    }

}