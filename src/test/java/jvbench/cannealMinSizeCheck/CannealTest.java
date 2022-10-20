package jvbench.cannealMinSizeCheck;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class CannealTest {

    Canneal scalarCanneal;
    Canneal vectorCanneal;

    @BeforeEach
    public void setUp() {
        scalarCanneal = new Canneal(
                1,
                15000,
                2000,
                "/canneal/input/100.nets",
                128
        );

        vectorCanneal = new Canneal(
                1,
                15000,
                2000,
                "/canneal/input/100.nets",
                128
        );

        scalarCanneal.init();
        vectorCanneal.init();

        scalarCanneal.scalar();
        vectorCanneal.vector();
    }


    @Test
    public void minSizeCheck() {
        assertEquals(scalarCanneal.getRoutingCost(), vectorCanneal.getRoutingCost());
    }

}