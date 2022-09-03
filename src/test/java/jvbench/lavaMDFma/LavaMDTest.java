package jvbench.lavaMDFma;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class LavaMDTest {

    @BeforeEach
    public void setup() {
        LavaMD.init("src/main/resources/lavaMD/input/lavaMD_96.input");
    }

    @Test
    void fma() {
        LavaMD.serial();
        FourVector[] serial = LavaMD.fvCPU;
        LavaMD.vector();
        FourVector[] vector = LavaMD.fvCPUV;

        for (int i = 0; i < vector.length; i++) {
            assertTrue(serial[i].v - 0.001 <= vector[i].v &&  vector[i].v <= serial[i].v + 0.001);
            assertTrue(serial[i].x - 0.001 <= vector[i].x &&  vector[i].x <= serial[i].x + 0.001);
            assertTrue(serial[i].y - 0.001 <= vector[i].y &&  vector[i].y <= serial[i].y + 0.001);
            assertTrue(serial[i].z - 0.001 <= vector[i].z &&  vector[i].z <= serial[i].z + 0.001);
        }
    }
}