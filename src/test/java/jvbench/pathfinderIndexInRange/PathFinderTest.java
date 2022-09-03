package jvbench.pathfinderIndexInRange;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

class PathFinderTest {

    @BeforeEach
    void setUp() {
        PathFinder.init(System.getProperty("input", "src/main/resources/pathfinder/input/pathfinder_32_32.input"));
    }

    @Test
    void indexInRange() {
        System.out.println(Arrays.toString(PathFinder.scalar()));
        System.out.println(Arrays.toString(PathFinder.vector()));
        assertArrayEquals(PathFinder.scalar(), PathFinder.vector());
    }
}
