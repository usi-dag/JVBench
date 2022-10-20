package jvbench.pathfinderConditionInsideLoop;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class PathFinderTest {

    @BeforeEach
    void setUp() {
        PathFinder.init(System.getProperty("input", "/pathfinder/input/pathfinder_32_32.input"));
    }

    @Test
    void conditionInsideLoop() {
        assertArrayEquals(PathFinder.scalar(), PathFinder.vector());
    }
}
