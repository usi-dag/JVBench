package jvbench.pathfinder;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorSpecies;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.Random;

public class PathFinder {

    static private final VectorSpecies<Integer> SPECIES = IntVector.SPECIES_MAX;
    static private final int SPECIES_LENGTH = SPECIES.length();

    private static int rows;
    private static int cols;
    private static  int [] wall   ;
    private static  int [] result ;

    private static final int NUM_RUNS = 100;

    private static final Random rand = new Random();

    private String inputFileName;
    private String outFileName;

    public static void init(String inFilename) {

        try(BufferedReader reader = new BufferedReader(new InputStreamReader(PathFinder.class.getResourceAsStream(inFilename)))) {

            String [] params = reader.readLine().split(" ");
            rows = Integer.parseInt(params[0]);
            cols = Integer.parseInt(params[1]);

            wall = new int[rows * cols];
            result = new int[cols];

            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    wall[i * cols + j] = Integer.parseInt(reader.readLine());
                }
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }


    }


    public static int [] scalar() {
        int min;
        int [] dst = new int[cols];
        int [] src, temp;

        for (int j = 0; j < NUM_RUNS; j++) {
            src = new int[cols];
            for (int x = 0; x < cols; x++) {
                result[x] = wall[x];
            }

            dst = result;
            for (int t = 0; t < rows-1; t++) {
                temp = src;
                src = dst;
                dst = temp;
                for (int n = 0; n < cols; n++) {
                    min = src[n];
                    if (n > 0) {
//                    min = Math.min(min, src[n - 1]); // @IntrinsicCandidate
                        min = Math.min(min, src[n - 1]);
                    }
                    if (n < cols-1) {
                        min = Math.min(min, src[n + 1]);
                    }

                    dst[n] = wall[(t+1)*cols + n] + min;
                }
            }
        }

        return dst;
    }

    public static int [] vector() {
        int [] src;
        int [] dst = new int[cols];
        int [] temp;
        int limit = SPECIES.loopBound(cols - 1);

        for (int j=0; j<NUM_RUNS; j++) {
            src = new int[cols];
            for (int x = 0; x < cols; x++) {
                result[x] = wall[x];
            }
            dst = result;

            IntVector xSrcSlideUp;
            IntVector xSrcSlideDown;
            IntVector xSrc;
            IntVector xNextRow;

            for (int t = 0; t < rows-1; t++) {
                int n;

                temp = src;
                src = dst;
                dst = temp;

                if (src.length > 2) {
                    dst[0] = wall[(t + 1) * cols] + Math.min(src[0], src[1]);
                }

                for (n = 1; n < limit; n += SPECIES_LENGTH) {
                    xNextRow = IntVector.fromArray(SPECIES, src, n);
                    xSrc = xNextRow;
                    xSrcSlideUp = IntVector.fromArray(SPECIES, src, n+1);


                    xSrcSlideDown = IntVector.fromArray(SPECIES, src, n - 1);

                    xSrc = xSrc.min(xSrcSlideUp);
                    xSrc = xSrc.min(xSrcSlideDown);

                    xNextRow = IntVector.fromArray(SPECIES, wall, (t+1)*cols + n);
                    xNextRow = xNextRow.add(xSrc);
                    xNextRow.intoArray(dst, n);
                }


                for (; n < cols; n++) {
                    int min = src[n];
                    if (n > 0) {
                        min = Math.min(min, src[n - 1]);
                    }
                    if (n < cols-1) {
                        min = Math.min(min, src[n + 1]);
                    }

                    dst[n] = wall[(t+1)*cols + n] + min;
                }


            }
        }

        return dst;
    }
}
