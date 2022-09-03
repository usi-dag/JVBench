package jvbench.streamcluster;

import org.junit.jupiter.api.Test;

class StreamClusterTest {

    String input = "src/main/resources/streamcluster/input/streamcluster_128_128.input";
    long kMin = 3;
    long kMax = 10;
    int dim =128;

    int chunkSize =   128;
    int clusterSize = 10;
    PStream stream;

    
    @Test
    void standard() {
        stream = StreamCluster.init(input);
        StreamCluster.streamCluster(stream, kMin, kMax, dim, chunkSize, clusterSize, "", false);
        int [] is_a_median = new int[(int) StreamCluster.centers.getNum()];
//        for( int i =0 ; i< StreamCluster.centers.getNum(); i++ ) {
//            is_a_median[StreamCluster.centers.getPoints()[i].getAssign()] = 1;
//        }
        float [] scalarWeights = new float[StreamCluster.centers.getNum()];
        float [][] scalarCoords = new float[StreamCluster.centers.getNum()][StreamCluster.centers.getDim()];
        for (int i = 0; i < StreamCluster.centers.getNum(); i++) {
//            if( is_a_median[i] != 0 ) {
                scalarWeights[i] = StreamCluster.centers.getPoints()[i].getWeight();

                for (int j = 0; j < StreamCluster.centers.getDim(); j++) {
                    scalarCoords[i][j] = StreamCluster.centers.getPoints()[i].getCoord(j, false);
                }
//            }
        }

        stream = StreamCluster.init(input);
        StreamCluster.streamCluster(stream, kMin, kMax, dim, chunkSize, clusterSize, "", true);

        float [] vectorWeights = new float[StreamCluster.centers.getNum()];
        float [][] vectorCoords = new float[StreamCluster.centers.getNum()][StreamCluster.centers.getDim()];
        for (int i = 0; i < StreamCluster.centers.getNum(); i++) {
            vectorWeights[i] = StreamCluster.centers.getPoints()[i].getWeight();

            for (int  j = 0; j < StreamCluster.centers.getDim(); j++) {
                vectorCoords[i][j] = StreamCluster.centers.getPoints()[i].getCoord(j, false);
            }
        }

//        System.out.println(Arrays.toString(scalarCoords[0]));
//        System.out.println(Arrays.toString(vectorCoords[0]));

//        assertArrayEquals(scalarWeights, vectorWeights);
//        assertArrayEquals(scalarCoords, vectorCoords);
    }

}