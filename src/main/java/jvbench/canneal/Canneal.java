package jvbench.canneal;


public class Canneal {

    // 1 15000 2000 input/400000.nets 128
    private int numThreads = 1;
    private int swapPerTemp = 10000;
    private int initialTemperature = 2000;
    private String netlistFilename = "/canneal/input/2500000.nets";
    private int numTempStep = 300;
    private static NetList netlist;
    private static AnnealerThread aT;
    private static boolean USE_VECTOR_API =true;


    public static int numberSwapCostOperation = 0;
    public static int meanLenOfArray = 0;

    public static int numberOfZeroVector = 0;


    public static void main(String[] args) {
        Canneal canneal = new Canneal();
        canneal.init();
        aT.run();
        System.out.println("final routing is: " + netlist.totalRoutingCost());

        System.out.println("number of swap: " + numberSwapCostOperation);
        System.out.println("mean len array " + meanLenOfArray / (numberSwapCostOperation * 2));
        System.out.println("# zero vector " + numberOfZeroVector);
    }

    public Canneal() {

    }

    public Canneal(int numThreads, int swapPerTemp, int initialTemperature, String netlistFilename, int numTempStep) {
        this.numThreads = numThreads;
        this.swapPerTemp = swapPerTemp;
        this.initialTemperature = initialTemperature;
        this.netlistFilename = netlistFilename;
        this.numTempStep = numTempStep;
    }


    public void init() {
        netlist = new NetList(netlistFilename);
        aT = new AnnealerThread(netlist, numThreads, swapPerTemp, initialTemperature, numTempStep);
    }

    public double getRoutingCost() {
        return netlist.totalRoutingCost();
    }

    public void scalar() {
        USE_VECTOR_API = false;
        aT.run();
    }

    public void vector() {
        USE_VECTOR_API = true;
        aT.run();
    }

    public static boolean useVectorAPI() {
        return USE_VECTOR_API;
    }


}
