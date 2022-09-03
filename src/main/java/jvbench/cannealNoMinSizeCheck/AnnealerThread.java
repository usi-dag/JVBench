package jvbench.cannealNoMinSizeCheck;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorSpecies;

public class AnnealerThread implements Runnable {

    private final NetList netlist;
    private boolean keepGoingGlobalFlag;
    private final int movesPerThreadTemp;
    private final int startTemp;
    private final int numberTempSteps;
    private boolean [] mask;

    private final static VectorSpecies<Integer> SPECIES = IntVector.SPECIES_MAX;
    private final static int SPECIES_LENGTH = SPECIES.length();

    public AnnealerThread(NetList netlist, int nThreads, int swapPerTemp, int startTemp, int numberTempSteps) {
        this.netlist = netlist;
        keepGoingGlobalFlag = true;
        movesPerThreadTemp = swapPerTemp/nThreads;
        this.startTemp = startTemp;
        this.numberTempSteps = numberTempSteps;
    }

    @Override
    public void run() {
        int acceptedGoodMoves = 0;
        int acceptedBadMoves = -1;
        double T = startTemp;
        Rng rng = new Rng();

        NetListElem a = netlist.getRandomElement(Elem.NO_MATCHING_ELEMENT.code, rng);
        NetListElem b = netlist.getRandomElement(Elem.NO_MATCHING_ELEMENT.code, rng);

        int tempStepsCompleted = 0;

        if (Canneal.useVectorAPI()) {
            mask = new boolean[SPECIES.length()];
            for (int i=0; i< SPECIES_LENGTH; i+=2) {mask[i] = true; mask[i+1] = false;}
        }

        while (keepGoing(tempStepsCompleted, acceptedGoodMoves, acceptedBadMoves)) {
            T = T / 1.5;
            acceptedGoodMoves = 0;
            acceptedBadMoves = 0;

            for (int i = 0; i < movesPerThreadTemp; i++) {
                a = b;
                b = netlist.getRandomElement(a.getId(), rng);

                double deltaCost;
                if (Canneal.useVectorAPI()) {
                    deltaCost = calculateDeltaRoutingCostVector(a, b);
                } else {
                    deltaCost = calculateDeltaRoutingCost(a, b);
                }

                MoveDecision isGoodMove = acceptMove(deltaCost, T, rng);
                if (isGoodMove == MoveDecision.ACCEPTED_BAD) {
                    acceptedBadMoves++;
                    netlist.swapLocations(a, b);
                } else if (isGoodMove == MoveDecision.ACCEPTED_GOOD) {
                    acceptedGoodMoves++;
                    netlist.swapLocations(a, b);
                } else if (isGoodMove == MoveDecision.REJECTED) {}
            }
            tempStepsCompleted++;
        }
    }

    private MoveDecision acceptMove(double deltaCost, double t, Rng rng) {
        if (deltaCost < 0) {
            return MoveDecision.ACCEPTED_GOOD;
        } else {
            double randomValue = rng.drand();
            double boltzman = Math.exp(- deltaCost/t);
            if (boltzman > randomValue) {
                return MoveDecision.ACCEPTED_BAD;
            } else {
                return MoveDecision.REJECTED;
            }
        }
    }

    private double calculateDeltaRoutingCostVector(NetListElem a, NetListElem b) {
        double deltaCost = 0.0;

       Location aLoc = a.getPresentLoc();
       Location bLoc = b.getPresentLoc();
        deltaCost = a.swapCostVector(aLoc, bLoc, netlist.getLocationX(), netlist.getLocationY());
        deltaCost = deltaCost + b.swapCostVector(bLoc, aLoc, netlist.getLocationX(), netlist.getLocationY());
        return deltaCost;
    }

    private double calculateDeltaRoutingCost(NetListElem a, NetListElem b) {
        double deltaCost = 0.0;

        Location aLoc = a.getPresentLoc();
        Location bLoc = b.getPresentLoc();
        deltaCost = a.swapCost(aLoc, bLoc);
        deltaCost = deltaCost + b.swapCost(bLoc, aLoc);

        return deltaCost;
    }

    private boolean keepGoing(int tempStepsCompleted, int acceptedGoodMoves, int acceptedBadMoves) {
        boolean rv;

        if (numberTempSteps == -1) {
            rv = keepGoingGlobalFlag && (acceptedGoodMoves > acceptedBadMoves);
            if (!rv) {
                keepGoingGlobalFlag = false;
            }
        } else {
            rv = tempStepsCompleted < numberTempSteps;
        }

        return rv;
    }
}
