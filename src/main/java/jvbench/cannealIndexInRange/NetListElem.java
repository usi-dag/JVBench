package jvbench.cannealIndexInRange;

import jdk.incubator.vector.IntVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

import java.util.Arrays;

public class NetListElem {

    private int id;
    private String itemName;
    private NetListElem [] fanin =  new NetListElem[16];
    private NetListElem [] fanout =   new NetListElem[16];
    private int [] faninId;
    private int [] fanoutId;
    private Location presentLoc;
    private int faninIndex = 0;
    private int fanoutIndex = 0;
    protected int faninIdIndex = 0;
    protected int fanoutIdIndex = 0;



    private final static VectorSpecies<Integer> INT_SPECIES = IntVector.SPECIES_MAX;
    private final static int INT_SPECIES_LENGTH = INT_SPECIES.length();

    private static final IntVector faninZeroV = IntVector.zero(INT_SPECIES);
    private static final IntVector fanoutZeroV = IntVector.zero(INT_SPECIES);




    public NetListElem(int xLen, int yLen) {
        faninId = new int[INT_SPECIES_LENGTH];
        fanoutId = new int[INT_SPECIES_LENGTH];

    }

    /**
     * Calculates the routing cost using the manhatten distance
     * @param loc the location to compute the cost
     * @return the cost of the location
     */
    public double routingCostGivenLoc(Location loc) {
        double faninCost = 0.0;
        double fanoutCost = 0.0;

        for (int i = 0; i < faninIndex; ++i) {
            Location faninLoc = fanin[i].getPresentLoc();
            faninCost += Math.abs(loc.getX() - faninLoc.getX());
            faninCost += Math.abs(loc.getY() - faninLoc.getY());
        }


        for (int i = 0; i < fanoutIndex; ++i) {
            Location fanoutLoc = fanout[i].getPresentLoc();
            fanoutCost += Math.abs(loc.getX() - fanoutLoc.getX());
            fanoutCost += Math.abs(loc.getY() - fanoutLoc.getY());
        }

        return faninCost + fanoutCost;
    }

    public double swapCostVector(Location oldLoc, Location newLoc, int [] locationX, int [] locationY) {

//        Canneal.numberSwapCostOperation +=1;


        int faninSize = faninIdIndex;
        int fanoutSize = fanoutIdIndex;


        double noSwap = 0;
        double yesSwap = 0;


        IntVector noSwapVector = faninZeroV;
        IntVector yesSwapVector = fanoutZeroV;
        IntVector oldLocX = IntVector.broadcast(INT_SPECIES, oldLoc.getX());
        IntVector oldLocY = IntVector.broadcast(INT_SPECIES, oldLoc.getY());
        IntVector newLocX = IntVector.broadcast(INT_SPECIES, newLoc.getX());
        IntVector newLocY = IntVector.broadcast(INT_SPECIES, newLoc.getY());
//        VectorMask<Integer> faninMask = INT_SPECIES.indexInRange(0, faninIdIndex);
//        VectorMask<Integer> fanoutMask = INT_SPECIES.indexInRange(0, fanoutIdIndex);

        int i = 0;
        for (; i < faninSize; i+= INT_SPECIES_LENGTH) {
            VectorMask<Integer> mask = INT_SPECIES.indexInRange(i, faninSize);
            IntVector faninLocVX = IntVector.fromArray(INT_SPECIES, locationX, 0, faninId, i, mask);
            IntVector faninLocVY = IntVector.fromArray(INT_SPECIES, locationY, 0, faninId, i, mask);

            noSwap += (oldLocX.sub(faninLocVX)).abs().add((oldLocY.sub(faninLocVY)).abs()).reduceLanes(VectorOperators.ADD, mask);

            yesSwap += (newLocX.sub(faninLocVX)).abs().add((newLocY.sub(faninLocVY)).abs()).reduceLanes(VectorOperators.ADD, mask);

        }


//        noSwap += noSwapVector.reduceLanes(VectorOperators.ADD);
//        yesSwap += yesSwapVector.reduceLanes(VectorOperators.ADD);

//        noSwapVector = IntVector.zero(INT_SPECIES);
//        yesSwapVector = IntVector.zero(INT_SPECIES);

        i=0;
        for (; i < fanoutSize; i+= INT_SPECIES_LENGTH) {
            VectorMask<Integer> mask = INT_SPECIES.indexInRange(i, faninSize);
            IntVector fanoutLocVX = IntVector.fromArray(INT_SPECIES,  locationX, 0, fanoutId, i, mask);
            IntVector fanoutLocVY = IntVector.fromArray(INT_SPECIES, locationY, 0, fanoutId, i, mask);

            noSwap += (oldLocX.sub(fanoutLocVX)).abs().add((oldLocY.sub(fanoutLocVY)).abs()).reduceLanes(VectorOperators.ADD, mask);

            yesSwap += (newLocX.sub(fanoutLocVX)).abs().add((newLocY.sub(fanoutLocVY)).abs()).reduceLanes(VectorOperators.ADD, mask);
        }

        return yesSwap - noSwap;
    }

    public double swapCost(Location oldLoc, Location newLoc) {

        int faninSize = faninIndex;
        int fanoutSize = fanoutIndex;
        double noSwap = 0;
        double yesSwap = 0;

        for (int i = 0; i < faninSize; i++) {
            Location faninLoc = fanin[i].getPresentLoc();

            noSwap += Math.abs(oldLoc.getX() - faninLoc.getX());
            noSwap += Math.abs(oldLoc.getY() - faninLoc.getY());

            yesSwap += Math.abs(newLoc.getX() - faninLoc.getX());
            yesSwap += Math.abs(newLoc.getY() - faninLoc.getY());

        }

        for (int i = 0; i < fanoutSize; i++) {
            Location fanoutLoc = fanout[i].getPresentLoc();

            noSwap += Math.abs(oldLoc.getX() - fanoutLoc.getX());
            noSwap += Math.abs(oldLoc.getY() - fanoutLoc.getY());

            yesSwap += Math.abs(newLoc.getX() - fanoutLoc.getX());
            yesSwap += Math.abs(newLoc.getY() - fanoutLoc.getY());

        }
        return yesSwap - noSwap;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public String getItemName() {
        return itemName;
    }

    public void setItemName(String itemName) {
        this.itemName = itemName;
    }

    public NetListElem [] getFanin() {
        return fanin;
    }

    public void setFanin(NetListElem [] fanin) {
        this.fanin = fanin;
    }

    public NetListElem [] getFanout() {
        return fanout;
    }

    public void setFanout(NetListElem [] fanout) {
        this.fanout = fanout;
    }

    public Location getPresentLoc() {
        return presentLoc;
    }

    public void setPresentLoc(Location presentLoc) {
        this.presentLoc = presentLoc;
    }

    public void addFaninElement(NetListElem element) {
        if (faninIndex == this.fanin.length) {
            this.fanin = Arrays.copyOf(this.fanin, faninIndex * 2);
        }

        this.fanin[faninIndex++] = element;
    }

    public void addFanoutElement(NetListElem element) {
        if (fanoutIndex == this.fanout.length) {
            this.fanout = Arrays.copyOf(this.fanout, fanoutIndex * 2);
        }

        this.fanout[fanoutIndex++] = element;
    }

    public void addFaninIdElement(int id) {
        if (faninIdIndex == this.faninId.length) {
            this.faninId = Arrays.copyOf(this.faninId, faninIdIndex * 2);
        }

        this.faninId[faninIdIndex++] = id;
    }

    public void addFanoutIdElement(int id) {
        if (fanoutIdIndex == this.fanoutId.length) {
            this.fanoutId = Arrays.copyOf(this.fanoutId, fanoutIdIndex * 2);
        }

        this.fanoutId[fanoutIdIndex++] = id;
    }

}
