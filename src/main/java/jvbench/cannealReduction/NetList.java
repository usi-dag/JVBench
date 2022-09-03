package jvbench.cannealReduction;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;


public class NetList {

    private int maxX;
    private int maxY;
    private int chipSize;
    private int unusedElem = 0;
    private NetListElem[] elements; //store the actual elements here
    private Location[][] locations; //store the actual locations here
    private final Map<String, NetListElem> elemNames = new HashMap<>();
    private int [] locationX;
    private int [] locationY;


    public NetList(String fileName) {
        try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
            String [] options = reader.readLine().split("\\s+");

            if (options.length != 3) {
                throw new IllegalArgumentException("Options element should be 3: NUM_ELEMENTS, MAX_X and, MAX_Y");
            }
            int numElements = Integer.parseInt(options[0]);
            maxX = Integer.parseInt(options[1]);
            maxY = Integer.parseInt(options[2]);
            chipSize = maxX * maxY;

            if (!(numElements < chipSize)) {
                throw new IllegalArgumentException("NUM_ELEMENTS > (maxX * maxY)");
            }

            elements = new NetListElem[chipSize];

            Location[] yVec = new Location[maxY];

            locations = new Location[maxX][maxY];
            for (int i = 0; i < maxX; i++) {
                locations[i] = yVec;
            }

            locationX = new int[chipSize];
            locationY = new int[chipSize];

            int iElem = 0;

            for (int x = 0; x < maxX; x++) {
                for (int y = 0; y < maxY; y++) {
                    Location loc = new Location();
                    loc.setX(x);
                    loc.setY(y);
                    locationX[iElem] = x;
                    locationY[iElem] = y;
                    NetListElem elem = new NetListElem(10, 10);
                    elem.setId(iElem);
                    elem.setPresentLoc(loc);
                    elements[iElem] = elem;
                    iElem++;
                }
            }


            int i = 0;
            for (String loc = reader.readLine(); loc != null; loc = reader.readLine()) {
                i++;
                String [] location = loc.split("\\s+");
                String name = location[0];

                NetListElem presentElem = createElemIfNecessary(name);

                presentElem.setItemName(name);

                int type = Integer.parseInt(location[1]);

                String faninName;

                for (int j = 2; j < location.length; j++) {
                    faninName = location[j];

                    if (Objects.equals(faninName, "END")) {
                        break;
                    }

                    NetListElem faninElem = createElemIfNecessary(faninName);
                    presentElem.addFaninElement(faninElem);
                    faninElem.addFanoutElement(presentElem);

                    // for vectorization
                    presentElem.addFaninIdElement(faninElem.getId());
                    faninElem.addFanoutIdElement(presentElem.getId());
                }
            }

        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
    }

    public RandomPair<NetListElem> getRandomPair(Rng rng) {
        long idA = rng.rand(chipSize);
        NetListElem elemA = elements[(int) idA];

        long idB = rng.rand(chipSize);
        NetListElem elemB = elements[(int) idB];

        while (idB == idA) {
            idB = rng.rand(chipSize);
            elemB = elements[(int) idB];
        }

        return new RandomPair<>(elemA, elemB);
    }

    public void swapLocations(NetListElem a, NetListElem b) {
        int idA = a.getId();
        int idB = b.getId();

        this.locationX[idA] = b.getPresentLoc().getX();
        this.locationY[idA] = b.getPresentLoc().getY();
        this.locationX[idB] = a.getPresentLoc().getX();
        this.locationY[idB] = a.getPresentLoc().getY();

        Location locA = a.getPresentLoc();
        Location locB = b.getPresentLoc();
        a.setPresentLoc(locB);
        b.setPresentLoc(locA);
    }

    public void shuffle(Rng rng) {
        for (int i  = 0; i < maxX * maxY * 1000; i++) {
            RandomPair<NetListElem> randomPair = getRandomPair(rng);
            NetListElem a =  randomPair.getFirst();
            NetListElem b = randomPair.getSecond();
            swapLocations(a, b);
        }
    }

    public NetListElem netlistElemFromName(String name) {
        return elemNames.get(name);
    }

    public double totalRoutingCost() {
        double rval = 0;
        for (NetListElem elem : elemNames.values()) {
            rval += elem.routingCostGivenLoc(elem.getPresentLoc());
        }
        return rval / 2;
    }

    public NetListElem getRandomElement(int differentFrom, Rng rng) {
        int id = (int) rng.rand(chipSize);
        NetListElem elem = elements[id];

        while (id == differentFrom) {
            id = (int) rng.rand(chipSize);
            elem = elements[id];
        }

        return elem;
    }


    protected NetListElem createElemIfNecessary(String name) {
        NetListElem rval;

        if (elemNames.containsKey(name)) {
            rval = elemNames.get(name);
        } else {
            rval = elements[unusedElem];
            elemNames.put(name, rval);
            unusedElem++;
        }
        return rval;
    }

    public NetListElem[] getElements() {
        return elements;
    }

    public int [] getLocationX() {
        return this.locationX;
    }

    public int [] getLocationY() {
        return this.locationY;
    }

}
