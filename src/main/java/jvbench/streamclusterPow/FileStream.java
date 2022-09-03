package jvbench.streamclusterPow;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class FileStream implements PStream {

    private String filename;
    private boolean isEOF = false;

    public FileStream(String filename) {
        this.filename = filename;
    }

    @Override
    public int read(float[] dest, int dim, int num) {

        try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
            int count = 0;
            String line;
            while ((line = reader.readLine()) != null && count < num) {

                String[] coordinate = line.split(" ");
                for (int k = 0; k < coordinate.length && k < dim; k++) {
                    dest[count * dim + k] = Float.parseFloat(coordinate[k]);
                }


                count++;
            }

            if (line == null) {
                isEOF = true;
            }

            return count;

        } catch (FileNotFoundException e) {
            System.err.println("File " + filename + " not found. " +  e.getMessage() );
        } catch (IOException e) {
            System.err.println("Unable to read file " + filename + ". " +  e.getMessage() );

        }
        return 0;
    }

    @Override
    public boolean feof() {
        return isEOF;
    }
}
