import edu.princeton.cs.algs4.Picture;

import java.util.ArrayList;


public class SeamCarver {
    private final static boolean FOR_VERTICAL_SEAM = false;
    private final static boolean FOR_HORISONTAL_SEAM = true;
    private Picture myPicture;
    private Matrix matrix; // has you

    // create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        if (picture == null) throw new IllegalArgumentException("1");
        myPicture = new Picture(picture);
        matrix = new Matrix(FOR_VERTICAL_SEAM);
    }

    // constructor for test purposes getting 16 bit rgb list
    /*private SeamCarver(String fileName, int h, int w) {
        In in = new In(fileName);
        myPicture = new Picture(w, h);
        matrix = new Matrix(FOR_VERTICAL_SEAM);
        // init matr with rgb
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                matrix.rgb[y][x] = Integer.parseInt(in.readString(), 16);
            }
        }
    }*/

    // current picture
    public Picture picture() {
        if (matrix.transposed) {
            if (width() != matrix.lenghtY() || height() != matrix.lenghtX()) updatePicture();
        }
        else {
            if (width() != matrix.lenghtX() || height() != matrix.lenghtY()) updatePicture();
        }
        updatePicture();
        return myPicture;
    }

    // width of current picture
    public int width() {
        return myPicture.width();
    }

    // height of current picture
    public int height() {
        return myPicture.height();
    }

    // hlpr return chanel int
    private int chanel(int x, int y, byte colour) {
        int chanel = -1;
        int rgbInt;
        rgbInt = matrix.rgb[y][x];


        switch (colour) {
            case 0:
                chanel = (rgbInt >> 16) & 0xFF;
                break;
            case 1:
                chanel = (rgbInt >> 8) & 0xFF;
                break;
            case 2:
                chanel = (rgbInt >> 0) & 0xFF;
        }
        return chanel;
    }

    // hlpr calc central difference of some chanel by specified coordinate
    private int diff(int x, int y, byte colour, byte coord) {
        int diff = -1;
        switch (coord) {
            case 0:
                diff = chanel(x + 1, y, colour) - chanel(x - 1, y, colour);
                break;
            case 1:
                diff = chanel(x, y + 1, colour) - chanel(x, y - 1, colour);
        }
        return diff;
    }

    // hlpr calc square of coord-gradient
    private double d2(int x, int y, byte coord) {
        byte RED = 0;
        byte GREEN = 1;
        byte BLUE = 2;
        double dx2 = (Math.pow(diff(x, y, RED, coord), 2) + Math
                .pow(diff(x, y, GREEN, coord), 2)
                + Math
                .pow(diff(x, y, BLUE, coord), 2));
        return dx2;
    }


    // energy of pixel at column x and row y
    public double energy(int x, int y) {
        if (matrix.transposed == FOR_HORISONTAL_SEAM) {
            int tempX = x;
            x = y;
            y = tempX;
        }
        if (x < 0 || y < 0 || x >= matrix.lenghtX() || y >= matrix.lenghtY())
            throw new IllegalArgumentException("1");
        if (x == 0 || y == 0 || x == matrix.lenghtX() - 1 || y == matrix.lenghtY() - 1)
            return 1000;
        byte X = 0;
        byte Y = 1;

        return Math.sqrt(d2(x, y, X) + d2(x, y, Y));
    }

    // private analog of energy method for mandatory recount of energy
    private double energyMand(int x, int y) {
        if (x < 0 || y < 0 || x >= matrix.lenghtX() || y >= matrix.lenghtY())
            throw new IllegalArgumentException("1");

        // if (matrix.matrE[y][x].e != -1) return matrix.matrE[y][x].e;
        if (x == 0 || y == 0 || x == matrix.lenghtX() - 1 || y == matrix.lenghtY() - 1)
            return 1000;
        byte X = 0;
        byte Y = 1;
        return Math.sqrt(d2(x, y, X) + d2(x, y, Y));
    }

    // energy matrix class
    private class Matrix {
        private int[][] rgb; // pixel color at RGB format (1 integer)
        private boolean transposed;

        public Matrix(boolean transposition) {
            rgb = new int[height()][width()];
            // init matr with empty energy
            for (int y = 0; y < height(); y++) {
                for (int x = 0; x < width(); x++) {
                    rgb[y][x] = myPicture.getRGB(x, y);
                }
            }

            if (transposition == FOR_HORISONTAL_SEAM)
                this.transpose();
            else transposed = FOR_VERTICAL_SEAM;

        }

        private void transpose(boolean orientation) {
            if (orientation != transposed) transpose();
        }

        private void transpose() {
            int oldX = rgb[0].length;
            int oldY = rgb.length;
            int[][] tempRgb = new int[oldX][oldY];

            for (int x = 0; x < oldX; x++) {
                for (int y = 0; y < oldY; y++) {
                    tempRgb[x][y] = rgb[y][x];
                }
            }
            rgb = tempRgb;
            if (transposed == FOR_VERTICAL_SEAM) transposed = FOR_HORISONTAL_SEAM;
            else transposed = FOR_VERTICAL_SEAM;
        }

        private int lenghtX() {
            return rgb[lenghtY() - 1].length;
        }


        private int lenghtY() {
            return rgb.length;
        }

        // removing spec pixel from 1 row
        private void removePixel(int x, int y) {
            int[] tempArray = new int[lenghtX() - 1]; // create array of less lenght on 1
            if (x > 0)
                System.arraycopy(rgb[y], 0, tempArray, 0, x); // copy old positions from left to x
            // copy right positions into new array begining from x+1
            if (x < lenghtX() - 1) {
                System.arraycopy(rgb[y], x + 1, tempArray, x, lenghtX() - (x + 1));
            }
            rgb[y] = tempArray;
        }
    }

    private class Sp {
        private int[][] verticeTo; // previous vertice
        private double[][] distTo; // distance to current pixel (sum energy)
        private final int xMax; // matrix X MAX
        private final int yMax; // matrix y MAX

        public Sp() {
            yMax = matrix.lenghtY();
            xMax = matrix.lenghtX();
            verticeTo = new int[yMax
                    + 1][xMax]; // for simplicity of further calculations 1st row is always = 0
            distTo = new double[yMax][xMax];
            for (int i = 0; i < xMax; i++) verticeTo[1][i] = i;
            for (int y = 0; y < yMax; y++) {
                for (int x = 0; x < xMax; x++) distTo[y][x] = Double.POSITIVE_INFINITY;
            }
            // relax all vert
            for (int y = 0; y < yMax; y++) {
                for (int x = 0; x < xMax; x++)
                    relax(x, y);
            }
        }

        // relax vertice
        private void relax(int x, int y) {
            if (y == 0) distTo[0][x] = energyMand(x, y);
            else {
                ArrayList<Integer> adjMap = adj(x);
                for (int adj : adjMap) {
                    // StdOut.println("x is..." + x);
                    if (distTo[y][x] > distTo[y - 1][adj] + energyMand(x, y)) {
                        distTo[y][x] = distTo[y - 1][adj] + energyMand(x, y);
                        verticeTo[y][x] = adj;
                    }
                }
            }
        }


        // find adjacent energies [Energy, Coordinate]
        private ArrayList<Integer> adj(int x) {
            ArrayList<Integer> energyMaps = new ArrayList<Integer>(3);
            if (x != 0) energyMaps.add(x - 1);
            energyMaps.add(x);
            if (x != xMax - 1) energyMaps.add(x + 1);
            return energyMaps;
        }

        // return function
        private int[] returnSeam() {
            int[] returnSeam = new int[yMax];
            double minDist = Double.POSITIVE_INFINITY;
            int firstX = -1;
            for (int i = 0; i < xMax; i++) {
                if (minDist > distTo[yMax - 1][i]) {
                    minDist = distTo[yMax - 1][i];
                    firstX = i;
                }
            }
            // last pixel is... (copy of before last)
            returnSeam[yMax - 1] = firstX;

            if (yMax > 1) {
                for (int i = yMax - 2; i >= 0; i--) {
                    returnSeam[i] = verticeTo[i + 1][returnSeam[i + 1]];
                }
            }
            return returnSeam;
        }
    }

    // sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        matrix.transpose(FOR_HORISONTAL_SEAM);
        Sp sp = new Sp();
        return sp.returnSeam();
    }


    // sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        matrix.transpose(FOR_VERTICAL_SEAM);
        Sp sp = new Sp();
        return sp.returnSeam();
    }

    // remove seam from matrix. Matrix should be oriented such that any seam will be vertical
    private void shrinkMatr(int[] seam) {
        for (int i = 0; i < matrix.lenghtY(); i++) {
            matrix.removePixel(seam[i], i);
        }
    }

    // update picture with new matrix
    private void updatePicture() {
        matrix.transpose(FOR_VERTICAL_SEAM);
        Picture tempPicture = new Picture(matrix.lenghtX(), matrix.lenghtY());
        for (int y = 0; y < matrix.lenghtY(); y++) {
            for (int x = 0; x < matrix.lenghtX(); x++) {
                tempPicture.setRGB(x, y, matrix.rgb[y][x]);
            }
        }
        myPicture = tempPicture;
    }

    // validate seam
    private boolean validSeam(int[] seam, int size) {

        // check on <0 or >size
        for (int s : seam) {
            if (s < 0 || s >= size) return true;
        }
        // chek on adjacent coordinates
        for (int i = 1; i < seam.length; i++) {
            if (seam[i] != seam[i - 1] && seam[i] != seam[i - 1] - 1
                    && seam[i] != seam[i - 1] + 1)
                return true;
        }
        return false;
    }

    // remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (seam == null) throw new IllegalArgumentException("1");
        if (seam.length != width()) throw new IllegalArgumentException("3");
        if (validSeam(seam, height())) throw new IllegalArgumentException("3.5");
        if (height() < 2) throw new IllegalArgumentException("4");
        matrix.transpose(FOR_HORISONTAL_SEAM);
        shrinkMatr(seam);
        updatePicture();
    }

    // remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (seam == null) throw new IllegalArgumentException("1");
        if (seam.length != height()) throw new IllegalArgumentException("3");
        if (validSeam(seam, width())) throw new IllegalArgumentException("3.5");
        if (width() < 2) throw new IllegalArgumentException("4");
        matrix.transpose(FOR_VERTICAL_SEAM);
        shrinkMatr(seam);
        updatePicture();
    }

    // wait method
    private void waitSec(double time) {
        double startTime = System.currentTimeMillis();
        while (true) {
            if ((System.currentTimeMillis() - startTime) > (time * 1000))
                return;
        }
    }


    public static void main(String[] args) {

        Picture picture = new Picture(args[0]);
        SeamCarver seam = new SeamCarver(picture);
        for (int i = 1; i < 200; i++) {
            seam.removeHorizontalSeam(seam.findHorizontalSeam());
            if (i % 10 == 0) {
                seam.picture().show();
                seam.waitSec(0.1);
            }
        }
    }
}
