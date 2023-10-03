package de.labathome;

public class Hypercube {

    int dim;
    double[] centers;
    double[] halfwidths;
    double volume;

    public static Hypercube initFromCenterHalfwidths(double[] centers, double[] halfwidths) {
        Hypercube h = new Hypercube();
        h.dim = centers.length;
        h.centers = centers.clone();
        h.halfwidths = halfwidths.clone();
        h.computeVolume();
        return h;
    }

    public static Hypercube initFromRanges(double[] xmin, double[] xmax) {
        Hypercube h = new Hypercube();
        h.dim = xmin.length;
        h.centers = new double[h.dim];
        h.halfwidths = new double[h.dim];
        for (int i = 0; i < h.dim; ++i) {
            h.centers[i] = (xmax[i] + xmin[i]) / 2.0;
            h.halfwidths[i] = (xmax[i] - xmin[i]) / 2.0;
        }
        h.computeVolume();
        return h;
    }

    public Hypercube clone() {
        return Hypercube.initFromCenterHalfwidths(centers, halfwidths);
    }

    private void computeVolume() {
        volume = 1.0;
        for (int i = 0; i < dim; ++i) {
            volume *= halfwidths[i] * 2.0;
        }
    }
}