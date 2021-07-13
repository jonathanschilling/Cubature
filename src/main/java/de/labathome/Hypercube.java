package de.labathome;

public class Hypercube {

	int dim;
	double[] centers;
	double[] halfwidths;
	double volume;

	public Hypercube initFromCenterHalfwidths(double[] centers, double[] halfwidths) {
		dim = centers.length;
		this.centers = centers.clone();
		this.halfwidths = halfwidths.clone();
		computeVolume();
		return this;
	}

	public Hypercube initFromRanges(double[] xmin, double[] xmax) {
		dim = xmin.length;
		centers = new double[dim];
		halfwidths = new double[dim];
		for(int i=0; i<dim; ++i) {
			centers[i]    = (xmax[i]+xmin[i])/2.0;
			halfwidths[i] = (xmax[i]-xmin[i])/2.0;
		}
		computeVolume();
		return this;
	}

	public Hypercube clone() {
		return new Hypercube().initFromCenterHalfwidths(centers, halfwidths);
	}

	private void computeVolume() {
		volume = 1.0;
		for (int i=0; i<dim; ++i) {
			volume *= halfwidths[i]*2.0;
		}
	}
}