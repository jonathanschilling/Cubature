package de.labathome;

import de.labathome.Cubature.Region;

public class TestCubature {

	public static void main(String[] args) {
		testHyperCube();
		testCutRegion();
	}

	public static void testHyperCube() {
		int dim = 2;
		double[] center = new double[] { 1.0, 3.0 };
		double[] halfwidth = new double[] { 0.5, 0.5 };

		double expVol = 1.0;
		for (int i=0; i<dim; ++i) {
			expVol *= 2.0 * halfwidth[i];
		}

		Cubature.HyperCube h = new Cubature.HyperCube();
		h.initFromCenterAndHalfwidth(dim, center, halfwidth);

		System.out.println("volume: "+h.vol + " (expected "+expVol+")");
	}

	public static void testCutRegion() {
		int dim = 2;
		double[] center = new double[] { 1.0, 3.0 };
		double[] halfwidth = new double[] { 0.5, 0.5 };

		double expVol = 1.0;
		for (int i=0; i<dim; ++i) {
			expVol *= 2.0 * halfwidth[i];
		}

		Cubature.HyperCube h = new Cubature.HyperCube();
		h.initFromCenterAndHalfwidth(dim, center, halfwidth);

		int fdim = 1;

		Region R = new Region();
		R.init(h, fdim);
		R.splitDim = 1;
		
		Region R2 = R.clone();

		Cubature.cutRegion(R, R2);
		
		System.out.println("after splitting:\nvol 1 = "+R.h.vol+"\nvol 2 = "+R2.h.vol+"\n expected: 2*"+expVol/2.0);
	}

}
