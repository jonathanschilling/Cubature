package de.labathome;

import java.util.Locale;

public class ThreeDimGaussianMemberExample {
	
	private double sigma;
	
	public ThreeDimGaussianMemberExample(double newSigma) {
		sigma = newSigma;
	}
	
	public double[][] gaussianNd(double[][] x, Object fdata) {
		int dim = x.length;
		int nPoints = x[0].length;
		double[][] fval = new double[1][nPoints];
		for (int i = 0; i < nPoints; ++i) {
			double sum = 0.0;
			for (int d = 0; d < dim; ++d) {
				sum += x[d][i] * x[d][i];
			}
			fval[0][i] = Math.exp(-sigma * sum);
		}
		return fval;
	}

	public static void ex_ThreeDimGaussian() {
		double[] xmin = { -2.0, -2.0, -2.0 };
		double[] xmax = { 2.0, 2.0, 2.0 };

		double sigma = 0.5;
		ThreeDimGaussianMemberExample exampleInstance = new ThreeDimGaussianMemberExample(sigma);
		
		double[][] val_err = Cubature.integrate(exampleInstance, "gaussianNd",
				xmin, xmax,
				1.0e-4, 0.0, Cubature.Error.INDIVIDUAL,
				0,
				sigma);
		
		System.out.println(String.format(Locale.ENGLISH,
				"Computed integral = %.8f +/- %g", val_err[0][0], val_err[1][0]));
	}

	public static void main(String[] args) {
		ex_ThreeDimGaussian();
	}
}
