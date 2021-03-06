package tests;

import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import org.junit.jupiter.api.Test;

import de.labathome.Cubature;
import de.labathome.CubatureError;

public class TestCubature {

	public static final double K_2_SQRTPI = 1.12837916709551257390;

	/**
	 * starting point for porting the original test integrands from test.c of stevenj/cubature into here
	 * @param x [dim][nPoints] evaluation points
	 * @param fdata additional data for integrand
	 * @return [fdim][nPoints] function values at x
	 */
	public static double[][] cubatureIntegrand(double[][] x, Object fdata) {
		int[] integrands = (int[])fdata;
		int dim = x.length;
		int nPoints = x[0].length;
		double[][] ret = new double[integrands.length][nPoints];

		for (int fdim=1; fdim<=integrands.length; ++fdim) {
			for (int p=0; p<nPoints; ++p) {
				double val = Double.NaN;

				if (integrands[fdim-1] == 0) {
					/* simple smooth (separable) objective: prod. cos(x[i]). */
					val = 1.0;
					for (int i = 0; i<dim; ++i) {
						val *= Math.cos(x[i][p]);
					}
				} else if (integrands[fdim-1] == 1) {
					/* integral of exp(-x^2), rescaled to (0,infinity) limits */
					double scale = 1.0;
					val = 0;
					for (int i = 0; i < dim; ++i) {
						if (x[i][p] > 0) {
							double z = (1 - x[i][p]) / x[i][p];
							val += z * z;
							scale *= K_2_SQRTPI / (x[i][p] * x[i][p]);
						} else {
							scale = 0;
							break;
						}
					}
					val = Math.exp(-val) * scale;
				}

				ret[fdim][p] = val;
			}
		}

		return ret;
	}

	public static class F0_class {
		int dim;
		public F0_class(int mydim) {
			dim = mydim;
		}

		public double[][] eval(double[][] x, Object fdata) {
			double[] ret = new double[x[0].length];
			for (int d=0; d<dim; ++d) {
				for (int i=0; i<x[0].length; ++i) {
					ret[i] = 2.0*x[d][i];
				}
			}
			return new double[][] { ret };
		}
	}

	public static class F1_class {
		double a = 0.1;
		public double[][] eval(double[][] x, Object fdata) {
			double[] ret = new double[x[0].length];
			for (int i=0; i<x[0].length; ++i) {
				double dx = x[0][i] - 0.5;
				double sum = dx * dx;
				ret[i] = K_2_SQRTPI / (2. * a) * Math.exp(-sum / (a * a));
			}
			return new double[][] { ret };
		}
	}

	public static class CosClass {
		int dim;
		public CosClass(int mydim) {
			dim = mydim;
		}
		public double[][] eval(double[][] x, Object fdata) {
			double[] ret = new double[x[0].length];
			for (int i=0; i<x[0].length; ++i) {
				ret[i] = 1.0;
				for (int d=0; d<dim; ++d) {
					ret[i] *= Math.cos(x[d][i]);
				}
			}
			return new double[][] { ret };
		}
	}

	public static class ExpClass {
		public double tau;
		public ExpClass(double tau) {
			this.tau = tau;
		}
		public double[][] tauExp(double[][] x) {
			double[] ret = new double[x[0].length];
			for (int i=0; i<x[0].length; ++i) {
				ret[i] = Math.exp(x[0][i]/tau);
			}
			return new double[][] { ret };
		}
	}

	@Test
	public void demoIntegration() {

//		ExpClass c1 = new ExpClass(1.0);
		CosClass cc = new CosClass(2);
//		F0_class f0 = new F0_class(2);
//		F1_class f1 = new F1_class();

		double tol = 1.0e-9;

		double[][] F = Cubature.integrate(cc, "eval",
				new double[] { 0.0, 0.0 }, // lower integration limit
				new double[] { 1.0, 1.0 }, // upper integration limit
				tol, // relative tolerance
				0.0, // no absolute tolerance requirement
				CubatureError.INDIVIDUAL,
				100000, // max. number of function evaluations
				null);  // no extra info

		if (F != null) {
			System.out.println("integration result: "+F[0][0]+" +/- " + F[1][0]);
			assertTrue(F[1][0] < tol);
		} else {
			System.out.println("no integration performed");
			fail("no integration performed");
		}
	}
}
