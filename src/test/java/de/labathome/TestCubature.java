package de.labathome;

import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import org.junit.jupiter.api.Test;

public class TestCubature {

	public static final double K_2_SQRTPI = 1.12837916709551257390;

	
	public static double[][] kbf_static(double[][] x, Object fdata) {
		double[] ret = new double[x[0].length];
		for (int i=0; i<x[0].length; ++i) {
			//double p1 = Math.cos(Math.pow(x[0][i], x[1][i]));
			double p1 = Math.cos(x[0][i] + x[1][i]);
			double p2 = Math.sin(x[0][i]*x[0][i]);
			double p3 = Math.exp(- x[0][i]*x[0][i] - x[1][i]);
			ret[i] = p1*p1 * p2 * p3;
		}
		return new double[][] { ret };
	}

	public static double[][] f_Zw(double[][] param, Object fdata) {
		int nPoints = param[0].length;
		double[] ret = new double[nPoints];
		for (int i=0; i<nPoints; ++i) {
			double x     = param[0][i];
			double y     = param[1][i];
			double theta = param[1][i];

			ret[i] = (x*x+y)*(x*x+y) * Math.sinh(theta) 
					* (Math.cos(theta) + Math.sin(theta*theta))
					* (Math.cos(theta) + Math.sin(theta*theta))
					* (Math.cos(theta) + Math.sin(theta*theta))
					/ ( theta * (Math.exp(-x*x/(y*y)) + Math.exp(x/y)) );
		}
		return new double[][] { ret };
	}

	public static double[][] f_Ef(double[][] param, Object fdata) {
		int nPoints = param[0].length;
		double[] ret = new double[nPoints];
		for (int i=0; i<nPoints; ++i) {
			double x     = param[0][i];
			ret[i] = Math.pow(x, Math.pow(x, x));
		}
		return new double[][] { ret };
	}

	public static double[][] f_Jkngk(double[][] param, Object fdata) {
		int nPoints = param[0].length;
		double[][] ret = new double[5][nPoints];
		for (int i=0; i<nPoints; ++i) {
			double alpha = param[0][i];
			double beta  = param[1][i];

			ret[0][i] = alpha + beta*beta*Math.cos(alpha);
			ret[1][i] = Math.sin(alpha-beta) * Math.cos(beta - 2.0*alpha);
			ret[2][i] = Math.exp(beta - alpha) * alpha;
			ret[3][i] = 1.0;
			//ret[4][i] = (alpha-beta)/ (alpha + beta*beta) * (Math.exp(-alpha) + Math.exp(-beta));
			ret[4][i] = Math.exp(-alpha) + Math.exp(-beta);
		}
		return ret;
	}

	//	
	//	@Test
	//	public void testZw() {
	//		double tol = 1.0e-9;
	//
	//		double[][] F = Cubature.integrate(TestCubature.class, "f_Zw",
	//				new double[] { -Math.PI/2.0, -1.0,    0.0  }, // lower integration limit
	//				new double[] {  Math.PI/2.0,  1.0, Math.PI }, // upper integration limit
	//				tol, // relative tolerance
	//				0.0, // no absolute tolerance requirement
	//				Cubature.Error.INDIVIDUAL,
	//				100000, // max. number of function evaluations
	//				null);  // no extra info
	//		if (F != null) {
	//			System.out.println("integration result: "+F[0][0]+" +/- " + F[1][0]);
	//			assertTrue(F[1][0] < tol);
	//		} else {
	//			System.out.println("no integration performed");
	//			fail("no integration performed");
	//		}
	//	}

	@Test
	public void testEf() {
		double tol = 1.0e-12;

		double[][] F = Cubature.integrate(TestCubature.class, "f_Ef",
				new double[] { 1.0 }, // lower integration limit
				new double[] { 2.0 }, // upper integration limit
				tol, // relative tolerance
				0.0, // no absolute tolerance requirement
				Cubature.Error.INDIVIDUAL,
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

	@Test
	public void testJkngk() {
		double tol = 1.0e-9;

		double[][] F = Cubature.integrate(TestCubature.class, "f_Jkngk",
				new double[] { 0.0, -10.0 }, // lower integration limit
				new double[] { 1.0,  10.0 }, // upper integration limit
				tol, // relative tolerance
				tol, // no absolute tolerance requirement
				Cubature.Error.INDIVIDUAL,
				1000000, // max. number of function evaluations
				null);  // no extra info
		if (F != null) {
			System.out.println("integration result: "
					+                    "["+F[0][0]+","+F[0][1]+","+F[0][2]+","+F[0][3]+","+F[0][4]+"]\n"
					+"                +/- ["+F[1][0]+","+F[1][1]+","+F[1][2]+","+F[1][3]+","+F[1][4]+"]");
			assertTrue(F[1][0] < tol);
			assertTrue(F[1][1] < tol);
			assertTrue(F[1][2] < tol);
			assertTrue(F[1][3] < tol);
			assertTrue(F[1][4] < tol);
		} else {
			System.out.println("no integration performed");
			fail("no integration performed");
		}
	}

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

	@Test
	public void testCubature() {
		int fdimsToTest = 3;

		for (int maxFDim=1; maxFDim<=fdimsToTest; ++maxFDim) {
			int[] integrands = new int[maxFDim];
			for (int fdim=1; fdim<=maxFDim; ++fdim) {
				integrands[fdim-1] = fdim-1;
			}

			
			

		}





	}




	@Test
	public void demoIntegration() {

		//		class ExpClass {
		//			public double tau;
		//			public ExpClass(double tau) {
		//				this.tau = tau;
		//			}
		//			@SuppressWarnings("unused")
		//			public double[][] tauExp(double[][] x) {
		//				double[] ret = new double[x[0].length];
		//				for (int i=0; i<x[0].length; ++i) {
		//					ret[i] = Math.exp(x[0][i]/tau);
		//				}
		//				return new double[][] { ret };
		//			}
		//		}
		//
		//		ExpClass c1 = new ExpClass(1.0);

		class CosClass {
			int dim;
			public CosClass(int mydim) {
				dim = mydim;
			}
			@SuppressWarnings("unused")
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
		CosClass cc = new CosClass(2);

		class F0_class {
			int dim;
			public F0_class(int mydim) {
				dim = mydim;
			}

			@SuppressWarnings("unused")
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
		F0_class f0 = new F0_class(2);

		class F1_class {
			double a = 0.1;
			@SuppressWarnings("unused")
			public double[][] eval(double[][] x, Object fdata) {
				double[] ret = new double[x[0].length];
				for (int i=0; i<x[0].length; ++i) {
					double dx = x[0][i] - 0.5;
					double sum = dx * dx;
					ret[i] = 1.12837916709551257390 / (2. * a) * Math.exp(-sum / (a * a));
				}
				return new double[][] { ret };
			}
		}
		F1_class f1 = new F1_class();



		class KlarasBadFunction {
			@SuppressWarnings("unused")
			public double[][] eval(double[][] x, Object fdata) {
				double[] ret = new double[x[0].length];
				for (int i=0; i<x[0].length; ++i) {
					//double p1 = Math.cos(Math.pow(x[0][i], x[1][i]));
					double p1 = Math.cos(x[0][i] + x[1][i]);
					double p2 = Math.sin(x[0][i]*x[0][i]);
					double p3 = Math.exp(- x[0][i]*x[0][i] - x[1][i]);
					ret[i] = p1*p1 * p2 * p3;
				}
				return new double[][] { ret };
			}
		}
		KlarasBadFunction kbf = new KlarasBadFunction();

		double tol = 1.0e-9;

		double[][] F = Cubature.integrate(cc, "eval",
				//double[][] F = integrate(f0, "eval",
				//double[][] F = integrate(f1, "eval",
				//double[][] F = integrate(kbf, "eval",
				//double[][] F = integrate(MDAI.class, "kbf_static",
				new double[] { 0.0, 0.0 }, // lower integration limit
				new double[] { 1.0, 1.0 }, // upper integration limit
				//				new double[] {   4.0, 4.1 }, // lower integration limit
				//				new double[] { -10.0, Math.PI }, // upper integration limit

				tol, // relative tolerance
				0.0, // no absolute tolerance requirement
				//Double.POSITIVE_INFINITY,
				Cubature.Error.INDIVIDUAL,
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
