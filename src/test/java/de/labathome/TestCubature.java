package de.labathome;

public class TestCubature {

	

	public static double[][] kbf_static(double[][] x) {
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

	public static void demoIntegration() {

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
			public double[][] eval(double[][] x) {
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
			public double[][] eval(double[][] x) {
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
			public double[][] eval(double[][] x) {
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
			public double[][] eval(double[][] x) {
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


		double[][] F = Cubature.integrate(cc, "eval",
				//double[][] F = integrate(f0, "eval",
				//double[][] F = integrate(f1, "eval",
				//double[][] F = integrate(kbf, "eval",
				//double[][] F = integrate(MDAI.class, "kbf_static",
				new double[] { 0.0, 0.0 }, // lower integration limit
				new double[] { 1.0, 1.0 }, // upper integration limit
				//				new double[] {   4.0, 4.1 }, // lower integration limit
				//				new double[] { -10.0, Math.PI }, // upper integration limit

				1.0e-9, // relative tolerance
				0.0, // no absolute tolerance requirement
				//Double.POSITIVE_INFINITY,
				Cubature.Error.INDIVIDUAL,
				100000); // max. number of function evaluations

		if (F != null) {
			System.out.println("integration result: "+F[0][0]+" +/- " + F[1][0]);
		} else {
			System.out.println("no integration performed");
		}
	}

	public static void main(String[] args) {

		//		demoStaticEval();
		//		demoMemberEval();

		demoIntegration();
	}


}
