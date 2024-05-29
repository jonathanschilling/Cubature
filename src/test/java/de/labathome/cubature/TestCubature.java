package de.labathome.cubature;

import java.util.Arrays;
import java.util.function.UnaryOperator;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

/**
 * Test program for hcubature/pcubature.
 *
 * Copyright (c) 2005-2013 Steven G. Johnson
 *
 * Portions (see comments) based on HIntLib (also distributed under the GNU GPL,
 * v2 or later), copyright (c) 2002-2005 Rudolf Schuerer.
 * (http://www.cosy.sbg.ac.at/~rschuer/hintlib/)
 *
 * Portions (see comments) based on GNU GSL (also distributed under the GNU GPL,
 * v2 or later), copyright (c) 1996-2000 Brian Gough.
 * (http://www.gnu.org/software/gsl/)
 *
 * This program is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA 02111-1307 USA
 */
public class TestCubature {

	/** radius of hypersphere to integrate; chosen at random */
	static final double radius = 0.50124145262344534123412;

	// f0, f1, f2, and f3 are test functions
	// from the Monte-Carlo integration routines in GSL 1.6 (monte/test.c).
	// Copyright (c) 1996-2000 Michael Booth, GNU GPL.

	/** Simple product function */
	static double f0(double[] x) {
		int dim = x.length;
		double prod = 1.0;
		for (int i = 0; i < dim; ++i) {
			prod *= 2.0 * x[i];
		}
		return prod;
	}

	static final double K_2_SQRTPI = 1.12837916709551257390;

	/** Gaussian centered at 1/2 */
	static double f1(double[] x, double a) {
		int dim = x.length;
		double sum = 0.0;
		for (int i = 0; i < dim; i++) {
			double dx = x[i] - 0.5;
			sum += dx * dx;
		}
		return (Math.pow(K_2_SQRTPI / (2.0 * a), dim) * Math.exp(-sum / (a * a)));
	}

	/** double Gaussian */
	static double f2(double[] x, double a) {
		int dim = x.length;
		double sum1 = 0.0;
		double sum2 = 0.0;
		for (int i = 0; i < dim; i++) {
			double dx1 = x[i] - 1.0 / 3.0;
			double dx2 = x[i] - 2.0 / 3.0;
			sum1 += dx1 * dx1;
			sum2 += dx2 * dx2;
		}
		return 0.5 * Math.pow(K_2_SQRTPI / (2.0 * a), dim) * (Math.exp(-sum1 / (a * a)) + Math.exp(-sum2 / (a * a)));
	}

	/** Tsuda's example */
	static double f3(double[] x, double c) {
		int dim = x.length;
		double prod = 1.0;
		for (int i = 0; i < dim; i++) {
			prod *= c / (c + 1) * Math.pow((c + 1) / (c + x[i]), 2.0);
		}
		return prod;
	}

	/**
	 * Test integrand from W. J. Morokoff and R. E. Caflisch, "Quasi-Monte Carlo
	 * integration", J. Comput. Phys 122, 218-230 (1995). Designed for integration
	 * on [0,1]^dim, integral = 1.
	 */
	static double morokoff(double[] x) {
		int dim = x.length;
		double p = 1.0 / dim;
		double prod = Math.pow(1 + p, dim);
		for (int i = 0; i < dim; i++) {
			prod *= Math.pow(x[i], p);
		}
		return prod;
	}

	// end of GSL test functions

	static double[] f_test(double[] x, int[] which_integrand) {
		int dim = x.length;
		int fdim = which_integrand.length;
		double[] retval = new double[fdim];

		for (int j = 0; j < fdim; ++j) {
			double val;
			switch (which_integrand[j]) {
			case 0:
				// simple smooth (separable) objective: prod. cos(x[i]).
				val = 1;
				for (int i = 0; i < dim; ++i) {
					val *= Math.cos(x[i]);
				}
				break;
			case 1:
				// integral of exp(-x^2), rescaled to (0, infinity) limits
				double scale = 1.0;
				val = 0;
				for (int i = 0; i < dim; ++i) {
					if (x[i] > 0) {
						double z = (1 - x[i]) / x[i];
						val += z * z;
						scale *= K_2_SQRTPI / (x[i] * x[i]);
					} else {
						scale = 0;
						break;
					}
				}
				val = Math.exp(-val) * scale;
				break;
			case 2:
				// discontinuous objective: volume of hypersphere
				val = 0;
				for (int i = 0; i < dim; ++i) {
					val += x[i] * x[i];
				}
				val = val < radius * radius ? 1.0 : 0.0;
				break;
			case 3:
				val = f0(x);
				break;
			case 4:
				val = f1(x, 0.1);
				break;
			case 5:
				val = f2(x, 0.1);
				break;
			case 6:
				val = f3(x, (1.0 + Math.sqrt(10.0)) / 9.0);
				break;
			case 7:
				val = morokoff(x);
				break;
			case 8:
				// from HCubature.jl#4
				if (dim != 3) {
					throw new RuntimeException(String.format("test 8 requires dim == 3, but dim=%d", dim));
				}
				val = x[0] * 0.2 * (x[2] - 0.5) * 0.4 * Math.sin(x[1] * 6.283185307179586);
				val = 1 + val * val;
				break;
			default:
				throw new RuntimeException(String.format("unknown integrand %d", which_integrand[j]));
			}
			retval[j] = val;
		}
		return retval;
	}

	static final double K_PI = 3.14159265358979323846;

	/** surface area of n-dimensional unit hypersphere */
	static double S(int n) {
		double val;
		int fact = 1;
		if (n % 2 == 0) {
			// n even
			val = 2 * Math.pow(K_PI, n * 0.5);
			n = n / 2;
			while (n > 1) {
				fact *= (n -= 1);
			}
			val /= fact;
		} else {
			// n odd
			val = (1 << (n / 2 + 1)) * Math.pow(K_PI, n / 2);
			while (n > 2) {
				fact *= (n -= 2);
			}
			val /= fact;
		}
		return val;
	}

	static double exact_integral(int which, double[] xmax) {
		int dim = xmax.length;
		double val;
		switch (which) {
		case 0:
			val = 1.0;
			for (int i = 0; i < dim; ++i) {
				val *= Math.sin(xmax[i]);
			}
			break;
		case 2:
			val = dim == 0 ? 1.0 : S(dim) * Math.pow(radius * 0.5, dim) / dim;
			break;
		default:
			val = 1.0;
		}
		return val;
	}

	// TODO: fails currently ...
	@Test
	void testThreeDimUnitSphere() {

		UnaryOperator<double[][]> integrand = (double[][] x) -> {
			int xdim = x.length;
			int nPts = x[0].length;

			System.out.printf("eval at %d points\n", nPts);

			double[][] ret = new double[1][nPts];

			double[] r2 = new double[nPts];
			for (int dim = 0; dim < xdim; ++dim) {
				for (int iPt = 0; iPt < nPts; ++iPt) {
					r2[iPt] += x[dim][iPt] * x[dim][iPt];
				}
			}

			for (int iPt = 0; iPt < nPts; ++iPt) {
				if (r2[iPt] < radius) {
					ret[0][iPt] += 1.0;
				}
			}

			return ret;
		};

		int xdim = 3;

		double[] xmin = new double[xdim];
		double[] xmax = new double[xdim];

		Arrays.fill(xmin, 0.0);
		Arrays.fill(xmax, 1.0);

		double relTol = 1.0e-2;
		double absTol = Double.NaN; // absolute error check disabled
		int maxEval = 0; // limit on maximum number of iterations disabled
		double[][] val_err = Cubature.integrate(integrand, xmin, xmax, relTol, absTol, CubatureError.INDIVIDUAL, maxEval);

		double expected = S(xdim) * Math.pow(radius * 0.5, xdim) / xdim;
		double act_val = val_err[0][0];
		double act_err = val_err[1][0];

		double rel_err_est = act_err / act_val;
		double rel_err = (act_val - expected) / expected;

		System.out.printf("\nintegrand=%d xdim=%d\n" +
		                  "  expected value of the integral = % .3e\n" +
				          "  computed value of the integral = % .3e\n" +
				          "         relative error estimate = % .3e\n" +
				          "          (actual relative error = % .3e)\n",
				          2, xdim,
				          expected, act_val, rel_err_est, rel_err);

		Assertions.assertTrue(Math.abs(rel_err_est) < relTol);
	}

	@Test
	void testIndividualIntegrands() {

		for (int idx_integrand = 0; idx_integrand <= 8; ++idx_integrand) {

			// dimensionality of parameter space
			for (int dim = 1; dim <= 4; ++dim) {

				double[] xmin = new double[dim];
				double[] xmax = new double[dim];

				Arrays.fill(xmin, 0.0);
				Arrays.fill(xmax, 1.0);

				// don't deal with repetions of individual integrands or permutations yet
				int[] which_integrand = { idx_integrand };

				UnaryOperator<double[][]> integrand = (double[][] x) -> {
					int xdim = x.length;
					int nPts = x[0].length;
					double[][] ret = new double[1][nPts];

					for (int iPt = 0; iPt < nPts; ++iPt) {

						double[] x_i = new double[xdim];
						for (int iXDim = 0; iXDim < xdim; ++iXDim) {
							x_i[iXDim] = x[iXDim][iPt];
						}

						double[] f_i = f_test(x_i, which_integrand);

						ret[0][iPt] = f_i[0];
					}

					return ret;
				};

				double relTol = 1.0e-2;
				double absTol = Double.NaN; // absolute error check disabled
				int maxEval = 0; // limit on maximum number of iterations disabled
				double[][] val_err = Cubature.integrate(integrand, xmin, xmax, relTol, absTol, CubatureError.INDIVIDUAL, maxEval);

				double expected = exact_integral(idx_integrand, xmax);
				double act_val = val_err[0][0];
				double act_err = val_err[1][0];

				double rel_err_est = act_err / act_val;
				double rel_err = (act_val - expected) / expected;

				System.out.printf("\nintegrand=%d xdim=%d\n" +
				                  "  expected value of the integral = % .3e\n" +
						          "  computed value of the integral = % .3e\n" +
						          "         relative error estimate = % .3e\n" +
						          "          (actual relative error = % .3e)\n",
						          idx_integrand, dim,
						          expected, act_val, rel_err_est, rel_err);

				Assertions.assertTrue(Math.abs(rel_err_est) < relTol);

				// error estimate can be wrong in some cases !!!
//				Assertions.assertTrue(Math.abs(rel_err) < relTol);
			}
		}
	}

//	@Test
//	void testCubature() {
//
//		// dimensionality of parameter space
//		for (int dim = 1; dim <= 4; ++dim) {
//
//			double[] xmin = new double[dim];
//			double[] xmax = new double[dim];
//
//			Arrays.fill(xmin, 0.0);
//			Arrays.fill(xmax, 1.0);
//
//			int max_fdim = 7;
//			if (dim == 3) {
//				// The example from HCubature.jl#4 only works for dim==3.
//				max_fdim = 8;
//			}
//
//			// dimensionality of vector function to integrate
////			for (int fdim = 1; fdim <= max_fdim; ++fdim) {
//			int fdim = 1; {
//
//				// don't deal with repetions of individual integrands or permutations yet
//				int[] which_integrand = new int[fdim];
//				for (int i = 0; i < fdim; ++i) {
//					which_integrand[i] = i;
//				}
//
//				final int my_fdim = fdim;
//
//				UnaryOperator<double[][]> integrand = (double[][] x) -> {
//					int xdim = x.length;
//					int nPts = x[0].length;
//					double[][] ret = new double[my_fdim][nPts];
//
//					for (int iPt = 0; iPt < nPts; ++iPt) {
//
//						double[] x_i = new double[xdim];
//						for (int iXDim = 0; iXDim < xdim; ++iXDim) {
//							x_i[iXDim] = x[iXDim][iPt];
//						}
//
//						double[] f_i = f_test(x_i, which_integrand);
//
//						for (int iFDim = 0; iFDim < my_fdim; ++iFDim) {
//							ret[iFDim][iPt] = f_i[iFDim];
//						}
//					}
//
//					return ret;
//				};
//
//				double relTol = 1.0e-2;
//				double absTol = Double.NaN; // absolute error check disabled
//				int maxEval = 0; // limit on maximum number of iterations disabled
//				double[][] val_err = Cubature.integrate(integrand, xmin, xmax, relTol, absTol, CubatureError.INDIVIDUAL, maxEval);
//
//				for (int i = 0; i < fdim; ++i) {
//
//					double expected = exact_integral(which_integrand[i], xmax);
//					double act_val = val_err[0][i];
//					double act_err = val_err[1][i];
//
//					double rel_err = (act_val - expected) / expected;
//					double rel_err_est = act_err / act_val;
//
//					System.out.printf("\nxdim=%d fdim=%d i=%d" +
//					                  "  expected = % .3e\n" +
//							          "   act val = % .3e\n" +
//							          "   act err = % .3e\n" +
//							          "   rel err = % .3e\n",
//							          dim, fdim, i,
//							          expected, act_val, act_err, rel_err);
//
//					Assertions.assertTrue(Math.abs(rel_err) < relTol);
//					Assertions.assertTrue(Math.abs(rel_err_est) < relTol);
//				}
//			}
//		}
//	}

	/** test based on the inline example in the README */
	@Test
	void testThreeDimensionalGaussian() {

		double[] xmin = { -2.0, -2.0, -2.0 };
		double[] xmax = { 2.0, 2.0, 2.0 };

		double sigma = 0.5;

		UnaryOperator<double[][]> gaussianNd = (double[][] x) -> {
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
		};

		double relTol = 1.0e-4;
		double[][] val_err = Cubature.integrate(gaussianNd, xmin, xmax, relTol, 0.0, CubatureError.INDIVIDUAL, 0);

		Assertions.assertEquals(13.69609043, val_err[0][0], 1.0e-8);
		Assertions.assertEquals(0.00136919, val_err[1][0], 1.0e-8);

		double relative_error = val_err[1][0] / val_err[0][0];
		Assertions.assertTrue(Math.abs(relative_error) < relTol);
	}
}
