package de.labathome;

import java.lang.reflect.Method;

/**
 * Based on rule75genzmalik.cpp in HIntLib-0.0.10:
 * An embedded cubature rule of degree 7 (embedded rule degree 5)
 * due to A. C. Genz and A. A. Malik.
 *
 * See:
 * A. C. Genz and A. A. Malik,
 * "An imbedded [sic] family of fully symmetric numerical integration rules",
 * SIAM J. Numer. Anal. 20 (3), 580-588 (1983).
*/
public class Rule75GenzMalik extends Rule {

	public static int num0_0(    int dim) { return 1; }
	public static int numR0_0fs( int dim) { return 2 * dim; }
	public static int numRR0_0fs(int dim) { return 2 * dim * (dim-1); }
	public static int numR_Rfs(  int dim) { return (1 << dim); }

	/* Based on orbitrule.cpp in HIntLib-0.0.10 */
	/*
	 * ls0 returns the least-significant 0 bit of n (e.g. it returns 0 if the LSB is
	 * 0, it returns 1 if the 2 LSBs are 01, etc.).
	 */
	public static int ls0(int n) {
		final byte[] bits = {
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,

				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7,

				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,

				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
				0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8 };
		int bit = 0;
		while ((n & 0xff) == 0xff) {
			n >>= 8;
			bit += 8;
		}
		return bit + bits[n & 0xff];
	}

	/**
	 * Evaluate the integration points for all 2^n points (+/-r,...+/-r)
	 *
	 * A Gray-code ordering is used to minimize the number of coordinate updates in
	 * p, although this doesn't matter as much now that we are saving all pts.
	 *
	 * @param pts_offset index offset in pts
	 * @param dim number of integration dimensions
	 * @param p re-usable array for points source
	 * @param c center
	 * @param r radius
	 */
	public void evalR_Rfs(int pts_offset, int dim, double[] p, final double[] c, final double[] r) {
		int i;

		/** 0/1 bit = +/- for corresponding element of r[] */
		int signs = 0;

		/**
		 * We start with the point where r is ADDed in every coordinate (this implies
		 * signs=0).
		 */
		for (i = 0; i < dim; ++i) {
			p[i] = c[i] + r[i];
		}

		/** Loop through the points in Gray-code ordering */
		for (i = 0;; ++i) {
			int mask, d;

			for (int d2=0; d2<dim; ++d2) { pts[d2][pts_offset] = p[d2]; }
			pts_offset++;

			/** which coordinate to flip */
			d = ls0(i);
			if (d >= dim) {
				break;
			}

			/** flip the d-th bit and add/subtract r[d] */
			mask = 1 << d;
			signs ^= mask;
			p[d] = (signs & mask) != 0 ? c[d] - r[d] : c[d] + r[d];
		}
	}

	public void evalRR0_0fs(int pts_offset, int dim, double[] p, final double[] c, final double[] r) {
		int i, j;

		for (i = 0; i < dim - 1; ++i) {
			p[i] = c[i] - r[i];
			for (j = i + 1; j < dim; ++j) {
				p[j] = c[j] - r[j];
				for (int d=0; d<dim; ++d) { pts[d][pts_offset] = p[d]; }
				pts_offset++;

				p[i] = c[i] + r[i];
				for (int d=0; d<dim; ++d) { pts[d][pts_offset] = p[d]; }
				pts_offset++;

				p[j] = c[j] + r[j];
				for (int d=0; d<dim; ++d) { pts[d][pts_offset] = p[d]; }
				pts_offset++;

				p[i] = c[i] - r[i];
				for (int d=0; d<dim; ++d) { pts[d][pts_offset] = p[d]; }
				pts_offset++;

				p[j] = c[j]; /* Done with j -> Restore p[j] */
			}
			p[i] = c[i]; /* Done with i -> Restore p[i] */
		}
	}

	public void evalR0_0fs4d(int pts_offset, int dim, double[] p, final double[] c, final double[] r1,
			final double[] r2) {
		int i;

		// center
		for (int d=0; d<dim; ++d) { pts[d][pts_offset] = p[d]; }
		pts_offset++;

		for (i = 0; i < dim; i++) {
			p[i] = c[i] - r1[i];
			for (int d=0; d<dim; ++d) { pts[d][pts_offset] = p[d]; }
			pts_offset++;

			p[i] = c[i] + r1[i];
			for (int d=0; d<dim; ++d) { pts[d][pts_offset] = p[d]; }
			pts_offset++;

			p[i] = c[i] - r2[i];
			for (int d=0; d<dim; ++d) { pts[d][pts_offset] = p[d]; }
			pts_offset++;

			p[i] = c[i] + r2[i];
			for (int d=0; d<dim; ++d) { pts[d][pts_offset] = p[d]; }
			pts_offset++;

			p[i] = c[i];
		}
	}

	/* temporary arrays of length dim */
	double[] widthLambda, widthLambda2, p;

	/* constants */
	final double weight2, weight4;
	final double weightE2, weightE4;

	/* dimension-dependent constants */
	double weight1, weight3, weight5;
	double weightE1, weightE3;

	public Rule75GenzMalik(int dim, int fdim, int num_points) {
		super(dim, fdim, num_points);
		if (dim<2) {
			/** this rule does not support 1d integrals */
			throw new RuntimeException("75GenzMalik only support integrals of dim>=2");
		}

		if (dim >= Integer.SIZE) {
			/** Because of the use of a bit-field in evalR_Rfs, we are limited
			 to be < 32 dimensions (or however many bits are in unsigned).
			 This is not a practical limitation...long before you reach
			 32 dimensions, the Genz-Malik cubature becomes excruciatingly
			 slow and is superseded by other methods (e.g. Monte-Carlo). */
			throw new RuntimeException("75GenzMalik only support integrals of dim<"+Integer.SIZE);
		}

		weight1  = (12824 - 9120*dim + 400*dim*dim) / 19683.0;
		weight3  = (1820 - 400*dim) / 19683.0;
		weight5  = 6859.0 / (19683.0*(1 << dim));
		weightE1 = (729 - 950*dim + 50*dim*dim)	/ 729.0;
		weightE3 = (265 - 100*dim) / 1458.0;

		weight2  = 980.0 / 6561.0;
		weight4  = 200.0 / 19683.0;
		weightE2 = 245.0 / 486.0;
		weightE4 =  25.0 / 729.0;

		p = new double[dim];
		widthLambda = new double[dim];
		widthLambda2 = new double[dim];
	}

	@Override
	public void evalError(Object o, Method m, Region[] R, int nR, Object fdata) {

		alloc_rule_pts(nR);

		/* lambda2 = sqrt(9/70), lambda4 = sqrt(9/10), lambda5 = sqrt(9/19) */
		final double lambda2 = 0.3585685828003180919906451539079374954541;
		final double lambda4 = 0.9486832980505137995996680633298155601160;
		final double lambda5 = 0.6882472016116852977216287342936235251269;

		final double ratio = (lambda2 * lambda2) / (lambda4 * lambda4);

		int npts = 0, iR, i, j;

		for (iR = 0; iR < nR; ++iR) {
			double[] center = R[iR].h.centers;
			double[] halfwidth = R[iR].h.halfwidths;

			for (i = 0; i < dim; ++i) {
				p[i] = center[i];
			}

			for (i = 0; i < dim; ++i) {
				widthLambda2[i] = halfwidth[i] * lambda2;
			}
			for (i = 0; i < dim; ++i) {
				widthLambda[i] = halfwidth[i] * lambda4;
			}

			/* Evaluate points in the center, in (lambda2,0,...,0) and
			 (lambda3=lambda4, 0,...,0).  */
			evalR0_0fs4d(npts, dim, p, center, widthLambda2, widthLambda);
			npts += num0_0(dim) + 2 * numR0_0fs(dim);

			/* Calculate points for (lambda4, lambda4, 0, ...,0) */
			evalRR0_0fs(npts, dim, p, center, widthLambda);
			npts += numRR0_0fs(dim);

			/* Calculate points for (lambda5, lambda5, ..., lambda5) */
			for (i = 0; i < dim; ++i) {
				widthLambda[i] = halfwidth[i] * lambda5;
			}
			evalR_Rfs(npts, dim, p, center, widthLambda);
			npts += numR_Rfs(dim);
		}

		try {
			// evaluate function
			vals = (double[][]) m.invoke(o, (Object)(pts), fdata); // [fdim][15]
		} catch (Exception e) {
			throw new RuntimeException(e);
		}

		/* we are done with the points, and so we can re-use the pts
		 array to store the maximum difference diff[i] in each dimension
		 for each hypercube */
		double[][] diff = pts;
		for (i = 0; i < dim /* * nR*/; ++i)
			for (j=0; j<nR; ++j)
				diff[i][j] = 0.0;

		for (j = 0; j < fdim; ++j) {
			double[] v = new double[num_points];

			for (iR = 0; iR < nR; ++iR) {
				System.arraycopy(vals[j], iR*num_points, v, 0, num_points);

				double result, res5th;
				double val0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0;
				int k, k0 = 0;
				/* accumulate j-th function values into j-th integrals
				 NOTE: this relies on the ordering of the eval functions
				 above, as well as on the internal structure of
				 the evalR0_0fs4d function */

				val0 = v[0]; /* central point */
				k0 += 1;

				for (k = 0; k < dim; ++k) {
					double v0 = v[(k0 + 4 * k) + 0];
					double v1 = v[(k0 + 4 * k) + 1];
					double v2 = v[(k0 + 4 * k) + 2];
					double v3 = v[(k0 + 4 * k) + 3];

					sum2 += v0 + v1;
					sum3 += v2 + v3;

					diff[k][iR] += Math.abs(
							           v0 + v1 - 2 * val0
							- ratio * (v2 + v3 - 2 * val0));
				}
				k0 += 4 * k;

				for (k = 0; k < numRR0_0fs(dim); ++k) {
					sum4 += v[k0 + k];
				}
				k0 += k;

				for (k = 0; k < numR_Rfs(dim); ++k) {
					sum5 += v[k0 + k];
				}

				/* Calculate fifth and seventh order results */
				result = R[iR].h.volume * (weight1  * val0 + weight2  * sum2 + weight3  * sum3 + weight4  * sum4 + weight5 * sum5);
				res5th = R[iR].h.volume	* (weightE1 * val0 + weightE2 * sum2 + weightE3 * sum3 + weightE4 * sum4);

				R[iR].ee[j].val = result;
				R[iR].ee[j].err = Math.abs(res5th - result);
			}
		}

		/* figure out dimension to split: */
		for (iR = 0; iR < nR; ++iR) {
			double maxdiff = 0;
			int dimDiffMax = 0;

			for (i = 0; i < dim; ++i) {
				if (diff[i][iR] > maxdiff) {
					maxdiff = diff[i][iR];
					dimDiffMax = i;
				}
			}
			R[iR].splitDim = dimDiffMax;
		}
	}
}