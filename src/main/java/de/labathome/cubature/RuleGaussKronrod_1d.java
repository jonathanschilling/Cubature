package de.labathome.cubature;

import java.util.function.UnaryOperator;

/** 1d 15-point Gaussian quadrature rule, based on qk15.c and qk.c in
GNU GSL (which in turn is based on QUADPACK). */
public class RuleGaussKronrod_1d extends Rule {

	public static final double DBL_EPSILON = Math.ulp(1.0);

	/* Gauss quadrature weights and kronrod quadrature abscissae and
	 weights as evaluated with 80 decimal digit arithmetic by
	 L. W. Fullerton, Bell Labs, Nov. 1981. */

	final int n = 8;

	/* abscissae of the 15-point kronrod rule */
	/* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule.
	 xgk[0], xgk[2], ... to optimally extend the 7-point gauss rule */
	final double xgk[] = {
			0.991455371120812639206854697526329,
			0.949107912342758524526189684047851,
			0.864864423359769072789712788640926,
			0.741531185599394439863864773280788,
			0.586087235467691130294144838258730,
			0.405845151377397166906606412076961,
			0.207784955007898467600689403773245,
			0.000000000000000000000000000000000
	};

	/* weights of the 7-point gauss rule */
	final double wg[] = {
			0.129484966168869693270611432679082,
			0.279705391489276667901467771423780,
			0.381830050505118944950369775488975,
			0.417959183673469387755102040816327
	};

	/* weights of the 15-point kronrod rule */
	final double wgk[] = {
			0.022935322010529224963732008058970,
			0.063092092629978553290700663189204,
			0.104790010322250183839876322541518,
			0.140653259715525918745189590510238,
			0.169004726639267902826583426598550,
			0.190350578064785409913256402421014,
			0.204432940075298892414161999234649,
			0.209482141084727828012999174891714
	};

	public RuleGaussKronrod_1d(int dim, int fdim, int num_points) {
		super(dim, fdim, num_points);
	}

	@Override
	public void evalError(UnaryOperator<double[][]> integrand, Region[] R, int nR) {

		alloc_rule_pts(nR);

		for (int iR = 0; iR < nR; ++iR) {
			double center = R[iR].h.centers[0];
			double halfwidth = R[iR].h.halfwidths[0];

			int j;

			int npts = 0;
			pts[0][iR*15 + npts++] = center;

			// 7-point Gauss-Legendre
			for (j = 1; j < n-1; j+=2) {
				double w = halfwidth * xgk[j];
				pts[0][iR*15 + npts++] = center - w;
				pts[0][iR*15 + npts++] = center + w;
			}

			// extended to 15-point by Kronrod
			for (j = 0; j < n; j+=2) {
				double w = halfwidth * xgk[j];
				pts[0][iR*15 + npts++] = center - w;
				pts[0][iR*15 + npts++] = center + w;
			}

			R[iR].splitDim = 0; /* no choice but to divide 0th dimension */
		}

		try {
			// evaluate function
			vals = integrand.apply(pts); // [fdim][nR*15]
		} catch (Exception e) {
			throw new RuntimeException(e);
		}

		for (int k=0; k<fdim; ++k) {
			double[] vk = vals[k]; // all 15 points in the given fdim
			for (int iR = 0; iR < nR; ++iR) {
				double halfwidth = R[iR].h.halfwidths[0];

				double result_gauss = vk[iR*15+0] * wg[n / 2 - 1];
				double result_kronrod = vk[iR*15+0] * wgk[n - 1];
				double result_abs = Math.abs(result_kronrod);
				double result_asc, mean, err;

				/* accumulate integrals */
				int npts = 1;
				for (int j = 0; j < (n - 1) / 2; ++j) {
					int j2 = 2 * j + 1;
					double v = vk[iR*15+npts] + vk[iR*15+npts+1]; // +/- w are next to each other
					result_gauss   += wg[j]   * v;
					result_kronrod += wgk[j2] * v;
					result_abs     += wgk[j2] * (Math.abs(vk[iR*15+npts]) + Math.abs(vk[iR*15+npts+1]));
					npts += 2;
				}
				for (int j = 0; j < n / 2; ++j) {
					int j2 = 2 * j;
					result_kronrod += wgk[j2] * (vk[iR*15+npts] + vk[iR*15+npts+1]); // +/- w are next to each other
					result_abs     += wgk[j2] * (Math.abs(vk[iR*15+npts]) + Math.abs(vk[iR*15+npts+1]));
					npts += 2;
				}

				/* integration result */
				R[iR].ee[k].val = result_kronrod * halfwidth;

				/**
				 * error estimate (from GSL, probably dates back to QUADPACK)
				 * ... not completely clear to me why we don't just use
				 * Math.abs(result_kronrod - result_gauss) * halfwidth
				 */
				mean = result_kronrod * 0.5;
				result_asc = wgk[n - 1] * Math.abs(vk[iR*15+0] - mean);
				npts = 1;
				for (int j = 0; j < (n - 1) / 2; ++j) {
					int j2 = 2 * j + 1;
					result_asc += wgk[j2] * (Math.abs(vk[iR*15+npts] - mean) + Math.abs(vk[iR*15+npts+1] - mean));
					npts += 2;
				}
				for (int j = 0; j < n / 2; ++j) {
					int j2 = 2 * j;
					result_asc += wgk[j2] * (Math.abs(vk[iR*15+npts] - mean) + Math.abs(vk[iR*15+npts+1] - mean));
					npts += 2;
				}

				err = Math.abs(result_kronrod - result_gauss) * halfwidth;

				result_abs *= halfwidth;
				result_asc *= halfwidth;
				if (result_asc != 0 && err != 0) {
					double scale = Math.pow(200 * err / result_asc, 1.5);
					err = (scale < 1) ? result_asc * scale : result_asc;
				}
				if (result_abs > Double.MIN_NORMAL / (50 * DBL_EPSILON)) {
					double min_err = 50 * DBL_EPSILON * result_abs;
					if (min_err > err) {
						err = min_err;
					}
				}

				R[iR].ee[k].err = err;
			}
		}

		for (int iR = 0; iR < nR; ++iR) {
			R[iR].errmax = errMax(fdim, R[iR].ee);
		}
	}
}
