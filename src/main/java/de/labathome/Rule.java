package de.labathome;

import java.lang.reflect.Method;

/** adaptive integration, analogous to adaptintegrator.cpp in HIntLib */
public abstract class Rule {
	int dim, fdim; /* the dimensionality & number of functions */
	int num_points; /* number of evaluation points */
	int num_regions; /* max number of regions evaluated at once */


	/** Set this to true in order to get some debugging output messages */
	public static boolean _debugMessages = false;

	/** points to eval: [dim][num_regions*num_points] */
	double[][] pts;

	/** [fdim][num_regions*num_points] */
	double[][] vals;

	public Rule(int dim, int fdim, int num_points) {
		this.dim = dim;
		this.fdim = fdim;
		this.num_points = num_points;
		this.num_regions = 0;
	}

	public abstract void evalError(Object o, Method m, Region[] R, int nR, Object fdata);

	public void alloc_rule_pts(int _num_regions) {
		if (_num_regions > num_regions) {
			/* allocate extra so that
				 repeatedly calling alloc_rule_pts with
				 growing num_regions only needs
				 a logarithmic number of allocations */
			num_regions = _num_regions*2;
			pts = new double[dim][num_regions*num_points];
			// allocation has to be be done by integrand
			//vals = new double[fdim][num_regions*num_points];
		}
	}

	public void cubature(Object o, Method m, int maxEval, double relTol, double absTol, double[] val, double[] err, Hypercube h, CubatureError norm, Object fdata) {

		int numEval = 0;
		RegionHeap regions = new RegionHeap(fdim);

		int i, j;
		Region[] R = null; /* array of regions to evaluate */
		int nR_alloc = 0;

		nR_alloc = 2;
		R = new Region[nR_alloc];
		R[0] = new Region().init(h, fdim);

		evalError(o, m, new Region[] { R[0] } , 1, fdata);
		R[0].errmax = errMax(R[0].fdim, R[0].ee);

		regions.add(R[0]);

		numEval += num_points;

		boolean converged = false;
		while (numEval < maxEval || maxEval==0) {
			if (converged(fdim, regions.ee, absTol, relTol, norm)) {
				if (_debugMessages) System.out.println("converged after "+numEval+" function evaluations");
				converged = true;
				break;
			}

			//				boolean parallel = true;
			//
			//				/** maximize potential parallelism */
			//				if (parallel) {
			/** adapted from I. Gladwell, "Vectorization of one
								 dimensional quadrature codes," pp. 230--238 in
								 _Numerical Integration. Recent Developments,
								 Software and Applications_, G. Fairweather and
								 P. M. Keast, eds., NATO ASI Series C203, Dordrecht
								 (1987), as described in J. M. Bull and
								 T. L. Freeman, "Parallel Globally Adaptive
								 Algorithms for Multi-dimensional Integration,"
								 http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.6638
								 (1994).

								 Basically, this evaluates in one shot all regions
								 that *must* be evaluated in order to reduce the
								 error to the requested bound: the minimum set of
								 largest-error regions whose errors push the total
								 error over the bound.

								 [Note: Bull and Freeman claim that the Gladwell
								 approach is intrinsically inefficent because it
								 "requires sorting", and propose an alternative
								 algorithm that "only" requires three passes over the
								 entire set of regions.  Apparently, they didn't
								 realize that one could use a heap data structure, in
								 which case the time to pop K biggest-error regions
								 out of N is only O(K log N), much better than the
								 O(N) cost of the Bull and Freeman algorithm if K <<
								 N, and it is also much simpler.] */
			int nR = 0;
			EstErr[] ee = new EstErr[fdim];
			for (j = 0; j < fdim; ++j) {
				ee[j] = new EstErr();
				ee[j].err = regions.ee[j].err;
				ee[j].val = regions.ee[j].val;
			}

			do {
				if (nR + 2 > nR_alloc) {
					nR_alloc = (nR + 2) * 2;

					Region[] R_old = R.clone();
					R = new Region[nR_alloc];
					for (i=0; i<R_old.length; ++i) { R[i] = R_old[i]; }
				}
				R[nR] = regions.poll();
				for (j = 0; j < fdim; ++j) {
					ee[j].err -= R[nR].ee[j].err;
				}
				R[nR+1] = R[nR].cut();

				numEval += num_points * 2;
				nR += 2;
				if (converged(fdim, ee, absTol, relTol, norm)) {
					break; /* other regions have small errs */
				}
			} while (regions.size() > 0 && (numEval < maxEval || maxEval==0));
			evalError(o,m,R, nR, fdata);
			for (i=0; i<nR; ++i) { R[i].errmax = errMax(R[i].fdim, R[i].ee); regions.add(R[i]); }
			//				} else { /** minimize number of function evaluations per call to 2 */
			//					/** get worst region */
			//					R[0] = regions.poll();
			//
			//					R[1] = R[0].cut();
			//
			//					evalError(o, m, R, 2);
			//					R[0].errmax = errMax(R[0].fdim, R[0].ee);
			//					R[1].errmax = errMax(R[1].fdim, R[1].ee);
			//
			//					regions.add(R[0]);
			//					regions.add(R[1]);
			//
			//					numEval += num_points * 2;
			//				}
		}

		if (!converged) {
			System.out.println("Cubature did not converge after "+numEval+" function evaluations!");
		}

		/** re-sum integral and errors */
		for (j = 0; j < fdim; ++j)
			val[j] = err[j] = 0;
		Region[] _regions = regions.toArray(new Region[regions.size()]);
		for (i = 0; i < _regions.length; ++i) {
			for (j = 0; j < fdim; ++j) {
				val[j] += _regions[i].ee[j].val;
				err[j] += _regions[i].ee[j].err;
			}
		}
	}

	private static boolean converged(int fdim, EstErr[] ee, double absTol,
			double relTol, CubatureError norm) {
		if (Double.isNaN(relTol) && Double.isNaN(absTol)) {
			throw new RuntimeException("Either relTol or absTol or both have to be not NaN in order to define a valid convergence criterion");
		}

		boolean converged = false;
		int j;
		switch (norm) {
		case INDIVIDUAL: {
			for (j = 0; j < fdim; ++j) {
				if (!Double.isNaN(absTol)) {
					converged |= ee[j].err < absTol;
				}
				if (!Double.isNaN(relTol)) {
					converged |= ee[j].err < Math.abs(ee[j].val)*relTol;
				}
			}
		}
		break;
		case PAIRED: {
			for (j = 0; j+1 < fdim; j += 2) {
				double maxerr, serr, err, maxval, sval, val;
				/* scale to avoid overflow/underflow */
				maxerr = ee[j].err > ee[j+1].err ? ee[j].err : ee[j+1].err;
				maxval = ee[j].val > ee[j+1].val ? ee[j].val : ee[j+1].val;
				serr = maxerr > 0 ? 1/maxerr : 1;
				sval = maxval > 0 ? 1/maxval : 1;
				err = Math.sqrt((ee[j].err*serr)*(ee[j].err*serr) + (ee[j+1].err*serr)*(ee[j+1].err*serr)) * maxerr;
				val = Math.sqrt((ee[j].val*sval)*(ee[j].val*sval) + (ee[j+1].val*sval)*(ee[j+1].val*sval)) * maxval;
				if (!Double.isNaN(absTol)) {
					converged |= err < absTol;
				}
				if (!Double.isNaN(relTol)) {
					converged |= err < Math.abs(val)*relTol;
				}
			}
			if (j < fdim) {
				/* fdim is odd, do last dimension individually */
				if (!Double.isNaN(absTol)) {
					converged |= ee[j].err < absTol;
				}
				if (!Double.isNaN(relTol)) {
					converged |= ee[j].err < Math.abs(ee[j].val)*relTol;
				}
			}
		}
		break;
		case L1: {
			double err = 0, val = 0;
			for (j = 0; j < fdim; ++j) {
				err += ee[j].err;
				val += Math.abs(ee[j].val);
			}
			if (!Double.isNaN(absTol)) {
				converged |= err < absTol;
			}
			if (!Double.isNaN(relTol)) {
				converged |= err < val*relTol;
			}
		}
		break;
		case LINF: {
			double err = 0, val = 0;
			for (j = 0; j < fdim; ++j) {
				double absval = Math.abs(ee[j].val);
				if (ee[j].err > err) err = ee[j].err;
				if (absval > val) val = absval;
			}
			if (!Double.isNaN(absTol)) {
				converged |= err < absTol;
			}
			if (!Double.isNaN(relTol)) {
				converged |= err < val*relTol;
			}
		}
		break;
		case L2: {
			double maxerr = 0, maxval = 0, serr, sval, err = 0, val = 0;
			/* scale values by 1/max to avoid overflow/underflow */
			for (j = 0; j < fdim; ++j) {
				double absval = Math.abs(ee[j].val);
				if (ee[j].err > maxerr) maxerr = ee[j].err;
				if (absval > maxval) maxval = absval;
			}
			serr = maxerr > 0 ? 1/maxerr : 1;
			sval = maxval > 0 ? 1/maxval : 1;
			for (j = 0; j < fdim; ++j) {
				err += (ee[j].err * serr)*(ee[j].err * serr);
				val += (Math.abs(ee[j].val) * sval)*(Math.abs(ee[j].val) * sval);
			}
			err = Math.sqrt(err) * maxerr;
			val = Math.sqrt(val) * maxval;
			if (!Double.isNaN(absTol)) {
				converged |= err < absTol;
			}
			if (!Double.isNaN(relTol)) {
				converged |= err < val*relTol;
			}
		}
		break;
		default: {
			converged = false;
		}
		break;
		}
		return converged;
	}


	protected static double errMax(int fdim, EstErr[] ee) {
		double errmax = 0;
		int k;
		for (k = 0; k < fdim; ++k)
			if (ee[k].err > errmax)
				errmax = ee[k].err;
		return errmax;
	}
}