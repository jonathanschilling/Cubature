package de.labathome;

import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.PriorityQueue;

/* Adaptive multidimensional integration of a vector of integrands.
*
* Copyright (c) 2005-2013 Steven G. Johnson
* Ported to Java in 2019 by Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
*
* Portions (see comments) based on HIntLib (also distributed under
* the GNU GPL, v2 or later), copyright (c) 2002-2005 Rudolf Schuerer.
*     (http://www.cosy.sbg.ac.at/~rschuer/hintlib/)
*
* Portions (see comments) based on GNU GSL (also distributed under
* the GNU GPL, v2 or later), copyright (c) 1996-2000 Brian Gough.
*     (http://www.gnu.org/software/gsl/)
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*
*/

/** Adaptive multidimensional integration on hypercubes (or, really,
hyper-rectangles) using cubature rules.

A cubature rule takes a function and a hypercube and evaluates
the function at a small number of points, returning an estimate
of the integral as well as an estimate of the error, and also
a suggested dimension of the hypercube to subdivide.

Given such a rule, the adaptive integration is simple:

1) Evaluate the cubature rule on the hypercube(s).
Stop if converged.

2) Pick the hypercube with the largest estimated error,
and divide it in two along the suggested dimension.

3) Goto (1).

The basic algorithm is based on the adaptive cubature described in

A. C. Genz and A. A. Malik, "An adaptive algorithm for numeric
integration over an N-dimensional rectangular region,"
J. Comput. Appl. Math. 6 (4), 295-302 (1980).

and subsequently extended to integrating a vector of integrands in

J. Berntsen, T. O. Espelid, and A. Genz, "An adaptive algorithm
for the approximate calculation of multiple integrals,"
ACM Trans. Math. Soft. 17 (4), 437-451 (1991).

Note, however, that we do not use any of code from the above authors
(in part because their code is Fortran 77, but mostly because it is
under the restrictive ACM copyright license).  I did make use of some
GPL code from Rudolf Schuerer's HIntLib and from the GNU Scientific
Library as listed in the copyright notice above, on the other hand.

I am also grateful to Dmitry Turbiner (dturbiner@alum.mit.edu), who
implemented an initial prototype of the "vectorized" functionality
for evaluating multiple points in a single call (as opposed to
multiple functions in a single call).  (Although Dmitry implemented
a working version, I ended up re-implementing this feature from
scratch as part of a larger code-cleanup, and in order to have
a single code path for the vectorized and non-vectorized APIs.  I
subsequently implemented the algorithm by Gladwell to extract
even more parallelism by evalutating many hypercubes at once.)

TODO:

* For high-dimensional integrals, it would be nice to implement
a sparse-grid cubature scheme using Clenshaw-Curtis quadrature.
Currently, for more than 7 dimensions or so, quasi Monte Carlo methods win.

* Berntsen et. al also describe a "two-level" error estimation scheme
that they claim makes the algorithm more robust.  It might be
nice to implement this, at least as an option (although I seem
to remember trying it once and it made the number of evaluations
substantially worse for my test integrands).

*/
public class Cubature {
	
	public static final double DBL_MIN     = 2.22507385850720138309023271733240406e-308;
	public static final double DBL_EPSILON = 2.22044604925031308084726333618164062e-16;

	/**
	 * Different ways of measuring the absolute and relative error when we have
	 * multiple integrands, given a vector e of error estimates in the individual
	 * components of a vector v of integrands. These are all equivalent when there
	 * is only a single integrand.
	 */
	public enum Error {
		/** individual relerr criteria in each component */
		INDIVIDUAL,

		/**
		 * paired L2 norms of errors in each component, mainly for integrating vectors
		 * of complex numbers; this assumes that the real component is at j and the corresponding
		 * imaginary part is in j+1 for j the index in fdim
		 */
		PAIRED,

		/** abserr is L_2 norm |e|, and relerr is |e|/|v| */
		L2,

		/** abserr is L_1 norm |e|, and relerr is |e|/|v| */
		L1,

		/** abserr is L_\infty norm |e|, and relerr is |e|/|v| */
		LINF
	}
	
	/**
	 * Integrate the given function (o.method) in the range [xmin:xmax] to a relative tolerance relTol or absolute tolerance absTol or until maxEval function evaluations were used.
	 * @param o Object or Class in which {@code method} is defined
	 * @param method integrand method; must have same call signature as {@code double[][] eval(double[][] x, Object fdata)} 
	 * @param xmin vector of lower integration bounds
	 * @param xmax vector of upper integration bounds
	 * @param relTol relative tolerance on function values
	 * @param absTol absolute tolerance on function values
	 * @param norm Error norm for vector-valued integrands
	 * @param maxEval absolute tolerance on function values
	 * @param fdata any Object that shall be passed directly to the integrand; can be used to specify additional info/parameter/...
	 * @return [0:val, 1:err][fdim]
	 */
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public static double[][] integrate(Object o, String method, double[] xmin, double[] xmax, double relTol, double absTol, Error norm, int maxEval, Object fdata) {
		double[][] ret = null;

		// dimensionality of parameter space
		final int dim;
		if (xmin==null || xmin.length==0) {
			throw new RuntimeException("xmin cannot be null or empty");
		} else {
			dim = xmin.length;
		}

		if (xmax==null || xmax.length != dim) {
			throw new RuntimeException("xmin and xmax must have the same length but have "+xmin.length+" and "+xmax.length);
		}

		// dimensionality of the integrand
		final int fdim;

		// determine dimensionality of integrand by evaluating it once in the center of the given integration interval
		double[][] _center = new double[dim][1];
		for (int i=0; i<dim; ++i) {
			_center[i][0] = (xmin[i]+xmax[i])/2.0;
		}
		try {
			Method m = null;
			if (o instanceof Class) {
				//System.out.println("static member method");
				m = ((Class)o).getDeclaredMethod(method, double[][].class, Object.class);
			} else {
				//System.out.println("method of object instance");
				m = o.getClass().getDeclaredMethod(method, double[][].class, Object.class);
			}

			double[][] f = (double[][]) m.invoke(o, (Object)_center, fdata);
			if (f == null || f.length==0 || f[0]==null || f[0].length==0) {
				throw new RuntimeException("Evaluation of given method at interval center failed");
			} else {
				fdim=f.length;
			}
			
			Rule r = null;
			if (dim==1) {
				r = new RuleGaussKronrod_1d(dim, fdim, 15);
			} else {
				// 5-7 Genz-Malik for dim>1
				int numPoints = Rule75GenzMalik.num0_0(dim)
						+   2 * Rule75GenzMalik.numR0_0fs(dim)
						+       Rule75GenzMalik.numRR0_0fs(dim)
						+       Rule75GenzMalik.numR_Rfs(dim);

				//System.out.println("dim="+dim+" fdim="+fdim);
				
				r = new Rule75GenzMalik(dim, fdim, numPoints);
			}

			double[] val = new double[fdim];
			double[] err = new double[fdim];
			Arrays.fill(err, Double.POSITIVE_INFINITY);

			Hypercube h = new Hypercube().initFromRanges(xmin, xmax);

			r.cubature(o, m, maxEval, relTol, absTol, val, err, h, norm, fdata);

			ret = new double[][] { val, err };

		} catch (Exception e) {
			throw new RuntimeException(e);
		}

		return ret;
	}
	
	/**
	 * Integrand interface for Cubature.
	 */
	public static interface Integrand {
		/**
		 * Evaluate the nFDim-dimensional integrand at nPoints nDim-dimensional locations.
		 * @param x [nDim][nPoints] locations where to evaluate the integrand
		 * @param fdata some arbitrary data that is passed from the integrate() call into the integrand
		 * @return [nFDim][nPoints] function values at x
		 */
		public abstract double[][] eval(final double[][] x, Object fdata);
	}

	/**
	 * Set this to true in order to get some debugging output messages
	 */
	public static boolean _debugMessages = false;
	
	private static class Hypercube {

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

	private static class EstErr {
		double val;
		double err;
	}

	private static double errMax(int fdim, EstErr[] ee) {
		double errmax = 0;
		int k;
		for (k = 0; k < fdim; ++k)
			if (ee[k].err > errmax)
				errmax = ee[k].err;
		return errmax;
	}

	private static class Region implements Comparable<Region> {
		Hypercube h;
		int splitDim;
		int fdim; /* dimensionality of vector integrand */
		EstErr[] ee; /* array of length fdim */
		double errmax; /* max ee[k].err */

		public Region init(Hypercube h, int fdim) {
			this.h = h.clone();
			this.splitDim = 0;
			this.fdim = fdim;
			this.ee = new EstErr[fdim];
			for (int i=0; i<fdim; ++i) {
				ee[i] = new EstErr();
			}
			this.errmax = Double.POSITIVE_INFINITY;
			return this;
		}

		// return new Region and modify current one
		public Region cut() {
			h.halfwidths[splitDim] *= 0.5;
			h.volume *= 0.5;
			Region R2 = new Region().init(h.clone(), fdim);
			h.centers[splitDim] -= h.halfwidths[splitDim];
			R2.h.centers[splitDim] += h.halfwidths[splitDim];
			return R2;
		}

		@Override
		public int compareTo(Region o) {
			if (this.errmax == o.errmax) return 0;

			/**
			 * Java implements the PriorityQueue type as a MinHeap, which means that the root
			 * is always the smallest element as judged by the compareTo method of the defining
			 * type. The integration algorithm however expects a MaxHeap (where the root is the 
			 * largest element) in order to work on the worst regions first.
			 * Hence, we need to reverse the order of the Java PriorityQueue by reversing the sign
			 * of the compareTo method output. 
			 */
			return (this.errmax < o.errmax) ? 1 : -1;
		}
	}


	private static boolean converged(int fdim, EstErr[] ee, double absTol,
			double relTol, Error norm) {
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

	private static class RegionHeap extends PriorityQueue<Region> {
		private static final long serialVersionUID = 3467825258360722116L;

		int fdim;
		EstErr[] ee;

		public RegionHeap(int _fdim) {
			super();
			fdim = _fdim;
			ee = new EstErr[fdim];
			for (int i=0; i<fdim; ++i) {
				ee[i] = new EstErr();
			}
		}

		public Region poll() {
			Region ret = super.poll();
			for (int i = 0; i < fdim; ++i) {
				ee[i].val -= ret.ee[i].val;
				ee[i].err -= ret.ee[i].err;
			}
			return ret;
		}

		public boolean add(Region e) {
			for (int i = 0; i < fdim; ++i) {
				ee[i].val += e.ee[i].val;
				ee[i].err += e.ee[i].err;
			}
			return super.add(e);
		}
	}

	/** adaptive integration, analogous to adaptintegrator.cpp in HIntLib */
	private static abstract class Rule {
		int dim, fdim; /* the dimensionality & number of functions */
		int num_points; /* number of evaluation points */
		int num_regions; /* max number of regions evaluated at once */

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

		public void cubature(Object o, Method m, int maxEval, double relTol, double absTol, double[] val, double[] err, Hypercube h, Error norm, Object fdata) {

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
	}

	/** 1d 15-point Gaussian quadrature rule, based on qk15.c and qk.c in
	 GNU GSL (which in turn is based on QUADPACK). */
	private static class RuleGaussKronrod_1d extends Rule {

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
		public void evalError(Object o, Method m, Region[] R, int nR, Object fdata) {

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
				vals = (double[][]) m.invoke(o, (Object)(pts), fdata); // [fdim][15]
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

					/* error estimate
						 (from GSL, probably dates back to QUADPACK
						 ... not completely clear to me why we don't just use
						 Math.abs(result_kronrod - result_gauss) * halfwidth */
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
						double scale = Math.pow((200 * err / result_asc), 1.5);
						err = (scale < 1) ? result_asc * scale : result_asc;
					}
					if (result_abs > DBL_MIN / (50 * DBL_EPSILON)) {
						double min_err = 50 * DBL_EPSILON * result_abs;
						if (min_err > err)
							err = min_err;
					}

					R[iR].ee[k].err = err;
				}
			}

			for (int iR = 0; iR < nR; ++iR) {
				R[iR].errmax = errMax(fdim, R[iR].ee);
			}
		}
	}

	/** Based on rule75genzmalik.cpp in HIntLib-0.0.10: An embedded
	 cubature rule of degree 7 (embedded rule degree 5) due to A. C. Genz
	 and A. A. Malik.  See:

	 A. C. Genz and A. A. Malik, "An imbedded [sic] family of fully
	 symmetric numerical integration rules," SIAM
	 J. Numer. Anal. 20 (3), 580-588 (1983).
	 */
	private static class Rule75GenzMalik extends Rule {

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
			final byte[] bits = { 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1,
					0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0,
					1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2,
					0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7, 0, 1, 0, 2, 0, 1, 0,
					3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1,
					0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0,
					1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
					0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8, };
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
			for (i = 0; i < dim; ++i)
				p[i] = c[i] + r[i];

			/** Loop through the points in Gray-code ordering */
			for (i = 0;; ++i) {
				int mask, d;

				for (int d2=0; d2<dim; ++d2) { pts[d2][pts_offset] = p[d2]; }
				pts_offset++;

				/** which coordinate to flip */
				d = ls0(i); 
				if (d >= dim)
					break;

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

				for (i = 0; i < dim; ++i)
					p[i] = center[i];

				for (i = 0; i < dim; ++i)
					widthLambda2[i] = halfwidth[i] * lambda2;
				for (i = 0; i < dim; ++i)
					widthLambda[i] = halfwidth[i] * lambda4;

				/* Evaluate points in the center, in (lambda2,0,...,0) and
				 (lambda3=lambda4, 0,...,0).  */
				evalR0_0fs4d(npts, dim, p, center, widthLambda2, widthLambda);
				npts += num0_0(dim) + 2 * numR0_0fs(dim);

				/* Calculate points for (lambda4, lambda4, 0, ...,0) */
				evalRR0_0fs(npts, dim, p, center, widthLambda);
				npts += numRR0_0fs(dim);

				/* Calculate points for (lambda5, lambda5, ..., lambda5) */
				for (i = 0; i < dim; ++i)
					widthLambda[i] = halfwidth[i] * lambda5;
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
						double v0 = v[ k0 + 4 * k];
						double v1 = v[(k0 + 4 * k) + 1];
						double v2 = v[(k0 + 4 * k) + 2];
						double v3 = v[(k0 + 4 * k) + 3];

						sum2 += v0 + v1;
						sum3 += v2 + v3;

						diff[k][iR] += Math.abs(
								v0 + v1 - 2 * val0 - ratio * (v2 + v3 - 2 * val0));
					}
					k0 += 4 * k;

					for (k = 0; k < numRR0_0fs(dim); ++k)
						sum4 += v[k0 + k];
					k0 += k;

					for (k = 0; k < numR_Rfs(dim); ++k)
						sum5 += v[k0 + k];

					/* Calculate fifth and seventh order results */
					result = R[iR].h.volume
							* (weight1 * val0 + weight2 * sum2 + weight3 * sum3
									+ weight4 * sum4 + weight5 * sum5);
					res5th = R[iR].h.volume
							* (weightE1 * val0 + weightE2 * sum2 + weightE3 * sum3
									+ weightE4 * sum4);

					R[iR].ee[j].val = result;
					R[iR].ee[j].err = Math.abs(res5th - result);
				}
			}

			/* figure out dimension to split: */
			for (iR = 0; iR < nR; ++iR) {
				double maxdiff = 0;
				int dimDiffMax = 0;

				for (i = 0; i < dim; ++i)
					if (diff[i][iR] > maxdiff) {
						maxdiff = diff[i][iR];
						dimDiffMax = i;
					}
				R[iR].splitDim = dimDiffMax;
			}
		}
	}
}
