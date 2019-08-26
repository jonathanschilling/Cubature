package de.labathome;

import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.PriorityQueue;

import aliceinnets.python.jyplot.JyPlot;

public class MDAI {

	public static final double DBL_MIN     = 2.22507385850720138309023271733240406e-308;
	public static final double DBL_EPSILON = 2.22044604925031308084726333618164062e-16;

	/**
	 * 1d sine
	 * @param x [dim=1][nPoints]
	 * @return [fdim=1][nPoints]
	 */
	public static double[][] oD_oD_sin(double[][] x) {
		double[] sin_x = new double[x[0].length];
		for (int i=0; i<x[0].length; ++i) {
			sin_x[i] = Math.sin(x[0][i]);
		}
		return new double[][] { sin_x };
	}

	/**
	 * 2d: sine and cosine
	 * @param x [dim=1][nPoints]
	 * @return [fdim=2: sin, cos][nPoints]
	 */
	public static double[][] oD_tD_sinCos(double[][] x) {
		double[] sin_x = new double[x[0].length];
		double[] cos_x = new double[x[0].length];
		for (int i=0; i<x[0].length; ++i) {
			sin_x[i] = Math.sin(x[0][i]);
			cos_x[i] = Math.cos(x[0][i]);
		}
		return new double[][] { sin_x, cos_x };
	}

	/**
	 * two-dimensional Gauss
	 * @param xy [dim=2: x,y][nPoints]
	 * @return [fdim=1][nPoints]
	 */
	public static double[][] tD_oD_gauss(double[][] xy) {
		double[] gauss = new double[xy[0].length];
		for (int i=0; i<xy[0].length; ++i) {
			gauss[i] = Math.exp(-(xy[0][i]*xy[0][i]+xy[1][i]*xy[1][i])/2.0);
		}
		return new double[][] { gauss };
	}

	/**
	 * two-dimensional test function: sin(x), cos(y)
	 * @param xy [dim=2: x,y][nPoints]
	 * @return [fdim=2: sin(x), cos(y)][nPoints]
	 */
	public static double[][] tD_tD_sinCos(double[][] xy) {
		double[] sin_x = new double[xy[0].length];
		double[] cos_y = new double[xy[1].length];
		for (int i=0; i<xy[0].length; ++i) {
			sin_x[i] = Math.sin(xy[0][i]);
		}
		for (int i=0; i<xy[1].length; ++i) {
			cos_y[i] = Math.cos(xy[1][i]);
		}
		return new double[][] { sin_x, cos_y };
	}


	/**
	 * Evaluate a method of a given class in the interval [xmin:xmax] at nPoints
	 * @param o object in which the given method is defined
	 * @param method name of a method which takes double[dim][nPoints] as argument and returns double[fdim][nPoints]
	 * @param xmin[dim] left interval borders
	 * @param xmax[dim] right interval borders
	 * @param nPoints number of points at which to evaluate the given method
	 * @return [fdim][nPoints] evaluated function values
	 */
	public static double[][] eval(Object o, String method, double[] xmin, double[] xmax, int nPoints) {
		double[][] x = new double[xmin.length][nPoints];
		for (int j=0; j<xmin.length; ++j) {
			for (int i=0; i<nPoints; ++i) {
				x[j][i] = xmin[j] + (xmax[j]-xmin[j])*i/(nPoints-1.0);
			}
		}

		try {
			Method m = o.getClass().getDeclaredMethod(method, double[][].class);
			return (double[][]) m.invoke(o, (Object)x);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	/**
	 * demonstrate the eval() method: evaluate a VectorFunctionND in a given interval at a given number of points
	 */
	public static void demoStaticEval() {
		double[] xmin = new double[] { 0.0        , 0.0         };
		double[] xmax = new double[] { 2.0*Math.PI, 2.0*Math.PI };
		int nPoints = 100;

		double[][] x = new double[xmin.length][nPoints];
		for (int j=0; j<xmin.length; ++j) {
			for (int i=0; i<nPoints; ++i) {
				x[j][i] = xmin[j] + (xmax[j]-xmin[j])*i/(nPoints-1.0);
			}
		}

		double[][] sinCos_x = MDAI.eval(new MDAI(), "oD_tD_sinCos", xmin, xmax, nPoints);

		JyPlot plt = new JyPlot();

		plt.figure();
		plt.plot(x[0], sinCos_x[0], "r.-", "label='sin(x)'");
		plt.plot(x[1], sinCos_x[1], "b.-", "label='cos(x)'");
		plt.xlabel("x");
		plt.ylabel("'y=f(x)'");
		plt.legend();
		plt.grid(true);

		plt.show();
		plt.exec();
	}


	public static void demoMemberEval() {
		double[] xmin = new double[] { 0.0 };
		double[] xmax = new double[] { 1.0 };
		int nPoints = 10;

		double[][] x = new double[xmin.length][nPoints];
		for (int j=0; j<xmin.length; ++j) {
			for (int i=0; i<nPoints; ++i) {
				x[j][i] = xmin[j] + (xmax[j]-xmin[j])*i/(nPoints-1.0);
			}
		}


		class ExpClass {
			public double tau;
			public ExpClass(double tau) {
				this.tau = tau;
			}
			@SuppressWarnings("unused")
			public double[][] tauExp(double[][] x) {
				double[] ret = new double[x[0].length];
				for (int i=0; i<x[0].length; ++i) {
					ret[i] = Math.exp(x[0][i]/tau);
				}
				return new double[][] { ret };
			}
		}

		ExpClass c1 = new ExpClass(1.0);
		ExpClass c2 = new ExpClass(2.0);

		double[][] exp1_x = MDAI.eval(c1, "tauExp", xmin, xmax, nPoints);
		double[][] exp2_x = MDAI.eval(c2, "tauExp", xmin, xmax, nPoints);

		JyPlot plt = new JyPlot();

		plt.figure();
		plt.plot(x[0], exp1_x[0], "r.-", "label='exp(x)'");
		plt.plot(x[0], exp2_x[0], "b.-", "label='exp(x/2)'");
		plt.xlabel("x");
		plt.ylabel("'y=f(x)'");
		plt.legend();
		plt.grid(true);

		plt.show();
		plt.exec();
	}

	public static class Hypercube {

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

	public static class EstErr {
		double val;
		double err;
	}

	static double errMax(int fdim, EstErr[] ee) {
		double errmax = 0;
		int k;
		for (k = 0; k < fdim; ++k)
			if (ee[k].err > errmax)
				errmax = ee[k].err;
		return errmax;
	}

	public static class Region implements Comparable<Region> {
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


	public static boolean converged(int fdim, EstErr[] ee, double absTol,
			double relTol) {
		//	#define ERR(j) ee[j].err
		//	#define VAL(j) ee[j].val

		/* Body of convergence test, shared between hcubature.c and
	   pcubature.c.  We use an #include file because the two routines use
	   somewhat different data structures, and define macros ERR(j) and
	   VAL(j) to get the error and value estimates, respectively, for
	   integrand j. */

		int j;
		//	#    define SQR(x) ((x) * (x))
		//	     switch (norm) {
		//		 case ERROR_INDIVIDUAL:
		for (j = 0; j < fdim; ++j) {
			if (ee[j].err > absTol && ee[j].err > Math.abs(ee[j].val)*relTol) {
			//if (ee[j].err > absTol || ee[j].err > Math.abs(ee[j].val)*relTol) {
				return false;
			}
		}
		return true;

		//		 case ERROR_PAIRED:
		//		      for (j = 0; j+1 < fdim; j += 2) {
		//			   double maxerr, serr, err, maxval, sval, val;
		//			   /* scale to avoid overflow/underflow */
		//			   maxerr = ERR(j) > ERR(j+1) ? ERR(j) : ERR(j+1);
		//			   maxval = VAL(j) > VAL(j+1) ? VAL(j) : VAL(j+1);
		//			   serr = maxerr > 0 ? 1/maxerr : 1;
		//			   sval = maxval > 0 ? 1/maxval : 1;
		//			   err = sqrt(SQR(ERR(j)*serr) + SQR(ERR(j+1)*serr)) * maxerr;
		//			   val = sqrt(SQR(VAL(j)*sval) + SQR(VAL(j+1)*sval)) * maxval;
		//			   if (err > reqAbsError && err > val*reqRelError)
		//				return 0;
		//		      }
		//		      if (j < fdim) /* fdim is odd, do last dimension individually */
		//			   if (ERR(j) > reqAbsError && ERR(j) > fabs(VAL(j))*reqRelError)
		//				return 0;
		//		      return 1;
		//
		//		 case ERROR_L1: {
		//		      double err = 0, val = 0;
		//		      for (j = 0; j < fdim; ++j) {
		//			   err += ERR(j);
		//			   val += fabs(VAL(j));
		//		      }
		//		      return err <= reqAbsError || err <= val*reqRelError;
		//		 }
		//
		//		 case ERROR_LINF: {
		//		      double err = 0, val = 0;
		//		      for (j = 0; j < fdim; ++j) {
		//			   double absval = fabs(VAL(j));
		//			   if (ERR(j) > err) err = ERR(j);
		//			   if (absval > val) val = absval;
		//		      }
		//		      return err <= reqAbsError || err <= val*reqRelError;
		//		 }
		//
		//		 case ERROR_L2: {
		//		      double maxerr = 0, maxval = 0, serr, sval, err = 0, val = 0;
		//		      /* scale values by 1/max to avoid overflow/underflow */
		//		      for (j = 0; j < fdim; ++j) {
		//			   double absval = fabs(VAL(j));
		//			   if (ERR(j) > maxerr) maxerr = ERR(j);
		//			   if (absval > maxval) maxval = absval;
		//		      }
		//		      serr = maxerr > 0 ? 1/maxerr : 1;
		//		      sval = maxval > 0 ? 1/maxval : 1;
		//		      for (j = 0; j < fdim; ++j) {
		//			   err += SQR(ERR(j) * serr);
		//			   val += SQR(fabs(VAL(j)) * sval);
		//		      }
		//		      err = sqrt(err) * maxerr;
		//		      val = sqrt(val) * maxval;
		//		      return err <= reqAbsError || err <= val*reqRelError;
		//		 }
		//	     }
		//	     return 1; /* unreachable */


	}

	@SuppressWarnings("serial")
	public static class RegionHeap extends PriorityQueue<Region> {
		int fdim;
		EstErr[] ee;
		
		public RegionHeap(int _fdim) {
			super();
			fdim = _fdim;
			System.out.println("init RegionHeap with fdim="+fdim);
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
	
	public static abstract class Rule {
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

		public abstract void evalError(Object o, Method m, Region[] R);

		public void alloc_rule_pts(int _num_regions) {
			if (_num_regions > num_regions) {
				/* allocate extra so that
					 repeatedly calling alloc_rule_pts with
					 growing num_regions only needs
					 a logarithmic number of allocations */
				num_regions = _num_regions*2;
				pts = new double[dim][num_regions*num_points];
				//vals = new double[fdim][num_regions*num_points];
			}
		}


		public void cubature(Object o, Method m, int maxEval, double relTol, double absTol, double[] val, double[] err, Hypercube h) {

			int numEval = 0;
			RegionHeap regions = new RegionHeap(fdim);
			
			int i, j;
			Region[] R = null; /* array of regions to evaluate */
			int nR_alloc = 0;
			

			System.out.println("rulecubature");

			nR_alloc = 2;
			R = new Region[nR_alloc];
			R[0] = new Region().init(h, fdim);

			evalError(o, m, new Region[] { R[0] } );
			R[0].errmax = errMax(R[0].fdim, R[0].ee);

			regions.add(R[0]);

			numEval += num_points;

			while (numEval < maxEval || maxEval==0) {
				if (converged(fdim, regions.ee, absTol, relTol)) {
					System.out.println("converged after "+numEval+" function evaluations");
					break;
				}
					

				//				if (parallel) { /* maximize potential parallelism */
				//					/* adapted from I. Gladwell, "Vectorization of one
				//					 dimensional quadrature codes," pp. 230--238 in
				//					 _Numerical Integration. Recent Developments,
				//					 Software and Applications_, G. Fairweather and
				//					 P. M. Keast, eds., NATO ASI Series C203, Dordrecht
				//					 (1987), as described in J. M. Bull and
				//					 T. L. Freeman, "Parallel Globally Adaptive
				//					 Algorithms for Multi-dimensional Integration,"
				//					 http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.6638
				//					 (1994).
				//
				//					 Basically, this evaluates in one shot all regions
				//					 that *must* be evaluated in order to reduce the
				//					 error to the requested bound: the minimum set of
				//					 largest-error regions whose errors push the total
				//					 error over the bound.
				//
				//					 [Note: Bull and Freeman claim that the Gladwell
				//					 approach is intrinsically inefficent because it
				//					 "requires sorting", and propose an alternative
				//					 algorithm that "only" requires three passes over the
				//					 entire set of regions.  Apparently, they didn't
				//					 realize that one could use a heap data structure, in
				//					 which case the time to pop K biggest-error regions
				//					 out of N is only O(K log N), much better than the
				//					 O(N) cost of the Bull and Freeman algorithm if K <<
				//					 N, and it is also much simpler.] */
				//					size_t nR = 0;
				//					for (j = 0; j < fdim; ++j)
				//						ee[j] = regions.ee[j];
				//					do {
				//						if (nR + 2 > nR_alloc) {
				//							nR_alloc = (nR + 2) * 2;
				//							R = (region *) realloc(R, nR_alloc * sizeof(region));
				//							if (!R)
				//								goto bad;
				//						}
				//						R[nR] = heap_pop(&regions);
				//						for (j = 0; j < fdim; ++j)
				//							ee[j].err -= R[nR].ee[j].err;
				//						if (cut_region(R + nR, R + nR + 1))
				//							goto bad;
				//						numEval += r->num_points * 2;
				//						nR += 2;
				//						if (converged(fdim, ee, reqAbsError, reqRelError, norm))
				//							break; /* other regions have small errs */
				//					} while (regions.n > 0 && (numEval < maxEval || !maxEval));
				//					if (eval_regions(nR, R, f, fdata, r)
				//							|| heap_push_many(&regions, nR, R))
				//						goto bad;
				//				} else { /* minimize number of function evaluations */
				//R[0] = heap_pop(&regions); /* get worst region */
				R[0] = regions.poll();

				R[1] = R[0].cut();

				evalError(o, m, R);
				R[0].errmax = errMax(R[0].fdim, R[0].ee);
				R[1].errmax = errMax(R[1].fdim, R[1].ee);

				regions.add(R[0]);
				regions.add(R[1]);
				//					
				//					if (cut_region(R, R + 1) || eval_regions(2, R, f, fdata, r)
				//							|| heap_push_many(&regions, 2, R))
				//						goto bad;
				numEval += num_points * 2;
				//}
			}
			
			if (numEval > maxEval) {
				System.out.println("did not converge after "+numEval+"/"+maxEval+" function evaluations");
			}

			/* re-sum integral and errors */
			for (j = 0; j < fdim; ++j)
				val[j] = err[j] = 0;
			Region[] _regions = regions.toArray(new Region[regions.size()]);
			for (i = 0; i < _regions.length; ++i) {
				for (j = 0; j < fdim; ++j) {
					val[j] += _regions[i].ee[j].val;
					err[j] += _regions[i].ee[j].err;
				}
				//destroy_region(&regions.items[i]);
			}

		}
	}

	public static class RuleGaussKronrod_1d extends Rule {

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
		public void evalError(Object o, Method m, Region[] R) {
			int nR = R.length;

			//System.out.println("evalError for "+nR+" regions");

			alloc_rule_pts(nR);


			for (int iR = 0; iR < nR; ++iR) {
				double center = R[iR].h.centers[0];
				double halfwidth = R[iR].h.halfwidths[0];

				//System.out.println(String.format("region %d:", iR));

				int j;

				int npts = 0;
				pts[0][iR*15 + npts++] = center;
				//System.out.println("pts["+(npts-1)+"]="+pts[0][npts-1]);

				// 7-point Gauss-Legendre
				for (j = 1; j < n-1; j+=2) {
					double w = halfwidth * xgk[j];
					pts[0][iR*15 + npts++] = center - w;
					//System.out.println("pts["+(npts-1)+"]="+pts[0][npts-1]);
					pts[0][iR*15 + npts++] = center + w;
					//System.out.println("pts["+(npts-1)+"]="+pts[0][npts-1]);
				}

				// extended to 15-point by Kronrod
				for (j = 0; j < n; j+=2) {
					double w = halfwidth * xgk[j];
					pts[0][iR*15 + npts++] = center - w;
					//System.out.println("pts["+(npts-1)+"]="+pts[0][npts-1]);
					pts[0][iR*15 + npts++] = center + w;
					//System.out.println("pts["+(npts-1)+"]="+pts[0][npts-1]);
				}

				R[iR].splitDim = 0; /* no choice but to divide 0th dimension */
			}

			try {
				// evaluate function
				vals = (double[][]) m.invoke(o, (Object)(pts)); // [fdim][15]
			} catch (Exception e) {
				//e.printStackTrace();
				throw new RuntimeException(e);
			}

			for (int k=0; k<fdim; ++k) {
				double[] vk = vals[k]; // all 15 points in the given fdim
				//const double *vk = vals + k;
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
						//System.out.println("v="+v);
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
					//r.get(0).ee[k].val = result_kronrod * halfwidth;
					R[iR].ee[k].val = result_kronrod * halfwidth;
					//ret[0][k] = result_kronrod * halfwidth;
					//ret[k] = result_gauss * halfwidth;

//					System.out.println("  gauss="+result_gauss*halfwidth);
//					System.out.println("kronrod="+result_kronrod*halfwidth);

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

	/**
	 * Integrate the given function (o.method) in the range [xmin:xmax] to a relative tolerance relTol or until maxEval function evaluations were used.
	 * @param o
	 * @param method
	 * @param xmin
	 * @param xmax
	 * @param relTol
	 * @param maxEval
	 * @return [0:val, 1:err][fdim]
	 */
	public static double[][] integrate(Object o, String method, double[] xmin, double[] xmax, double relTol, double absTol, int maxEval) {
		double[][] ret = null;

		// dimensionality of parameter space
		final int dim;
		if (xmin==null || xmin.length==0) {
			throw new RuntimeException("xmin cannot be null or empty");
		} else {
			dim = xmin.length;
			System.out.println("dim="+dim);
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
			Method m = o.getClass().getDeclaredMethod(method, double[][].class);
			double[][] f = (double[][]) m.invoke(o, (Object)_center);
			if (f == null || f.length==0 || f[0]==null || f[0].length==0) {
				throw new RuntimeException("Evaluation of given method at interval center failed");
			} else {
				fdim=f.length;
			}
			System.out.println("fdim="+fdim);

			Rule r = null;
			if (dim==1) {
				r = new RuleGaussKronrod_1d(dim, fdim, 15);
			} else {
				// 5-7 Genz-Malik for dim>1
			}




			double[] val = new double[fdim];
			double[] err = new double[fdim];
			Arrays.fill(err, Double.POSITIVE_INFINITY);

			Hypercube h = new Hypercube().initFromRanges(xmin, xmax);

			r.cubature(o, m, maxEval, relTol, absTol, val, err, h);

			ret = new double[][] { val, err };

		} catch (Exception e) {
			throw new RuntimeException(e);
		}


		return ret;
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

//		class CosClass {
//			@SuppressWarnings("unused")
//			public double[][] eval(double[][] x) {
//				double[] ret = new double[x[0].length];
//				for (int i=0; i<x[0].length; ++i) {
//					ret[i] = Math.cos(x[0][i]);
//				}
//				return new double[][] { ret };
//			}
//		}
//		CosClass cc = new CosClass();
//
//		class F0_class {
//			@SuppressWarnings("unused")
//			public double[][] eval(double[][] x) {
//				double[] ret = new double[x[0].length];
//				for (int i=0; i<x[0].length; ++i) {
//					ret[i] = 2.0*x[0][i];
//				}
//				return new double[][] { ret };
//			}
//		}
//		F0_class f0 = new F0_class();

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
		
		
		//double[][] F = integrate(cc, "eval",
		//double[][] F = integrate(f0, "eval",
		double[][] F = integrate(f1, "eval",
				new double[] { 0.0    }, // lower integration limit
				new double[] { 1.0    }, // upper integration limit
				1.0e-6, // relative tolerance
				0.0, // no absolute tolerance requirement
				//Double.POSITIVE_INFINITY,
				110); // max. number of function evaluations

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
