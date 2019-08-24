package de.labathome;

import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

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
	
	public static class Region {
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
	}


	/**
	 * Integrate the given function (o.method) in the range [xmin:xmax] to a relative tolerance relTol or until maxEval function evaluations were used.
	 * @param o
	 * @param method
	 * @param xmin
	 * @param xmax
	 * @param relTol
	 * @param maxEval
	 * @return
	 */
	public static double[] integrate(Object o, String method, double[] xmin, double[] xmax, double[] relTol, double[] absTol, int maxEval) {
		double[] ret = null;

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

		if ((relTol==null || relTol.length==0) && (absTol==null || absTol.length==0)) {
			throw new RuntimeException("relTol and absTol cannot be both null or empty");
		}
		if (relTol != null && relTol.length != dim && relTol.length != 1) {
			throw new RuntimeException("relTol should have length of 1 (same tol for all dimensions) or same as dim, but has "+relTol.length);
		}
		if (absTol != null && absTol.length != dim && absTol.length != 1) {
			throw new RuntimeException("absTol should have length of 1 (same tol for all dimensions) or same as dim, but has "+absTol.length);
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
//
//			// dummy result
//			ret = new double[fdim];
//			for (int i=0; i<fdim; ++i) {
//				ret[i] = f[i][0];
//			}
	
			
			double[] val = new double[fdim];
			double[] _err = new double[fdim];
			Arrays.fill(_err, Double.POSITIVE_INFINITY);
	
			Hypercube h = new Hypercube().initFromRanges(xmin, xmax);
			int numEval = 0;
			
			if (dim==1) {
				// actually one-dimensional adaptive quadrature: 15-point adptive Gauss-Kronrod
				
				// evalError
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
	
				double center = h.centers[0];
				double halfwidth = h.halfwidths[0];
				
				double[][] pts = new double[1][15]; // 15*dim
				//double[][] vals = new double[fdim][15];
				
				int j;
				
				int npts = 0;
				pts[0][npts++] = center;
				System.out.println("pts["+(npts-1)+"]="+pts[0][npts-1]);
				
				// 7-point Gauss-Legendre
				for (j = 1; j < n-1; j+=2) {
					double w = halfwidth * xgk[j];
					pts[0][npts++] = center - w;
					System.out.println("pts["+(npts-1)+"]="+pts[0][npts-1]);
					pts[0][npts++] = center + w;
					System.out.println("pts["+(npts-1)+"]="+pts[0][npts-1]);
				}
				
				// extended to 15-point by Kronrod
				for (j = 0; j < n; j+=2) {
					double w = halfwidth * xgk[j];
					pts[0][npts++] = center - w;
					System.out.println("pts["+(npts-1)+"]="+pts[0][npts-1]);
					pts[0][npts++] = center + w;
					System.out.println("pts["+(npts-1)+"]="+pts[0][npts-1]);
				}
					
				// evaluate function
				double[][] vals = (double[][]) m.invoke(o, (Object)pts); // [fdim][15]

				ret = new double[fdim];

				/*for (*/int k=0;/* k<fdim; ++k) {*/
					double[] vk = vals[k]; // all 15 points in the given fdim
					//const double *vk = vals + k;
					//for (iR = 0; iR < nR; ++iR) {
						//const double halfwidth = R[iR].h.data[1];
						double result_gauss = vk[0] * wg[n / 2 - 1];
						double result_kronrod = vk[0] * wgk[n - 1];
						//double result_abs = Math.abs(result_kronrod);
						//double result_asc, mean, err;

						/* accumulate integrals */
						npts = 1;
						for (j = 0; j < (n - 1) / 2; ++j) {
							int j2 = 2 * j + 1;
							double v = vk[npts] + vk[npts+1]; // +/- w are next to each other
							System.out.println("v="+v);
							result_gauss += wg[j] * v;
							result_kronrod += wgk[j2] * v;
							//result_abs += wgk[j2] * (Math.abs(vk[fdim * npts]) + Math.abs(vk[fdim * npts + fdim]));
							npts += 2;
						}
						for (j = 0; j < n / 2; ++j) {
							int j2 = 2 * j;
							result_kronrod += wgk[j2] * (vk[npts] + vk[npts+1]); // +/- w are next to each other
							//result_abs += wgk[j2] * (Math.abs(vk[fdim * npts]) + Math.abs(vk[fdim * npts + fdim]));
							npts += 2;
						}

						/* integration result */
						//r.get(0).ee[k].val = result_kronrod * halfwidth;
						
						ret[k] = result_kronrod * halfwidth;
						//ret[k] = result_gauss * halfwidth;
						
						System.out.println("  gauss="+result_gauss*halfwidth);
						System.out.println("kronrod="+result_kronrod*halfwidth);
//
//						/* error estimate
//						 (from GSL, probably dates back to QUADPACK
//						 ... not completely clear to me why we don't just use
//						 Math.abs(result_kronrod - result_gauss) * halfwidth */
//						mean = result_kronrod * 0.5;
//						result_asc = wgk[n - 1] * Math.abs(vk[0] - mean);
//						npts = 1;
//						for (j = 0; j < (n - 1) / 2; ++j) {
//							int j2 = 2 * j + 1;
//							result_asc += wgk[j2]
//									* (Math.abs(vk[fdim * npts] - mean)
//											+ Math.abs(vk[fdim * npts + fdim] - mean));
//							npts += 2;
//						}
//						for (j = 0; j < n / 2; ++j) {
//							int j2 = 2 * j;
//							result_asc += wgk[j2]
//									* (Math.abs(vk[fdim * npts] - mean)
//											+ Math.abs(vk[fdim * npts + fdim] - mean));
//							npts += 2;
//						}
//						err = Math.abs(result_kronrod - result_gauss) * halfwidth;
//						result_abs *= halfwidth;
//						result_asc *= halfwidth;
//						if (result_asc != 0 && err != 0) {
//							double scale = Math.pow((200 * err / result_asc), 1.5);
//							err = (scale < 1) ? result_asc * scale : result_asc;
//						}
//						if (result_abs > DBL_MIN / (50 * DBL_EPSILON)) {
//							double min_err = 50 * DBL_EPSILON * result_abs;
//							if (min_err > err)
//								err = min_err;
//						}
//						r.get(0).ee[k].err = err;
//
//						/* increment vk to point to next batch of results */
//						//vk += 15 * fdim;
					/*}*/
				
				
	
	
			} else {
				// multi-dimensional adaptive cubature: 7th (5th) order Genz-Malik 
	
	
	
	
	
			}
	






















		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		

		return ret;
	}



	public static void demoIntegration() {

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

		class CosClass {
			@SuppressWarnings("unused")
			public double[][] eval(double[][] x) {
				double[] ret = new double[x[0].length];
				for (int i=0; i<x[0].length; ++i) {
					ret[i] = Math.cos(x[0][i]);
				}
				return new double[][] { ret };
			}
		}
		CosClass cc = new CosClass();



		double[] F = integrate(cc, "eval",
				new double[] { 0.0    }, // lower integration limit
				new double[] { 1.0    }, // upper integration limit
				new double[] { 1.0e-3 }, // relative tolerance
				null, // no absolute tolerance requirement
				100); // max. number of function evaluations

		if (F != null) {
			System.out.println("integration result: "+F[0]);
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
