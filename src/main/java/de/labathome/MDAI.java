package de.labathome;

import java.lang.reflect.Method;

import aliceinnets.python.jyplot.JyPlot;

public class MDAI {

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
	 * @param c class in which the given method is defined
	 * @param method name of a method which takes double[dim][nPoints] as argument and returns double[fdim][nPoints]
	 * @param xmin[dim] left interval borders
	 * @param xmax[dim] right interval borders
	 * @param nPoints number of points at which to evaluate the given method
	 * @return [fdim][nPoints] evaluated function values
	 */
	public static double[][] eval(Class<?> c, String method, double[] xmin, double[] xmax, int nPoints) {
		double[][] x = new double[xmin.length][nPoints];
		for (int j=0; j<xmin.length; ++j) {
			for (int i=0; i<nPoints; ++i) {
				x[j][i] = xmin[j] + (xmax[j]-xmin[j])*i/(nPoints-1.0);
			}
		}
		
		try {
			Method m = c.getDeclaredMethod(method, double[][].class);
			return (double[][]) m.invoke((Object)MDAI.class, (Object)x);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	public static void demoEval() {
		double[] xmin = new double[] { 0.0        , 0.0         };
		double[] xmax = new double[] { 2.0*Math.PI, 2.0*Math.PI };
		int nPoints = 100;
		
		double[][] x = new double[xmin.length][nPoints];
		for (int j=0; j<xmin.length; ++j) {
			for (int i=0; i<nPoints; ++i) {
				x[j][i] = xmin[j] + (xmax[j]-xmin[j])*i/(nPoints-1.0);
			}
		}
		
		double[][] sinCos_x = MDAI.eval(MDAI.class, "oD_tD_sinCos", xmin, xmax, nPoints);
		
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
	
	
	public static void main(String[] args) {
		
		demoEval();
		
	}

}
