package de.labathome;

/**
 * Integrand interface for Cubature.
 */
public interface Integrand {
	/**
	 * Evaluate the nFDim-dimensional integrand at nPoints nDim-dimensional locations.
	 * @param x [nDim][nPoints] locations where to evaluate the integrand
	 * @param fdata some arbitrary data that is passed from the integrate() call into the integrand
	 * @return [nFDim][nPoints] function values at x
	 */
	public abstract double[][] eval(final double[][] x, Object fdata);
}