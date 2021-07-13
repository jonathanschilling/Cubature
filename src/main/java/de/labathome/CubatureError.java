package de.labathome;

/**
 * Different ways of measuring the absolute and relative error when we have
 * multiple integrands, given a vector e of error estimates in the individual
 * components of a vector v of integrands. These are all equivalent when there
 * is only a single integrand.
 */
public enum CubatureError {
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
