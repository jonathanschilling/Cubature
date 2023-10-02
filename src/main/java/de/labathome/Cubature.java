package de.labathome;

import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.function.UnaryOperator;

/* Adaptive multidimensional integration of a vector of integrands.
*
* Copyright (c) 2005-2013 Steven G. Johnson
* Ported to Java in 2019 by Jonathan Schilling (jonathan.schilling@mail.de)
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

	/**
	 * Integrate the given function (o.method) in the range [xmin:xmax] to a relative tolerance relTol or absolute tolerance absTol or until maxEval function evaluations were used.
	 * @param integrand function to integrate; must map t[xdim][nEvalPoints]) -> f[fdim][nEvalPoints]
	 * @param xmin [xdim] vector of lower integration bounds
	 * @param xmax [xdim] vector of upper integration bounds
	 * @param relTol relative tolerance on function values
	 * @param absTol absolute tolerance on function values
	 * @param norm Error norm for vector-valued integrands
	 * @param maxEval absolute tolerance on function values
	 * @param fdata any Object that shall be passed directly to the integrand; can be used to specify additional info/parameter/...
	 * @return [0:val, 1:err][fdim] value and error of integrals
	 */
	public static double[][] integrate(UnaryOperator<double[][]> integrand, double[] xmin, double[] xmax, double relTol, double absTol, CubatureError norm, int maxEval) {

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

		// determine dimensionality of integrand by evaluating it once in the center of the given integration interval
		double[][] _center = new double[dim][1];
		for (int i=0; i<dim; ++i) {
			_center[i][0] = (xmin[i]+xmax[i])/2.0;
		}

		double[][] f = integrand.apply(_center);

		// dimensionality of the integrand
		final int fdim;
		if (f == null || f.length==0 || f[0]==null || f[0].length==0) {
			throw new RuntimeException("Evaluation of given method at interval center failed");
		} else {
			fdim=f.length;
		}

		double[] val = new double[fdim];
		double[] err = new double[fdim];
		Arrays.fill(err, Double.POSITIVE_INFINITY);

		Hypercube h = new Hypercube().initFromRanges(xmin, xmax);

		Rule r = null;
		if (dim==1) {
			r = new RuleGaussKronrod_1d(dim, fdim, 15);
		} else {
			// 5-7 Genz-Malik for dim>1
			int numPoints = Rule75GenzMalik.num0_0(dim)
					+   2 * Rule75GenzMalik.numR0_0fs(dim)
					+       Rule75GenzMalik.numRR0_0fs(dim)
					+       Rule75GenzMalik.numR_Rfs(dim);
			r = new Rule75GenzMalik(dim, fdim, numPoints);
		}

		r.cubature(integrand, maxEval, relTol, absTol, val, err, h, norm);

		return new double[][] { val, err };
	}
}
