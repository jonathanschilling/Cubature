package de.labathome;

public class Region implements Comparable<Region> {
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
