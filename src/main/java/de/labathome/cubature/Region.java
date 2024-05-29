package de.labathome.cubature;

public class Region implements Comparable<Region> {

    Hypercube h;
    int splitDim;

    /** dimensionality of vector integrand */
    int fdim;

    /** array of length fdim */
    ErrorEstimate[] ee;

    /** max ee[k].err */
    double errmax;

    public static Region init(Hypercube h, int fdim) {
        Region r = new Region();
        r.h = h.clone();
        r.splitDim = 0;
        r.fdim = fdim;
        r.ee = new ErrorEstimate[fdim];
        for (int i = 0; i < fdim; ++i) {
            r.ee[i] = new ErrorEstimate();
        }
        r.errmax = Double.POSITIVE_INFINITY;
        return r;
    }

    /** return new Region and modify current one */
    public Region cut() {
        h.halfwidths[splitDim] *= 0.5;
        h.volume *= 0.5;
        Region r2 = Region.init(h.clone(), fdim);
        h.centers[splitDim] -= h.halfwidths[splitDim];
        r2.h.centers[splitDim] += h.halfwidths[splitDim];
        return r2;
    }

    @Override
    public int compareTo(Region o) {
        if (this.errmax == o.errmax) {
            return 0;
        }

        /**
         * Java implements the PriorityQueue type as a MinHeap, which means that the
         * root is always the smallest element as judged by the compareTo method of the
         * defining type. The integration algorithm however expects a MaxHeap (where the
         * root is the largest element) in order to work on the worst regions first.
         * Hence, we need to reverse the order of the Java PriorityQueue by reversing
         * the sign of the compareTo method output.
         */
        return (this.errmax < o.errmax) ? 1 : -1;
    }
}
