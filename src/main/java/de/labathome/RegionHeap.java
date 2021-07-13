package de.labathome;

import java.util.PriorityQueue;

public class RegionHeap extends PriorityQueue<Region> {
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