package com.nian.firstproject.server.myHierarchical.dendrogram;

import com.nian.firstproject.shared.Branch;
import com.nian.firstproject.shared.Cluster;

//import com.nian.firstproject.shared.Pattern;

public class Dendrogram {

	public Cluster root;
	public static Branch rootBranch;

	// public Pattern pattern;

	/**
	 * public static void main(String[] args) { Cluster c1 = new Cluster("c1",
	 * 1); Cluster c2 = new Cluster("c2", 2); Cluster c3 = new Cluster("c3", 3);
	 * Cluster c4 = new Cluster("c4", 4); Cluster c5 = new Cluster("c5", 5);
	 * Cluster c23 = new Cluster(c2, c3, 0.3); Cluster c45 = new Cluster(c4, c5,
	 * 0.5); Cluster c123 = new Cluster(c1, c23, 0.7); Cluster root = new
	 * Cluster(c123, c45, 0.9); Dendrogram den = new Dendrogram(root); Branch br
	 * = den.rootBranch; br.print(); }
	 */

	public Dendrogram(Cluster root) {
		this.root = root;
		rootBranch = root.createBranch(0);
		System.out.println();
	}

	public void draw(Branch branch) {
		drawLine(branch);
		if (branch.getDown() != null) {
			draw(branch.getDown());
		}
		if (branch.getUp() != null) {
			draw(branch.getUp());
		}

	}

	public Branch getRootBranch() {
		return rootBranch;
	}

	private void drawLine(Branch branch) {

	}

}
