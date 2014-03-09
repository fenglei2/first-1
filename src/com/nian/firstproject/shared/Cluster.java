package com.nian.firstproject.shared;

import java.util.ArrayList;
import java.util.List;

import com.google.gwt.user.client.rpc.IsSerializable;
import com.nian.firstproject.shared.Branch;

public class Cluster implements IsSerializable {

	Cluster down = null;
	Cluster up = null;
	private List<Integer> patternIndexes;// combination id
	double width;
	Branch branch;
	String label;

	public Cluster() {

	}

	public Branch createBranch(double midline) {
		branch = new Branch(this, midline);
		if (down != null && down.size() > 1) {
			branch.setDown(down.createBranch(midline + down.up.size()));
		}
		if (up != null && up.size() > 1) {
			branch.setUp(up.createBranch(midline - up.down.size()));
		}
		return branch;

	}

	public int size() {
		return patternIndexes.size();

	}

	public String getLabel() {
		return label;
	}

	public Cluster(String label, int patternIndex) {
		patternIndexes = new ArrayList<Integer>();
		patternIndexes.add(patternIndex);
		this.label = label;
	}

	public Cluster(Cluster down, Cluster up, double width) {
		this.down = down;
		this.up = up;
		this.width = width;

		patternIndexes = new ArrayList<Integer>();
		for (Integer index : down.getPatternIndexes()) {
			patternIndexes.add(index);
		}
		for (Integer index : up.getPatternIndexes()) {
			patternIndexes.add(index);
		}
	}

	public List<Integer> getPatternIndexes() {
		return patternIndexes;
	}

	public Cluster getDown() {
		return down;
	}

	public Cluster getUp() {
		return up;
	}

	public double getWidth() {
		return width;
	}

	@Override
	public String toString() {
		return label;
	}
}
