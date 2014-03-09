package com.nian.firstproject.shared;

import java.util.List;

//import com.google.gwt.i18n.client.NumberFormat;
import com.google.gwt.user.client.rpc.IsSerializable;

public class Branch implements IsSerializable {

	// NumberFormat de=NumberFormat.getFormat("0.00");

	Branch down = null;// left
	Branch up = null;// right
	Point downLeft;// topleft
	Point upLeft;// topright
	Point downRight;// (bottom left)
	Point upRight;// (bottom right)
	double width = -1;
	int size;

	double midline;
	// height between two clusters which will combine

	double gap = 0.2;

	Cluster cluster;

	public Cluster getCluster() {
		return cluster;
	}

	public Branch() {

	}

	public Branch(Cluster cluster, double midline) {
		this.cluster = cluster;

		// where is height???????
		width = cluster.getWidth();

		// use where : size????
		size = cluster.size();

		this.midline = midline;

		// distance to mid-line

		// + - no gap x=y; down=left; up=right

		double xdown = midline;

		double xup = midline;

		double ydown = midline - gap;

		double yup = midline + gap;

		// double ydown = midline;

		// double yup = midline;

		if (cluster.getDown() != null && cluster.getDown().size() > 1) {
			ydown += (cluster.getDown().getUp().size());
			xdown = cluster.getDown().getWidth();
		} else {
			ydown += 1;
			// yleft=height-LEAF_HEIGHT;
		}
		if (cluster.getUp() != null && cluster.getUp().size() > 1) {
			yup -= cluster.getUp().getDown().size();
			xup = cluster.getUp().getWidth();
		} else {
			yup -= 1;
			// yright=height-LEAF_HEIGHT;
		}

		downLeft = new Point(width, ydown);
		downRight = new Point(xdown, ydown);

		upLeft = new Point(width, yup);
		upRight = new Point(xup, yup);

	}

	public void print() {
		System.out.println(this);
		if (down != null) {
			down.print();
		}
		if (up != null) {
			up.print();
		}
	}

	public void getAllBranches(List<Branch> list) {
		list.add(this);
		if (down != null) {
			down.getAllBranches(list);
		}
		if (up != null) {
			up.getAllBranches(list);
		}
	}

	public Branch getDown() {
		return down;
	}

	public void setDown(Branch down) {
		this.down = down;
	}

	public Branch getUp() {
		return up;
	}

	public void setUp(Branch up) {
		this.up = up;
	}

	public Point getUpLeft() {
		return upLeft;
	}

	public void setUpLeft(Point upLeft) {
		this.upLeft = upLeft;
	}

	public Point getUpRight() {
		return upRight;
	}

	public void setUpRight(Point upRight) {
		this.upRight = upRight;
	}

	public Point getDownLeft() {
		return downLeft;
	}

	public void setDownLeft(Point downLeft) {
		this.downLeft = downLeft;
	}

	public Point getDownRight() {
		return downRight;
	}

	public void setDownRight(Point downRight) {
		this.downRight = downRight;
	}

	public double getWidth() {
		return width;
	}

	public void setWidth(double width) {
		this.width = width;
	}

	@Override
	public String toString() {
		return "[" + upLeft + " " + upRight + " " + downLeft + " " + downRight
				+ "] width=" + width + " size=" + size;
	}

}
