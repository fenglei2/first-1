package com.nian.firstproject.shared;

import com.google.gwt.user.client.rpc.IsSerializable;

//use to draw (x,y)
public class Point implements IsSerializable {
	double x;
	double y;

	public Point() {

	}

	public double getX() {
		return x;
	}

	public void setX(double x) {
		this.x = x;
	}

	public double getY() {
		return y;
	}

	public void setY(double y) {
		this.y = y;
	}

	public String toString() {
		return "(" + x + "," + y + ")";
	}

	public Point(double x, double y) {
		this.x = x;
		this.y = y;
	}
}
