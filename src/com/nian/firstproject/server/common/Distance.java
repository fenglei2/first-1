package com.nian.firstproject.server.common;

import com.nian.firstproject.shared.Pattern;

import Jama.Matrix;

public class Distance {
	private Pattern destPattern;
	private double distance;

	public double getDistance() {
		return this.distance;
	}

	public Pattern getDestination() {
		return this.destPattern;
	}

	public static double calculateDissimlarity(Pattern pattern1,
			Pattern pattern2, Matrix disMatrix) {
		double distance;
		int size = disMatrix.getColumnDimension();
		if (pattern1.getIndex() >= size || pattern2.getIndex() >= size) {
			System.err.println(" pattern index error ");
			return -1;
		}

		if (pattern1.getIndex() > pattern2.getIndex()) {
			distance = disMatrix.get(pattern1.getIndex(), pattern2.getIndex());
		} else {
			distance = disMatrix.get(pattern2.getIndex(), pattern1.getIndex());
		}
		return distance;
	}

	public static double calculateDissimlarity(int a, int b, Matrix disMatrix) {
		double distance = -1;
		if (a > b) {
			distance = disMatrix.get(a, b);
		} else {
			distance = disMatrix.get(b, a);
		}
		return distance;
	}

	public static double calculateDistance(Matrix srcVector, Matrix destVector) {
		Matrix differenceVector = destVector.minus(srcVector);
		double patternDistance = differenceVector.normF();
		return patternDistance;
	}
}
