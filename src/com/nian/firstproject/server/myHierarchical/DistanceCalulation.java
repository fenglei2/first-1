package com.nian.firstproject.server.myHierarchical;

import java.util.ArrayList;
import java.util.List;

import com.nian.firstproject.shared.Cluster;

import Jama.Matrix;

public class DistanceCalulation {

	List<List<Double>> rows;
	List<List<Double>> matrix;// original rows
	List<Cluster> clusters;
	Matrix copheneticMatrix;

	private boolean process = false;

	public DistanceCalulation(Matrix disMatrix) {

		int size = disMatrix.getColumnDimension();
		copheneticMatrix = new Matrix(size, size);

		rows = new ArrayList<List<Double>>();
		matrix = new ArrayList<List<Double>>();
		for (int i = 0; i < size; i++) {
			rows.add(new ArrayList<Double>());
			matrix.add(new ArrayList<Double>());
		}

		for (int row = 0; row < size; row++) {
			for (int col = 0; col <= row; col++) {
				double dist = disMatrix.get(row, col);
				rows.get(row).add(dist);
				matrix.get(row).add(dist);
			}
		}

	}

	public List<Cluster> mergeTwoClosestClusters(DistanceInfo distanceInfo,
			List<Cluster> partition) {

		int col = distanceInfo.column;
		int row = distanceInfo.row;
		double dis = distanceInfo.distance;

		if (process)
			System.out.println("closcs.umbc pair: " + row + "-" + col + " dis="
					+ dis);
		this.clusters = partition;
		// distanceInfo is the shortest
		merge(distanceInfo.row, distanceInfo.column);
		// after merge, the size of rows decrease 1

		Cluster newCluster = new Cluster(clusters.get(distanceInfo.row),
				clusters.get(distanceInfo.column), distanceInfo.distance);

		setcopheneticMatrix(row, col, dis);

		clusters.remove(col);
		clusters.remove(row);
		clusters.add(row, newCluster);
		return clusters;

	}

	private void setcopheneticMatrix(int row, int col, double dis) {

		for (int r : clusters.get(row).getPatternIndexes()) {

			for (int c : clusters.get(col).getPatternIndexes()) {
				if (r > c) {
					copheneticMatrix.set(r, c, dis);
				} else {
					copheneticMatrix.set(c, r, dis);
				}
			}
		}

	}

	// setCell(min, i, distance_k_rs);
	private void setCell(int row, int column, double value) {
		if (row >= column)
			rows.get(row).set(column, value);
		else
			rows.get(column).set(row, value);
	}

	public void merge(int row, int column) {
		int min = Math.min(row, column);
		int max = Math.max(row, column);
		// 0 1 (2 3) 4 5 combine 2,3
		// need to update 02 03, 12 13, 42 43, 52 53
		for (int i = 0; i < rows.size(); i++) {
			if (i != row && i != column) {
				double distance = getDistance(row, column, i);
				if (process)
					System.out.println("merge " + row + "," + column
							+ " newdis to " + i + "=" + distance);

				setCell(min, i, distance);// update one object's distance,
											// the other ready to be removed

			}
		}
		rows.remove(max);// remove another one, so it looks like merge 2
							// objects,
		// treat these 2 objects as 1 object

		for (int i = max; i < rows.size(); i++) {
			rows.get(i).remove(max);
		}
		// correct the indexes of the sorted distances
	}

	// update (i,comb1)(i,comb2), comb1 comb2 are the values need to be
	// combined,
	// so all other values dis(i,comb1) dis(i,comb2)(i=0~n,except comb1,comb2)
	// need to be updated
	private double getDistance(int comb1, int comb2, int i) {
		List<Integer> cluster1 = clusters.get(comb1).getPatternIndexes();
		List<Integer> cluster2 = clusters.get(comb2).getPatternIndexes();
		List<Integer> cluster = clusters.get(i).getPatternIndexes();
		List<Integer> clusterCombine = combine(cluster1, cluster2);

		return clusterDistance(clusterCombine, cluster);
	}

	private double clusterDistance(List<Integer> clusterCombine,
			List<Integer> cluster) {
		if (clusterCombine.size() == 0 || cluster.size() == 0) {
			System.err.println("objects in cluster error");
			return Double.MIN_VALUE;
		}
		// System.out.println("combine: "+clusterCombine+" cluster: "+cluster);
		double sum = 0;

		for (int x : clusterCombine) {
			for (int y : cluster) {
				double dis = 0;
				if (x > y) {
					dis = matrix.get(x).get(y);
				} else {
					dis = matrix.get(y).get(x);
				}
				sum += dis;
			}
		}

		return sum / (clusterCombine.size() * cluster.size());
	}

	private List<Integer> combine(List<Integer> cluster1, List<Integer> cluster2) {
		List<Integer> clusterCombine = new ArrayList<Integer>();
		clusterCombine.addAll(cluster1);
		clusterCombine.addAll(cluster2);

		return clusterCombine;
	}

	// get the minimal dissimilarity
	public DistanceInfo getClosestPair() {
		if (rows.size() == 1)
			return null;

		double min = Double.MAX_VALUE;
		int row = -1;
		int column = -1;
		for (int i = 0; i < rows.size(); i++) {
			for (int j = 0; j < i; j++) {
				if (rows.get(i).get(j) < min) {
					min = rows.get(i).get(j);
					row = i;
					column = j;
				}
			}
		}

		// System.out.println("min dis: "+min+" i="+row+" j="+column);

		return new DistanceInfo(min, row, column);
	}

	public Matrix getCopheneticMatrix() {
		return copheneticMatrix;
	}

	public class DistanceInfo implements Comparable<DistanceInfo> {
		double distance;
		int row;
		int column;

		DistanceInfo(double distance, int row, int column) {
			this.distance = distance;
			this.row = Math.min(row, column);
			this.column = Math.max(row, column);
		}

		@Override
		public int compareTo(DistanceInfo o) {
			return Double.compare(distance, o.distance);
		}

		@Override
		public String toString() {
			return "<" + row + "," + column + "> = " + distance;
		}

	}

}
