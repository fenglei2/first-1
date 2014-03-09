package com.nian.firstproject.server.myHierarchical;

import java.util.ArrayList;
import java.util.List;

import com.nian.firstproject.shared.Cluster;

import com.nian.firstproject.server.common.DataModel1;
import com.nian.firstproject.server.myHierarchical.DistanceCalulation.DistanceInfo;

import Jama.Matrix;

public class HierarchicalCluster {

	private boolean PRINT_CLUSTER_RESULT;
	private boolean process = false;

	Matrix dissimilarityMatrix;
	Matrix copheneticMatrix;

	private List<List<Cluster>> allClusters;

	// List<Cluster> partition;

	public HierarchicalCluster(Matrix disMatrix, boolean print) {
		this.PRINT_CLUSTER_RESULT = print;
		this.dissimilarityMatrix = disMatrix;
		// this.partition.size=partition.size();
	}

	// ----ok---------

	// we need to merge clusters distanceInfo.row and
	// distanceInfo.column
	// into one cluster
	// this happens in two steps:
	// 1. merging the distances
	// 2. merging the clusters themselves that are stored in clusters
	// list
	public List<Cluster> partition() {
		// create a cluster for each pattern. pattern label for every,
		// combine become cluster
		List<Cluster> currPartition = new ArrayList<Cluster>();

		for (int i = 0; i < dissimilarityMatrix.getColumnDimension(); i++) {

			String label = DataModel1.getInstance().getDataCombinations()
					.get(i).getLabel();
			currPartition.add(new Cluster(label, i));
		}

		DistanceCalulation cal = new DistanceCalulation(dissimilarityMatrix);

		@SuppressWarnings("unused")
		int iterationCount = 0;
		allClusters = new ArrayList<List<Cluster>>();
		printClusterResult(currPartition);
		for (DistanceInfo shortestDistInfo = cal.getClosestPair(); shortestDistInfo != null; shortestDistInfo = cal
				.getClosestPair()) {
			allClusters.add(copy(currPartition));

			currPartition = cal.mergeTwoClosestClusters(shortestDistInfo,
					currPartition);
			printClusterResult(currPartition);
			iterationCount++;
		}
		allClusters.add(currPartition);

		copheneticMatrix = cal.getCopheneticMatrix();

		return currPartition;
	}

	private List<Cluster> copy(List<Cluster> currPartition) {
		List<Cluster> copy = new ArrayList<Cluster>();
		copy.addAll(currPartition);
		return copy;
	}

	private void printClusterResult(List<Cluster> partition) {
		if (PRINT_CLUSTER_RESULT) {
			System.out.print("cluster size = " + partition.size() + ": ");
			for (int i = 0; i < partition.size(); i++) {
				List<Integer> print = partition.get(i).getPatternIndexes();
				System.out.print("[");
				for (int ii = 0; ii < print.size() - 1; ii++) {
					System.out.print(print.get(ii) + 1 + ",");
				}
				System.out.print(print.get(print.size() - 1) + 1 + "] ");
			}
			System.out.println();
		}
		if (process)
			System.out.println();

	}

	public Matrix getDissimilarityMatrix() {
		return dissimilarityMatrix;
	}

	public Matrix getCopheneticMatrix() {
		return copheneticMatrix;
	}

	public List<List<Cluster>> getAllClusters() {
		return allClusters;
	}

}
