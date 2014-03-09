package com.nian.firstproject.server.PAM;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import com.nian.firstproject.shared.Pattern;
import com.nian.firstproject.server.common.Distance;

import Jama.Matrix;

/*
 * This is the one I am using 
 */
public class MyKMedoidsDis {
	public int k;
	int n;
	private Pattern[] initMedoids;// once set, never change

	List<Integer> medoids;// medoids for each iteration
	List<Integer> nonmedoids;
	private List<Pattern> dataSet;

	private double sse;
	private double minsse = Double.MAX_VALUE;
	List<Integer[]> comb = new ArrayList<Integer[]>();
	private boolean process = true;
	HashMap<Integer, List<Pattern>> cluster_result = new HashMap<Integer, List<Pattern>>();

	private Matrix disMatrix;

	public MyKMedoidsDis() {
	}

	public MyKMedoidsDis(Matrix matrix, int k) {
		initial(matrix, k);
	}

	public void initial(Matrix matrix, int k) {
		setDataSet(matrix.getColumnDimension());
		this.disMatrix = matrix;
		this.k = k;
		initMedoids = chooseInitialMedoids();
		setInitialMedoids();

		if (process)
			printMedoids();

	}

	public void cluster_random() {

		cluster_simple();

		if (process) {
			System.out.print("k=" + k + " [");
			for (int me : medoids) {
				System.out.print(me + " ");
			}
			System.out.print("] ");
			System.out.println("sse: " + minsse);
		}

	}

	private void cluster_simple() {
		HashMap<Integer, List<Pattern>> result = cluster();
		cluster_result = result;
		minsse = sse;
	}

	// add 0->n-1 patterns to dataset
	private void setDataSet(int size) {
		dataSet = new ArrayList<Pattern>();
		for (int i = 0; i < size; i++) {
			Pattern pat = new Pattern(i);
			dataSet.add(pat);
		}

	}

	public void cluster_swap() {

		HashMap<Integer, List<Pattern>> result = cluster();
		cluster_result = result;
		minsse = sse;

		/*
		 * if (process) { System.out.println("build sse: " + sse);
		 * System.out.println(" --------   build result ------- ");
		 * printmap(cluster_result);
		 * System.out.println("\n  ---   start swap ---\n"); }
		 * 
		 * for (int i = 0; i < initMedoids.length; i++) { for (int h = 0; h <
		 * nonmedoids.size(); h++) { // swap medoid i and nonmedoid h swap(i,
		 * h);
		 * 
		 * if (process) { System.out.print("round -> i: " + i + " h: " + h +
		 * "   "); printMedoidNonMedoidInt(); }
		 * 
		 * result = cluster(); if (process) System.out.println("SSE: " + sse);
		 * 
		 * if (sse < minsse) { cluster_result = result; minsse = sse; } else {
		 * swap(i, h); }
		 * 
		 * if (process == true) System.out.println("----" + i + "--" + h +
		 * "------------------------------ \n");
		 * 
		 * }
		 * 
		 * }
		 * 
		 * if (process) {
		 * System.out.println(" -----------PAM DONE --------------");
		 * System.out.println("final sse in cluster_swap: " + minsse); }
		 */

	}

	private void setInitialMedoids() {
		if (medoids == null) {
			medoids = new ArrayList<Integer>();
		} else {
			medoids.clear();
		}
		for (Pattern pat : initMedoids) {
			medoids.add(pat.getIndex());
		}
		setNonMedoids();
	}

	// swap medoid and nonmedoid---------------
	@SuppressWarnings("unused")
	private void swap(int i, int h) {
		if (process) {
			System.out
					.println("<before swap>: " + medoids + " | " + nonmedoids);
		}
		int value = medoids.get(i);
		medoids.set(i, nonmedoids.get(h));
		nonmedoids.set(h, value);
		if (process) {
			System.out.println("<After swap>: " + medoids + " | " + nonmedoids);
			System.out.println();
		}

	}

	private void setNonMedoids() {
		nonmedoids = new ArrayList<Integer>();

		for (int i = 0; i < disMatrix.getColumnDimension(); i++) {
			if (!medoids.contains(i)) {
				nonmedoids.add(i);
			}
		}

	}

	@SuppressWarnings("unused")
	private void printMedoidNonMedoidInt() {
		System.out.print("meoids: ");
		for (int medoid : medoids) {
			System.out.print(medoid + " ");
		}
		System.out.print(" | ");
		System.out.print("nonmeoids: ");
		for (int nonmedoid : nonmedoids) {
			System.out.print(nonmedoid + " ");
		}
		System.out.println();

	}

	public HashMap<Integer, Integer> getClusterAssignmt() {
		HashMap<Integer, Integer> clusterAssignmt = new HashMap<Integer, Integer>();

		for (int key : cluster_result.keySet()) {
			for (Pattern pat : cluster_result.get(key)) {
				clusterAssignmt.put(pat.getIndex(), key);
			}
		}

		return clusterAssignmt;
	}

	public HashMap<Integer, List<Pattern>> cluster() {
		HashMap<Integer, List<Pattern>> clusteredPatterns = null;

		if (process) {
			System.out.print("new round ->");
			printMedoids();
		}

		clusteredPatterns = new HashMap<Integer, List<Pattern>>();
		for (int i = 0; i < k; i++) {
			clusteredPatterns.put(i, new ArrayList<Pattern>());
		}

		List<Integer> medoidsInt = new ArrayList<Integer>();
		for (int i = 0; i < medoids.size(); i++) {
			int index = medoids.get(i);
			Pattern medo = new Pattern(index, i);
			clusteredPatterns.get(i).add(medo);
			medoidsInt.add(index);
		}

		// for every point in whole dataset, find the nearest medoid to it
		for (Pattern p : dataSet) {

			if (medoidsInt.contains(p.getIndex())) {
				continue;
			}

			double minDistance = Double.MAX_VALUE;

			for (int i = 0; i < medoids.size(); i++) {

				double distance = Distance.calculateDissimlarity(p.getIndex(),
						medoids.get(i), disMatrix);

				if (distance < minDistance) {
					minDistance = distance;
					// each node associate to the closest medoid
					p.setClusterNo(i);
				}

			}// inner for

			// p.getClusterNo is the nearest medoid's cluster no
			clusteredPatterns.get(p.getClusterNo()).add(p);

		}// outer for

		calculateSSE(clusteredPatterns);

		return clusteredPatterns;
	}

	@SuppressWarnings("unused")
	private void printmap(HashMap<Integer, List<Pattern>> clusteredPatterns) {
		for (int key : clusteredPatterns.keySet()) {
			List<Pattern> pats = clusteredPatterns.get(key);
			System.out.print("cluster" + key + ": ");
			for (Pattern p : pats) {
				System.out.print(p.getIndex() + 1 + " ");
			}
			System.out.println();
		}

	}

	private void calculateSSE(HashMap<Integer, List<Pattern>> clusteredPatterns) {
		sse = 0;
		for (int key : clusteredPatterns.keySet()) {
			List<Pattern> pats = clusteredPatterns.get(key);
			// medoid=medoids.get(key);
			for (Pattern p : pats) {
				double distance = Distance.calculateDissimlarity(p.getIndex(),
						medoids.get(key), disMatrix);
				sse += distance;
			}
		}
		sse /= disMatrix.getColumnDimension();

	}

	private Pattern[] chooseInitialMedoids() {
		// the first chosen medoid
		Pattern medoidPattern = getFirstMedoid();
		List<Pattern> medoids = new ArrayList<Pattern>();
		medoidPattern.setClusterNo(0);
		// put the first selected medoid into medoids list
		medoids.add(medoidPattern);

		// select other medoids from 1 to k-1
		for (int i = 1; i < k; i++) {
			medoidPattern = getNextMedoid(medoids);
			// System.out.println("new medoids: "+medoidPattern.getIndex());

			// ClusterNo to find its size and survival curve
			medoidPattern.setClusterNo(i);
			medoids.add(medoidPattern);
		}

		// put list to array
		Pattern[] medoidPatterns = new Pattern[k];
		for (int i = 0; i < k; i++) {
			Pattern p = medoids.get(i);
			medoidPatterns[i] = p;
		}

		return medoidPatterns;
	}

	private Pattern getFirstMedoid() {
		Pattern firstMedoid = null;
		double minTotalDistance = Double.MAX_VALUE;

		for (Pattern candidatePattern : dataSet) {
			double totalDistance = 0; // total distance for each candidate
			for (Pattern mido : dataSet) {

				double distance = Distance.calculateDissimlarity(
						candidatePattern, mido, disMatrix);
				totalDistance += distance;

			}
			if (totalDistance < minTotalDistance) {

				minTotalDistance = totalDistance;
				firstMedoid = candidatePattern;

			}

		}

		return firstMedoid;
	}

	private Pattern getNextMedoid(List<Pattern> medoidPatterns) {

		// candidatePatterns = all dataset - selectedMedoids
		// for the first one: allPatterns = fromPatterns = dataSet, skip this
		// function, return alldata
		List<Pattern> candidates = getCandidates(medoidPatterns);

		// select maximum sum of dissimilarity to other unselected nodes
		Pattern nextMedoid = maxDisSumToNonMedoids(candidates, medoidPatterns);

		return nextMedoid;
	}

	private Pattern maxDisSumToNonMedoids(List<Pattern> candidates,
			List<Pattern> medoidPatterns) {
		Pattern nextMedoid = null;
		double maxTotalDistance = -1;
		for (Pattern candidatePattern : candidates) {
			double sumOfCij = 0; // total distance for each candidate
			for (Pattern unselected : candidates) {

				double dij = Distance.calculateDissimlarity(candidatePattern,
						unselected, disMatrix);
				double distToMedoid = selectClosestMedoid(medoidPatterns,
						unselected) - dij;
				double distance = Math.max(distToMedoid - dij, 0);
				sumOfCij += distance;
				// System.out.println("dis: "+distance);

			}
			if (sumOfCij > maxTotalDistance) {

				maxTotalDistance = sumOfCij;
				nextMedoid = candidatePattern;

			}

		}
		return nextMedoid;
	}

	private double selectClosestMedoid(List<Pattern> medoidPatterns,
			Pattern candidate) {
		// Pattern closestMedoid=null;

		double minDis = Double.MAX_VALUE;
		for (Pattern mido : medoidPatterns) {
			double distance = Distance.calculateDissimlarity(candidate, mido,
					disMatrix);
			// System.out.println("distance: "+distance);
			if (distance < minDis) {
				minDis = distance;
			}

		}

		// System.out.println(" min : "+minDis+"\n");

		return minDis;
	}

	private List<Pattern> getCandidates(List<Pattern> medoidPatterns) {

		List<Pattern> candidates = new ArrayList<Pattern>();
		candidates.addAll(dataSet);
		if (medoidPatterns.size() != dataSet.size()) {
			for (int i = 0; i < medoidPatterns.size(); i++) {
				candidates.remove(medoidPatterns.get(i));
			}
		}

		return candidates;
	}

	// print medoids list
	private void printMedoids() {
		System.out.print("The medoid list: ");
		for (int pat : medoids) {
			System.out.print(pat + 1 + " ");

		}
		System.out.println();
	}

	public double getMinsse() {
		return minsse;
	}

	public void setMinsse(double minsse) {
		this.minsse = minsse;
	}

	public HashMap<Integer, List<Pattern>> getCluster_result() {
		return cluster_result;
	}

	public void setCluster_result(HashMap<Integer, List<Pattern>> cluster_result) {
		this.cluster_result = cluster_result;
	}

	public List<Pattern> getDataSet() {
		return dataSet;
	}

	public void setDataSet(List<Pattern> dataSet) {
		this.dataSet = dataSet;
	}

	public Matrix getDisMatrix() {
		return disMatrix;
	}

	public void setDisMatrix(Matrix disMatrix) {
		this.disMatrix = disMatrix;
	}

	public void printarr(Pattern[] newMedoids) {
		System.out.print("new medoids ");
		for (Pattern pat : newMedoids) {
			System.out.print(pat.getIndex() + " ");
		}
		System.out.println();
	}

	public void setInitialMedoids(Matrix matrix, int[] initMedoids, int k) {
		n = matrix.getColumnDimension();
		setDataSet(n);
		this.disMatrix = matrix;
		this.k = k;
		setRandomInitialMedoids(initMedoids);

	}

	private void setRandomInitialMedoids(int[] init) {
		if (medoids == null) {
			medoids = new ArrayList<Integer>();
		} else {
			medoids.clear();
		}
		for (int in : init) {
			medoids.add(in);
		}
		setNonMedoids();
	}

}
