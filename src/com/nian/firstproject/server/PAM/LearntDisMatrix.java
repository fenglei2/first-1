package com.nian.firstproject.server.PAM;

import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import com.nian.firstproject.shared.Pattern;
import com.nian.firstproject.server.common.DataModel1;

import Jama.Matrix;

public class LearntDisMatrix {

	private Matrix dissimilarityMatrix;
	private int rounds;
	private HashMap<Integer, HashMap<Integer, Double>> learntDissimilarity;

	private boolean process = false;


	private int size;

	int[] count;

	public LearntDisMatrix(int rounds) {

	
		this.dissimilarityMatrix = DataModel1.getInstance()
				.getLogRankTestMatrix();

		this.rounds = rounds;

	
		this.size = dissimilarityMatrix.getColumnDimension();

		count = new int[size];

	}

	public void run(boolean show) {
		int size = dissimilarityMatrix.getColumnDimension();

		initialLeartDisMatrix(size);

		Random random = new Random();
		HashMap<Integer, Integer> clusterAssign = null;
		int k;
		for (int i = 0; i < rounds; i++) {
			if (process)
				System.out.println("rounds: " + i);

			k = random.nextInt(size - 2) + 2;// k picked from 2 to 79(if
												// size=80)

			clusterAssign = usePamRandom(k);

			setValueForLearntDisMatrix(clusterAssign);

		}

		takeAverageValue(show);
	}

	private int[] randomComb(int n, int k) {
		Random random = new Random();
		int[] set = new int[k];
		int j = 0;
		for (int i = 0; i < n; i++) {
			int r = random.nextInt(n - i);
			if (r < k) {
				set[j++] = i;
				k--;
			}
		}
		return set;
	}

	private HashMap<Integer, Integer> usePamRandom(int k) {
		MyKMedoidsDis pam = new MyKMedoidsDis();
		int[] initMedoids = randomComb(size, k);
		for (int me : initMedoids) {
			count[me]++;
		}
		pam.setInitialMedoids(dissimilarityMatrix, initMedoids, k);

		pam.cluster_random();

		HashMap<Integer, Integer> clusterAssign = pam.getClusterAssignmt();

		if (process) {
			System.out.println(clusterAssign);
			printResult(pam.getCluster_result());
			System.out.println("--------------PAM--------------------");
		}

		return clusterAssign;

	}

	@SuppressWarnings("unused")
	private HashMap<Integer, Integer> usePam(int k) {
		MyKMedoidsDis pam = new MyKMedoidsDis();

		pam.initial(dissimilarityMatrix, k);
		pam.cluster_swap();

		HashMap<Integer, Integer> clusterAssign = pam.getClusterAssignmt();

		if (process) {
			System.out.println(clusterAssign);
			printResult(pam.getCluster_result());
			System.out.println("--------------PAM--------------------");
		}
		return clusterAssign;

	}

	private void initialLeartDisMatrix(int size) {
		learntDissimilarity = new HashMap<Integer, HashMap<Integer, Double>>();
		for (int i = 0; i < size; i++) {
			learntDissimilarity.put(i, new HashMap<Integer, Double>());
			for (int j = 0; j <= i; j++) {
				learntDissimilarity.get(i).put(j, (double) 0);
			}
		}
	}

	private void printResult(HashMap<Integer, List<Pattern>> resultMap) {
		System.out.println("\nresult map: ");

		for (int i : resultMap.keySet()) {
			System.out.print("cluster " + i + ": ");
			for (Pattern pat : resultMap.get(i)) {
				System.out.print(pat.getIndex() + " ");
			}
			System.out.println();
		}

	}

	private void takeAverageValue(boolean show) {
		int size = learntDissimilarity.keySet().size();
		System.out.println("learnt dis_matrix :  size = " + size + "\n");
		if (show) {
			System.out.print("\t");
			for (int i = 0; i < size; i++) {
				System.out.print(i + "\t");
			}
			System.out.println();
		}
		for (int i = 0; i < size; i++) {
			if (show)
				System.out.print(i + "\t");
			for (int j = 0; j <= i; j++) {
				HashMap<Integer, Double> row = learntDissimilarity.get(i);
				double newValue = row.get(j) / rounds;
				DecimalFormat df = new DecimalFormat("#.###");
				if (show)
					System.out.print(df.format(newValue) + "\t");
				row.put(j, newValue);

			}
			if (show)
				System.out.println();
		}
		// Util.printWholeMatrix(learntDissimilarity, size);

	}

	private void setValueForLearntDisMatrix(
			HashMap<Integer, Integer> clusterAssignmt) {
		int size = learntDissimilarity.keySet().size();

		for (int i = 0; i < size; i++) {
			for (int j = 0; j <= i; j++) {

				HashMap<Integer, Double> row = learntDissimilarity.get(i);

				double oldValue = row.get(j);
				double newValue = inSameCluster(i, j, clusterAssignmt);

				// System.out.println("i="+i+" j="+j+"  newvalue: "+newValue);
				row.put(j, oldValue + newValue);
			}
		}
	}

	private int inSameCluster(int a, int b,
			HashMap<Integer, Integer> clusterAssignmt) {
		if (clusterAssignmt.get(a) == clusterAssignmt.get(b)) {
			return 0;
		} else {
			return 1;
		}

	}

	public Matrix getLearntDissimilarity() {
		int size = this.dissimilarityMatrix.getColumnDimension();
		return mapToMatrix(size, learntDissimilarity);
	}

	public static Matrix mapToMatrix(int size,
			HashMap<Integer, HashMap<Integer, Double>> map) {
		Matrix matrix = new Matrix(size, size);
		// HashMap<Integer, HashMap<Integer, Double>>
		for (Integer key : map.keySet()) {
			HashMap<Integer, Double> row = map.get(key);
			for (Integer col : row.keySet()) {
				matrix.set((int) key, (int) col, (double) row.get(col));
			}
		}
		return matrix;
	}

}
