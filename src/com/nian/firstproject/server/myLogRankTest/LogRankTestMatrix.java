package com.nian.firstproject.server.myLogRankTest;

import java.util.List;

import com.nian.firstproject.shared.Pattern;
import com.nian.firstproject.server.common.DataModel1;

import Jama.Matrix;

public class LogRankTestMatrix {

	private List<Pattern> dataSet;
	private boolean calcuateProcess = true;
	Matrix matrix;

	public LogRankTestMatrix() {
		dataSet = DataModel1.getInstance().getDataCombinations();
	}

	public void calTestStaticsMatrix() {
		int n = dataSet.size();
		double[][] testStaticsMatrix = new double[n][n];
		Pattern pat1 = null;
		Pattern pat2 = null;
		for (int i = 0; i < n; i++) {
			pat1 = dataSet.get(i);
			for (int j = 0; j <= i; j++) {

				pat2 = dataSet.get(j);
				if (calcuateProcess) {
					System.out.println(i + " -- " + j + " size1:"
							+ pat1.getCensor().size() + " size2:"
							+ pat2.getCensor().size());
				}

				testStaticsMatrix[i][j] = compare(pat1, pat2);
			}
		}
		matrix = new Matrix(testStaticsMatrix);
	}

	public double compare(Pattern pattern1, Pattern pattern2) {
		PatternLogRankTest patternLogRanktest = new PatternLogRankTest(
				pattern1, pattern2);
		return Math.abs(patternLogRanktest.testStatistic);
	}

	public Matrix getMatrix() {
		return matrix;
	}

}
