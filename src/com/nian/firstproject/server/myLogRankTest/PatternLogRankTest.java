package com.nian.firstproject.server.myLogRankTest;

import java.util.List;

import com.nian.firstproject.shared.Pattern;

public class PatternLogRankTest extends LogRankTest {

	public PatternLogRankTest(Pattern pattern1, Pattern pattern2) {

		censor1 = getDoubleArray(pattern1.getCensor());
		censor2 = getDoubleArray(pattern2.getCensor());
		time1 = getDoubleArray(pattern1.getSurvival());
		time2 = getDoubleArray(pattern2.getSurvival());
		testStatistic();

	}

	private double[] getDoubleArray(List<Double> list) {
		Object[] arr = list.toArray();
		double[] time = new double[arr.length];
		for (int i = 0; i < time.length; i++) {
			time[i] = Double.parseDouble(arr[i].toString());
		}
		return time;
	}

	public PatternLogRankTest(double[] time1, double[] censor1, double[] time2,
			double[] censor2) {
		super(time1, censor1, time2, censor2);
		// TODO Auto-generated constructor stub
	}

}
