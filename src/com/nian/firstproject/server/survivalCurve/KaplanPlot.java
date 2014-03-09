package com.nian.firstproject.server.survivalCurve;

import java.util.Arrays;
import java.util.HashMap;
import java.util.SortedSet;
import java.util.TreeSet;

public class KaplanPlot {
	public double[] time;// x
	public double[] censor;
	public double[] probabilityOfSurvival;// y
	private boolean process = false;

	public KaplanPlot() {

	}

	public KaplanPlot(Double[] time2, Double[] censor2) {
		time = new double[time2.length];
		censor = new double[censor2.length];
		// censor = new double[time2.length];
		HashMap<Double, Integer> map = new HashMap<Double, Integer>();

		for (int i = 0; i < time2.length; i++) {
			double tm = time2[i];
			while (map.keySet().contains(tm)) {
				tm += 0.0000000001;
			}
			map.put(tm, i);
			time[i] = tm;
		}
		Arrays.sort(this.time);

		SortedSet<Double> keys = new TreeSet<Double>(map.keySet());
		int i = 0;
		for (Double key : keys) {
			int index = map.get(key);
			censor[i++] = censor2[index];
		}
		if (process) {
			System.out.println("--------------Estimating------------");
		}
		KaplanMeierEstimate kaplan = new KaplanMeierEstimate(this.time,
				this.censor);
		// if(process)
		// System.out.println("--------------Kaplan estimate done------------");
		probabilityOfSurvival = kaplan.estimate;
	}

	public void plot() {
		System.out.println("---------probability--------------");
		for (double num : probabilityOfSurvival) {
			System.out.println(num);
		}
		System.out.println("---------time--------------");
		for (double num : time) {
			System.out.println(num);
		}

	}

}
