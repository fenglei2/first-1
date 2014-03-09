package com.nian.firstproject.server.common;

//not used now, for mysql testing, success; please ignore

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import com.nian.firstproject.shared.Branch;
import com.nian.firstproject.shared.Cluster;
import com.nian.firstproject.shared.Curve;
import com.nian.firstproject.shared.Pattern;
import com.nian.firstproject.client.blobs.DataFilter;
import com.nian.firstproject.server.MySQLConnection;
import com.nian.firstproject.server.PAM.LearntDisMatrix;
import com.nian.firstproject.server.myHierarchical.HierarchicalCluster;
import com.nian.firstproject.server.myHierarchical.dendrogram.Dendrogram;
import com.nian.firstproject.server.myLogRankTest.LogRankTestMatrix;
import com.nian.firstproject.server.survivalCurve.KaplanPlot;
import Jama.Matrix;

public class DataModel {

	private static volatile DataModel instance = null;
	private static double maxSurvivalTime;

	// number of combinations
	public static int SIZE;

	List<Pattern> dataCombinations;
	Map<String, Pattern> allData;
	double cluster[];
	double prob[];

	String[] featureStr;
	Matrix logRankTestMatrix;// initial dissimilarity matrix
	Matrix learntMatrix;// learnt dissimilarity matrix
	List<List<Cluster>> allClusters;// all hierarchical clusters
	public static KaplanPlot kap;

	private Cluster root;
	List<Curve> curves;
	public Branch rootBranch;

	public static DataModel getInstance() {
		if (instance == null) {
			synchronized (DataModel.class) {
				if (instance == null) {
					instance = new DataModel();
				}
			}
		}
		return instance;
	}

	private DataModel() {
	}

	public void start(DataFilter filter) {
		// read data file
		dataPreprocess(filter);
		logRankTest();
		// calculate learnt dissimilarity
		standardizeRound();
		hierarchicalCluster();
		dendrogram();
		KaplanMeierCurve();
	}

	private void dendrogram() {
		System.out.println("\n =========  Dendrogram ==========");
		Dendrogram den = new Dendrogram(root);
		rootBranch = den.getRootBranch();
		if (Config.SHOW_DENDROGRAM.bool) {
			rootBranch.print();
		}
		System.out.println();
	}

	public void KaplanMeierCurve() {
		List<Cluster> partition = allClusters.get(0);

		List<Double[]> times = new ArrayList<Double[]>();
		List<Double[]> censors = new ArrayList<Double[]>();

		// List<Double> timeList = new ArrayList<Double>();

		int ct = 0;

		for (Cluster c : partition) {

			List<Integer> patternIndexes = c.getPatternIndexes();
			List<Double> timeList = new ArrayList<Double>();
			List<Double> censorList = new ArrayList<Double>();
			for (int patternIdx : patternIndexes) {

				List<Double> time = getPattern(patternIdx).getSurvival();
				timeList.addAll(time);
				List<Double> censor = getPattern(patternIdx).getCensor();
				censorList.addAll(censor);
			}

			// timeList[ct].size();

			if (Config.SHOW_CLUSTERSIZE.bool) {

				System.out.println("cluster" + ct++ + " size: "
						+ timeList.size());

			}

			times.add(doubleArr(timeList));
			censors.add(doubleArr(censorList));
		}

		curves = new ArrayList<Curve>();

		System.out
				.println("==============  Kaplan estimating cluster ==============");
		for (int i = 0; i < times.size(); i++) {
			Double[] time = times.get(i);// one cluster, every patient's time
			kap = new KaplanPlot(time, censors.get(i));
			curves.add(new Curve(kap.probabilityOfSurvival, kap.time,
					dataCombinations.get(i).getLabel()));

			setMaxSurvivalTime(Math.max(maxSurvivalTime,
					kap.time[time.length - 1]));

		}

		System.out.println("maxtime: " + maxSurvivalTime);

		// GraphingData.plot(list);
	}

	public void hierarchicalCluster() {
		System.out
				.println("\n=========   Hierarchical clustering  =========\n");
		HierarchicalCluster hi = new HierarchicalCluster(learntMatrix,
				Config.SHOW_CLUSTERING_RESULT.bool);// show
		root = hi.partition().get(0);
		allClusters = hi.getAllClusters();
	}

	public void standardizeRound() {

		if (!Config.READ_LEARNT_MATRIX_RESULT.bool) {
			System.out
					.println("\n=========== Run Standardize rounds ========================\n");
			LearntDisMatrix learnt = new LearntDisMatrix(Config.ROUND_TIMES.num);
			learnt.run(Config.SHOW_LEARNT_MATRIX.bool);
			learntMatrix = learnt.getLearntDissimilarity();
		} else {
			// System.out.println("\n======= Read Standardize rounds result from file =======\n");
			// learntMatrix = disMatrixFromFile(Config.LEARNT_PATH.name, SIZE,
			// Config.CLUSTER_ALGORITHM.alg);
		}
	}

	public void logRankTest() {
		if (dataCombinations == null || dataCombinations.size() == 0) {
			throw new IllegalArgumentException("Combination is empty");
		}

		if (!Config.READ_LEARNT_MATRIX_RESULT.bool) {
			if (!Config.READ_TEST_RESULT.bool) {
				System.out
						.println("\n=======Run Log rank cs.umbc============\n");
				LogRankTestMatrix logRankTest = new LogRankTestMatrix();
				logRankTest.calTestStaticsMatrix();
				logRankTestMatrix = logRankTest.getMatrix();
				// Util.writeResult2File(logRankTestMatrix,
				// Config.TEST_OUTPUT_PATH.name);
			} else {
				// System.out.println("\n=======Get Log rank test result from file============\n");
				// logRankTestMatrix =
				// Util.disMatrixFromFile(Config.TEST_PATH.name, SIZE,
				// Config.CLUSTER_ALGORITHM.alg);
			}
		}

		if (Config.SHOW_TEST_MATRIX.bool) {
			logRankTestMatrix.print(8, 3);
		}
	}

	private void dataPreprocess(DataFilter filter) {
		System.out.println("\n==============  Pre-processing  ==============");

		int minPatientsNo = filter.getMinPatientsNo();

		String[] factorStr = filter.getfactorStr();

		if (factorStr == null) {
			throw new IllegalArgumentException("Data stream is empty");
		}
		// all combinations, including the combination with small number of
		// patients
		// key: combination
		// value: double[]
		allData = new HashMap<String, Pattern>();

		// after remove the combinations with small number of patients
		dataCombinations = new ArrayList<Pattern>();
		try {

			// if we do not read file from top to bottom, we select certain
			// columns
			// very slow in local testing environment,
			for (int j = 0; j < MySQLConnection.getreturntime().size(); j++) {

				StringBuffer combinationLabel = new StringBuffer();

				int len = factorStr.length;
				double[] features = new double[len];

				featureStr = new String[len];

				for (int i = 0; i < factorStr.length; i++) {

					features[i] = MySQLConnection.getReturnfactors().get(i)
							.get(j);

					featureStr[i] = translate(factorStr, i, features);

					combinationLabel.append(featureStr[i]);

				}// for i

				String label = combinationLabel.toString();

				if (allData.get(label) == null) {
					Pattern pattern = new Pattern(features, label);

					allData.put(label, pattern);

				}

				Pattern pattern = allData.get(label);

				pattern.setSurvival(MySQLConnection.getreturntime().get(j));
				pattern.setCensor(MySQLConnection.getreturndelta().get(j));

				if (allData.get(label) == null) {
					Pattern pattern3 = new Pattern(label);
					allData.put(label, pattern3);
				}

			}// while

			System.out.println("Total dataset: "
					+ MySQLConnection.getreturntime().size());// good
			System.out.println("Combination size: " + allData.size());

			// Remove the combinations with small number of patients
			Iterator<String> classItr = allData.keySet().iterator();
			int index = 0;
			while (classItr.hasNext()) {

				String com = classItr.next();
				Pattern pat = allData.get(com);

				if (pat.getSurvival().size() >= minPatientsNo) {
					pat.setIndex(index++);

					dataCombinations.add(pat);

				}

			}// while

			SIZE = dataCombinations.size();

			System.out.println("After removing small sampleSet(<"
					+ minPatientsNo + "): " + SIZE);

		}

		catch (Exception e) {
			System.out.println("Exception while reading csv file: " + e);
		}

	}

	public List<Pattern> getDataCombinations() {
		return dataCombinations;
	}

	public Matrix getLearntMatrix() {
		return learntMatrix;
	}

	public Matrix getLogRankTestMatrix() {
		return logRankTestMatrix;
	}

	public Map<String, Pattern> getAllData() {
		return allData;
	}

	public static int getSIZE() {
		return SIZE;
	}

	public double getMaxSurvivalTime() {
		return maxSurvivalTime;
	}

	public static void setMaxSurvivalTime(double maxSurvivalTime) {
		DataModel.maxSurvivalTime = maxSurvivalTime;
	}

	public List<Curve> getCurves() {
		return curves;
	}

	public void setCurves(List<Curve> curves) {
		this.curves = curves;
	}

	public Branch getRootBranch() {
		return rootBranch;
	}

	public void setRootBranch(Branch rootBranch) {
		this.rootBranch = rootBranch;
	}

	public Pattern getPattern(int index) {
		return dataCombinations.get(index);
	}

	public static Double[] doubleArr(List<Double> survival) {

		Double[] arr = new Double[survival.size()];
		for (int i = 0; i < arr.length; i++) {
			arr[i] = survival.get(i);
		}
		return arr;
	}

	private String translate(String[] factorStr, int i, double[] features) {

		if ("Grade".equalsIgnoreCase(factorStr[i])) {
			featureStr[i] = "G" + (int) (features[i]);
		}

		// String[] featureStr;
		else if ("size".equalsIgnoreCase(factorStr[i])) {
			featureStr[i] = "T" + (int) (features[i]);
		}

		else if ("nodes".equalsIgnoreCase(factorStr[i])) {
			featureStr[i] = "N" + (int) (features[i] - 1);
		}

		else if ("ER".equalsIgnoreCase(factorStr[i])) {

			if ((int) features[i] == 1)

			{
				featureStr[i] = "ER+";
			}

			else if ((int) features[i] == 2)

			{
				featureStr[i] = "ER-";
			}

		}

		else if ("agedx".equalsIgnoreCase(factorStr[i])) {
			featureStr[i] = "age" + (int) (features[i]);
		}

		else if ("race".equalsIgnoreCase(factorStr[i])) {
			featureStr[i] = "race" + (int) (features[i]);
		}

		else

		{
			featureStr[i] = "#" + (int) features[i];
		}

		return featureStr[i];

	}

}