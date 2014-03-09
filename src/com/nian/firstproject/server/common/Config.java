package com.nian.firstproject.server.common;

public enum Config {

	//DATA_PATH("Breast Cancer Dataset.csv"),// raw dataset input
	READ_TEST_RESULT(false),// read log rank cs.umbc result from file
	READ_LEARNT_MATRIX_RESULT(false),// read learnt matrix from file
	
	//show the matrix result in console
	SHOW_TEST_MATRIX(true),
	SHOW_LEARNT_MATRIX(true),
	SHOW_CLUSTERING_RESULT(true),
	SHOW_DENDROGRAM(true),
	SHOW_CORRELATION(true),
	SHOW_CLUSTERSIZE(true),
	
	//TEST_PATH("chen_logrankcs.umbc.txt"),// if READ_TEST_RESULT = true, read log rank cs.umbc result from file directly, without running log rank cs.umbc code
	//TEST_OUTPUT_PATH("cs.umbc_output.txt"),// log rank cs.umbc result output to file
	//LEARNT_PATH("learnt.txt"),// if READ_LEARNT_MATRIX_RESULT = true, learnt dissimilarity result from file directly without running code
	//LEARNT_OUTPUT_PATH("learnt_output.txt"),//learnt dissimilarity result output to file
	
	MIN_SAMPLE(100),// the minimum number in a combination
	ROUND_TIMES(1000); //number of runs : PAM
	

	String name;
	boolean bool;
	int num;

	Config(int num) {
		this.num = num;
	}

	Config(boolean bool) {
		this.bool = bool;
	}

	Config(String name) {
		this.name = name;
	}

}
