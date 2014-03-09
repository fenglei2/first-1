package com.nian.firstproject.shared;

import java.util.ArrayList;
import java.util.List;

import com.google.gwt.user.client.rpc.IsSerializable;

/**
 * Pattern class includes all information of a combination
 */

public class Pattern implements IsSerializable {

	protected double[] features;// N-dimension input data
	private String label;// the label of one combination
	private int clusterNo; // 0,1...k-1 (k is the number of clusters), belongs
							// to which clusters???how many clusters???
	private List<Double> survival = new ArrayList<Double>();// survival time
															// list
	private List<Double> censor = new ArrayList<Double>();// death:1 or alive:0

	private int index;// combination id

	public Pattern(Pattern p) {
		this.features = p.features;
		this.label = p.label;
		this.clusterNo = p.clusterNo;
		this.index = p.index;
	}

	public Pattern(int index, int clusterNo) {
		this.index = index;
		this.clusterNo = clusterNo;
	}

	public Pattern(int index) {
		this.index = index;
	}

	public Pattern(String label) {
		this.label = label;
	}

	public int getClusterNo() {
		return clusterNo;
	}

	public void setClusterNo(int clusterNo) {
		this.clusterNo = clusterNo;
	}

	public int getIndex() {
		return index;
	}

	public void setIndex(int index) {
		this.index = index;
	}

	public List<Double> getCensor() {
		return censor;
	}

	public void setCensor(Double censor) {
		this.censor.add(censor);
	}

	public List<Double> getSurvival() {
		return survival;
	}

	public void setSurvival(Double survival) {
		this.survival.add(survival);
	}

	public Pattern(double[] features, String label) {
		this.features = features;
		this.label = label;
	}

	public Pattern(double[] features, int index) {
		this.features = features;
		this.index = index;
	}

	/**
	 * Get the feature vector of the pattern
	 * 
	 * @return the <b>row feature vector</b> of the pattern
	 */

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	public double[] getFeatures() {
		return features;
	}

	public void setFeatures(double[] features) {
		this.features = features;
	}
}
