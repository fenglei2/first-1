package com.nian.firstproject.shared;

import com.google.gwt.user.client.rpc.IsSerializable;

public class Curve implements IsSerializable {
	public double[] probabilityOfSurvival;// y
	public double[] time;// x
	String label;

	public Curve(double[] prob, double[] time, String label) {
		this.probabilityOfSurvival = prob;
		this.time = time;
		this.label = label;
		// this.tenNo=tenNo;

	}

	public Curve(String label) {
		this.label = label;
	}

	public double[] getProbabilityOfSurvival() {
		return probabilityOfSurvival;
	}

	public void setProbabilityOfSurvival(double[] probabilityOfSurvival) {
		this.probabilityOfSurvival = probabilityOfSurvival;
	}

	public double[] getTime() {
		return time;
	}

	public void setTime(double[] time) {
		this.time = time;
	}

	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}

	public Curve() {

	}

}
