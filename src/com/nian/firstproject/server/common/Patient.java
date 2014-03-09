package com.nian.firstproject.server.common;

public class Patient {

	public double time;
	public double censor;
	public int group;

	public Patient(double time, double censor, int group) {
		this.time = time;
		this.censor = censor;
		this.group = group;
	}

	public Patient(double time, double censor) {
		this.time = time;
		this.censor = censor;
	}

	public String toString() {
		return time + "-" + censor + " ";

	}
}
