package com.nian.firstproject.shared;

//not used now, 
//if we use mysql database(or google cloud instance), we may need it

import com.google.gwt.user.client.rpc.IsSerializable;

/**
 * A class for the database table
 */
public class BreastDataset implements IsSerializable {
	// private double id;

	private double Grade;
	private double size;
	private double nodes;
	private double ER;
	private double agedx;
	private double race;
	// private double vs;

	private double time;
	private double delta;

	public BreastDataset() {
		// just here because GWT wants it.
	}

	/**
	 * A user of the system
	 * 
	 * @param id
	 * @param username
	 * @param password
	 */
	public BreastDataset(double Grade, double size, double nodes, double ER,
			double agedx, double race) {
		// this.setId(id);
		this.setGrade(Grade);
		this.setSize(size);
		this.setNodes(nodes);
		this.setER(ER);
		this.setAgedx(agedx);
		this.setRace(race);
		// this.setVs(vs);
		// this.setTime(time);
		// this.setDelta(delta);

	}

	public BreastDataset(double time, double delta) {
		// this.setId(id);

		this.setTime(time);
		this.setDelta(delta);

	}

	/**
	 * @param username
	 *            the username to set
	 */
	public void setGrade(double Grade) {
		this.Grade = Grade;
	}

	public double getGrade() {
		return Grade;
	}

	public void setSize(double size) {
		this.size = size;
	}

	public double getSize() {
		return size;
	}

	public void setNodes(double nodes) {
		this.nodes = nodes;
	}

	public double getNodes() {
		return nodes;
	}

	public void setER(double ER) {
		this.ER = ER;
	}

	public double getER() {
		return ER;
	}

	public void setAgedx(double agedx) {
		this.agedx = agedx;
	}

	public double getAgedx() {
		return agedx;
	}

	public void setRace(double race) {
		this.race = race;
	}

	public double getRace() {
		return race;
	}

	public void setTime(double time) {
		this.time = time;
	}

	public double getTime() {
		return time;
	}

	public void setDelta(double delta) {
		this.delta = delta;
	}

	public double getDelta() {
		return delta;
	}

}