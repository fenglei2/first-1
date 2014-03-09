package com.nian.firstproject.server.myLogRankTest;

public class Event {

	// unique among sorted array, current event's days
	public double time;

	// current alive, includes current survival time's death, exclude the
	// paitents already counted
	// with survival days less than current event's days
	public double currTotal;

	// current event's death number, eg: how many 1s with survival days=time
	public double ndeath;

	// risk of death, death
	public double risk;

	// number of alive patients from group2
	public double naliveFromGroup2;

	public int group;
	public double noPatients;

	public double E;

	public double proportionSurvival;

	// public double tenNo;
	public Event(double time, double currtotal, double ndeath,
			double naliveFromGroup2, int group) {
		this.time = time;
		this.currTotal = currtotal;
		this.ndeath = ndeath;
		this.naliveFromGroup2 = naliveFromGroup2;
		this.group = group;
		this.risk = ndeath / currTotal;
		this.E = naliveFromGroup2 * risk;
		// this.

	}

	public double getCurrTotal(int group) {

		return currTotal;
	}

	public Event(double time, double currtotal, double ndeath) {
		this.time = time;
		this.currTotal = currtotal;
		this.ndeath = ndeath;
		this.proportionSurvival = (currTotal - ndeath) / currTotal;

		// System.out.println("dd: "+currTotal+" "+ndeath+" "+proportionSurvival);
	}

	public String toString() {
		return currTotal + " " + risk + " " + ndeath + " " + naliveFromGroup2
				+ " " + group + "\n";
	}

}
