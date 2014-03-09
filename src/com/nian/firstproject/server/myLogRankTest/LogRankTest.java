package com.nian.firstproject.server.myLogRankTest;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.nian.firstproject.server.common.Patient;
import com.nian.firstproject.server.common.PatientComparator;

public class LogRankTest {
	/**
	 * The test statistic.
	 */

	public double testStatistic;

	/**
	 * The survival times of the patients in group 1.
	 */

	public double[] time1;

	/**
	 * The censor indicators for the patients in group 1, <br>
	 * censor1[i]=1.0: death; <br>
	 * censor1[i]=0.0: censored.
	 */

	public double[] censor1;

	/**
	 * The survival times of the patients in group 2.
	 */

	public double[] time2;

	/**
	 * The censor indicators for the patients in group 2, <br>
	 * censor2[i]=1.0: death; <br>
	 * censor2[i]=0.0: censored.
	 */

	public double[] censor2;

	public int totalDeath;
	public int O1;// number of death from group 1, O1
	public int O2;

	public LogRankTest() {

	}

	public LogRankTest(double[] time1, double[] censor1, double[] time2,
			double[] censor2) {
		this.time1 = time1;
		this.censor1 = censor1;
		this.time2 = time2;
		this.censor2 = censor2;
		testStatistic();
	}

	public double testStatistic() {

		if (time1.length != censor1.length) {
			throw new IllegalArgumentException(
					"The time vector and the censor vector in group 1 must "
							+ "have the same length.");
		}
		if (time2.length != censor2.length) {
			throw new IllegalArgumentException(
					"The time vector and the censor vector in group 2 must "
							+ "have the same length.");
		}
		testStatistic = 0.0;
		Patient[] group1 = getGroup(time1, censor1, 1);
		Patient[] group2 = getGroup(time2, censor2, 2);

		List<Patient> allPatients = merge2Groups(group1, group2);// return E2
		Map<Double, Event> events = createEvents(allPatients);
		double E2 = calculateE2(events);
		double E1 = totalDeath - E2;

		testStatistic = (O1 - E1) * (O1 - E1) / E1 + (O2 - E2) * (O2 - E2) / E2;

		return testStatistic;
	}

	private List<Patient> merge2Groups(Patient[] group1, Patient[] group2) {
		int start1 = 0, start2 = 0, n1 = group1.length, n2 = group2.length;

		Arrays.sort(group1, new PatientComparator());
		Arrays.sort(group2, new PatientComparator());

		List<Patient> sorted = new ArrayList<Patient>();
		while (start1 < n1 && start2 < n2) {
			if (group1[start1].time < group2[start2].time) {
				sorted.add(group1[start1++]);
			} else {
				sorted.add(group2[start2++]);
			}
		}

		while (start1 < n1) {
			sorted.add(group1[start1++]);
		}
		while (start2 < n2) {
			sorted.add(group2[start2++]);
		}
		return sorted;
	}

	private double calculateE2(Map<Double, Event> allEvents) {
		double e2 = 0;
		for (double key : allEvents.keySet()) {
			Event event = allEvents.get(key);
			e2 += event.E;
		}
		return e2;
	}

	private Map<Double, Event> createEvents(List<Patient> allPatients) {
		// change to unique survival days array
		Map<Double, Event> map = new HashMap<Double, Event>();
		totalDeath = 0;
		O1 = 0;
		O2 = 0;
		for (int i = 0; i < allPatients.size(); i++) {
			Patient patient = allPatients.get(i);
			double key = patient.time;
			double censor = patient.censor;

			Event event;

			totalDeath += censor;
			if (patient.group == 1) {
				O1 += censor;
			} else {
				O2 += censor;
			}

			if (!map.containsKey(key)) {// if survival day exist in map
				event = new Event(key, allPatients.size() - i, censor,
						gruop2AliveNo(i, allPatients), patient.group);

			} else {
				event = map.get(key);
				event.ndeath += censor;
				event.risk = event.ndeath / event.currTotal;
				event.E = event.naliveFromGroup2 * event.risk;
			}
			// System.out.println(event);
			map.put(key, event);// update map
		}

		// System.out.println("o1:"+O1+" o2:"+O2+" totaldeath: "+totalDeath);//ok

		return map;
	}

	private int gruop2AliveNo(int curr, List<Patient> allPatients) {
		int ct = 0;
		for (int i = curr; i < allPatients.size(); i++) {
			Patient p = allPatients.get(i);
			if (p.group == 2) {
				ct++;
			}
		}
		return ct;
	}

	private Patient[] getGroup(double[] time, double[] censor, int no) {
		int n = time.length;
		Patient[] group = new Patient[n];
		for (int i = 0; i < n; i++) {
			Patient p = new Patient(time[i], censor[i], no);
			group[i] = p;

		}
		return group;
	}

	// public static void main(String[] args) {
	//
	//
	// double[] time1={8,12,15,25,37,55,72};
	// double[] censor1={1,1,0,0,1,1,0};
	// double[] time2={1,1,4,5,6,9,9,22};
	// double[] censor2={1,1,1,1,0,0,1,1};
	//
	// LogRankTest log = new LogRankTest(time1,censor1,time2,censor2);
	// double r=log.testStatistic(time1, censor1, time2, censor2);
	// System.out.println(r);
	//
	// Patient p0 = new Patient(8, 1, 1);
	// Patient p1 = new Patient(12, 1, 1);
	// Patient p2 = new Patient(15, 0, 1);
	// Patient p3 = new Patient(25, 0, 1);
	// Patient p4 = new Patient(37, 1, 1);
	// Patient p5 = new Patient(55, 1, 1);
	// Patient p6 = new Patient(72, 0, 1);
	// Patient[] arr1 = new Patient[7];
	// arr1[0] = p0;
	// arr1[1] = p1;
	// arr1[2] = p2;
	// arr1[3] = p3;
	// arr1[4] = p4;
	// arr1[5] = p5;
	// arr1[6] = p6;
	//
	// Patient p10 = new Patient(1, 1, 2);
	// Patient p11 = new Patient(1, 1, 2);
	// Patient p12 = new Patient(4, 1, 2);
	// Patient p13 = new Patient(5, 1, 2);
	// Patient p14 = new Patient(6, 0, 2);
	// Patient p15 = new Patient(9, 0, 2);
	// Patient p16 = new Patient(9, 1, 2);
	// Patient p17 = new Patient(22, 1, 2);
	// Patient[] arr2 = new Patient[8];
	// arr2[0] = p10;
	// arr2[1] = p11;
	// arr2[2] = p12;
	// arr2[3] = p13;
	// arr2[4] = p14;
	// arr2[5] = p15;
	// arr2[6] = p16;
	// arr2[7] = p17;
	//
	//
	//
	// Map<Double, Event> events = log.createEvents(log.merge2Groups(arr1,
	// arr2));
	// System.out.println(events);
	//
	// double E2 = log.calculateE2(events);
	// System.out.println("e2: " + E2);
	//
	// double E1 = log.totalDeath - E2;
	//
	// log.testStatistic = (log.O1 - E1) * (log.O1 - E1) / E1 + (log.O2 - E2) *
	// (log.O2 - E2) / E2;
	//
	// System.out.println("result: " + log.testStatistic);
	//
	// }

}
