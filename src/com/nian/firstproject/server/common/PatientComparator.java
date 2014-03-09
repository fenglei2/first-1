package com.nian.firstproject.server.common;

import java.util.Comparator;

public class PatientComparator implements Comparator<Patient> {

	@Override
	public int compare(Patient p0, Patient p1) {

		if (p0.time - p1.time > 0) {
			return 1;
		} else if (p0.time - p1.time < 0) {
			return -1;
		} else {
			return 0;
		}

	}

}
