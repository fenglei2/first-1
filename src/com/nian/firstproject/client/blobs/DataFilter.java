package com.nian.firstproject.client.blobs;

//used in mysql connection, same function as BlobDataFilter, 

import com.google.gwt.user.client.rpc.IsSerializable;

public class DataFilter implements IsSerializable {

	private int minPatientsNo;
	private String[] factorStr;

	public int getMinPatientsNo() {
		return minPatientsNo;
	}

	public void setMinPatientsNo(int minPatientsNo) {
		this.minPatientsNo = minPatientsNo;
	}

	public void setfactorStr(String[] factorStr) {

		this.factorStr = factorStr;
	}

	public String[] getfactorStr() {
		return factorStr;
	}

}
