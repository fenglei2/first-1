package com.nian.firstproject.client.blobs;

import com.google.gwt.user.client.rpc.IsSerializable;

public class BlobDataFilter implements IsSerializable {

	private long start = 0;
	private long limit = 50;

	private String key;
	private int timeCol;
	private int censorCol;
	private int minPatientsNo;
	private int[] selectItems;

	private double cuttingHeight;

	public void setLimit(long start, long limit) {
		this.start = start;
		this.limit = limit;
	}

	public double getCuttingHeight() {
		return cuttingHeight;
	}

	public void setCuttingHeight(double cuttingHeight) {
		this.cuttingHeight = cuttingHeight;
	}

	public int getTimeCol() {
		return timeCol;
	}

	public void setTimeCol(int timeCol) {
		this.timeCol = timeCol;
	}

	public int getCensorCol() {
		return censorCol;
	}

	public void setCensorCol(int censorCol) {
		this.censorCol = censorCol;
	}

	public int getMinPatientsNo() {
		return minPatientsNo;
	}

	public void setMinPatientsNo(int minPatientsNo) {
		this.minPatientsNo = minPatientsNo;
	}

	/**
	 * get start of range
	 * 
	 * @return
	 */
	public long getRangeStart() {
		return start;
	}

	/**
	 * get end of range finish - this is a work around for the offset
	 * 
	 * @return
	 */
	public long getRangeFinish() {
		long finish = start + limit;
		return finish;
	}

	public void setBlobKey(String key) {
		this.key = key;
	}

	public String getBlobKey() {
		return key;
	}

	public void setSelectItems(int[] selectItems) {

		this.selectItems = selectItems;
	}

	public int[] getSelectItems() {
		return selectItems;
	}

}
