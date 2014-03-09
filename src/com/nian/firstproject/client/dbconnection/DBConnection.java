package com.nian.firstproject.client.dbconnection;

import com.google.gwt.user.client.rpc.RemoteService;
import com.nian.firstproject.client.blobs.DataFilter;
import com.nian.firstproject.client.blobs.PlotData;

import com.nian.firstproject.shared.UserInfo;

public interface DBConnection extends RemoteService {

	public String upload(String fpath);

	public UserInfo authenticateUser(String username, String password)
			throws Exception;

	public PlotData getPlotData(DataFilter filter) throws Exception;

}
