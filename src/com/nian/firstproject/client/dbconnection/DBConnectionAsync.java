package com.nian.firstproject.client.dbconnection;

import com.google.gwt.user.client.rpc.AsyncCallback;
import com.nian.firstproject.client.blobs.DataFilter;
import com.nian.firstproject.client.blobs.PlotData;
import com.nian.firstproject.shared.UserInfo;

public interface DBConnectionAsync {

	public void upload(String fpath, AsyncCallback<String> callback);

	public void authenticateUser(String username, String password,
			AsyncCallback<UserInfo> callback);

	public void getPlotData(DataFilter filter, AsyncCallback<PlotData> callback);

}
