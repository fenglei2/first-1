package com.nian.firstproject.client.dbconnection;

import com.google.gwt.core.client.GWT;
import com.google.gwt.user.client.rpc.ServiceDefTarget;

public class DBCInit {

	public static DBConnectionAsync init() {
		DBConnectionAsync rpcService = (DBConnectionAsync) GWT
				.create(DBConnection.class);
		ServiceDefTarget target = (ServiceDefTarget) rpcService;
		String moduleRelativeURL = GWT.getModuleBaseURL() + "MySQLConnection";
		target.setServiceEntryPoint(moduleRelativeURL);
		// rpcService
		return rpcService;
	}

}