package com.nian.firstproject.client.rpc;

//this package collect data from client side to server side 
//to run algorithm
//to call data

import com.google.gwt.core.client.GWT;

public class RpcInit {

	public static RpcServiceAsync init() {
		RpcServiceAsync call = (RpcServiceAsync) GWT.create(RpcService.class);
		return call;
	}
	
}
