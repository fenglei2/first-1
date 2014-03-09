package com.nian.firstproject.client.rpc;


import com.google.gwt.user.client.rpc.RemoteService;
import com.google.gwt.user.client.rpc.RemoteServiceRelativePath;
import com.nian.firstproject.client.blobs.BlobData;
import com.nian.firstproject.client.blobs.BlobDataFilter;
import com.nian.firstproject.client.blobs.PlotData;

/**
 * The client side stub for the RPC service.
 */
@RemoteServiceRelativePath("rpc")
public interface RpcService extends RemoteService {
	
	public String getBlobStoreUrl();
	
	public BlobData[] getBlobs(BlobDataFilter filter);
	 
	public boolean deleteBlob(BlobDataFilter filter);
	
	public String getContent(BlobDataFilter filter);

	public PlotData getPlotData(BlobDataFilter filter);

	public String[] getPopulateData(BlobDataFilter filter);
	
	
}
