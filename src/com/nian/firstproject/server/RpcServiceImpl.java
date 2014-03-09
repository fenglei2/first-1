package com.nian.firstproject.server;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;

import com.google.appengine.api.blobstore.BlobKey;
import com.google.appengine.api.blobstore.BlobstoreInputStream;
import com.google.appengine.api.blobstore.BlobstoreService;
import com.google.appengine.api.blobstore.BlobstoreServiceFactory;
import com.google.gwt.user.server.rpc.RemoteServiceServlet;
import com.nian.firstproject.client.blobs.BlobData;
import com.nian.firstproject.client.blobs.BlobDataFilter;
import com.nian.firstproject.client.blobs.PlotData;
import com.nian.firstproject.client.rpc.RpcService;
//import com.nian.firstproject.server.common.DataModel;
import com.nian.firstproject.server.common.DataModel1;
import com.nian.firstproject.shared.Curve;
//import com.nian.firstproject.shared.Pattern;
//import com.nian.firstproject.shared.Pattern;

/**
 * The server side implementation of the RPC service.
 */
@SuppressWarnings("serial")
public class RpcServiceImpl extends RemoteServiceServlet implements RpcService {

	public static Object[] factors;
	private BlobstoreService blobstoreService = BlobstoreServiceFactory
			.getBlobstoreService();

	/**
	 * init a blobstore url
	 */
	public String getBlobStoreUrl() {
		String url = blobstoreService.createUploadUrl("/upload");
		return url;
	}

	/**
	 * get blob info list
	 * 
	 * @param filter
	 * @return
	 */
	public BlobData[] getBlobs(BlobDataFilter filter) {

		BlobInfoJdo db = new BlobInfoJdo();
		BlobData[] r = db.getBlobs(filter);

		return r;
	}

	@Override
	public String getContent(BlobDataFilter filter) {

		String blobKey = filter.getBlobKey();

		if (blobKey == null) {
			return "";
		}

		StringBuffer buf = new StringBuffer();
		try {
			BlobstoreInputStream is = new BlobstoreInputStream(new BlobKey(
					blobKey));
			InputStreamReader reader = new InputStreamReader(is);
			@SuppressWarnings("resource")
			BufferedReader in = new BufferedReader(reader);
			String line;

			while ((line = in.readLine()) != null) {
				System.out.println(line);
				buf.append(line);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return buf.toString();
	}

	public boolean deleteBlob(BlobDataFilter filter) {

		String blobKey = filter.getBlobKey();

		if (blobKey == null) {
			return false;
		}

		BlobKey blobKeys = new BlobKey(blobKey);
		blobstoreService.delete(blobKeys);

		return true;
	}

	@SuppressWarnings("static-access")
	@Override
	public PlotData getPlotData(BlobDataFilter filter) {

		String blobKey = filter.getBlobKey();

		if (blobKey == null) {
			return null;
		}
		PlotData plotData = null;
		try {
			BlobstoreInputStream is = new BlobstoreInputStream(new BlobKey(
					blobKey));
			InputStreamReader reader = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(reader);
			// int index=0;
			DataModel1 model = DataModel1.getInstance();
			model.start(br, filter);
			plotData = new PlotData();
			plotData.setCurves(convert(model.getCurves()));
			plotData.setMaxSurvivalTime(model.getMaxSurvivalTime());
			plotData.setRoot(model.getRootBranch());

			// plotData.setPatterns(convert1(model.getPatterns()));

		} catch (IOException e) {
			e.printStackTrace();
		}

		return plotData;
	}

	private Curve[] convert(List<Curve> curves) {
		Curve[] arr = new Curve[curves.size()];
		for (int i = 0; i < arr.length; i++) {
			arr[i] = curves.get(i);
		}
		return arr;
	}

	@Override
	public String[] getPopulateData(BlobDataFilter filter) {
		String blobKey = filter.getBlobKey();
		if (blobKey == null) {
			return null;
		}
		String[] factors = null;
		try {
			BlobstoreInputStream is = new BlobstoreInputStream(new BlobKey(
					blobKey));
			InputStreamReader reader = new InputStreamReader(is);
			BufferedReader br = new BufferedReader(reader);

			factors = getFactors(br);

		} catch (IOException e) {
			e.printStackTrace();
		}

		return factors;
	}

	private String[] getFactors(BufferedReader br) {
		if (br == null) {
			throw new IllegalArgumentException("Data stream is empty");
		}
		String[] factors = null;

		try {
			String strLine = br.readLine();
			factors = strLine.split("\t");

		} catch (Exception e) {
			System.out.println("Exception while reading  file: " + e);
		}
		for (String s : factors) {
			System.out.print(s + " | ");
		}
		System.out.println();
		return factors;
	}

}
