
//test github


package com.nian.firstproject.client;

import com.nian.firstproject.client.blobs.BlobData;
import com.nian.firstproject.client.blobs.BlobDataFilter;
import com.nian.firstproject.client.rpc.RpcInit;
import com.nian.firstproject.client.rpc.RpcServiceAsync;
import com.nian.firstproject.client.widges.DeleteWidget;
import com.nian.firstproject.client.widges.SelectWidget;
import com.google.gwt.event.dom.client.ChangeEvent;
import com.google.gwt.event.dom.client.ChangeHandler;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.rpc.AsyncCallback;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.Grid;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.VerticalPanel;

//The data processing in Google App server

public class Blobs extends Composite {

	private RpcServiceAsync rpc;
	private VerticalPanel pWidget = new VerticalPanel();
	private BlobDataFilter filter = new BlobDataFilter();

	public Blobs() {

		initWidget(pWidget);
		rpc = RpcInit.init();
	}

	public void draw() {
		getBlobsList();
	}

	private void process(BlobData[] blobData) {

		pWidget.clear();

		if (blobData == null) {
			return;
		}

		int r = blobData.length + 1;
		int c = 7;
		Grid grid = new Grid(r, c);
		pWidget.add(grid);

		grid.setWidget(0, 0, new HTML("Filename"));
		grid.setWidget(0, 2, new HTML("Size"));
		grid.setWidget(0, 5, new HTML("Delete"));
		grid.setWidget(0, 6, new HTML("Select"));

		// retrieve data from google server from by key
		for (int i = 0; i < blobData.length; i++) {

			String link = "<a target=\"_blank\" href=\"/serve?blob-key="
					+ blobData[i].getKey() + "\">" + blobData[i].getFilename()
					+ "</a>";
			// na: filename, si: size
			HTML na = new HTML(link);
			HTML si = new HTML(Long.toString(blobData[i].getSize()));
			DeleteWidget dw = new DeleteWidget(blobData[i].getKey());
			SelectWidget sw = new SelectWidget(blobData[i].getKey());

			grid.setWidget(i + 1, 0, na);
			grid.setWidget(i + 1, 2, si);
			grid.setWidget(i + 1, 5, dw);
			grid.setWidget(i + 1, 6, sw);

			dw.addChangeHandler(new ChangeHandler() {
				public void onChange(ChangeEvent event) {
					getBlobsList();
				}
			});
		}

	}

	public void getBlobsList() {

		rpc.getBlobs(filter, new AsyncCallback<BlobData[]>() {
			public void onSuccess(BlobData[] blobData) {
				process(blobData);

			}

			public void onFailure(Throwable caught) {
				Window.alert("couldn't get blob list, rpc failed");
			}
		});

	}

}
