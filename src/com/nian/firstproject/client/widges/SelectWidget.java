package com.nian.firstproject.client.widges;

import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.rpc.AsyncCallback;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.PushButton;
import com.google.gwt.user.client.ui.Widget;
import com.nian.firstproject.client.FileUploader;
import com.nian.firstproject.client.GetCurveRoll;
import com.nian.firstproject.client.MoreCancer;
import com.nian.firstproject.client.blobs.BlobDataFilter;
import com.nian.firstproject.client.blobs.PlotData;
import com.nian.firstproject.client.rpc.RpcInit;
import com.nian.firstproject.client.rpc.RpcServiceAsync;
import com.nian.firstproject.shared.Curve;

public class SelectWidget extends Composite implements ClickHandler {

	private RpcServiceAsync rpc;
	private PushButton selectButton = new PushButton("Select");
	private String blobKey;

	public SelectWidget(String key) {
		this.blobKey = key;
		initWidget(selectButton);
		selectButton.addClickHandler(this);

		rpc = RpcInit.init();

	}

	@Override
	public void onClick(ClickEvent event) {
		selectButton.setEnabled(false);
		selectButton.setHTML("Waiting");
		Widget sender = (Widget) event.getSource();
		if (sender == selectButton) {

			GetCurveRoll.pWidget.setStyleName(blobKey);

			populateData();
			
			MoreCancer mc=new MoreCancer();
			FileUploader.pWidget.add(mc);

		}
	}

	private void populateData() {
		BlobDataFilter filter = new BlobDataFilter();
		filter.setBlobKey(blobKey);

		rpc.getPopulateData(filter, new AsyncCallback<String[]>() {

			@Override
			public void onSuccess(String[] result) {
				// result has all factors,
				System.out.print("factor returned: ");

				for (String s : result) {
					System.out.print(s + ",");
				}
				System.out.println();

				selectButton.setEnabled(true);
				selectButton.setHTML("Select");

			}

			@Override
			public void onFailure(Throwable caught) {
				System.out.println("POPULATE DATA ERROR: " + caught);
			}

		});

	}

	@SuppressWarnings("unused")
	private void checkresult(PlotData result) {

		System.out.println("maxtimeeee: " + result.getMaxSurvivalTime());

		Curve[] curves = result.getCurves();
		for (Curve c : curves) {
			System.out.print("label: " + c.getLabel());
			for (int i = 0; i < c.getTime().length / 10; i++) {
				System.out.print(" " + c.getTime()[i] + " "
						+ c.getProbabilityOfSurvival()[i]);
			}

			System.out.println();

		}

	}

}
