package com.nian.firstproject.client;

import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.rpc.AsyncCallback;

import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.HTML;

import com.google.gwt.user.client.ui.PushButton;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.widgetideas.graphics.client.Color;
import com.nian.firstproject.client.blobs.BlobDataFilter;

import com.nian.firstproject.client.blobs.PlotData;

import com.nian.firstproject.client.rpc.RpcInit;
import com.nian.firstproject.client.rpc.RpcServiceAsync;

import com.nian.firstproject.client.widges.DendgWidge;

public class Cutting extends Composite {

	private RpcServiceAsync rpc;
	private VerticalPanel pWidget;
	private TextBox text2;
	private PushButton submitButton2;
	private PushButton submitButton3;
	public static double cuttingHeight;

	public Cutting() {
		initWidget(getTheWidget());
		rpc = RpcInit.init();
	}

	private VerticalPanel getTheWidget() {
		if (pWidget == null) {

			pWidget = new VerticalPanel();

			pWidget.setSpacing(10);

			pWidget.add(new HTML("Select a height to cutting"));// 0
			pWidget.add(cuttingH());// 1

			pWidget.add(new HTML());// 2
			pWidget.add(submitBtn2());// 3

			pWidget.add(new HTML());// 4
			pWidget.add(submitBtn3());// 5

		}
		return pWidget;
	}

	private PushButton submitBtn2() {

		submitButton2 = new PushButton("Cut");
		submitButton2.setEnabled(false);

		submitButton2.addClickHandler(new ClickHandler() {

			@Override
			public void onClick(ClickEvent event) {

				// draw red boxes in dendrogram

				String txt2 = text2.getText();

				if (txt2.isEmpty() || !isDouble(txt2)) {
					Window.alert("Enter a height(0.0-1.0) to cut!");
				}

				Double cuttingHeight = Double.valueOf(text2.getText());

				if (cuttingHeight < 0.0 || cuttingHeight > 1.0)

				{
					Window.alert("height should between 0.0 and 1.0 !");
				}

				showBox(cuttingHeight);

				// two lines should in showBox method

			}

		});

		return submitButton2;
	}

	public boolean isDouble(String input) {
		try {
			Double.parseDouble(input);
			return true;
		} catch (Exception e) {
			return false;
		}
	}

	private void showBox(final Double cuttingHeight) {

		BlobDataFilter filter = new BlobDataFilter();
		filter.setCuttingHeight(cuttingHeight);

		rpc.getPlotData(filter, new AsyncCallback<PlotData>() {

			@Override
			public void onSuccess(PlotData result) {

				DendgWidge.canvas.setLineWidth(0.5);
				DendgWidge.canvas.setStrokeStyle(Color.RED);

				DendgWidge.canvas.saveContext();

				DendgWidge.canvas
						.strokeRect

						(DendgWidge.centerX - 600
								* (cuttingHeight / DendgWidge.MAX_WIDTH),
								(DendgWidge.centerY + 50 * (DendgWidge.upmost - 1 + DendgWidge.DIST_PAD)),
								600 * ((cuttingHeight + 0.05) / DendgWidge.MAX_WIDTH),
								50 * (DendgWidge.downmost - DendgWidge.upmost + 2));

				DendgWidge.canvas
						.moveTo(DendgWidge.centerX - 600
								* (cuttingHeight / DendgWidge.MAX_WIDTH),
								DendgWidge.centerY
										+ 50
										* (DendgWidge.linedownMost + 1 + DendgWidge.DIST_PAD));

				DendgWidge.canvas.lineTo(DendgWidge.centerX - 600
						* (-0.05 / DendgWidge.MAX_WIDTH), DendgWidge.centerY
						+ 50
						* (DendgWidge.linedownMost + 1 + DendgWidge.DIST_PAD));

				DendgWidge.canvas.restoreContext();

				DendgWidge.canvas.stroke();

				submitButton2.setEnabled(true);
				submitButton2.setHTML("Cut");

				submitButton3.setEnabled(true);
				submitButton3.setHTML("Clear cut");

			}

			@Override
			public void onFailure(Throwable caught) {
				System.out.println("getPlotData error: " + caught);

			}

		});

	}

	private TextBox cuttingH() {
		text2 = new TextBox();
		text2.setText("0.7");
		text2.setWidth("100px");
		return text2;
	}

	public VerticalPanel getpWidget() {
		return pWidget;
	}

	private PushButton submitBtn3() {

		submitButton3 = new PushButton("Clear cut");
		submitButton3.setEnabled(false);
		submitButton3.addClickHandler(new ClickHandler() {

			@Override
			public void onClick(ClickEvent event) {

				BlobDataFilter filter = new BlobDataFilter();
				filter.setCuttingHeight(cuttingHeight);

				rpc.getPlotData(filter, new AsyncCallback<PlotData>() {

					@Override
					public void onSuccess(PlotData result) {

						// VerticalPanel charts = CommonWidge.getCharts();
						// charts.clear();

					}

					@Override
					public void onFailure(Throwable caught) {
						System.out.println("getPlotData error: " + caught);

					}
				});

			}

		});

		return submitButton3;

	}

}