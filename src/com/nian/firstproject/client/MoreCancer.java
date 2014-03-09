package com.nian.firstproject.client;

import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.TextBox;
//import com.nian.firstproject.client.dbconnection.DBCInit;
//import com.nian.firstproject.client.dbconnection.DBConnectionAsync;
//import com.nian.firstproject.shared.BreastDataset;
import com.google.gwt.user.client.ui.HTML;
import com.google.gwt.user.client.ui.PushButton;

public class MoreCancer extends Composite {

	private VerticalPanel pWidget;
	private HorizontalPanel hPanel;
	private TextBox datasetBox;
	private PushButton factorButton;
	private TextBox factorBox;
	private PushButton levelButton;
	private TextBox levelBox;
	private TextBox levelRangeBox;
	private PushButton submitButton;

	// private DBConnectionAsync a1;

	public MoreCancer() {
		initWidget(getTheWidget());
		// a1 = DBCInit.init();
	}

	public VerticalPanel getTheWidget() {
		if (pWidget == null) {

			pWidget = new VerticalPanel();

			pWidget.setSpacing(10);

			pWidget.add(submitButton());

			hPanel = new HorizontalPanel();
			hPanel.add(new HTML("Name of your dataset:"));// 4
			hPanel.add(datasetBox());// 5
			pWidget.add(hPanel);

			pWidget.add(factorButton());// 5

		}
		return pWidget;
	}

	private TextBox datasetBox() {
		datasetBox = new TextBox();

		// ??
		datasetBox.setText("");
		datasetBox.setWidth("100px");
		return datasetBox;
	}

	private TextBox factorBox() {
		factorBox = new TextBox();

		// ??
		factorBox.setText("");
		factorBox.setWidth("100px");
		return factorBox;
	}

	private TextBox levelBox() {
		levelBox = new TextBox();

		// ??
		levelBox.setText("");
		levelBox.setWidth("100px");
		return levelBox;
	}

	private TextBox levelRangeBox() {
		levelRangeBox = new TextBox();

		// ??
		levelRangeBox.setText("");
		levelRangeBox.setWidth("100px");
		return levelRangeBox;
	}

	private PushButton factorButton() {

		factorButton = new PushButton("Add Factor");
		factorButton.setEnabled(true);

		factorButton.addClickHandler(new ClickHandler() {

			@Override
			public void onClick(ClickEvent event) {

				Window.alert("onclick");

				hPanel = new HorizontalPanel();
				hPanel.add(new HTML("Factor name:"));// 4
				hPanel.add(factorBox());// 5
				pWidget.add(hPanel);
				pWidget.add(levelButton());
			}

		});

		return factorButton;
	}

	private PushButton levelButton() {

		levelButton = new PushButton("Add Level");
		levelButton.setEnabled(true);

		levelButton.addClickHandler(new ClickHandler() {

			@Override
			public void onClick(ClickEvent event) {

				Window.alert("onclick");

				hPanel = new HorizontalPanel();
				hPanel.add(new HTML("Level name:"));// 4
				hPanel.add(levelBox());// 5
				pWidget.add(hPanel);
				hPanel = new HorizontalPanel();
				hPanel.add(new HTML("from:"));// 4
				hPanel.add(levelRangeBox());// 5
				hPanel.add(new HTML("to:"));// 4
				hPanel.add(levelRangeBox());// 5
				pWidget.add(hPanel);

			}

		});

		return levelButton;
	}

	public VerticalPanel getpWidget() {
		return pWidget;
	}

	private PushButton submitButton() {

		submitButton = new PushButton("Submit");
		submitButton.setEnabled(true);

		submitButton.addClickHandler(new ClickHandler() {

			@Override
			public void onClick(ClickEvent event) {

			}

		});

		return submitButton;

	}

}