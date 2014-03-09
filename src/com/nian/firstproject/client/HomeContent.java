package com.nian.firstproject.client;

//website description here

import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Label;

public class HomeContent extends Composite {

	HorizontalPanel horizontalPanel;

	public HomeContent() {
		getTheWidget();
	}

	public HorizontalPanel getTheWidget() {
		if (horizontalPanel == null) {

			horizontalPanel = new HorizontalPanel();
			initWidget(horizontalPanel);

			horizontalPanel.setSpacing(40);
			Label lblNewLabel = new Label(
					"Go to Predict part of our web, "
							+ "with your patient's information, you can predict their survival rate and hazard rate.");
			horizontalPanel.add(lblNewLabel);

			Label lblNewLabel_1 = new Label(
					"Still, go to the Predict part to see "
							+ "the survival curves, hazard curves, dendrogram using our dataset. ");
			horizontalPanel.add(lblNewLabel_1);

			Label lblNewLabel_2 = new Label(
					"Also, go to the Filter part you can upload you own"
							+ " dataset to generate your curves and dendrogram.");
			horizontalPanel.add(lblNewLabel_2);
		}
		return horizontalPanel;

	}

}
