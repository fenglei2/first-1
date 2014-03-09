package com.nian.firstproject.client;

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.event.logical.shared.SelectionEvent;
import com.google.gwt.event.logical.shared.SelectionHandler;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.TabPanel;
import com.nian.firstproject.client.widges.CommonWidge;

/**
 * Entry point classes define <code>onModuleLoad()</code>.
 */
public class HomePanel implements EntryPoint {

	// static TabPanel tb;

	public void onModuleLoad() {

		TabPanel tabPanel = new TabPanel();
		// 1

		tabPanel.add(new HomeContent(), "EMOH");// 0
		// tabPanel.add(new UploadC(), "Filter");// 1

		tabPanel.add(new GetCurveRoll(), "Predict");// 1
		tabPanel.add(new FileUploader(), "Filter");// 2

		tabPanel.add(new Login(), "nigol");// 3
		tabPanel.setStyleName("customizedTabPanel");
		tabPanel.setVisible(true);

		tabPanel.selectTab(0);
		tabPanel.addSelectionHandler(new SelectionHandler<Integer>() {

			@Override
			public void onSelection(SelectionEvent<Integer> event) {
				if (event.getSelectedItem() == 1) {

					CommonWidge.getCurveRoll();

				}

				else if (event.getSelectedItem() == 2) {
					CommonWidge.getFileUploader();

				}

				else if (event.getSelectedItem() == 3) {

					CommonWidge.getLogin();

				}

				else if (event.getSelectedItem() == 0) {

					CommonWidge.getHomeContent();

				}

			}
		});
		RootPanel.get().add(tabPanel);

	}

}
