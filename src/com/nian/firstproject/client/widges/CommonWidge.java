package com.nian.firstproject.client.widges;

import com.nian.firstproject.client.FileUploader;
import com.nian.firstproject.client.GetCurveRoll;
import com.nian.firstproject.client.HomeContent;
import com.nian.firstproject.client.Login;
import com.nian.firstproject.client.MoreCancer;
import com.nian.firstproject.client.UploadC;
import com.google.gwt.user.client.ui.DockPanel;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.TabPanel;
import com.google.gwt.user.client.ui.VerticalPanel;

public class CommonWidge {
	public static VerticalPanel getCurveRoll() {

		TabPanel tabPanel = (TabPanel) RootPanel.get().getWidget(0);
		GetCurveRoll gcr = (GetCurveRoll) tabPanel.getWidget(1);

		return gcr.getTheWidget();
	}

	public static HorizontalPanel getHomeContent() {

		TabPanel tabPanel = (TabPanel) RootPanel.get().getWidget(0);

		HomeContent mc = (HomeContent) tabPanel.getWidget(0);

		return mc.getTheWidget();

	}

	public static VerticalPanel getMoreCancer() {

		TabPanel tabPanel = (TabPanel) RootPanel.get().getWidget(0);

		MoreCancer mc = (MoreCancer) tabPanel.getWidget(2);

		return mc.getTheWidget();
	}

	public static VerticalPanel getUploadC() {

		TabPanel tabPanel = (TabPanel) RootPanel.get().getWidget(0);

		UploadC uc = (UploadC) tabPanel.getWidget(2);

		return uc.getpWidget();
	}

	public static VerticalPanel getFileUploader() {

		TabPanel tabPanel = (TabPanel) RootPanel.get().getWidget(0);

		FileUploader fu = (FileUploader) tabPanel.getWidget(2);

		return fu.getTheWidget();
	}

	public static DockPanel getLogin() {

		TabPanel tabPanel = (TabPanel) RootPanel.get().getWidget(0);

		Login l = (Login) tabPanel.getWidget(3);

		return l.getMain();
	}

	public static VerticalPanel getCutting() {
		return null;

	}

}
