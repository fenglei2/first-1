package com.nian.firstproject.client;

//ignore this class, for mysql only, not used now

import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.rpc.AsyncCallback;

import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.FileUpload;
import com.google.gwt.user.client.ui.HTML;

import com.google.gwt.user.client.ui.PushButton;

import com.google.gwt.user.client.ui.VerticalPanel;

import com.nian.firstproject.client.dbconnection.DBCInit;

import com.nian.firstproject.client.dbconnection.DBConnectionAsync;

public class UploadC extends Composite {

	private VerticalPanel pWidget;

	private DBConnectionAsync a1;

	FileUpload fileUpload;
	private PushButton uploadButton;
	private PushButton defineButton;

	public UploadC() {
		initWidget(getTheWidget());
		a1 = DBCInit.init();
	}

	public VerticalPanel getTheWidget() {
		if (pWidget == null) {

			pWidget = new VerticalPanel();

			pWidget.setSpacing(10);

			fileUpload = new FileUpload();// 0
			pWidget.add(fileUpload);// 1

			pWidget.add(new HTML());// 2
			pWidget.add(submitUpload());// 3

			pWidget.add(new HTML());// 4
			pWidget.add(submitDefine());// 5

		}
		return pWidget;
	}

	private PushButton submitUpload() {

		uploadButton = new PushButton("Upload");

		uploadButton.setEnabled(true);

		uploadButton.addClickHandler(new ClickHandler() {

			@Override
			public void onClick(ClickEvent event) {

				Window.alert("onclick");

				String fpath = fileUpload.getFilename();

				a1.upload(fpath, new AsyncCallback<String>() {

					@Override
					public void onSuccess(String result) {
						// TODO Auto-generated method stub
						String f = "yes";
						System.out.println(result);
						// TODO Auto-generated method stub
						if (result.equals(f)) {
							System.out.println("onsuccess");
							Window.alert("File Upload");

						} else {
							Window.alert("Filis not uploaded");
						}

					}

					@Override
					public void onFailure(Throwable caught) {
						// TODO Auto-generated method stub
						Window.alert("Failure");
					}
				});

			}

		});

		return uploadButton;
	}

	public VerticalPanel getpWidget() {
		return pWidget;
	}

	private PushButton submitDefine() {

		defineButton = new PushButton("Define");

		defineButton.setEnabled(true);

		defineButton.addClickHandler(new ClickHandler() {

			@Override
			public void onClick(ClickEvent event) {

				MoreCancer mc = new MoreCancer();
				pWidget.add(mc);

			}

		});

		return defineButton;

	}

}