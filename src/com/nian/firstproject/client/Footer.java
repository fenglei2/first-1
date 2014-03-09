package com.nian.firstproject.client;

import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.Label;

public class Footer extends Composite {
	public Footer() {

		Label lblNewLabel = new Label(
				" Copyright &copy; 2013. All right reserved.");
		initWidget(lblNewLabel);
	}

}
