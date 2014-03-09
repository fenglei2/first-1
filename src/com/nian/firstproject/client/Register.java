package com.nian.firstproject.client;

import java.sql.DriverManager;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.RootPanel;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Label;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;

public class Register extends Composite {
	public Register() {

		VerticalPanel verticalPanel = new VerticalPanel();
		initWidget(verticalPanel);

		HorizontalPanel horizontalPanel_6 = new HorizontalPanel();
		verticalPanel.add(horizontalPanel_6);

		HorizontalPanel horizontalPanel = new HorizontalPanel();
		verticalPanel.add(horizontalPanel);

		Label lblBecomeAMember = new Label("Become a member?");
		horizontalPanel.add(lblBecomeAMember);

		HorizontalPanel horizontalPanel_1 = new HorizontalPanel();
		verticalPanel.add(horizontalPanel_1);

		HorizontalPanel horizontalPanel_2 = new HorizontalPanel();
		verticalPanel.add(horizontalPanel_2);

		Label lblUsername = new Label("Create Username");
		horizontalPanel_2.add(lblUsername);

		final TextBox textBox_1 = new TextBox();
		horizontalPanel_2.add(textBox_1);

		HorizontalPanel horizontalPanel_3 = new HorizontalPanel();
		verticalPanel.add(horizontalPanel_3);

		Label lblCreatePassword = new Label("Create Password:");
		horizontalPanel_3.add(lblCreatePassword);

		final TextBox textBox_2 = new TextBox();
		horizontalPanel_3.add(textBox_2);

		HorizontalPanel horizontalPanel_4 = new HorizontalPanel();
		verticalPanel.add(horizontalPanel_4);

		Label lblConfirmPassword = new Label("Confirm Password: ");
		horizontalPanel_4.add(lblConfirmPassword);

		TextBox textBox_3 = new TextBox();
		horizontalPanel_4.add(textBox_3);

		HorizontalPanel horizontalPanel_5 = new HorizontalPanel();
		verticalPanel.add(horizontalPanel_5);

		Button btnRegister = new Button("Register");
		horizontalPanel_5.add(btnRegister);
		btnRegister.addClickHandler(new ClickHandler() {

			@Override
			public void onClick(ClickEvent event) {

				String username = textBox_1.getText();
				String password = textBox_2.getText();

				String sql = "insert into userinfo(username,password) values(textBox_1.getText(),textBox_2.getText())";
				try {

					java.sql.Connection conn = DriverManager
							.getConnection(
									"jdbc:google:rdbms://Local instance MySQL55/prognosticsys",
									"root", "root");
					java.sql.PreparedStatement stmt = conn
							.prepareStatement(sql);
					stmt.setString(1, username);
					stmt.setString(2, password);
					stmt.executeUpdate();
				} catch (Exception e1) {
					e1.printStackTrace();
				}

				Login login = new Login();
				RootPanel.get().getWidget(0).setVisible(false);
				RootPanel.get().add(login);

			}

		});

	}

}
