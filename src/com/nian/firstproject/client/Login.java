package com.nian.firstproject.client;

//ignore this class, for mysql only

import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.rpc.AsyncCallback;
import com.google.gwt.user.client.rpc.ServiceDefTarget;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.DockPanel;
import com.google.gwt.user.client.ui.Label;
import com.google.gwt.user.client.ui.HasHorizontalAlignment;
import com.google.gwt.user.client.ui.Image;
import com.google.gwt.user.client.ui.HasVerticalAlignment;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.TextBox;
import com.google.gwt.core.client.GWT;
import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.ui.Button;
import com.nian.firstproject.client.dbconnection.DBConnection;
import com.nian.firstproject.client.dbconnection.DBConnectionAsync;
import com.nian.firstproject.shared.UserInfo;

public class Login extends Composite {

	private DockPanel dockPanel;

	public Login() {

		initWidget(getMain());
	}

	public DockPanel getMain() {

		if (dockPanel == null) {

			dockPanel = new DockPanel();
			dockPanel.setVerticalAlignment(HasVerticalAlignment.ALIGN_MIDDLE);

			dockPanel.setSize("600", "422px");

			Label lblPrognosticSystemFor = new Label(
					"Please login: do we need this part if users want to store their dataset here for their own use?");
			lblPrognosticSystemFor.setStyleName("h2");
			lblPrognosticSystemFor
					.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_RIGHT);
			dockPanel.add(lblPrognosticSystemFor, DockPanel.NORTH);
			lblPrognosticSystemFor.setSize("500px", "111px");

			Image image = new Image("WEB-INF/classes/images/images1.jpeg");
			dockPanel.add(image, DockPanel.WEST);
			image.setSize("200", "200");
			dockPanel.setCellVerticalAlignment(image,
					HasVerticalAlignment.ALIGN_MIDDLE);

			VerticalPanel verticalPanel = new VerticalPanel();
			dockPanel.add(verticalPanel, DockPanel.EAST);

			HorizontalPanel horizontalPanel_2 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_2);

			Button btnRegisterToBe = new Button("Register?");

			btnRegisterToBe.addClickHandler(new ClickHandler() {

				@Override
				public void onClick(ClickEvent event) {

					// Register r=new Register();

				}

			});
			horizontalPanel_2.add(btnRegisterToBe);

			HorizontalPanel horizontalPanel = new HorizontalPanel();
			verticalPanel.add(horizontalPanel);

			Label lblNewLabel = new Label("Username");
			lblNewLabel.setStyleName("generalword");
			horizontalPanel.add(lblNewLabel);

			final TextBox textBox = new TextBox();
			textBox.setName("username");
			horizontalPanel.add(textBox);

			HorizontalPanel horizontalPanel_1 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_1);

			Label lblNewLabel_1 = new Label("Password");
			lblNewLabel_1.setStyleName("generalword");
			horizontalPanel_1.add(lblNewLabel_1);

			final TextBox textBox_1 = new TextBox();
			textBox_1.setName("password");
			horizontalPanel_1.add(textBox_1);

			Button btnLogin = new Button("Login");
			btnLogin.addClickHandler(new ClickHandler() {

				@Override
				public void onClick(ClickEvent event) {

					DBConnectionAsync rpcService = (DBConnectionAsync) GWT
							.create(DBConnection.class);
					ServiceDefTarget target = (ServiceDefTarget) rpcService;
					String moduleRelativeURL = GWT.getModuleBaseURL()
							+ "MySQLConnection";
					target.setServiceEntryPoint(moduleRelativeURL);
					rpcService.authenticateUser(textBox.getText(),
							textBox_1.getText(), new AsyncCallback<UserInfo>() {
								public void onFailure(Throwable caught) {
									Window.alert("You got to help me. I don't know what to do. I can't make decisions. I'm a president!");
								}

								public void onSuccess(UserInfo result) {
									Window.alert("Hey I'm a user with id "
											+ result.getId());
								}

							});

				}

			});

			verticalPanel.add(btnLogin);

		}
		return dockPanel;

	}

}
