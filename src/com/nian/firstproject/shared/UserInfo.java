package com.nian.firstproject.shared;

import com.google.gwt.user.client.rpc.IsSerializable;

// not used now
//if we need login, we may use it

/**
 * A class for the user database table
 */
public class UserInfo implements IsSerializable {
	private int id;

	private String username;
	private String password;

	@SuppressWarnings("unused")
	private UserInfo() {
		// just here because GWT wants it.
	}

	/**
	 * A user of the system
	 * 
	 * @param id
	 * @param username
	 * @param password
	 */
	public UserInfo(int id, String username, String password) {
		this.setId(id);
		this.setUsername(username);
		this.setPassword(password);

	}

	public UserInfo(String username, String password) {

		this.setUsername(username);
		this.setPassword(password);

	}

	/**
	 * @param id
	 *            the id to set
	 */
	public void setId(int id) {
		this.id = id;
	}

	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}

	public void setUsername(String username) {
		this.username = username;
	}

	/**
	 * @return the id
	 */
	public String getUsername() {
		return username;
	}

	public void setPassword(String password) {
		this.password = password;
	}

	/**
	 * @return the id
	 */
	public String getPassword() {
		return password;
	}

}