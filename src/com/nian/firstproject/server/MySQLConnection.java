package com.nian.firstproject.server;

//ignore this class, not used now

//(client) use server to query mysql database

import java.sql.Connection;

import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

import java.util.ArrayList;
import java.util.List;

import java.io.File;
import java.io.FileInputStream;

import java.io.IOException;

import com.google.gwt.user.server.rpc.RemoteServiceServlet;
import com.nian.firstproject.client.blobs.DataFilter;
import com.nian.firstproject.client.blobs.PlotData;
import com.nian.firstproject.client.dbconnection.DBConnection;
import com.nian.firstproject.server.common.DataModel;

import com.nian.firstproject.shared.Curve;
import com.nian.firstproject.shared.UserInfo;

public class MySQLConnection extends RemoteServiceServlet implements
		DBConnection {

	private static final long serialVersionUID = 1L;

	private String _url = "jdbc:mysql://localhost:3306/prognosticsys?user=root&password=root";

	static List<List<Double>> returnfactors;

	static List<Double> returnonecolumn;

	static List<Double> returntime;
	static List<Double> returndelta;

	PlotData plotData;

	/**
	 * Constructor
	 */
	public MySQLConnection() {
	}

	/**
	 * Gets the connection for all of our commands
	 * 
	 * @return
	 * @throws Exception
	 */
	private Connection getConnection() throws Exception {

		Class.forName("com.mysql.jdbc.Driver").newInstance();
		Connection conn = DriverManager.getConnection(_url);
		return conn;
	}

	@Override
	public String upload(String fpath) {
		// TODO Auto-generated method stub
		// call();
		String ss = "no";
		int k = 0;
		FileInputStream fis = null;
		PreparedStatement pstmt = null;
		Connection conn = null;

		System.out.print(fpath);
		try {
			conn = getConnection();
			System.out.print("Hello try");
			conn.setAutoCommit(false);
			File file = new File(fpath);
			fis = new FileInputStream(file);
			pstmt = conn
					.prepareStatement("insert into uploadfile(fpath, fbody) values (?,?)");

			pstmt.setString(1, fpath);
			pstmt.setAsciiStream(2, fis, (int) file.length());
			k = pstmt.executeUpdate();
			System.out.print("Bye");
			if (k == 1) {
				ss = "yes";

			}
			conn.commit();
		} catch (Exception e) {
			System.err.print("Error: " + e.getMessage());
			e.printStackTrace();
		} finally {
			try {
				pstmt.close();
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			try {
				fis.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			try {
				conn.close();
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		return ss;
	}

	public UserInfo authenticateUser(String username, String password)
			throws Exception {
		UserInfo returnuser = null;
		Connection conn = null;
		PreparedStatement pstmt = null;
		ResultSet result = null;

		try {
			conn = getConnection();
			pstmt = conn
					.prepareStatement("SELECT * FROM userinfo WHERE username = ? AND password = ?");
			pstmt.setString(1, username);
			pstmt.setString(2, password);
			result = pstmt.executeQuery();
			while (result.next()) {
				returnuser = new UserInfo(result.getInt("id"),
						result.getString("username"),
						result.getString("password"));
			}
		} catch (SQLException sqle) {
			sqle.printStackTrace();
		} finally {
			// Cleanup
			result.close();
			pstmt.close();
			conn.close();
		}

		return returnuser;
	}

	protected void setreturnonecolumn(List<Double> returnonecolumn) {
		MySQLConnection.returnonecolumn = returnonecolumn;

	}

	public static List<Double> getreturnonecolumn() {
		return returnonecolumn;
	}

	protected void setreturntime(List<Double> returntime) {
		MySQLConnection.returntime = returntime;

	}

	public static List<Double> getreturntime() {
		return returntime;
	}

	protected void setreturndelta(List<Double> returndelta) {
		MySQLConnection.returndelta = returndelta;

	}

	public static List<Double> getreturndelta() {
		return returndelta;
	}

	protected void setReturnfactors(List<List<Double>> returnfactors) {
		MySQLConnection.returnfactors = returnfactors;

	}

	public static List<List<Double>> getReturnfactors() {
		return returnfactors;
	}

	public PlotData getPlotData(DataFilter filter) throws SQLException {

		String[] factorStr = filter.getfactorStr();

		returnfactors = new ArrayList<List<Double>>();

		returnonecolumn = new ArrayList<Double>();

		returntime = new ArrayList<Double>();
		returndelta = new ArrayList<Double>();

		Connection conn = null;
		PreparedStatement pstmt = null;
		ResultSet result = null;
		PreparedStatement pstmt1 = null;
		ResultSet result1 = null;

		try {

			conn = getConnection();
			// pstmt1=conn.prepareStatement()
			pstmt = conn
					.prepareStatement("SELECT time, delta FROM breastdataset");
			//
			result = pstmt.executeQuery();
			while (result.next()) {

				returntime.add(new Double(result.getDouble("time")));
				returndelta.add(new Double(result.getDouble("delta")));

			}

			for (int i = 0; i < factorStr.length; i++) {

				if ("Grade".equals(factorStr[i])) {

					pstmt1 = conn
							.prepareStatement("SELECT Grade FROM breastdataset");
					result1 = pstmt1.executeQuery();
					while (result1.next()) {

						returnonecolumn.add(result1.getDouble("Grade"));
					}

					returnfactors.add(returnonecolumn);

				}

				else if ("size".equals(factorStr[i])) {

					pstmt1 = conn
							.prepareStatement("SELECT size FROM breastdataset");
					result1 = pstmt1.executeQuery();
					while (result1.next()) {

						returnonecolumn.add(result1.getDouble("size"));
					}

					returnfactors.add(returnonecolumn);
				}

				else if ("nodes".equals(factorStr[i])) {

					pstmt1 = conn
							.prepareStatement("SELECT nodes FROM breastdataset");
					result1 = pstmt1.executeQuery();
					while (result1.next()) {

						returnonecolumn.add(result1.getDouble("nodes"));
					}

					returnfactors.add(returnonecolumn);
				}

				else if ("ER".equals(factorStr[i])) {

					pstmt1 = conn
							.prepareStatement("SELECT ER FROM breastdataset");
					result1 = pstmt1.executeQuery();
					while (result1.next()) {

						returnonecolumn.add(result1.getDouble("ER"));
					}

					returnfactors.add(returnonecolumn);
				}

				else if ("agedx".equals(factorStr[i])) {

					pstmt1 = conn
							.prepareStatement("SELECT agedx FROM breastdataset");
					result1 = pstmt1.executeQuery();
					while (result1.next()) {

						returnonecolumn.add(result1.getDouble("agedx"));
					}

					returnfactors.add(returnonecolumn);
				}

				else if ("race".equals(factorStr[i])) {

					pstmt1 = conn
							.prepareStatement("SELECT race FROM breastdataset");
					result1 = pstmt1.executeQuery();
					while (result1.next()) {

						returnonecolumn.add(result1.getDouble("race"));
					}
					returnfactors.add(returnonecolumn);

				}

				else {
				}

			}

		} catch (Exception sqle) {
			sqle.printStackTrace();
		} finally {
			// Cleanup
			result.close();
			pstmt.close();
			result1.close();
			pstmt1.close();
			conn.close();
		}

		try {

			DataModel model = DataModel.getInstance();
			model.start(filter);
			plotData = new PlotData();
			plotData.setCurves(convert(model.getCurves()));
			plotData.setMaxSurvivalTime(model.getMaxSurvivalTime());
			plotData.setRoot(model.getRootBranch());

		} catch (Exception e) {
			e.printStackTrace();
		}

		return plotData;
	}

	private Curve[] convert(List<Curve> curves) {
		Curve[] arr = new Curve[curves.size()];
		for (int i = 0; i < arr.length; i++) {
			arr[i] = curves.get(i);
		}
		return arr;
	}

}