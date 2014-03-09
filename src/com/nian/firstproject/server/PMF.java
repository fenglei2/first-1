package com.nian.firstproject.server;

import java.util.logging.Logger;

import javax.jdo.JDOHelper;
import javax.jdo.PersistenceManagerFactory;

public final class PMF {

	@SuppressWarnings("unused")
	private static final Logger log = Logger.getLogger(PMF.class.getName());

	// GAE xml
	private static final PersistenceManagerFactory pmfInstance = JDOHelper
			.getPersistenceManagerFactory("transactions-optional");

	// tomcat6 txt file settings
	// private static final PersistenceManagerFactory pmfInstance =
	// JDOHelper.getPersistenceManagerFactory("datanucleus.properties");

	private PMF() {
	}

	public static PersistenceManagerFactory get() {
		return pmfInstance;
	}

}
