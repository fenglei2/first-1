package com.nian.firstproject.client;

import ca.nanometrics.gflot.client.DataPoint;
import ca.nanometrics.gflot.client.SeriesHandler;

import com.google.gwt.event.dom.client.ClickEvent;
import com.google.gwt.event.dom.client.ClickHandler;
import com.google.gwt.user.client.Window;
import com.google.gwt.user.client.rpc.AsyncCallback;
import com.google.gwt.user.client.ui.Button;
import com.google.gwt.user.client.ui.Composite;
import com.google.gwt.user.client.ui.PushButton;
import com.google.gwt.user.client.ui.VerticalPanel;
import com.google.gwt.user.client.ui.HorizontalPanel;
import com.google.gwt.user.client.ui.Label;
import com.google.gwt.user.client.ui.TextBox;
import com.nian.firstproject.client.blobs.BlobDataFilter;
import com.nian.firstproject.client.blobs.PlotData;
import com.nian.firstproject.client.rpc.RpcInit;
import com.nian.firstproject.client.rpc.RpcServiceAsync;
import com.nian.firstproject.client.widges.DendgWidge;
import com.nian.firstproject.client.widges.SurvCurveWidge;
import com.nian.firstproject.client.widges.SurvCurveWidgePre;
import com.nian.firstproject.shared.Curve;
import com.swtoolbox.canvasfont.client.SWTBCanvasText;

import java.util.ArrayList;
import java.util.List;

import com.google.gwt.user.client.ui.RadioButton;
import com.google.gwt.user.client.ui.HasHorizontalAlignment;

public class GetCurveRoll extends Composite {
	
	//www.cap.org; home: select dataset;excel;csv,txt format;description;

	static double two;
	static double one;
	static int four;
	static int size;

	// private DBConnectionAsync rpc;
	private RpcServiceAsync rpc;
	int[] selectItems;// l
	static String[] factorStr;// m
	String[] factorLabels;// n
	static Curve[] curves;
	// static Curve[] finalcurves;//p
	static List<Curve> finalcurves;
	String label;
	ArrayList<String> factorStrAL;
	ArrayList<String> factorLabelsAL;
	ArrayList<Integer> selectItemsAL;

	PushButton submitButton;
	public static VerticalPanel pWidget;
	VerticalPanel verticalPanel;

	TextBox textBox;
	private String blobKey;

	public static SurvCurveWidge survivalP;
	public static SWTBCanvasText dendgP;
	public static SurvCurveWidgePre plot;
	// public static SurvCurveWidge plot;
	static HorizontalPanel horizontalPanel_1;
	static HorizontalPanel hp1;

	RadioButton rdbtnBreastcancer;
	RadioButton rdbtnLungCancer;

	RadioButton radioButton;

	RadioButton rdbtnNewRadioButton;

	RadioButton radioButton_1;

	RadioButton radioButton_2;

	RadioButton rdbtnT;

	RadioButton rdbtnT_1;

	RadioButton rdbtnT_2;

	RadioButton rdbtnT_3;

	RadioButton rdbtnN;

	RadioButton rdbtnN_1;

	RadioButton rdbtnN_2;

	RadioButton rdbtnN_3;

	RadioButton rdbtnN_4;

	RadioButton rdbtnPos_1;

	RadioButton rdbtnNeg_1;

	RadioButton rdbtnWhite;

	RadioButton rdbtnBlack;

	RadioButton rdbtnAsian;

	RadioButton rdbtnMale;
	RadioButton rdbtnFemale;

	RadioButton rdbtnYes;
	RadioButton rdbtnNo;

	RadioButton rdbtnYes_1;
	RadioButton rdbtnNo_1;

	RadioButton rdbtnRight;
	RadioButton rdbtnLeft;

	RadioButton rdbtnPos;
	RadioButton rdbtnNeg;

	RadioButton rdbtnPos_2;
	RadioButton rdbtnNeg_2;

	RadioButton rdbtnPos_3;
	RadioButton rdbtnNeg_3;

	RadioButton rdbtnPos_4;
	RadioButton rdbtnNeg_4;

	RadioButton rdbtnPos_5;

	RadioButton rdbtnNeg_5;

	RadioButton rdbtnPos_6;

	RadioButton rdbtnNeg_6;
	RadioButton radioButton_3;
	RadioButton radioButton_4;

	public static List<Curve> getfinalCurves() {
		return finalcurves;
	}

	public static void setfinalCurves(List<Curve> finalcurves) {
		GetCurveRoll.finalcurves = finalcurves;
	}

	public GetCurveRoll() {

		// factorStr = new String[6];
		// factorLabels = new String[6];
		factorStrAL = new ArrayList<String>();
		factorLabelsAL = new ArrayList<String>();

		// finalcurves = new ArrayList<Curve>();

		selectItemsAL = new ArrayList<Integer>();

		finalcurves = new ArrayList<Curve>();

		// selectItems = new int[6];

		plot = SurvCurveWidgePre.initial();
		// m=0;

		initWidget(getTheWidget());
		rpc = RpcInit.init();

	}

	public VerticalPanel getTheWidget() {
		if (pWidget == null) {
			pWidget = new VerticalPanel();
			pWidget.setSize("", "");
			pWidget.setStyleName("fontsize");
			pWidget.setSpacing(50);

			// VerticalPanel
			verticalPanel = new VerticalPanel();
			pWidget.add(verticalPanel);
			verticalPanel.setSize("", "");

			HorizontalPanel horizontalPanel_3 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_3);
			horizontalPanel_3.setWidth("1000");
			horizontalPanel_3.setSpacing(10);

			Label lblSelectCancerCategory = new Label("Cancer category is: ");
			lblSelectCancerCategory.setStyleName("customlabel");
			horizontalPanel_3.add(lblSelectCancerCategory);
			lblSelectCancerCategory.setWidth("303");

			rdbtnBreastcancer = new RadioButton("one", "Breast Cancer");
			rdbtnBreastcancer.setHTML("Breast Cancer");
			horizontalPanel_3.add(rdbtnBreastcancer);

			rdbtnLungCancer = new RadioButton("one", "Lung Cancer");
			horizontalPanel_3.add(rdbtnLungCancer);

			RadioButton rdbtnCancel = new RadioButton("new name", "Cancel");
			horizontalPanel_3.add(rdbtnCancel);

			HorizontalPanel horizontalPanel_2 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_2);
			horizontalPanel_2.setWidth("1000");
			horizontalPanel_2.setSpacing(10);

			Label lblGetPatientsSurvival = new Label("Survival and Hazard");
			lblGetPatientsSurvival
					.setHorizontalAlignment(HasHorizontalAlignment.ALIGN_CENTER);
			horizontalPanel_2.add(lblGetPatientsSurvival);

			HorizontalPanel horizontalPanel_4 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_4);
			horizontalPanel_4.setWidth("1000");
			horizontalPanel_4.setSpacing(10);

			Label lblAge = new Label("Age: ");
			lblAge.setStyleName("customlabel");
			horizontalPanel_4.add(lblAge);
			lblAge.setWidth("303");

			radioButton_3 = new RadioButton("new name", "<=50");
			horizontalPanel_4.add(radioButton_3);

			radioButton_4 = new RadioButton("new name", ">50");
			horizontalPanel_4.add(radioButton_4);

			RadioButton rdbtnNa_13 = new RadioButton("new name", "NA");
			horizontalPanel_4.add(rdbtnNa_13);

			HorizontalPanel horizontalPanel_5 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_5);
			horizontalPanel_5.setWidth("1000");
			horizontalPanel_5.setSpacing(10);

			Label lblGender = new Label("Gender: ");
			lblGender.setStyleName("customlabel");
			horizontalPanel_5.add(lblGender);
			lblGender.setWidth("303");

			rdbtnMale = new RadioButton("two", "Male");
			horizontalPanel_5.add(rdbtnMale);

			rdbtnFemale = new RadioButton("two", "Female");
			horizontalPanel_5.add(rdbtnFemale);

			RadioButton rdbtnNa_14 = new RadioButton("new name", "NA");
			horizontalPanel_5.add(rdbtnNa_14);

			HorizontalPanel horizontalPanel_6 = new HorizontalPanel();
			horizontalPanel_6.setStyleName("gwt-RadioButton");
			verticalPanel.add(horizontalPanel_6);
			horizontalPanel_6.setWidth("1000");
			horizontalPanel_6.setSpacing(10);

			Label lblRacialethnicGroup = new Label("Racial/Ethnic Group: ");
			lblRacialethnicGroup.setStyleName("customlabel");
			horizontalPanel_6.add(lblRacialethnicGroup);
			lblRacialethnicGroup.setWidth("303");

			rdbtnWhite = new RadioButton("three", "White");
			horizontalPanel_6.add(rdbtnWhite);
			rdbtnWhite.setWidth("114px");

			rdbtnBlack = new RadioButton("three", "Black");
			horizontalPanel_6.add(rdbtnBlack);

			rdbtnAsian = new RadioButton("three", "Asian");
			horizontalPanel_6.add(rdbtnAsian);

			RadioButton rdbtnOther = new RadioButton("three", "Other");
			horizontalPanel_6.add(rdbtnOther);

			HorizontalPanel horizontalPanel_7 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_7);
			horizontalPanel_7.setWidth("1000");
			horizontalPanel_7.setSpacing(10);

			Label lblHistologicalTUmor = new Label("Histological Tumor Type:");
			lblHistologicalTUmor.setStyleName("customlabel");
			horizontalPanel_7.add(lblHistologicalTUmor);

			HorizontalPanel horizontalPanel_8 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_8);
			horizontalPanel_8.setWidth("1000");
			horizontalPanel_8.setSpacing(10);

			Label lblHistologicalGrade = new Label("Histological Grade: ");
			lblHistologicalGrade.setStyleName("customlabel");
			horizontalPanel_8.add(lblHistologicalGrade);
			lblHistologicalGrade.setWidth("303");

			radioButton = new RadioButton("four", "1");
			horizontalPanel_8.add(radioButton);

			rdbtnNewRadioButton = new RadioButton("four", "2");
			horizontalPanel_8.add(rdbtnNewRadioButton);

			radioButton_1 = new RadioButton("four", "3");
			horizontalPanel_8.add(radioButton_1);

			radioButton_2 = new RadioButton("four", "4");
			horizontalPanel_8.add(radioButton_2);

			RadioButton rdbtnNa = new RadioButton("four", "NA");
			horizontalPanel_8.add(rdbtnNa);

			HorizontalPanel horizontalPanel_11_1 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_11_1);
			horizontalPanel_11_1.setWidth("1000");
			horizontalPanel_11_1.setSpacing(10);

			Label lblTumorSize = new Label("Tumor Size: ");
			lblTumorSize.setStyleName("customlabel");
			horizontalPanel_11_1.add(lblTumorSize);
			lblTumorSize.setWidth("303");

			rdbtnT = new RadioButton("five", "T1");
			horizontalPanel_11_1.add(rdbtnT);

			rdbtnT_1 = new RadioButton("five", "T2");
			horizontalPanel_11_1.add(rdbtnT_1);

			rdbtnT_2 = new RadioButton("five", "T3");
			horizontalPanel_11_1.add(rdbtnT_2);

			rdbtnT_3 = new RadioButton("five", "T4");
			horizontalPanel_11_1.add(rdbtnT_3);

			RadioButton rdbtnNa_1 = new RadioButton("five", "NA");
			horizontalPanel_11_1.add(rdbtnNa_1);

			HorizontalPanel horizontalPanel_12_1 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_12_1);
			horizontalPanel_12_1.setWidth("1000");
			horizontalPanel_12_1.setSpacing(10);

			Label lblLymphNodeStatus = new Label("Lymph Node Status: ");
			lblLymphNodeStatus.setStyleName("customlabel");
			horizontalPanel_12_1.add(lblLymphNodeStatus);
			lblLymphNodeStatus.setWidth("303");

			rdbtnN = new RadioButton("six", "N0");
			horizontalPanel_12_1.add(rdbtnN);

			rdbtnN_1 = new RadioButton("six", "N1");
			horizontalPanel_12_1.add(rdbtnN_1);

			rdbtnN_2 = new RadioButton("six", "N2");
			horizontalPanel_12_1.add(rdbtnN_2);

			rdbtnN_3 = new RadioButton("six", "N3");
			horizontalPanel_12_1.add(rdbtnN_3);

			rdbtnN_4 = new RadioButton("six", "N4");
			horizontalPanel_12_1.add(rdbtnN_4);

			RadioButton rdbtnNa_2 = new RadioButton("six", "NA");
			horizontalPanel_12_1.add(rdbtnNa_2);

			HorizontalPanel horizontalPanel_13 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_13);
			horizontalPanel_13.setWidth("1000");
			horizontalPanel_13.setSpacing(10);

			Label lblMetastaticDisease = new Label("Metastatic Disease: ");
			lblMetastaticDisease.setStyleName("customlabel");
			horizontalPanel_13.add(lblMetastaticDisease);

			rdbtnYes = new RadioButton("seven", "Yes");
			horizontalPanel_13.add(rdbtnYes);

			rdbtnNo = new RadioButton("seven", "No");
			horizontalPanel_13.add(rdbtnNo);

			RadioButton rdbtnNa_3 = new RadioButton("seven", "NA");
			horizontalPanel_13.add(rdbtnNa_3);

			HorizontalPanel horizontalPanel_9 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_9);
			horizontalPanel_9.setWidth("1000");
			horizontalPanel_9.setSpacing(10);

			Label lblVascularInvasion = new Label("Vascular Invasion: ");
			lblVascularInvasion.setStyleName("customlabel");
			horizontalPanel_9.add(lblVascularInvasion);
			lblVascularInvasion.setWidth("303");

			rdbtnYes_1 = new RadioButton("eight", "Yes");
			horizontalPanel_9.add(rdbtnYes_1);

			rdbtnNo_1 = new RadioButton("eight", "No");
			horizontalPanel_9.add(rdbtnNo_1);

			RadioButton rdbtnNa_4 = new RadioButton("eight", "NA");
			horizontalPanel_9.add(rdbtnNa_4);

			HorizontalPanel horizontalPanel_10 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_10);
			horizontalPanel_10.setWidth("1000");
			horizontalPanel_10.setSpacing(10);

			Label lblLaterality = new Label("Laterality: ");
			lblLaterality.setStyleName("customlabel");
			horizontalPanel_10.add(lblLaterality);
			lblLaterality.setWidth("303");

			rdbtnRight = new RadioButton("nine", "Right");
			horizontalPanel_10.add(rdbtnRight);

			rdbtnLeft = new RadioButton("nine", "Left");
			horizontalPanel_10.add(rdbtnLeft);

			RadioButton rdbtnNa_5 = new RadioButton("nine", "NA");
			horizontalPanel_10.add(rdbtnNa_5);

			HorizontalPanel horizontalPanel_14 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_14);
			horizontalPanel_14.setWidth("1000");
			horizontalPanel_14.setSpacing(10);

			Label lblMargins = new Label("Surgical Margins: ");
			lblMargins.setStyleName("customlabel");
			horizontalPanel_14.add(lblMargins);
			lblMargins.setWidth("303");

			rdbtnPos = new RadioButton("ten", "Pos");
			horizontalPanel_14.add(rdbtnPos);

			rdbtnNeg = new RadioButton("ten", "Neg");
			horizontalPanel_14.add(rdbtnNeg);

			RadioButton rdbtnNa_6 = new RadioButton("ten", "NA");
			horizontalPanel_14.add(rdbtnNa_6);

			HorizontalPanel horizontalPanel_15 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_15);
			horizontalPanel_15.setWidth("1000");
			horizontalPanel_15.setSpacing(10);

			Label lblEstrogenReceptor = new Label("Estrogen Receptor: ");
			lblEstrogenReceptor.setStyleName("customlabel");
			horizontalPanel_15.add(lblEstrogenReceptor);
			lblEstrogenReceptor.setWidth("303");

			rdbtnPos_1 = new RadioButton("eleven", "Pos");
			horizontalPanel_15.add(rdbtnPos_1);

			rdbtnNeg_1 = new RadioButton("eleven", "Neg");
			horizontalPanel_15.add(rdbtnNeg_1);

			RadioButton rdbtnNa_7 = new RadioButton("eleven", "NA");
			horizontalPanel_15.add(rdbtnNa_7);

			HorizontalPanel horizontalPanel_16 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_16);
			horizontalPanel_16.setWidth("1000");
			horizontalPanel_16.setSpacing(10);

			Label lblProgesteroneReceptor = new Label("Progesterone Receptor: ");
			lblProgesteroneReceptor.setStyleName("customlabel");
			horizontalPanel_16.add(lblProgesteroneReceptor);
			lblProgesteroneReceptor.setWidth("303");

			rdbtnPos_2 = new RadioButton("twelve", "Pos");
			horizontalPanel_16.add(rdbtnPos_2);

			rdbtnNeg_2 = new RadioButton("twelve", "Neg");
			horizontalPanel_16.add(rdbtnNeg_2);

			RadioButton rdbtnNa_8 = new RadioButton("twelve", "NA");
			horizontalPanel_16.add(rdbtnNa_8);

			HorizontalPanel horizontalPanel_17 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_17);
			horizontalPanel_17.setWidth("1000");
			horizontalPanel_17.setSpacing(10);

			Label lblHers = new Label("Her/2:");
			lblHers.setStyleName("customlabel");
			horizontalPanel_17.add(lblHers);
			lblHers.setWidth("303");

			rdbtnPos_3 = new RadioButton("thirteen", "Pos");
			horizontalPanel_17.add(rdbtnPos_3);

			rdbtnNeg_3 = new RadioButton("thirteen", "Neg");
			horizontalPanel_17.add(rdbtnNeg_3);

			RadioButton rdbtnNa_9 = new RadioButton("thirteen", "NA");
			horizontalPanel_17.add(rdbtnNa_9);

			HorizontalPanel horizontalPanel_19 = new HorizontalPanel();
			horizontalPanel_19.setSpacing(10);
			verticalPanel.add(horizontalPanel_19);
			horizontalPanel_19.setWidth("1000");

			Label lblCauseOfDeath = new Label("Cause of Death: Breast Cancer");
			lblCauseOfDeath.setStyleName("customlabel");
			horizontalPanel_19.add(lblCauseOfDeath);
			lblCauseOfDeath.setWidth("303");

			HorizontalPanel horizontalPanel_18 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_18);
			horizontalPanel_18.setWidth("1000");
			horizontalPanel_18.setSpacing(10);

			Label lblMarker = new Label("Marker 1: ");
			lblMarker.setStyleName("customlabel");
			horizontalPanel_18.add(lblMarker);
			lblMarker.setWidth("303");

			rdbtnPos_4 = new RadioButton("fourteen", "Pos");
			horizontalPanel_18.add(rdbtnPos_4);

			rdbtnNeg_4 = new RadioButton("fourteen", "Neg");
			horizontalPanel_18.add(rdbtnNeg_4);

			RadioButton rdbtnNa_10 = new RadioButton("fourteen", "NA");
			horizontalPanel_18.add(rdbtnNa_10);

			HorizontalPanel horizontalPanel = new HorizontalPanel();
			verticalPanel.add(horizontalPanel);
			horizontalPanel.setWidth("1000");

			Label lblMarker_1 = new Label("Marker 2: ");
			lblMarker_1.setStyleName("customlabel");
			horizontalPanel.add(lblMarker_1);
			lblMarker_1.setWidth("303");

			rdbtnPos_5 = new RadioButton("fifteen", "Pos");
			horizontalPanel.add(rdbtnPos_5);

			rdbtnNeg_5 = new RadioButton("fifteen", "Neg");
			horizontalPanel.add(rdbtnNeg_5);
			horizontalPanel.setSpacing(10);

			RadioButton rdbtnNa_11 = new RadioButton("fifteen", "NA");
			horizontalPanel.add(rdbtnNa_11);

			HorizontalPanel horizontalPanel_20 = new HorizontalPanel();
			verticalPanel.add(horizontalPanel_20);
			horizontalPanel_20.setWidth("1000");

			Label lblMaeker = new Label("Marker 3: ");
			lblMaeker.setStyleName("customlabel");
			horizontalPanel_20.add(lblMaeker);
			lblMaeker.setWidth("303");

			rdbtnPos_6 = new RadioButton("sixteen", "Pos");
			horizontalPanel_20.add(rdbtnPos_6);

			rdbtnNeg_6 = new RadioButton("sixteen", "Neg");
			horizontalPanel_20.add(rdbtnNeg_6);

			RadioButton rdbtnNa_12 = new RadioButton("sixteen", "NA");
			horizontalPanel_20.add(rdbtnNa_12);

			horizontalPanel_20.setSpacing(10);

			verticalPanel.add(getButton());

			horizontalPanel_1 = new HorizontalPanel();
			pWidget.add(horizontalPanel_1);
			horizontalPanel_1.setSize("", "");
			horizontalPanel_1.setSpacing(5);
			horizontalPanel_1.setVisible(false);

			Button btnNewButton = new Button("New button");
			btnNewButton.setText("back");
			horizontalPanel_1.add(btnNewButton);
			btnNewButton.addClickHandler(new ClickHandler() {

				@Override
				public void onClick(ClickEvent event) {

					// if (pWidget.getWidget(2) != null) {
					// pWidget.getWidget(2).removeFromParent();

					// }
					pWidget.getWidget(2).setVisible(false);
					pWidget.getWidget(1).setVisible(false);

					verticalPanel.setVisible(true);

				}

			});

			Button btnNewButton_3 = new Button("New button");
			btnNewButton_3.setText("PatientSurvivalCurve");
			horizontalPanel_1.add(btnNewButton_3);

			btnNewButton_3.addClickHandler(new ClickHandler() {

				@Override
				public void onClick(ClickEvent event) {
					if (pWidget.getWidget(2) != null) {
						pWidget.getWidget(2).removeFromParent();

					}

					pWidget.add(plot);
					pWidget.getWidget(2).setVisible(true);
				}

			});

			Button btnNewButton_4 = new Button("New button");
			btnNewButton_4.setText("PatientHazardCurve");
			horizontalPanel_1.add(btnNewButton_4);

			btnNewButton_4.addClickHandler(new ClickHandler() {

				@Override
				public void onClick(ClickEvent event) {
					if (pWidget.getWidget(2) != null) {
						pWidget.getWidget(2).removeFromParent();

					}
					// add hazard curves here?
					pWidget.add(SurvCurveWidge.createEmptyPlot());
					pWidget.getWidget(2).setVisible(true);
				}

			});

			Button btnNewButton_1 = new Button("New button");
			btnNewButton_1.setText("SurvivalCurves");
			horizontalPanel_1.add(btnNewButton_1);

			btnNewButton_1.addClickHandler(new ClickHandler() {

				@Override
				public void onClick(ClickEvent event) {
					if (pWidget.getWidget(2) != null) {
						pWidget.getWidget(2).removeFromParent();

					}

					pWidget.add(survivalP);
					pWidget.getWidget(2).setVisible(true);
				}

			});

			Button btnH = new Button("H");
			btnH.setText("HazardCurves");
			horizontalPanel_1.add(btnH);
			btnH.addClickHandler(new ClickHandler() {

				@Override
				public void onClick(ClickEvent event) {
					if (pWidget.getWidget(2) != null) {
						pWidget.getWidget(2).removeFromParent();

					}
					// add hazard curves here?
					pWidget.add(SurvCurveWidge.createEmptyPlot());
					pWidget.getWidget(2).setVisible(true);
				}

			});

			Button btnNewButton_2 = new Button("New button");
			btnNewButton_2.setText("Dendrogram");
			horizontalPanel_1.add(btnNewButton_2);
			btnNewButton_2.addClickHandler(new ClickHandler() {

				@Override
				public void onClick(ClickEvent event) {

					if (pWidget.getWidget(2) != null) {
						pWidget.getWidget(2).removeFromParent();

					}

					pWidget.add(dendgP);
					pWidget.getWidget(2).setVisible(true);

				}

			});

			pWidget.add(SurvCurveWidgePre.createEmptyPlot());// 2
			pWidget.getWidget(2).setVisible(false);

		}

		return pWidget;
	}

	private PushButton getButton() {

		submitButton = new PushButton();
		submitButton.setHTML("Submit");
		submitButton.addClickHandler(new ClickHandler() {

			@SuppressWarnings("deprecation")
			@Override
			public void onClick(ClickEvent event) {

				factorStrAL = new ArrayList<String>();
				factorLabelsAL = new ArrayList<String>();

				// finalcurves = new ArrayList<Curve>();

				selectItemsAL = new ArrayList<Integer>();

				plot = SurvCurveWidgePre.initial();

				if (rdbtnBreastcancer.isChecked()) {

					// GetCurveRoll.pWidget.setStyleName(blobKey);

					// populateData();

					if (radioButton.isChecked()) {
						selectItemsAL.add(0);
						factorStrAL.add("Grade");
						factorLabelsAL.add("G1");
					}

					else if (rdbtnNewRadioButton.isChecked()) {
						selectItemsAL.add(0);
						factorStrAL.add("Grade");
						factorLabelsAL.add("G2");
					}

					else if (radioButton_1.isChecked()) {
						selectItemsAL.add(0);
						factorStrAL.add("Grade");
						factorLabelsAL.add("G3");
					}

					else if (radioButton_2.isChecked()) {
					//	selectItemsAL.add(0);
					//	factorStrAL.add("Grade");
					//	factorLabelsAL.add("G4");

					} else {
					}

					if (rdbtnT.isChecked()) {
						selectItemsAL.add(1);
						factorStrAL.add("size");
						factorLabelsAL.add("T1");
					}

					else if (rdbtnT_1.isChecked()) {
						selectItemsAL.add(1);
						factorStrAL.add("size");
						factorLabelsAL.add("T2");
					}

					else if (rdbtnT_2.isChecked()) {
						selectItemsAL.add(1);
						factorStrAL.add("size");
						factorLabelsAL.add("T3");
					}

					else if (rdbtnT_3.isChecked()) {
						selectItemsAL.add(1);
						factorStrAL.add("size");
						factorLabelsAL.add("T4");
					} else {
					}

					if (rdbtnN.isChecked()) {
						selectItemsAL.add(2);
						factorStrAL.add("nodes");
						factorLabelsAL.add("N0");
					}

					else if (rdbtnN_1.isChecked()) {
						selectItemsAL.add(2);
						factorStrAL.add("nodes");
						factorLabelsAL.add("N1");
					}

					else if (rdbtnN_2.isChecked()) {
						selectItemsAL.add(2);
						factorStrAL.add("nodes");
						factorLabelsAL.add("N2");
					}

					else if (rdbtnN_3.isChecked()) {
						selectItemsAL.add(2);
						factorStrAL.add("nodes");
						factorLabelsAL.add("N3");
					}

					else if (rdbtnN_4.isChecked()) {
						selectItemsAL.add(2);
						factorStrAL.add("nodes");
						factorLabelsAL.add("N4");
					}

					else {
					}

					if (rdbtnPos_1.isChecked()) {
						selectItemsAL.add(3);
						factorStrAL.add("ER");
						factorLabelsAL.add("ER+");
					}

					else if (rdbtnNeg_1.isChecked()) {
						selectItemsAL.add(3);
						factorStrAL.add("ER");
						factorLabelsAL.add("ER-");
					}

					else {
					}

					// verticalPanel_1.add(SurvCurveWidgePre.createEmptyPlot());

					if (radioButton_3.isChecked()) {
						selectItemsAL.add(4);
						factorStrAL.add("agedx");
						factorLabelsAL.add("age1");
						// factorNum = factorNum + 1;
						// factorStr[factorNum - 1] = "agedx";
						// factorLabels[factorNum - 1] = "age1";

					}

					else if (radioButton_4.isChecked()) {
						selectItemsAL.add(4);
						factorStrAL.add("agedx");
						factorLabelsAL.add("age2");

					}

					else {

					}

					if (rdbtnWhite.isChecked()) {
						selectItemsAL.add(5);
						factorStrAL.add("race");
						factorLabelsAL.add("race1");
					}

					else if (rdbtnBlack.isChecked()) {
						selectItemsAL.add(5);
						factorStrAL.add("race");
						factorLabelsAL.add("race2");
					}

					else if (rdbtnAsian.isChecked()) {

						Window.alert("We do not have Asian Option in our dataset now, you can select other, or you can upload one for us.");

					}

					else {
					}

					/**
					 * if (rdbtnMale.isChecked() | rdbtnFemale.isChecked()) {
					 * Window.alert(
					 * "We do not have Gender factor in our dataset now, you can upload one for us."
					 * ); }
					 * 
					 * if (rdbtnYes.isChecked() | rdbtnNo.isChecked()) {
					 * Window.alert(
					 * "We do not have Metastatic Disease factor in our dataset now, you can select NA, or you can upload one for us."
					 * ); }
					 * 
					 * if (rdbtnYes_1.isChecked() | rdbtnNo_1.isChecked()) {
					 * Window.alert(
					 * "We do not have Vascular Invasion factor in our dataset now, you can select NA, or you can upload one for us."
					 * ); }
					 * 
					 * if (rdbtnRight.isChecked() | rdbtnLeft.isChecked()) {
					 * Window.alert(
					 * "We do not have Laterality factor in our dataset now, you can select NA, or you can upload one for us."
					 * ); }
					 * 
					 * if (rdbtnPos.isChecked() | rdbtnNeg.isChecked()) {
					 * Window.alert(
					 * "We do not have Margins factor in our dataset now, you can select NA, or you can upload one for us."
					 * ); }
					 * 
					 * if (rdbtnPos_2.isChecked() | rdbtnNeg_2.isChecked()) {
					 * Window.alert(
					 * "We do not have Progesterone Receptor factor in our dataset now, you can select NA, or you can upload one for us."
					 * ); }
					 * 
					 * if (rdbtnPos_3.isChecked() | rdbtnNeg_3.isChecked()) {
					 * Window.alert(
					 * "We do not have Her/2 factor in our dataset now, you can select NA, or you can upload one for us."
					 * ); }
					 * 
					 * if (rdbtnPos_4.isChecked() | rdbtnNeg_4.isChecked() |
					 * rdbtnPos_5.isChecked() | rdbtnNeg_5.isChecked() |
					 * rdbtnPos_6.isChecked() | rdbtnNeg_6.isChecked()) {
					 * Window.alert("What is this Marker for?"); }
					 */

				}// first breastdataset if over

				else {
					Window.alert("Select breast cancer;We do not have this dataset now, you can upload one for us.");
				}

				factorStr = new String[factorStrAL.size()];
				factorLabels = new String[factorLabelsAL.size()];

				selectItems = new int[selectItemsAL.size()];

				for (int i = 0; i < factorStr.length; i++) {

					factorStr[i] = factorStrAL.get(i);
					factorLabels[i] = factorLabelsAL.get(i);
					selectItems[i] = selectItemsAL.get(i);

				}

				verticalPanel.setVisible(false);// 0

				label = "";

				for (int i = 0; i < factorLabels.length; i++) {

					label += factorLabels[i];

				}

				horizontalPanel_1.setVisible(true);// 1
				// pWidget.add(plot);

				// hp4.setName(factorStr);
				// factorBox.setName(factors);

				// selectItems=new int[factorStrAL.size()];

				showPlots(selectItems);// 1

			}

		});

		return submitButton;

	}

	/**
	 * private static void drawCurve(String label) {
	 * 
	 * List<Curve> allCurves = SurvCurveWidge.getCurves(); boolean find = false;
	 * for (Curve curve : allCurves) {
	 * 
	 * // System.out.println("exist label: " + curve.getLabel()); // for(int
	 * i=0; i<allCurves.size();i++){
	 * 
	 * while (curve.getLabel().equals(label)) {
	 * 
	 * // hp4=new HorizontalPanel(); // hp4.setVisible(true);
	 * 
	 * pWidget.add(SurvCurveWidgePre.createPlot(curve, 120));//
	 * pWidget.getWidget(2); plot=SurvCurveWidgePre.createPlot(curve, 120);
	 * 
	 * 
	 * find = true;
	 * 
	 * break; }
	 * 
	 * if (!find) {
	 * 
	 * Window.alert("Sorry, no such combination!");
	 * 
	 * }
	 * 
	 * }
	 * 
	 * }
	 */

	public static boolean isInteger(String input) {
		try {
			Integer.parseInt(input);
			return true;
		} catch (Exception e) {
			return false;
		}
	}

	protected void setCurves(Curve[] curves) {
		GetCurveRoll.curves = curves;

	}

	public static Curve[] getCurves() {
		return curves;
	}

	private void showPlots(int[] selectItems) {

		blobKey = pWidget.getStyleName();

		BlobDataFilter filter = new BlobDataFilter();
		// DataFilter filter = new DataFilter();
		filter.setSelectItems(selectItems);
		filter.setBlobKey(blobKey);
		filter.setCensorCol(7);
		filter.setTimeCol(6);
		filter.setMinPatientsNo(100);

		rpc.getPlotData(filter, new AsyncCallback<PlotData>() {

			@Override
			public void onSuccess(PlotData result) {

				curves = result.getCurves();

				setCurves(curves);
				
				for(int i=0;i<curves.length;i++){			
					
				if (curves[i].getLabel().equals(label)) {
					finalcurves.add(curves[i]);

				}
				}
				setfinalCurves(finalcurves);

				for (int j = 0; j < finalcurves.size(); j++) {
					// System.out.println("label: "+curves.get(i).getLabel());
					size = finalcurves.get(j).getTime().length;
					double[] b = new double[size];

					for (int k = 0; k < size; k++) {

						b[j] = Math.abs(finalcurves.get(j).getTime()[k] - 120);

					}

					double c = getMin(b);

					for (int k = 0; k < size; k++) {

						if (finalcurves.get(j).getTime()[k] == 120 + c)

						{
							two = finalcurves.get(j).getProbabilityOfSurvival()[k];
						}
						if (finalcurves.get(j).getTime()[k] == 120 - c) {
							two = finalcurves.get(j).getProbabilityOfSurvival()[k];
						}
					}

					double three = two * 100;

					four = (int) (three);

					plot.getModel().addSeries(
							finalcurves.get(j).getLabel() + ":"
									+ String.valueOf(four));
				}

				for (int j = 0; j < plot.getModel().getHandlers().size(); j++) {
					SeriesHandler series = plot.getModel().getHandlers().get(j);

					Curve curve = finalcurves.get(j);

					// i=1
					for (int m = 0; m < curve.getTime().length; m++) {

						// System.out.println("---"+curve.getTime()[j]+" "+curve.getProbabilityOfSurvival()[j]);

						// System.out.println("---"+curve.getTime()[120]+" "+curve.getProbabilityOfSurvival()[120]);

						series.add(new DataPoint(curve.getTime()[m], curve
								.getProbabilityOfSurvival()[m]));
					}
				}

				if (pWidget.getWidget(2) != null) {
					pWidget.getWidget(2).removeFromParent();

					pWidget.add(plot);
					pWidget.getWidget(2).setVisible(true);
				}

				survivalP = SurvCurveWidge.createPlot(result.getCurves(),
						result.getMaxSurvivalTime());// 1
				dendgP = DendgWidge.createPlot(result.getRoot());// 2

			}

			@Override
			public void onFailure(Throwable caught) {
				System.out.println("getPlotData error: " + caught);

			}

		});

	}

	private static double getMin(double[] b) {
		double min, max;
		max = min = b[0];
		for (int i = 0; i < b.length; i++) {
			if (b[i] > max) {
				max = b[i];
			}
			if (b[i] < min) {
				min = b[i];
			}
		}

		return min;

	}

}
