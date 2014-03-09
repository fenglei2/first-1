package com.nian.firstproject.client.widges;

import ca.nanometrics.gflot.client.DataPoint;
import ca.nanometrics.gflot.client.PlotModel;
import ca.nanometrics.gflot.client.SeriesHandler;
import ca.nanometrics.gflot.client.SimplePlot;
import ca.nanometrics.gflot.client.options.AxisOptions;
import ca.nanometrics.gflot.client.options.GlobalSeriesOptions;
import ca.nanometrics.gflot.client.options.LegendOptions;
import ca.nanometrics.gflot.client.options.LegendOptions.LegendPosition;
import ca.nanometrics.gflot.client.options.LineSeriesOptions;
import ca.nanometrics.gflot.client.options.PlotOptions;

import com.google.gwt.uibinder.client.UiField;
import com.nian.firstproject.shared.Curve;

//Survival curve UI in "Predict survival time" tab

public class SurvCurveWidgePre extends SimplePlot {

	@UiField(provided = true)
	// public
	static SurvCurveWidgePre plot = null;
	private static Curve curve;

	public static Curve getCurves() {
		return curve;
	}

	public static void setCurves(Curve curve) {
		SurvCurveWidgePre.curve = curve;
	}

	private static double maxTime;
	private static int step;

	public SurvCurveWidgePre(PlotModel plotModel, PlotOptions plotOptions) {
		super(plotModel, plotOptions);
	}

	public static SurvCurveWidgePre createEmptyPlot() {
		initial();
		setData();
		return plot;
	}

	public static SurvCurveWidgePre createPlot(Curve oneCurve, double max) {

		setCurves(oneCurve);

		maxTime = round(max);
		System.out.println("maxtime " + maxTime);
		initial();
		drawCurves();
		return plot;
	}

	private static double round(double time) {
		int i = 0;
		int n = 10;
		// int step;

		if (time < 100) {
			while (time > n && n < 100) {
				n *= 10;
			}
			step = n / 10;
		} else if (time < 500) {
			step = 50;
		} else {
			step = 100;
		}
		// System.out.println("step: "+step+" time: "+time+" i: "+i);

		n = 10;
		while (time > n) {
			n *= 10;
			i++;
		}

		int power10 = (int) (Math.pow(10, i));
		int round = (int) time / power10 * power10;
		// System.out.println("n: "+n+" round: "+round);
		while (round < time) {
			round += step;
		}
		return round;
	}

	private static void drawCurves() {

		System.out.println("this label: " + curve.getLabel());

		plot.getModel().addSeries(curve.getLabel());// i ( size i)

		SeriesHandler series = plot.getModel().getHandlers().get(0);

		for (int j = 0; j < curve.getTime().length; j++) {

			series.add(new DataPoint(curve.getTime()[j], curve
					.getProbabilityOfSurvival()[j]));
		}

	}

	private static void setData() {

		plot.getModel().addSeries("");
		plot.getModel().getHandlers().get(0).add(new DataPoint(0, 0));
	}

	public static SurvCurveWidgePre initial() {
		PlotModel model = new PlotModel();
		// model

		PlotOptions plotOptions = new PlotOptions();

		plotOptions.setGlobalSeriesOptions(new GlobalSeriesOptions()
				.setLineSeriesOptions(new LineSeriesOptions().setLineWidth(2))
				.setShadowSize(0));

		plotOptions.setLegendOptions(new LegendOptions().setMargin(200, -100)
				.setNumOfColumns(5).setBackgroundOpacity(0)
				.setPosition(LegendPosition.SOUTH_EAST));
		plotOptions.addXAxisOptions(new AxisOptions()
				.setLabel("Survival time (month)").setMinimum(0)
				// .setMaximum(maxTime + 30));
				.setMaximum(120));
		plotOptions.addYAxisOptions(new AxisOptions()
				.setLabel("Probability of survival").setMinimum(0.5)
				.setMaximum(1));

		return plot = new SurvCurveWidgePre(model, plotOptions);
	}

}
