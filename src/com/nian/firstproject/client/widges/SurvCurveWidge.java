package com.nian.firstproject.client.widges;

import java.util.ArrayList;
import java.util.List;

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
//import com.google.gwt.user.client.Random;
import com.nian.firstproject.shared.Curve;

//Survival curve UI in "Group patients" tab

public class SurvCurveWidge extends SimplePlot {

	@UiField(provided = true)
	public static SurvCurveWidge plot = null;
	static List<Curve> curves;
	static double two;
	static double one;
	static int four;
	static int size;

	public static List<Curve> getCurves() {
		return curves;
	}

	public static void setCurves(List<Curve> curves) {
		SurvCurveWidge.curves = curves;
	}

	private static double maxTime;
	private static int step;

	public SurvCurveWidge(PlotModel plotModel, PlotOptions plotOptions) {
		super(plotModel, plotOptions);
	}

	public static SurvCurveWidge createEmptyPlot() {
		initial();
		// plot.setData();
		setData();
		return plot;
	}

	public static SurvCurveWidge createPlot(Curve[] list, double max) {

		setCurves(list);
		maxTime = round(max);
		System.out.println("maxtime " + maxTime);
		initial();
		drawCurves();

		// plot.generateRandomData();
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

		for (int i = 0; i < curves.size(); i++) {
			// System.out.println("label: "+curves.get(i).getLabel());

			size = curves.get(i).getTime().length;
			double[] b = new double[size];

			for (int j = 0; j < size; j++) {

				b[i] = Math.abs(curves.get(i).getTime()[j] - 120);

			}

			double c = getMin(b);

			for (int j = 0; j < size; j++) {

				if (curves.get(i).getTime()[j] == 120 + c)

				{
					two = curves.get(i).getProbabilityOfSurvival()[j];
				}
				if (curves.get(i).getTime()[j] == 120 - c) {
					two = curves.get(i).getProbabilityOfSurvival()[j];
				}
			}

			double three = two * 100;

			four = (int) (three);

			plot.getModel().addSeries(
					curves.get(i).getLabel() + ":" + String.valueOf(four));
		}

		for (int i = 0; i < plot.getModel().getHandlers().size(); i++) {
			SeriesHandler series = plot.getModel().getHandlers().get(i);

			Curve curve = curves.get(i);

			// i=1
			for (int j = 0; j < curve.getTime().length; j++) {

				// System.out.println("---"+curve.getTime()[j]+" "+curve.getProbabilityOfSurvival()[j]);

				// System.out.println("---"+curve.getTime()[120]+" "+curve.getProbabilityOfSurvival()[120]);

				series.add(new DataPoint(curve.getTime()[j], curve
						.getProbabilityOfSurvival()[j]));
			}
		}

	}

	private static void setCurves(Curve[] list) {
		curves = new ArrayList<Curve>();
		for (Curve cu : list) {
			curves.add(cu);
		}
	}

	private static void setData() {
		plot.getModel().addSeries("");
		plot.getModel().getHandlers().get(0).add(new DataPoint(0, 0));
	}

	public static void initial() {
		PlotModel model = new PlotModel();
		PlotOptions plotOptions = new PlotOptions();

		plotOptions.setGlobalSeriesOptions(new GlobalSeriesOptions()
				.setLineSeriesOptions(new LineSeriesOptions().setLineWidth(1))
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

		plot = new SurvCurveWidge(model, plotOptions);

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
