package com.nian.firstproject.client.widges;

import java.util.ArrayList;
import java.util.List;

import com.google.gwt.i18n.client.NumberFormat;
import com.google.gwt.widgetideas.graphics.client.Color;
import com.swtb.font.client.Simplex;
import com.swtoolbox.canvasfont.client.SWTBCanvasText;
import com.nian.firstproject.client.Cutting;
import com.nian.firstproject.client.GetCurveRoll;
import com.nian.firstproject.shared.Branch;
import com.nian.firstproject.shared.Curve;
import com.nian.firstproject.shared.Point;

//Dendrogram UI in "Group patients" tab

public class DendgWidge extends SWTBCanvasText {

	private static NumberFormat ROUND_TO_INT = NumberFormat.getFormat("0");
	private static NumberFormat ROUND_T0_ONEDEGIT = NumberFormat
			.getFormat("0.0");

	public static SWTBCanvasText canvas;

	public static int width = 1500, height = 2000;

	public static int centerX;//
	public static int centerY;//

	static List<Branch> allBranches;
	public static double MAX_WIDTH;// the distance between 2 last merging
									// static // clusters

	static double two;
	static double one;
	static int four;
	static int size;
	static double LEAF_WIDTH;

	public static int DIST_PAD;// the distance pad for a proper distance from
								// leftmost
	// node to scale

	public static int upmost;
	public static int downmost;
	public static int linedownMost;
	public static int lineupMost;

	private static Curve[] curves;

	private DendgWidge(int width, int height) {

		super(width, height);
	}

	public static SWTBCanvasText createPlot(Branch root) {
		createEmptyPlot();
		allBranches = new ArrayList<Branch>();
		root.getAllBranches(allBranches);

		MAX_WIDTH = Integer.valueOf(ROUND_TO_INT.format(root.getCluster()
				.getWidth()));
		System.out.println("max width: " + MAX_WIDTH + " "
				+ root.getCluster().getWidth());

		LEAF_WIDTH = MAX_WIDTH / 10;

		DIST_PAD = 2 - upMostPos(root);

		upmost = upMostPos(root);
		downmost = downMostPos(root);
		linedownMost = downMostPos2(root);
		lineupMost = upMostPos2(root);

		drawScale();
		draw(root);

		return canvas;
	}

	private static void draw(Branch root) {
		for (Branch branch : allBranches) {
			drawBranch(branch);
		}
	}

	private static void drawBranch(Branch branch) {
		int upSize = branch.getCluster().getUp().size();
		int downSize = branch.getCluster().getDown().size();

		if (upSize == 1) {
			branch.getUpRight().setX(branch.getUpLeft().getX() - LEAF_WIDTH);
			drawLeafLabel(branch.getCluster().getUp().getLabel(),
					branch.getUpRight());
		}

		if (downSize == 1) {
			branch.getDownRight()
					.setX(branch.getDownLeft().getX() - LEAF_WIDTH);
			drawLeafLabel(branch.getCluster().getDown().getLabel(),
					branch.getDownRight());
		}

		drawLine(branch.getUpLeft(), branch.getUpRight());
		drawLine(branch.getDownLeft(), branch.getDownRight());
		drawLine(branch.getUpLeft(), branch.getDownLeft());
	}

	private static void drawLine(Point p1, Point p2) {
		drawLine2(p1.getX(), p1.getY(), p2.getX(), p2.getY());
	}

	private static void drawLeafLabel(String label, Point point) {

		curves = GetCurveRoll.getCurves();

		for (Curve curve : curves) {

			if (curve.getLabel() == label) {

				size = curve.getTime().length;
				double[] b = new double[size];

				for (int i = 0; i < size; i++) {

					b[i] = Math.abs(curve.getTime()[i] - 120);

				}

				double c = getMin(b);

				for (int i = 0; i < size; i++) {

					if (curve.getTime()[i] == 120 + c)

					{
						two = curve.getProbabilityOfSurvival()[i];
					}
					if (curve.getTime()[i] == 120 - c) {
						two = curve.getProbabilityOfSurvival()[i];
					}
				}

			}

			double three = two * 100;

			four = (int) (three);

		}// end if

		label = label + ":" + String.valueOf(four);

		canvas.saveContext();
		canvas.setFont(new Simplex());
		canvas.beginPath();

		canvas.textTo(label, getX(point.getX()), getY2(point.getY()) - 4, 0,
				0.4);
		canvas.stroke();
		canvas.restoreContext();
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

	private static void drawLine(double x1, double y1, double x2, double y2) {
		canvas.beginPath();
		x1 = getX(x1);
		x2 = getX(x2);
		y1 = getY(y1);
		y2 = getY(y2);

		canvas.lineTo(x1, y1);
		canvas.lineTo(x2, y2);
		canvas.stroke();
	}

	private static void drawLine2(double x1, double y1, double x2, double y2) {
		canvas.beginPath();
		x1 = getX(x1);
		x2 = getX(x2);
		y1 = getY2(y1);
		y2 = getY2(y2);

		canvas.lineTo(x1, y1);
		canvas.lineTo(x2, y2);
		canvas.stroke();
	}

	private static double getX(double x) {

		// 1500 //2000
		return centerX - 600 * x;// 800
	}

	private static double getY(double y) {
		return centerY + 50 * y;// 50
	}

	private static double getY2(double y) {
		return centerY + 50 * (y + DIST_PAD);
	}

	private static void drawScale() {

		drawLine(0, 0, 1, 0);
		drawLine(0, -0.1, 0, 0);

		for (double h = 0; h <= 1; h += 0.2) {
			drawLine(h, -0.1, h, 0);
			// System.out.println("mh:" + ROUND_T0_ONEDEGIT.format(MAX_HEIGHT *
			// h));
			drawLabel(String.valueOf(ROUND_T0_ONEDEGIT.format(MAX_WIDTH * h)),
					h + 0.02, -0.5);
		}
	}

	private static void drawTitle() {
		canvas.saveContext();
		canvas.setFont(new Simplex());
		canvas.beginPath();
		canvas.textTo("Dendrogram", 30, height / 15, Math.toRadians(-90), 0.6);
		canvas.stroke();
		canvas.restoreContext();
	}

	private static void drawLabel(String label, double x, double y) {
		canvas.saveContext();
		canvas.setFont(new Simplex());
		canvas.beginPath();
		canvas.textTo(label, getX(x), getY(y), 0, 0.4);
		canvas.stroke();
		canvas.restoreContext();
	}

	public static SWTBCanvasText createEmptyPlot() {

		drawBorder();
		drawTitle();

		return canvas;
	}

	private static void drawBorder() {
		// 1500, 2000
		// 800
		canvas = new SWTBCanvasText(width, height);

		// ???????
		centerX = 800;

		// ????????
		centerY = 50;

		canvas.setLineWidth(1);
		canvas.setStrokeStyle(Color.BLACK);
		canvas.beginPath();
		canvas.moveTo(1, 1);
		canvas.lineTo(1, height - 1);
		canvas.lineTo(width - 1, height - 1);
		canvas.lineTo(width - 1, 1);
		canvas.closePath();

		canvas.stroke();

	}

	private static int upMostPos(Branch root) {
		Branch node = root;
		while (node.getUp() != null) {
			node = node.getUp();
		}
		return (int) node.getUpRight().getY();
	}

	private static int downMostPos(Branch root) {
		Branch node = root;
		while (node.getDown() != null) {
			node = node.getDown();
		}
		return (int) node.getDownRight().getY();
	}

	private static int downMostPos2(Branch root) {
		Branch node = root;
		while (node.getDown() != null
				& node.getDown().getDownLeft().getX() < Cutting.cuttingHeight) {
			node = node.getDown();
		}
		return (int) node.getDownRight().getY();
	}

	private static int upMostPos2(Branch root) {
		Branch node = root;
		while (node.getUp() != null
				& node.getUp().getUpLeft().getX() < Cutting.cuttingHeight) {
			node = node.getUp();
		}
		return (int) node.getUpRight().getY();
	}

}
