package examples;

import java.util.function.UnaryOperator;

import aliceinnets.python.jyplot.JyPlot;
import de.labathome.Cubature;
import de.labathome.CubatureError;

public class PlotEvalPositions {



	public static void plotEvalPositions() {

		double[] xmin = { -1.1, -1.1 };
		double[] xmax = {  1.1,  1.1 };

		JyPlot plt = new JyPlot();
		plt.figure();

		plt.write("circle1 = plt.Circle((0, 0), 1.0, color='r')");
		plt.write("fig = plt.gcf()");
		plt.write("ax = fig.gca()");
		plt.write("ax.add_patch(circle1)");

		UnaryOperator<double[][]> nDimUnitSphere = (double[][] x) -> {
			int nDim = x.length;
			int nPoints = x[0].length;
			double[][] fval = new double[1][nPoints];

			System.out.println("eval at "+nPoints+" pos");

			double radius;
			for (int i=0; i<nPoints; ++i) {

				if (nPoints < 1000) {
					plt.plot(x[0], x[1], "k.");
				}

				radius = 0.0;
				for (int dim=0; dim<nDim; ++dim) {
					radius += x[dim][i] * x[dim][i];
				}
				radius = Math.sqrt(radius);

				if (radius < 1.0) {
					fval[0][i] += 1.0;
				}

			}

			return fval;
		};

		double[][] val_err = Cubature.integrate(nDimUnitSphere,
				xmin, xmax,
				0.01, 0.0, CubatureError.INDIVIDUAL,
				0);

		double val = val_err[0][0];
		System.out.println("unit disk area: " +  val);

		plt.axis("equal");

		plt.show();
		plt.exec();
	}

	public static void main(String[] args) {
		plotEvalPositions();
	}

}
