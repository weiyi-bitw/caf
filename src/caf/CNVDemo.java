package caf;

import java.util.ArrayList;
import java.util.HashMap;

import obj.DataFile;
import obj.Genome;
import worker.ITComputer;

public class CNVDemo {

	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String datafile = args[0];
		String seed = args[1];
		String geneLocFile = args[2];
		final float pstart = 1f;
		final float pend = 6f;
		final float delp = 0.5f;
		final int wstart = 11;
		final int wend = 51;
		final int delw = 10;
		final int quantile = 5;
		final float precision = 1E-4f;
		final int maxIter = 100;
		
		System.out.println("==Correlation Attractor Finder============Wei-Yi Cheng==wc2302@columbia.edu===========\n");
		System.out.println("-- CAF DEMO");
		System.out.printf("%-25s%s\n", "Expression Data:", datafile);
		System.out.printf("%-25s%s\n", "Genome Data:", geneLocFile);
		System.out.printf("%-25s%s\n", "Seed:", seed);
		System.out.println("\n==DEMO default setting================================================================\n");
		System.out.printf("%-25s[%s - %s]\n", "Weight power range:",pstart, pend);
		System.out.printf("%-25s[%s - %s]\n", "Window size range:",wstart, wend);
		System.out.printf("%-25s%s\n", "Maximize MI ranking:",quantile);
		System.out.printf("%-25s%s\n", "Converge threshold:", precision);
		System.out.printf("%-25s%s\n", "Max Iterations:", maxIter);
		System.out.println("\n======================================================================================\n");
		
		System.out.println("Loading files...");
		DataFile ma = DataFile.parse(datafile);
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		HashMap<String, Integer> geneMap = ma.getRows();
		ArrayList<String> geneNames = ma.getProbes();
		
		Genome gn = Genome.parseGeneLocation(geneLocFile);
		gn.linkToDataFile(ma);
		
		ITComputer itc = new ITComputer(6, 3, 0, 1, true);
		
	}

}
