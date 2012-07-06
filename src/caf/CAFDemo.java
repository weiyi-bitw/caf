package caf;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import obj.DataFile;
import obj.DataFileD;
import obj.ValIdx;
import obj.ValIdxD;
import worker.ITComputer;

public class CAFDemo {

	private static double[] getWeightedMetaGene(double[][] data, double[] w, double power, int m, int n){
		double[] out = new double[n];
		double sum = 0;
		for(int i = 0; i < m; i++){
			if(w[i] > 0){
				double f = Math.exp(power*Math.log(w[i]));
				sum += f;
				for(int j = 0; j < n; j++){
					out[j] += data[i][j] * f;
				}
			}
		}
		for(int j = 0; j < n; j++){
			out[j] /= sum;
		}
		return out;
	}
	private static float calcMSE(float[] a, float[] b, int n){
		float err = 0;
		for(int i = 0; i < n; i++){
			err += (a[i] - b[i]) * (a[i] - b[i]);
		}
		return err / n;
	}
	private static double calcMSE(double[] a, double[] b, int n){
		double err = 0;
		for(int i = 0; i < n; i++){
			err += (a[i] - b[i]) * (a[i] - b[i]);
		}
		return err / n;
	}
	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		String datafile = args[0];
		String seed = "CENPA";
		double pow = Double.parseDouble(args[1]);
		int bins = Integer.parseInt(args[2]);
		int so = Integer.parseInt(args[3]);
		final double precision = 1E-4f;
		int maxIter = Integer.parseInt(args[4]);
		final int rank = 20;
		DecimalFormat df = new DecimalFormat("0.0000"); 
		System.out.println("Loading files...");
		
		DataFileD ma = DataFileD.parse(datafile);
		
		System.out.println("==Correlation Attractor Finder==========Wei-Yi Cheng==wc2302@columbia.edu======\n");
		System.out.println("-- CAF DEMO");
		System.out.printf("%-25s%s\n", "Expression Data:", datafile);
		System.out.println("\n==DEMO default setting=========================================================\n");
		System.out.printf("%-25s%s\n", "Weight Power:",pow);
		System.out.printf("%-25s%s\n", "Bins:",bins);
		System.out.printf("%-25s%s\n", "Spline order:",so);
		System.out.printf("%-25s%s\n", "Converge threshold:", precision);
		System.out.printf("%-25s%s\n", "Max Iterations:", maxIter);
		System.out.println("\n===============================================================================\n");
		
		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
		System.out.print("Please enter the seed gene (Default=CENPA):");
		String in = br.readLine();
		if(!in.equals("")){
			seed = in.toUpperCase();
		}
		System.out.println("\n===============================================================================\n");
		System.out.printf("%-25s%s\n", "Seed:", seed);
		System.out.println("\n===============================================================================\n");
		
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		double[][] data = ma.getData();
		HashMap<String, Integer> geneMap = ma.getRows();
		ArrayList<String> geneNames = ma.getProbes();
		
		ITComputer itc = new ITComputer(bins, so, 0, 1, true);
		itc.negateMI(true);
		
		double convergeTh = 5E-14;
		
		while(!seed.equals("")){
		
			if(geneNames.contains(seed)){
				
				int idx = geneMap.get(seed);
				double[] vec = data[idx];
				
				double[] wVec = itc.getAllDoubleMIWith(vec, data);
				double[] preWVec = new double[m];
				System.arraycopy(wVec, 0, preWVec, 0, m);
				
			// initial output
				ValIdxD[] out = new ValIdxD[m];
				for(int i = 0; i < m; i++){
					out[i] = new ValIdxD(i, wVec[i]);
				}
				Arrays.sort(out);
				int cnt = 0;
				System.out.println("Iteration " + cnt);
				System.out.printf("Rank\t%-15s\t%s\n", "Gene", "MI");
				for(int i = 0; i < rank; i++){
					System.out.printf("%s\t%-15s\t%s\n", (i+1), geneNames.get(out[i].idx) , df.format(out[i].val));
				}
				boolean converge = false;
				while(cnt < maxIter){
					double[] metaGene = getWeightedMetaGene(data, wVec, pow,  m, n);
					wVec = itc.getAllDoubleMIWith(metaGene, data);
				// output
					out = new ValIdxD[m];
					for(int i = 0; i < m; i++){
						out[i] = new ValIdxD(i, wVec[i]);
					}
					Arrays.sort(out);
					double err = calcMSE(wVec, preWVec, m);
					System.out.println("\nIteration " + (cnt+1) + "\tDelta: " + err);
					System.out.printf("Rank\t%-15s\t%s\n", "Gene", "MI");
					for(int i = 0; i < rank; i++){
						System.out.printf("%s\t%-15s\t%s\n", (i+1), geneNames.get(out[i].idx) , df.format(out[i].val));
					}
					
					if(err < convergeTh){
						System.out.println("Converged.");
						converge = true;
						break;
					}
					System.arraycopy(wVec, 0, preWVec, 0, m);
					cnt++;
				}
				
				if(converge){
					String outFile = seed + "_attractor.txt";
					System.out.println("\nAttractor was written to file " + outFile);
					PrintWriter pw = new PrintWriter(new FileWriter(outFile));
					pw.println("Rank\tGene\tMI");
					for(int i = 0; i < m; i++){
						pw.println((i+1) + "\t" + geneNames.get(out[i].idx) + "\t" + (out[i].val));
					}
					pw.close();	
				}else{
					System.out.println("Not converged.");
				}
				
				
			}else{
				System.out.println("Dataset does not contain seed gene " + seed + "!!");
			}
			System.out.print("\nPress < Enter > to exit, or enter another gene: ");
			in = br.readLine();
			seed = in;
			if(!seed.equals("")){
				System.out.println("\n===============================================================================\n");
				System.out.printf("%-25s%s\n", "Seed:", seed);
				System.out.println("\n===============================================================================\n");
			}
		}
		System.out.println("Thank you :)");
	}

}
