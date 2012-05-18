package caf;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import obj.DataFile;
import obj.ValIdx;
import worker.ITComputer;

public class PlayWithExp {
	private static int segment = 0;
	private static int numSegments = 1;
	private static long jobID;
	
	private static void SystemConfiguration() throws Exception{
		String numSegmentsProperty = System.getProperty("SGE_TASK_LAST");
		if (numSegmentsProperty != null && numSegmentsProperty.length() > 0) {
            try {
                numSegments = Integer.parseInt(numSegmentsProperty);
            } catch (NumberFormatException nfe) {
                System.out.println("WARNING: Couldn't parse number of segments property SGE_TASK_LAST=" + numSegmentsProperty + ", use 1");
            }
        }
		String jobIDProperty = System.getProperty("JOB_ID");
        jobID = System.currentTimeMillis();
        if (jobIDProperty != null && jobIDProperty.length() > 0) {
            try {
                jobID = Integer.parseInt(jobIDProperty);
            } catch (NumberFormatException nfe) {
                System.out.println("WARNING: Couldn't parse job id JOB_ID=" + jobIDProperty + ", use current time");
                jobID = System.currentTimeMillis();
            }
        }else if(numSegments != 1){
        	throw new RuntimeException("No job ID is assigned!!");
        }
        System.out.printf("%-25s%s\n", "Job ID: ",jobID);
		String segmentProperty = System.getProperty("SGE_TASK_ID");
        if (segmentProperty != null && segmentProperty.length() > 0) {
            try {
                segment = Integer.parseInt(segmentProperty) - 1;
            } catch (NumberFormatException nfe) {
                System.out.println("WARNING: Couldn't parse segment property SGE_TASK_ID=" + segmentProperty + ", use 0");
            }
        }
        System.out.println("Total Segments: " + numSegments + "\tThis Segment: " + segment);
        
		
	}
	private static float[] getWeightedMetaGene(float[][] data, double[] w, float power, int m, int n){
		float[] out = new float[n];
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
		SystemConfiguration();
		int pstart = 0;
		int pend = 20;
		float pdel = 0.1f;
		int tasks = (int) ((pend - pstart) / pdel);
		
		String seed = args[1].toUpperCase();
		
		int bins = 6;
		int so = 3;
		int maxIter = 100;
		
		DataFile ma = DataFile.parse(args[0]);
		
		System.out.println("DataFile: " + args[0]);
		
		int m = ma.getNumRows();
		int n = ma.getNumCols();
		float[][] data = ma.getData();
		HashMap<String, Integer> geneMap = ma.getRows();
		ArrayList<String> geneNames = ma.getProbes();
		if(!geneNames.contains(seed)){
			System.out.println("Dataset does not contain seed gene " + seed + "!!");
			System.exit(1);
		}
		ITComputer itc = new ITComputer(bins, so, 0, 1, true);
		
		double convergeTh = 5E-14;
		int start = tasks * segment / numSegments;
		int end = tasks * (segment+1)/numSegments;
		new File("output").mkdir();
		new File("output/" + jobID).mkdir();
		DecimalFormat df = new DecimalFormat("00.00");
		for(int cc = start; cc < end; cc++){
			float pow = pstart + pdel * cc;
			System.out.print("Power: " + pow + "...");
				
				int idx = geneMap.get(seed);
				float[] vec = data[idx];
				
				double[] wVec = itc.getAllDoubleMIWith(vec, data);
				double[] preWVec = new double[m];
				System.arraycopy(wVec, 0, preWVec, 0, m);
				
			// initial output
				ValIdx[] out = new ValIdx[m];
				for(int i = 0; i < m; i++){
					out[i] = new ValIdx(i, (float) wVec[i]);
				}
				Arrays.sort(out);
				int cnt = 0;
				boolean converge = false;
				while(cnt < maxIter){
					float[] metaGene = getWeightedMetaGene(data, wVec, pow,  m, n);
					wVec = itc.getAllDoubleMIWith(metaGene, data);
					double err = calcMSE(wVec, preWVec, m);
					if(err < convergeTh){
						System.out.println("Converged.");
						out = new ValIdx[m];
						for(int i = 0; i < m; i++){
							out[i] = new ValIdx(i, (float) wVec[i]);
						}
						converge = true;
						break;
					}
					System.arraycopy(wVec, 0, preWVec, 0, m);
					cnt++;
				}
				
				if(converge){
					Arrays.sort(out);
					String outFile = "output/" + jobID + "/" + seed + "_" + df.format(pow) + "_attractor.txt";
					PrintWriter pw = new PrintWriter(new FileWriter(outFile));
					pw.println("Rank\tGene\tMI");
					for(int i = 0; i < m; i++){
						pw.println((i+1) + "\t" + geneNames.get(out[i].idx) + "\t" + out[i].val);
					}
					pw.close();	
				}else{
					System.out.println("Not converged.");
				}
			
		}
		
	}

}
